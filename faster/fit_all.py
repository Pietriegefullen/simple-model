"""Run staged differential-evolution fits for every available sample and replica.

Each fit is independent.  A failure is recorded in ``fit_all.log`` and the
runner continues with the next sample, replica, fit mode, and model variant.
"""

import argparse
import hashlib
import json
import logging
import math
import traceback
from pathlib import Path

import data
import hashing
import model
import parameters
from fit_sample import fit


PATHWAYS = ['Hydrolysis', 'Fermentation', 'Hydro', 'Aceto', 'Homo', 'Fe3']
FIT_MODES = ('single', 'split')
STAGE_SURVIVORS = (10, 3)

# The three stages deliberately remain differential evolution.  The strategy,
# population size, mutation, and crossover move from broad exploration to a
# more concentrated search.  ``maxfun`` is filled in for each fit at runtime.
STAGE_SETTINGS = (
    {
        #'name': 'explore',
        'strategy': 'rand1bin',
        'popsize': 5,
        'mutation': (0.5, 1.0),
        'recombination': 0.7,
        'tol': 1e-3,
    },
    {
        #'name': 'refine_top_10',
        'strategy': 'best1bin',
        'popsize': 8,
        'mutation': (0.35, 0.7),
        'recombination': 0.8,
        'tol': 1e-5,
    },
    {
        #'name': 'polish_top_3',
        'strategy': 'best1bin',
        'popsize': 5,
        'mutation': (0.2, 0.45),
        'recombination': 0.9,
        'tol': 1e-7,
    },
)


def objective_config():
    """Return the objective used by the existing ``fit_sample.py`` CLI."""
    return {
        'loss_weight': {'CO2': 1., 'CH4': 1.},
        'reduction': {'CO2': 'mse', 'CH4': 'mse'},
        'transform': {
            'CO2': ['normalize'],
            'CH4': ['log', 'normalize'],
        },
    }


def variants():
    """The requested model variants, kept explicit in run metadata."""
    thermodynamics = {
        p.name: False
        for p in parameters.default_model_parameters()
        if p.name.endswith('_thermodynamics')
    }
    return (
        ('default', PATHWAYS, {}),
        ('no_thermodynamics', PATHWAYS, thermodynamics),
        ('without_hydro', [p for p in PATHWAYS if p != 'Hydro'], {}),
    )


def chosen_config(sample, replica, fit_mode, pathways, overrides):
    return {
        'sample': str(sample.sample_name),
        'validation_replica': str(replica.replica_number),
        't_start': None,
        't_end': None,
        'fit_mode': fit_mode,
        'pathways': list(pathways),
        'parameter_override': dict(overrides),
        'normalized_parameters': True,
        'algorithm': 'differential_evolution',
    }


def model_identifier(chosen):
    """Build the structural model ID used to filter restart checkpoints."""
    initial_parameters = parameters.ModelParameters({
        p.name: p for p in parameters.default_model_parameters()
    })
    for name, value in chosen['parameter_override'].items():
        initial_parameters[name].constant(value)
    pathway_model = model.configure_model(
        chosen['pathways'], chosen['normalized_parameters'], initial_parameters)
    return hashing.build_model_id(pathway_model.get_config(only_structure=True))


def stage_init_config(chosen, model_id, stage_index):
    """Use the repository's checkpoint loader for the two restart stages."""
    if stage_index == 0:
        return {'default': True}
    return {
        'best_N': STAGE_SURVIVORS[stage_index - 1],
        'sample': chosen['sample'],
        'validation_replica': chosen['validation_replica'],
        'model': model_id,
    }


def _seed_for(base_seed, *parts):
    if base_seed is None:
        return None
    digest = hashlib.blake2b(
        '|'.join(map(str, parts)).encode(), digest_size=4).digest()
    return (int(base_seed) + int.from_bytes(digest, 'big')) % (2**32 - 1)


def _objective_calls_for_model_budget(sample, fit_mode, model_calls):
    """One objective evaluation runs once per fitted replica."""
    replica_count = 1 if fit_mode == 'single' else max(1, len(sample.replicas) - 1)
    return int(math.ceil(model_calls / replica_count)), replica_count


def stage_algorithm_config(stage_index, model_calls, sample, fit_mode, workers,
                           seed, scipy_polish):
    settings = dict(STAGE_SETTINGS[stage_index])
    objective_calls, replicas_per_objective = _objective_calls_for_model_budget(
        sample, fit_mode, model_calls)
    settings.update({
        'workers': workers,
        'updating': 'deferred',
        'init': 'sobol',
        'polish': scipy_polish and stage_index == len(STAGE_SETTINGS) - 1,
        #'maxfun': objective_calls,
    })
    if seed is not None:
        settings['seed'] = seed
    return {'differential_evolution': settings}, objective_calls, replicas_per_objective

def run_one_fit(output, variant_name, pathways, overrides, sample, replica, fit_mode,
                stage_calls, workers, seed, scipy_polish, logger):
    chosen = chosen_config(sample, replica, fit_mode, pathways, overrides)
    model_id = model_identifier(chosen)
    stage_result = []

    for stage_index, model_calls in enumerate(stage_calls):
        stage_seed = _seed_for(
            seed, variant_name, sample.sample_name, replica.replica_number,
            fit_mode, stage_index)
        algorithm_config, objective_calls, replica_count = stage_algorithm_config(
            stage_index, model_calls, sample, fit_mode, workers, stage_seed,
            scipy_polish)
        logger.info(
            'starting variant=%s sample=%s replica=%s mode=%s stage=%s; '
            'budget=%d model evaluations (%d objective evaluations × %d replicas)',
            variant_name, sample.sample_name, replica.replica_number, fit_mode,
            stage_index + 1, model_calls, objective_calls, replica_count)
        
        fit(
            chosen,
            objective_config(),
            algorithm_config,
            stage_init_config(chosen, model_id, stage_index),
        )
        stage_result.append({
            'stage': stage_index + 1,
            'name': STAGE_SETTINGS[stage_index]['name'],
            'directory': None,
        })
        logger.info(
            'finished variant=%s sample=%s replica=%s mode=%s stage=%s',
            variant_name, sample.sample_name, replica.replica_number, fit_mode,
            stage_index + 1)
    return stage_result


def configure_logger(output):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    logger = logging.getLogger('fit_all')
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    handler = logging.FileHandler(output / 'fit_all.log')
    handler.setFormatter(logging.Formatter(
        '%(asctime)s %(levelname)s %(message)s'))
    logger.addHandler(handler)
    return logger


def selected_samples(dataset, requested_samples):
    sample_filter = None if not requested_samples else set(map(str, requested_samples))
    return [sample for sample in dataset.samples
            if sample.has_replicas() and
            (sample_filter is None or str(sample.sample_name) in sample_filter)]


def run_batch(args):
    output = Path(args.output)
    logger = configure_logger(output)
    dataset = data.get_data_before_day()
    samples = selected_samples(dataset, args.samples)
    replica_filter = None if not args.replicas else set(map(str, args.replicas))
    summary = {'succeeded': [], 'failed': [], 'skipped': []}

    for variant_name, variant_pathways, overrides in variants():
        for sample in samples:
            for replica in sample.replicas:
                if (replica_filter is not None and
                        str(replica.replica_number) not in replica_filter):
                    continue
                for fit_mode in args.fit_modes:
                    job = {
                        'variant': variant_name,
                        'sample': str(sample.sample_name),
                        'replica': str(replica.replica_number),
                        'fit_mode': fit_mode,
                    }
                    if args.dry_run:
                        summary['skipped'].append(job)
                        continue
                    try:
                        job['stages'] = run_one_fit(
                            output, variant_name, variant_pathways, overrides,
                            sample, replica, fit_mode, args.stage_calls,
                            args.workers, args.seed, args.scipy_polish, logger)
                        summary['succeeded'].append(job)
                    except Exception as ex:
                        job['traceback'] = traceback.format_exc()
                        summary['failed'].append(job)
                        logger.error(
                            'fit failed: variant=%s sample=%s replica=%s mode=%s\n%s',
                            variant_name, sample.sample_name, replica.replica_number,
                            fit_mode, job['traceback'])

                    # Persist progress after every fit so a long run retains a
                    # useful failure report even if it is later interrupted.
                    with (output / 'run_summary.json').open('w') as summary_file:
                        json.dump(summary, summary_file, indent=2)

    # This also makes --dry-run useful as an exact, machine-readable job list.
    with (output / 'run_summary.json').open('w') as summary_file:
        json.dump(summary, summary_file, indent=2)
    logger.info('batch complete: %d succeeded, %d failed, %d skipped',
                len(summary['succeeded']), len(summary['failed']),
                len(summary['skipped']))
    return summary


def parse_args():
    parser = argparse.ArgumentParser(
        description='Run staged differential-evolution fits for all samples and replicas.')
    default_output = Path(__file__).resolve().parent.parent / 'results'
    parser.add_argument('--output', default=str(default_output),
                        help='Directory for stage checkpoints, log, and summary.')
    parser.add_argument('--fit-modes', choices=FIT_MODES, nargs='+',
                        default=list(FIT_MODES))
    parser.add_argument('--samples', nargs='+',
                        help='Optional sample IDs to run; default is every sample.')
    parser.add_argument('--replicas', nargs='+',
                        help='Optional replica numbers to run; default is every replica.')
    parser.add_argument('--stage-calls', nargs=3, type=int,
                        metavar=('EXPLORE', 'REFINE', 'POLISH'),
                        default=(20000, 12000, 6000),
                        help='Per-stage ODE/model-evaluation budgets (default: 20000 12000 6000).')
    parser.add_argument('--workers', type=int, default=-1,
                        help='SciPy differential-evolution workers; use 1 for serial execution.')
    parser.add_argument('--seed', type=int,
                        help='Optional base seed; each fit/stage derives its own seed.')
    parser.add_argument('--scipy-polish', action='store_true',
                        help='After the final DE stage, also use SciPy\'s local polish step.')
    parser.add_argument('--dry-run', action='store_true',
                        help='List the selected jobs in run_summary.json without fitting.')
    args = parser.parse_args()
    if any(calls <= 0 for calls in args.stage_calls):
        parser.error('--stage-calls values must be positive.')
    return args


if __name__ == '__main__':
    run_batch(parse_args())
