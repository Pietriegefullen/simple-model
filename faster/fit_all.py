"""Run staged fits across selected samples, replicas, modes, and variants.

The stages mirror the interactive workflow used with ``fit_sample.py``:

1. Explore the full parameter space (``--default``).
2. Refine from the best 10 saved candidates.
3. Polish from the best 3 (or 2) saved candidates.

Each restart is an independent search run.  Exploration uses a Sobol initial
design by default; later stages use differential evolution.  Checkpoints are
shared with ``fit_sample.py``, so a later stage automatically sees candidates
from every earlier attempt for the same sample, replica, fit mode, and model.
"""

import argparse
import hashlib
import json
import logging
import math
import traceback
from pathlib import Path
import checkpoint
import data
import hashing
import model
import parameters
from fit_sample import ExistingFitGoodEnough, fit


PATHWAYS = ['Hydrolysis', 'Fermentation', 'Hydro', 'Aceto', 'Homo', 'Fe3']
FIT_MODES = ('single', 'split')
STAGE_SURVIVORS = (10, 3)
GOOD_ENOUGH_R2 = 0.9

STAGE_SETTINGS = (
    {
        'name': 'explore',
        'strategy': 'rand1bin',
        'popsize': 5,
        'mutation': (0.5, 1.0),
        'recombination': 0.7,
        'tol': 1e-3,
    },
    {
        'name': 'refine_top_10',
        'strategy': 'best1bin',
        'popsize': 8,
        'mutation': (0.35, 0.7),
        'recombination': 0.8,
        'tol': 1e-5,
    },
    {
        'name': 'polish_top_3',
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
        'loss_weight': {'CO2': 1.0, 'CH4': 1.0},
        'reduction': {'CO2': 'mse', 'CH4': 'mse'},
        'transform': {'CO2': ['normalize'], 'CH4': ['log', 'normalize']},
    }


def variants():
    """Return named model structures used by the batch runner."""
    thermodynamics = {
        parameter.name: False
        for parameter in parameters.default_model_parameters()
        if parameter.name.endswith('_thermodynamics')
    }
    return (
        ('default', PATHWAYS, {}),
        ('no_thermodynamics', PATHWAYS, thermodynamics),
        ('without_hydro', [pathway for pathway in PATHWAYS if pathway != 'Hydro'], {}),
    )


def chosen_config(sample, replica, fit_mode, pathways, overrides,
                  algorithm='differential_evolution'):
    return {
        'sample': str(sample.sample_name),
        'validation_replica': str(replica.replica_number),
        't_start': None,
        't_end': None,
        'fit_mode': fit_mode,
        'pathways': list(pathways),
        'parameter_override': dict(overrides),
        'normalized_parameters': True,
        'algorithm': algorithm,
    }


def model_identifier(chosen):
    """Build the structural model ID used to select restart checkpoints."""
    pathway_model = model.build_model(
        chosen['pathways'], chosen['normalized_parameters'])
    for name, value in chosen['parameter_override'].items():
        pathway_model.parameters()[name].constant(value)
    return hashing.build_model_id(pathway_model.get_config(only_structure=True))


def stage_init_config(chosen, model_id, stage_index, polish_best=None):
    """Use the repository's checkpoint loader for restart stages."""
    if stage_index == 0:
        return {'default': True}

    best_n = STAGE_SURVIVORS[stage_index - 1]
    if stage_index == 2 and polish_best is not None:
        best_n = polish_best
    return {
        'best_N': best_n,
        'sample': chosen['sample'],
        'validation_replica': chosen['validation_replica'],
        'model': model_id,
    }


def _seed_for(base_seed, *parts):
    if base_seed is None:
        return None
    digest = hashlib.blake2b(
        '|'.join(map(str, parts)).encode(), digest_size=4).digest()
    return (int(base_seed) + int.from_bytes(digest, 'big')) % 4_294_967_295


def _objective_calls_for_model_budget(sample, fit_mode, model_calls):
    """One objective evaluation runs once per fitted replica."""
    replica_count = (1 if fit_mode == 'single'
                     else max(1, len(sample.replicas) - 1))
    return int(math.ceil(model_calls / replica_count)), replica_count


def stage_algorithm_config(stage_index, model_calls, sample, fit_mode, workers,
                           seed, scipy_polish=False, initialization_only=False,
                           initialization='sobol', exploration_initializer=True,
                           exploration_call_limit=None,
                           exploration_stall_timeout=None):
    """Build the optimiser configuration for one stage.

    Exploration defaults to an exact-size space-filling design.  This avoids
    SciPy DE's Sobol population rounding and its subsequent mutation
    generations, while leaving refinement and polishing as DE stages.
    """
    settings = dict(STAGE_SETTINGS[stage_index])
    objective_calls, replicas_per_objective = _objective_calls_for_model_budget(
        sample, fit_mode, model_calls)
    settings.update({
        'workers': workers,
        'updating': 'deferred',
        'init': 'sobol',
        'polish': scipy_polish and stage_index == len(STAGE_SETTINGS) - 1,
        # The optimizer stops at the first completed generation at or beyond
        # this budget.  A generation can overshoot slightly when evaluated in
        # parallel, but a stalled run cannot consume the old unbounded SciPy
        # default of 1,000 generations.
        'max_evaluations': objective_calls,
    })
    if seed is not None:
        settings['seed'] = seed
    use_initializer = initialization_only or (
        stage_index == 0 and exploration_initializer)
    if use_initializer:
        # Do not queue a much larger initial design than the exploration call
        # limit.  This is especially important with multiprocessing, where a
        # pool can otherwise have already received the remainder of the design
        # when the interrupt is raised.
        samples = objective_calls
        if stage_index == 0 and exploration_call_limit is not None:
            samples = min(samples, exploration_call_limit)
        initial_sampling = {
            'init': initialization,
            'samples': samples,
            'seed': seed,
            'workers': workers,
        }
        if stage_index == 0 and exploration_stall_timeout is not None:
            initial_sampling['stall_timeout'] = exploration_stall_timeout
        return {'initial_sampling': initial_sampling}, samples, replicas_per_objective
    return {'differential_evolution': settings}, objective_calls, replicas_per_objective


def _serializable_r2(value):
    """Return an ordinary finite float for JSON summaries, or ``None``."""
    if value is None:
        return None
    value = float(value)
    return value if math.isfinite(value) else None


def _target_reached(r2, target):
    return target is not None and r2 is not None and r2 >= target


def run_one_fit(output, variant_name, pathways, overrides, sample, replica,
                fit_mode, stage_calls, workers, seed, scipy_polish, logger,
                initialization_only=False, initialization='sobol',
                stage_repeats=(3, 3, 1), r2_target=0.5, polish_best=3,
                continue_below_r2_target=False,
                exploration_initializer=True, explore_call_limit=490,
                explore_stall_timeout=10):
    """Run one fit through all stages, returning a JSON-ready stage summary.

    Exploration is retried until it reaches ``r2_target`` or exhausts its
    configured attempts.  Every exploration attempt is interrupted after
    ``explore_call_limit`` objective calls, before a full DE generation can
    stall.  It also terminates the initial-sampling worker pool if no result
    arrives within ``explore_stall_timeout`` seconds.  Exploration uses a
    scrambled Sobol initial design by default; refinement and polishing use
    DE.  Every attempt inherits retained checkpoints from earlier attempts.
    """
    del output  # Results are written through the checkpoint root used by fit_sample.
    chosen = chosen_config(
        sample, replica, fit_mode, pathways, overrides, 'differential_evolution')
    model_id = model_identifier(chosen)
    stage_results = []

    for stage_index, model_calls in enumerate(stage_calls):
        attempts = stage_repeats[stage_index]
        stage_name = STAGE_SETTINGS[stage_index]['name']
        reached_target = False
        stage_skipped = False
        call_limit_reached = False
        for attempt in range(1, attempts + 1):
            use_initializer = initialization_only or (
                stage_index == 0 and exploration_initializer)
            try:
                stage_seed = _seed_for(
                    seed, variant_name, sample.sample_name, replica.replica_number,
                    fit_mode, stage_index, attempt)
                algorithm_config, objective_calls, replica_count = stage_algorithm_config(
                    stage_index, model_calls, sample, fit_mode, workers, stage_seed,
                    scipy_polish, initialization_only, initialization,
                    exploration_initializer, explore_call_limit,
                    explore_stall_timeout)
                # Keep exploration bounded even when it is using an initial
                # design rather than DE.  Checkpoints from the sampled points
                # remain available to the following repeat or stage.
                call_limit = objective_calls if stage_index == 0 else None
                stage_chosen = dict(chosen)
                stage_chosen['algorithm'] = next(iter(algorithm_config))
                logger.info(
                    'starting variant=%s sample=%s replica=%s mode=%s stage=%s '
                    'attempt=%d/%d; budget=%d model evaluations '
                    '(%d objective evaluations x %d replicas)%s',
                    variant_name, sample.sample_name, replica.replica_number, fit_mode,
                    stage_name, attempt, attempts, objective_calls * replica_count,
                    objective_calls,
                    replica_count,
                    '; initialization=' + initialization if use_initializer else '')

                _, objective, _ = fit(
                    stage_chosen, objective_config(), algorithm_config,
                    stage_init_config(stage_chosen, model_id, stage_index, polish_best),
                    call_limit=call_limit,
                    skip_existing_r2=GOOD_ENOUGH_R2)
                best_r2 = _serializable_r2(objective.best_R2())
                result = {
                    'stage': stage_index + 1,
                    'name': stage_name,
                    'attempt': attempt,
                    'objective_calls': objective.call_count(),
                    'best_loss': float(objective.best_loss()),
                    'best_r2': best_r2,
                }
                stage_results.append(result)
                logger.info(
                    'finished variant=%s sample=%s replica=%s mode=%s stage=%s '
                    'attempt=%d; best_loss=%g best_r2=%s',
                    variant_name, sample.sample_name, replica.replica_number, fit_mode,
                    stage_name, attempt, result['best_loss'],
                    'n/a' if best_r2 is None else f'{best_r2:.3f}')

                if stage_index == 0 and _target_reached(best_r2, r2_target):
                    reached_target = True
                    break

            except checkpoint.StalledEvaluationInterrupt as error:
                call_limit_reached = True
                stage_results.append({
                    'stage': stage_index + 1,
                    'name': stage_name,
                    'attempt': attempt,
                    'status': 'stall_timeout',
                    'timeout_seconds': error.timeout,
                })
                logger.warning(
                    'exploration stalled for %.1f seconds: variant=%s sample=%s '
                    'replica=%s mode=%s stage=%s attempt=%d; workers terminated; '
                    'continuing with the next attempt',
                    error.timeout, variant_name, sample.sample_name,
                    replica.replica_number, fit_mode, stage_name, attempt)
                continue

            except checkpoint.CallLimitInterrupt:
                call_limit_reached = True
                stage_results.append({
                    'stage': stage_index + 1,
                    'name': stage_name,
                    'attempt': attempt,
                    'status': 'call_limit',
                })
                logger.info(
                    'call limit reached: variant=%s sample=%s replica=%s mode=%s '
                    'stage=%s attempt=%d; continuing with the next attempt',
                    variant_name, sample.sample_name, replica.replica_number, fit_mode,
                    stage_name, attempt)
                continue

            except KeyboardInterrupt:
                print()
                # Do not catch a second KeyboardInterrupt from input(): it is
                # deliberately allowed to reach run_batch and abort fit_all.
                logger.warning(
                    'interrupt during stage=%s attempt=%d for variant=%s sample=%s '
                    'replica=%s mode=%s. Enter y to skip this stage and continue '
                    'with the next one; Enter/n abandons only this attempt and '
                    'starts the next repeat; Ctrl-C now aborts the full batch.',
                    stage_name, attempt, variant_name, sample.sample_name,
                    replica.replica_number, fit_mode)
                try:
                    answer = input('Skip current stage? [y/N] ')
                except KeyboardInterrupt:
                    logger.warning(
                        'second interrupt while confirming; aborting the full batch.')
                    raise
                if answer.strip().lower() in {'y', 'yes'}:
                    stage_skipped = True
                    logger.info(
                        'variant=%s sample=%s replica=%s mode=%s stage=%s '
                        'attempt=%d skipped by user; continuing with the next stage',
                        variant_name, sample.sample_name, replica.replica_number, fit_mode,
                        stage_name, attempt)
                    stage_results.append({
                        'stage': stage_index + 1,
                        'name': stage_name,
                        'attempt': attempt,
                        'status': 'skipped_by_user',
                    })
                    break

                logger.info(
                    'variant=%s sample=%s replica=%s mode=%s stage=%s '
                    'attempt=%d interrupted; continuing with the next attempt',
                    variant_name, sample.sample_name, replica.replica_number, fit_mode,
                    stage_name, attempt)
                continue

        if stage_index == 0 and r2_target is not None and not reached_target \
                and not stage_skipped and not call_limit_reached:
            message = (
                f'exploration did not reach R2 >= {r2_target:g} after '
                f'{attempts} attempt(s) for variant={variant_name} '
                f'sample={sample.sample_name} replica={replica.replica_number} '
                f'mode={fit_mode}')
            if continue_below_r2_target:
                logger.warning('%s; continuing because requested.', message)
            else:
                raise RuntimeError(message)
        elif stage_index == 0 and not reached_target and call_limit_reached:
            logger.info(
                'exploration stopped at its call limit; continuing with the next stage')

    return stage_results


def configure_logger(output):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    logger = logging.getLogger('fit_all')
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    handler = logging.FileHandler(output / 'fit_all.log')
    handler.setFormatter(logging.Formatter('%(asctime)s %(levelname)s %(message)s'))
    logger.addHandler(handler)
    console = logging.StreamHandler()
    console.setFormatter(logging.Formatter('%(levelname)s %(message)s'))
    logger.addHandler(console)
    logger.propagate = False
    return logger


def selected_samples(dataset, requested_samples=None, excluded_samples=None):
    sample_filter = None if not requested_samples else set(map(str, requested_samples))
    excluded_filter = set(map(str, excluded_samples or ()))
    return [
        sample for sample in dataset.samples
        if sample.has_replicas()
        and (sample_filter is None or str(sample.sample_name) in sample_filter)
        and str(sample.sample_name) not in excluded_filter
    ]


def _write_summary(output, summary):
    with (output / 'run_summary.json').open('w') as summary_file:
        json.dump(summary, summary_file, indent=2)


def run_batch(args):
    output = Path(args.output)
    logger = configure_logger(output)
    dataset = data.get_data_before_day()
    samples = selected_samples(dataset, args.samples, args.exclude_samples)
    replica_filter = None if not args.replicas else set(map(str, args.replicas))
    variant_filter = set(args.variants)
    selected_variants = [item for item in variants() if item[0] in variant_filter]
    summary = {'succeeded': [], 'failed': [], 'skipped': []}

    for variant_name, variant_pathways, overrides in selected_variants:
        for sample in samples:
            for replica in sample.replicas:
                if (replica_filter is not None
                        and str(replica.replica_number) not in replica_filter):
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
                            args.workers, args.seed, args.scipy_polish, logger,
                            args.initialization_only, args.initialization,
                            args.stage_repeats, args.r2_target, args.polish_best,
                            args.continue_below_r2_target,
                            not args.explore_with_de, args.explore_call_limit,
                            args.explore_stall_timeout)
                        summary['succeeded'].append(job)
                    except ExistingFitGoodEnough as error:
                        job.update({
                            'reason': 'existing_fit_good_enough',
                            'best_r2': error.r2,
                        })
                        summary['skipped'].append(job)
                        logger.info(
                            'skipping existing good fit: variant=%s sample=%s '
                            'replica=%s mode=%s best_r2=%.3f',
                            variant_name, sample.sample_name, replica.replica_number,
                            fit_mode, error.r2)
                    except KeyboardInterrupt:
                        logger.warning('fit_all interrupted; stopping batch.')
                        return summary

                    except Exception:
                        job['traceback'] = traceback.format_exc()
                        summary['failed'].append(job)
                        logger.error(
                            'fit failed: variant=%s sample=%s replica=%s mode=%s\n%s',
                            variant_name, sample.sample_name, replica.replica_number,
                            fit_mode, job['traceback'])
                    finally:
                        _write_summary(output, summary)

    _write_summary(output, summary)
    logger.info('batch complete: %d succeeded, %d failed, %d skipped',
                len(summary['succeeded']), len(summary['failed']),
                len(summary['skipped']))
    return summary


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description='Run staged fits for selected samples, replicas, modes, and variants.')
    default_output = Path(__file__).resolve().parent.parent / 'results'
    parser.add_argument('--output', default=str(default_output),
                        help='Directory for fit_all.log and run_summary.json.')
    parser.add_argument('--fit-modes', choices=FIT_MODES, nargs='+',
                        default=list(FIT_MODES))
    parser.add_argument('--variants', choices=[variant[0] for variant in variants()],
                        nargs='+', default=[variant[0] for variant in variants()],
                        help='Model variants to run.')
    parser.add_argument('--samples', nargs='+',
                        help='Optional sample IDs; default is every sample.')
    parser.add_argument('--exclude-samples', nargs='+',
                        help='Sample IDs to omit, after applying --samples.')
    parser.add_argument('--replicas', nargs='+',
                        help='Optional replica numbers; default is every replica.')
    parser.add_argument('--stage-calls', nargs=3, type=int,
                        metavar=('EXPLORE', 'REFINE', 'POLISH'),
                        default=(1_500, 3_000, 1_500),
                        help=('Per-stage model-evaluation budgets (default: '
                              '1500 3000 1500).'))
    parser.add_argument('--stage-repeats', nargs=3, type=int,
                        metavar=('EXPLORE', 'REFINE', 'POLISH'),
                        default=(3, 3, 1),
                        help='Maximum independent search runs per stage.')
    parser.add_argument('--r2-target', type=float, default=0.5,
                        help='Exploration R2 required before refinement (default: 0.5).')
    parser.add_argument('--polish-best', type=int, choices=(2, 3), default=3,
                        help='Number of retained candidates used in stage 3.')
    parser.add_argument('--continue-below-r2-target', action='store_true',
                        help='Refine even when exploration exhausts retries below --r2-target.')
    parser.add_argument('--workers', type=int, default=-1,
                        help='Worker processes for candidate evaluation; use 1 for serial.')
    parser.add_argument('--seed', type=int,
                        help='Optional base seed; each fit/stage/attempt derives its own seed.')
    parser.add_argument('--scipy-polish', action='store_true',
                        help="Use SciPy's local polish step after the final DE stage.")
    parser.add_argument('--explore-with-de', action='store_true',
                        help=('Use DE for exploration instead of the default '
                              'exact-size Sobol initial design.'))
    parser.add_argument('--explore-call-limit', '--explore-de-call-limit',
                        dest='explore_call_limit', type=int, default=490,
                        help=('Interrupt every exploration attempt after this many '
                              'objective calls (default: 490).'))
    parser.add_argument('--explore-stall-timeout', type=float, default=10,
                        help=('Terminate exploration workers when no objective '
                              'evaluation finishes within this many seconds '
                              '(default: 10).'))
    parser.add_argument('--initialization-only', action='store_true',
                        help='Use the initial design, rather than DE, in every stage.')
    parser.add_argument('--initialization',
                        choices=('sobol', 'latinhypercube', 'halton', 'random'),
                        default='sobol',
                        help='Candidate design for initialization stages.')
    parser.add_argument('--dry-run', action='store_true',
                        help='List selected jobs in run_summary.json without fitting.')
    args = parser.parse_args(argv)

    if any(calls <= 0 for calls in args.stage_calls):
        parser.error('--stage-calls values must be positive.')
    if any(repeats <= 0 for repeats in args.stage_repeats):
        parser.error('--stage-repeats values must be positive.')
    if args.explore_call_limit <= 0:
        parser.error('--explore-call-limit must be positive.')
    if args.explore_stall_timeout <= 0:
        parser.error('--explore-stall-timeout must be positive.')
    if args.r2_target is not None and not -math.inf < args.r2_target <= 1:
        parser.error('--r2-target must be at most 1.')
    if args.scipy_polish and args.initialization_only:
        parser.error('--scipy-polish cannot be used with --initialization-only.')
    return args


if __name__ == '__main__':
    run_batch(parse_args())
