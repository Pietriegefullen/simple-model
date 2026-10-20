"""Refine local replica parameters after a joint multi-replica fit.

The multi-fit checkpoint supplies the chosen shared values and the local
incumbent for each replica. This command holds the shared values fixed and
runs a fresh, small-step CMA-ES search over only the remaining local
parameters. It never changes the source multi-fit or its inner checkpoints.

Examples
--------
python faster/polish_replicas.py JI7G
python faster/polish_replicas.py JI7G --replicas 13514 13526 --max-evaluations 2000
"""

import argparse
import copy
import json
from pathlib import Path

import fit_replicas
import fit_sample
import parameters
from run_checkpoint import checkpoint_id, resolve_checkpoint
from USER_VARIABLES import RESULTS_DIRECTORY


DEFAULT_SIGMA = 0.10
DEFAULT_MAX_EVALUATIONS = 3000


def _multi_metadata(checkpoint_data):
    """Validate and return the metadata required to polish a multi fit."""
    run_config = checkpoint_data.get('run_config', {})
    chosen = run_config.get('chosen', {})
    metadata = run_config.get('multi_replica')
    if chosen.get('fit_mode') != 'multi' or not isinstance(metadata, dict):
        raise ValueError('The requested checkpoint is not a multi-replica fit.')
    shared = metadata.get('shared')
    replicas = checkpoint_data.get('replicas')
    if not isinstance(shared, list) or not shared:
        raise ValueError('The multi-fit checkpoint does not declare shared parameters.')
    if not isinstance(replicas, list) or not replicas:
        raise ValueError('The multi-fit checkpoint does not contain replica incumbents.')
    if not isinstance(checkpoint_data.get('parameters'), dict):
        raise ValueError('The multi-fit checkpoint has no saved shared parameter values.')
    return run_config, metadata, shared, replicas


def select_replica_records(records, replica_ids=None):
    """Select source records by their unambiguous ``sample+replica`` IDs."""
    by_id = {}
    for record in records:
        replica_id = str(record.get('replica', ''))
        if not replica_id:
            raise ValueError('A multi-fit replica record has no replica ID.')
        if replica_id in by_id:
            raise ValueError(f'Multi-fit checkpoint repeats replica {replica_id}.')
        by_id[replica_id] = record

    if replica_ids is None:
        return list(by_id.values())
    requested = [str(replica_id).replace('/', '').replace(':', '')
                 for replica_id in replica_ids]
    missing = [replica_id for replica_id in requested if replica_id not in by_id]
    if missing:
        raise ValueError(
            'Replica(s) are not present in the selected multi fit: '
            + ', '.join(missing))
    return [by_id[replica_id] for replica_id in dict.fromkeys(requested)]


def _saved_inner_range(record):
    """Return the original inner search ranges when that checkpoint remains."""
    source = record.get('checkpoint')
    if not source:
        return None
    try:
        contents = fit_sample._requested_checkpoint(source)
    except (OSError, ValueError):
        return None
    range_config = contents.get('run_config', {}).get('range')
    return (parameters.ModelParameters(range_config)
            if isinstance(range_config, dict) else None)


def local_parameters_for_polish(source_data, record):
    """Recreate the local space, seed it, and lock the source shared values."""
    run_config, _, shared_names, _ = _multi_metadata(source_data)
    chosen = run_config['chosen']
    paths = chosen.get('pathways', fit_replicas.PATHWAYS)
    normalized = chosen.get('normalized_parameters', True)
    local = fit_replicas.build_parameter_template(paths, normalized).parameters()

    # Retain a narrowed local space from the original inner fit if it is
    # available. An archived/moved inner checkpoint is not fatal: the normal
    # declared model ranges still give a valid polishing search.
    saved_range = _saved_inner_range(record)
    if saved_range is not None:
        for saved in saved_range:
            if saved.name in local:
                local[saved.name].set(saved)

    local_values = record.get('parameters', {})
    if not isinstance(local_values, dict):
        raise ValueError(f"Replica {record.get('replica')!r} has no saved local parameters.")
    for name, value in local_values.items():
        if name in local:
            local[name].set(value)

    fixed = dict(chosen.get('parameter_override', {}))
    shared_values = source_data['parameters']
    missing_shared = [name for name in shared_names if name not in shared_values]
    if missing_shared:
        raise ValueError('Source checkpoint lacks shared values for: '
                         + ', '.join(missing_shared))
    for name, value in {**fixed, **{
            name: shared_values[name] for name in shared_names
    }}.items():
        if name not in local:
            raise ValueError(f'Fixed/shared parameter {name!r} is not in this model.')
        local[name].constant(value)
    return local


def polish_run_configuration(source_data, record):
    """Build the single-replica settings used by :func:`fit_sample.fit`."""
    source_config, _, shared_names, _ = _multi_metadata(source_data)
    source_chosen = source_config['chosen']
    shared_values = {name: source_data['parameters'][name] for name in shared_names}
    fixed = dict(source_chosen.get('parameter_override', {}))
    return {
        'sample': int(record['sample']),
        'validation_replica': int(record['replica_number']),
        't_start': source_chosen.get('t_start'),
        't_end': source_chosen.get('t_end'),
        'fit_mode': 'single',
        'pathways': list(source_chosen['pathways']),
        'parameter_override': {**fixed, **shared_values},
        'normalized_parameters': source_chosen.get('normalized_parameters', True),
        'algorithm': 'cma_es',
    }


def _polish_target(source_path, record, output_root=None):
    source_id = checkpoint_id(source_path) or Path(source_path).name
    root = Path(RESULTS_DIRECTORY if output_root is None else output_root)
    return root / 'joint_polish' / f'cp-{source_id}' / str(record['replica'])


def _write_manifest(path, source_path, selected_records, args, completed):
    path.mkdir(parents=True, exist_ok=True)
    document = {
        'source_checkpoint': str(Path(source_path).resolve()),
        'source_checkpoint_id': checkpoint_id(source_path),
        'replicas_requested': [str(record['replica']) for record in selected_records],
        'algorithm': 'cma_es',
        'sigma0': args.sigma,
        'max_evaluations': args.max_evaluations,
        'popsize': args.popsize,
        'seed': args.seed,
        'completed_replicas': completed,
    }
    with (path / 'manifest.json').open('w') as handle:
        json.dump(document, handle, indent=2, sort_keys=True)
        handle.write('\n')


def polish(source_checkpoint, replica_ids=None, sigma=DEFAULT_SIGMA,
           max_evaluations=DEFAULT_MAX_EVALUATIONS, popsize=None, seed=None,
           output_root=None, dry_run=False):
    """Polish selected source replicas and return their output directories."""
    source_path = resolve_checkpoint(source_checkpoint)
    source_data = fit_sample._requested_checkpoint(source_path)
    _, _, _, records = _multi_metadata(source_data)
    selected = select_replica_records(records, replica_ids)
    target_root = _polish_target(source_path, {'replica': '_'}, output_root).parent
    args = argparse.Namespace(sigma=sigma, max_evaluations=max_evaluations,
                              popsize=popsize, seed=seed)
    completed = []
    if not dry_run:
        _write_manifest(target_root, source_path, selected, args, completed)

    algorithm = {
        'sigma0': sigma,
        'max_evaluations': max_evaluations,
        'seed': seed,
        'workers': 1,
    }
    if popsize is not None:
        algorithm['popsize'] = popsize

    for record in selected:
        local = local_parameters_for_polish(source_data, record)
        chosen = polish_run_configuration(source_data, record)
        target = _polish_target(source_path, record, output_root)
        print(f"Polishing replica {record['replica']} with {len(local.variables())} local parameters")
        if not dry_run:
            fit_sample.fit(
                chosen, source_data['run_config']['objective'], {'cma_es': algorithm},
                {'range': local, 'parameters': copy.deepcopy(local)},
                checkpoint_keep_only_n=1, cp_target=str(target))
        completed.append(str(record['replica']))
        if not dry_run:
            _write_manifest(target_root, source_path, selected, args, completed)
    return [str(_polish_target(source_path, record, output_root)) for record in selected]


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description=('Polish local replica parameters from a saved multi-fit '
                     'checkpoint while freezing its shared values.'))
    parser.add_argument('checkpoint',
                        help='Multi-fit checkpoint path or unique full/partial cp-ID.')
    parser.add_argument('--replicas', nargs='+', metavar='SAMPLE_REPLICA',
                        help=('Full replica IDs to polish, e.g. 13514 or 1351/4; '
                              'default is all.'))
    parser.add_argument('--sigma', type=float, default=DEFAULT_SIGMA,
                        help='Initial CMA-ES step in unit-bound coordinates (default: 0.10).')
    parser.add_argument('--max-evaluations', type=int, default=DEFAULT_MAX_EVALUATIONS,
                        help='Completed CMA-ES evaluations per replica (default: 3000).')
    parser.add_argument('--popsize', type=int,
                        help="CMA-ES population size; default uses CMA-ES's dimension-aware value.")
    parser.add_argument('--seed', type=int, help='Random seed for CMA-ES.')
    parser.add_argument('--output-root',
                        help="Result root; default is this project's results directory.")
    parser.add_argument('--dry-run', action='store_true',
                        help='Validate the checkpoint and report selected replicas without fitting.')
    args = parser.parse_args(argv)
    if args.sigma <= 0:
        parser.error('--sigma must be positive.')
    if args.max_evaluations <= 0:
        parser.error('--max-evaluations must be positive.')
    if args.popsize is not None and args.popsize <= 0:
        parser.error('--popsize must be positive.')
    return args


def main(argv=None):
    args = parse_args(argv)
    outputs = polish(args.checkpoint, args.replicas, args.sigma,
                     args.max_evaluations, args.popsize, args.seed,
                     args.output_root, args.dry_run)
    print('Dry run: would write results to:' if args.dry_run else 'Polish results:')
    for output in outputs:
        print('  ' + output)


if __name__ == '__main__':
    main()
