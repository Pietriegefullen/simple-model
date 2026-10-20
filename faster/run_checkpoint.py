"""Replay an annotated checkpoint, with optional parameter and time overrides.

Unlike ``fit_sample.py --checkpoint``, this command does not start an
optimisation.  It reconstructs the checkpoint's saved model, objective,
replica split, parameter ranges, and optimiser-independent configuration, then
evaluates that parameter set once.

Examples
--------
python faster/run_checkpoint.py results/.../cp-ABCD-EFGH
python faster/run_checkpoint.py results/.../cp-ABCD-EFGH --p Hydro_v_max 0.02
python faster/run_checkpoint.py results/.../cp-ABCD-EFGH --t 5 40
"""

import argparse
import copy
from datetime import datetime, timezone
import json
import math
import re
from pathlib import Path

import parameters
from fit_sample import _requested_checkpoint, run
from user_input import parse, parse_override
from USER_VARIABLES import PROJECT_DIRECTORY, RESULTS_DIRECTORY


_NOT_SPECIFIED = object()
CHECKPOINT_ID_PATTERN = re.compile(r'(?:^|[_-])cp-([a-z0-9-]+)$', re.IGNORECASE)


def checkpoint_id(path):
    """Return the short checkpoint identifier embedded in a result filename."""
    match = CHECKPOINT_ID_PATTERN.search(Path(path).name)
    return None if match is None else match.group(1)


def resolve_checkpoint(identifier, results_directory=RESULTS_DIRECTORY):
    """Resolve a checkpoint file path or a unique full/partial checkpoint ID."""
    direct_path = Path(identifier)
    if direct_path.is_file():
        return direct_path

    root = Path(results_directory)
    if not root.is_dir():
        raise ValueError(f'Results directory does not exist: {root}')
    query = identifier.casefold().removeprefix('cp-')
    matches = [
        (found_id, path)
        for path in root.rglob('*') if path.is_file()
        for found_id in [checkpoint_id(path)]
        if found_id is not None and query in found_id.casefold()
    ]
    if not matches:
        raise ValueError(
            f'No checkpoint ID matching {identifier!r} in {results_directory}')
    if len(matches) > 1:
        suggestions = '\n'.join(
            f'  {found_id}: {path.relative_to(root)}'
            for found_id, path in matches[:12])
        remaining = '' if len(matches) <= 12 else (
            f'\n  ... and {len(matches) - 12} more')
        raise ValueError(
            f'Checkpoint ID {identifier!r} is ambiguous; use more characters:\n'
            f'{suggestions}{remaining}')
    return matches[0][1]


def replay_configuration(checkpoint_path, parameter_overrides=None,
                         t_start=_NOT_SPECIFIED, t_end=_NOT_SPECIFIED):
    """Return the saved run configuration and parameters with requested edits.

    The checkpoint itself is never changed.  Parameter overrides are added to
    the saved ``chosen['parameter_override']`` mapping, which is the same
    mechanism used by ``fit_sample`` to make a value constant for a run.
    """
    checkpoint_data = _requested_checkpoint(resolve_checkpoint(checkpoint_path))
    run_config = copy.deepcopy(checkpoint_data['run_config'])
    chosen = run_config.setdefault('chosen', {})
    saved_overrides = dict(chosen.get('parameter_override', {}))
    saved_overrides.update(parameter_overrides or {})
    chosen['parameter_override'] = saved_overrides

    if t_start is not _NOT_SPECIFIED:
        chosen['t_start'] = t_start
    if t_end is not _NOT_SPECIFIED:
        chosen['t_end'] = t_end

    return run_config, parameters.ModelParameters(checkpoint_data['parameters'])


def replay_checkpoint(checkpoint_path, parameter_overrides=None,
                      t_start=_NOT_SPECIFIED, t_end=_NOT_SPECIFIED):
    """Evaluate a checkpoint once and return its model, objective, and log."""
    run_config, initial_parameters = replay_configuration(
        checkpoint_path, parameter_overrides, t_start, t_end)
    return run(run_config, initial_parameters)


def prepare_output_directory(directory):
    """Create an empty report directory without overwriting prior outputs."""
    target = Path(directory)
    if not target.is_absolute():
        target = Path(PROJECT_DIRECTORY) / target
    if target.exists():
        if not target.is_dir():
            raise ValueError(f'Output path is not a directory: {target}')
        if any(target.iterdir()):
            raise ValueError(
                f'Output directory is not empty: {target}. '
                'Choose a new directory to preserve the existing report.')
    else:
        target.mkdir(parents=True)
    return target


def _json_metric(value):
    """Represent undefined or infinite scores as null in portable JSON."""
    if value is None or not math.isfinite(float(value)):
        return None
    return float(value)


def _json_safe(value):
    """Convert non-finite numeric values to JSON-valid, replayable strings."""
    if isinstance(value, dict):
        return {key: _json_safe(child) for key, child in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(child) for child in value]
    if isinstance(value, float) and not math.isfinite(value):
        if math.isnan(value):
            return 'NaN'
        return 'Infinity' if value > 0 else '-Infinity'
    return value


def save_replay_report(output_directory, source_path, run_config,
                       parameter_overrides, t_start, t_end, model_run):
    """Write collect_best-style figures and reproducibility metadata.

    ``collect_best`` owns the report styling and diagnostic figures.  This
    adapter only chooses names below the caller's report directory and records
    the effective configuration that produced them.
    """
    # Import lazily: a parameter-only replay configuration should not need the
    # plotting stack, and this keeps the reusable helpers lightweight.
    import matplotlib.pyplot as plt
    import collect_best
    import data

    output_directory = prepare_output_directory(output_directory)
    pathway_model, objective, run_log = model_run
    source_path = Path(source_path)
    source_data = _requested_checkpoint(source_path)
    effective_config = copy.deepcopy(run_config)
    effective_config['model'] = pathway_model.get_config(only_structure=True)
    source_id = checkpoint_id(source_path)
    report_checkpoint = collect_best.Checkpoint(
        parameters=pathway_model.parameters(),
        loss=float(objective.last_call()[-1]),
        run_config=effective_config,
        id=source_id or source_path.name,
        source=str(source_path),
    )

    dataset = data.get_data_before_carex()
    figure, _ = collect_best.plot_fit(
        report_checkpoint, dataset, model_run=model_run)
    fit_plot = output_directory / '00_fit.svg'
    figure.savefig(fit_plot, dpi=300)
    plt.close(figure)
    plot_paths = [fit_plot]

    sample = dataset[report_checkpoint.replica[0]]
    fitted_replicas = sample.get_split(
        report_checkpoint.replica[1], report_checkpoint.fit_mode)['fit']
    suptitle = collect_best.figure_suptitle(
        report_checkpoint, fitted_replicas, getattr(run_log, 'replica', None))

    def target(number, quantity):
        return output_directory / f'{number:02d}_{quantity}.svg'

    plot_paths.extend(collect_best._save_diagnostic_plots(
        report_checkpoint, model_run, str(output_directory), target, suptitle))

    result_r2 = _json_metric(objective.last_R2())
    metadata = {
        'report_type': 'checkpoint_replay',
        'created_at_utc': datetime.now(timezone.utc).isoformat(),
        'source_checkpoint': str(source_path.resolve()),
        'source_checkpoint_id': source_id,
        'source_checkpoint_parameters': source_data['parameters'],
        'effective_run_config': effective_config,
        'requested_parameter_overrides': parameter_overrides,
        'requested_time_overrides': {
            't_start': None if t_start is _NOT_SPECIFIED else t_start,
            't_end': None if t_end is _NOT_SPECIFIED else t_end,
        },
        'result': {
            'loss': float(objective.last_call()[-1]),
            'r2': result_r2,
        },
        'plots': [str(path.relative_to(output_directory)) for path in plot_paths],
    }
    metadata_path = output_directory / 'metadata.json'
    with metadata_path.open('w') as handle:
        json.dump(_json_safe(metadata), handle, indent=2, sort_keys=True,
                  allow_nan=False)
        handle.write('\n')
    return plot_paths, metadata_path


def replay_to_directory(checkpoint_path, output_directory,
                        parameter_overrides=None, t_start=_NOT_SPECIFIED,
                        t_end=_NOT_SPECIFIED):
    """Replay a checkpoint and write its standard plots and metadata."""
    source_path = resolve_checkpoint(checkpoint_path)
    run_config, initial_parameters = replay_configuration(
        source_path, parameter_overrides, t_start, t_end)
    model_run = run(run_config, initial_parameters)
    plot_paths, metadata_path = save_replay_report(
        output_directory, source_path, run_config, parameter_overrides or {},
        t_start, t_end, model_run)
    return (*model_run, plot_paths, metadata_path)


def _time_bound(value):
    parsed = parse(value)
    if not isinstance(parsed, (int, float)) or isinstance(parsed, bool):
        raise argparse.ArgumentTypeError('time bounds must be numbers or None')
    if parsed < 0:
        raise argparse.ArgumentTypeError('time bounds cannot be negative')
    return parsed


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description='Replay a saved annotated checkpoint once; do not optimise.')
    parser.add_argument('checkpoint', metavar='CHECKPOINT',
                        help='annotated checkpoint path or unique full/partial '
                             'checkpoint ID')
    parser.add_argument('--output-dir', metavar='DIRECTORY', required=True,
                        help='new or empty report directory, relative to the '
                             'project root unless absolute')
    parser.add_argument('--p', metavar=('NAME', 'VALUE'), nargs='+',
                        help='parameter/value pairs to hold fixed; accepts the '
                             'same syntax as fit_sample.py')
    parser.add_argument('--t', metavar=('T_START', 'T_END'), nargs=2,
                        type=_time_bound,
                        help='replace the checkpoint time window (use 0 for the '
                             'beginning of the experiment)')
    args = parser.parse_args(argv)

    if args.p is not None and len(args.p) % 2:
        parser.error('--p requires NAME VALUE pairs')
    args.parameter_overrides = parse_override(args.p)

    args.t_start = _NOT_SPECIFIED
    args.t_end = _NOT_SPECIFIED
    if args.t is not None:
        args.t_start, args.t_end = args.t
        # Keep the established fit_sample convention: zero represents an
        # unspecified lower bound, while an explicit upper bound of zero is
        # preserved and will naturally report no usable observations.
        if args.t_start == 0:
            args.t_start = None
        if (args.t_start is not None and args.t_end is not None
                and args.t_end <= args.t_start):
            parser.error('T_END must be greater than T_START')
    return args


def main(argv=None):
    args = parse_args(argv)
    _, objective, _, plot_paths, metadata_path = replay_to_directory(
        args.checkpoint, args.output_dir, args.parameter_overrides,
        args.t_start, args.t_end)
    loss = objective.last_call()[-1]
    r2 = objective.last_R2()
    r2_text = 'n/a' if r2 is None else f'{r2:.3f}'
    print(f'Replayed checkpoint: loss {loss:.6g}, R2 {r2_text}')
    print(f'Wrote {len(plot_paths)} plot(s) and {metadata_path}')


if __name__ == '__main__':
    main()
