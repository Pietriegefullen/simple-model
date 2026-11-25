"""Plot litre-normalised converted checkpoints with the current model v0.4.

Unlike :mod:`plot_checkpoint`, this command deliberately does not replay the
saved run configuration.  It uses only a converted checkpoint's ``parameters``
and its ``unit_conversion.replica`` provenance, builds the current v0.4
six-pathway model, then delegates rendering to the standard checkpoint-report
plotter.

Example::

    python faster/plot_converted_checkpoints.py \
      converted_checkpoints/model-3T3P-ICKT/*__liter__1351_4
"""

from __future__ import annotations

import argparse
import copy
import csv
import json
from pathlib import Path

import checkpoint_report as report
import converted_checkpoint
import data
import plot_checkpoint
from USER_VARIABLES import PLOT_DIRECTORY


def load_converted_checkpoint(path):
    """Load and minimally validate a converted checkpoint JSON document."""
    path = Path(path)
    with path.open() as handle:
        document = json.load(handle)
    metadata = converted_checkpoint.metadata_for_document(document)
    if not isinstance(document.get('parameters'), dict) or metadata is None:
        raise ValueError(f'{path} is not a converted checkpoint.')
    replica = metadata['unit_conversion'].get('replica')
    if not isinstance(replica, dict) or not {'sample', 'replica_number'} <= set(replica):
        raise ValueError(f'{path} has no converted-checkpoint replica provenance.')
    return document


def replica_for(document, dataset):
    """Find the one experimental replica associated with a converted file."""
    metadata = converted_checkpoint.metadata_for_document(document)
    return converted_checkpoint.replica_for_conversion(
        None if metadata is None else metadata['unit_conversion'], dataset)


def replay_current_v04(document, replica, *, gas_limits=True):
    """Evaluate a converted parameter set in the standard report units."""
    metadata = converted_checkpoint.metadata_for_document(document)
    if metadata is None:
        raise ValueError('Document is not a converted checkpoint.')
    version, pathways, normalized = converted_checkpoint.replay_spec(metadata)
    pathway_model = converted_checkpoint.build_model(version, pathways, normalized=normalized)
    # A converted archive can contain exploratory optimiser snapshots.  The
    # current model's replica-scaled gas guard turns an obviously runaway
    # trajectory into a short, plotable partial trajectory instead of waiting
    # for the integration watchdog.
    pathway_model.use_gas_limits = gas_limits
    pathway_model.parameters().set(document['parameters'])
    objective = converted_checkpoint.optimizer.build_objective_function(
        pathway_model, [replica], document['run_config']['objective'])
    run_log = pathway_model.predict_with_logging(
        replica, reference='gdw', microbe_unit='mg')
    return pathway_model, objective, run_log


def safe_name(path):
    """Turn an input filename into a stable directory component."""
    return ''.join(character if character.isalnum() or character in '-._' else '-'
                   for character in Path(path).name).strip('-')


def save_plots(source, document, dataset, output_root, *, overwrite=False,
               gas_limits=True):
    """Replay one source with the same report builder as ``plot_checkpoint``."""
    replica = replica_for(document, dataset)
    pathway_model, objective, run_log = replay_current_v04(
        document, replica, gas_limits=gas_limits)
    conversion = document['unit_conversion']
    metadata = converted_checkpoint.metadata_for_document(document)
    sample = conversion['replica']['sample']
    replica_number = conversion['replica']['replica_number']
    target = Path(output_root) / safe_name(source)
    target.mkdir(parents=True, exist_ok=True)

    # Converted multi-fit checkpoints have already been expanded to one
    # replica-specific parameter set.  Present every converted replay as a
    # single fit so the standard report can use its ordinary plot contract.
    run_config = copy.deepcopy(document['run_config'])
    run_config['model'] = pathway_model.get_config(only_structure=True)
    run_config['chosen'].update({
        'sample': sample,
        'validation_replica': replica_number,
        'fit_mode': 'single',
        'pathways': list(converted_checkpoint.replay_spec(metadata)[1]),
        'normalized_parameters': False,
    })
    checkpoint = report.Checkpoint(
        parameters=pathway_model.parameters(),
        loss=float(document.get('total_loss', float('nan'))),
        run_config=run_config,
        id='converted-' + safe_name(source),
        source=str(source),
        metadata=metadata,
    )
    result = plot_checkpoint.CheckpointPlot(
        checkpoint, (pathway_model, objective, run_log))
    return plot_checkpoint.save_single(
        result, target, dataset, prepare=False, overwrite=overwrite)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('checkpoints', nargs='+', type=Path,
                        help='One or more converted checkpoint files.')
    parser.add_argument('--output-dir', type=Path,
                        default=Path(PLOT_DIRECTORY) / 'converted_checkpoints',
                        help='Parent directory for plot folders.')
    parser.add_argument('--overwrite', action='store_true',
                        help='Replace existing PNGs.')
    parser.add_argument('--no-gas-limits', action='store_true',
                        help='Disable the current model’s runaway-gas guard.')
    args = parser.parse_args(argv)

    dataset = data.get_data_before_day()
    saved = []
    failures = []
    for source in args.checkpoints:
        document = load_converted_checkpoint(source)
        try:
            paths = save_plots(source, document, dataset, args.output_dir,
                               overwrite=args.overwrite,
                               gas_limits=not args.no_gas_limits)
        except Exception as error:
            failures.append({'checkpoint': str(source), 'error': str(error)})
            print(f'{source}: could not replay ({error})')
            continue
        saved.extend(paths)
        print(f'{source}: wrote {len(paths)} plot(s)')
    print(f'Wrote {len(saved)} plot(s) below {args.output_dir}')
    if failures:
        args.output_dir.mkdir(parents=True, exist_ok=True)
        failure_path = args.output_dir / 'failed_replays.csv'
        with failure_path.open('w', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=('checkpoint', 'error'))
            writer.writeheader()
            writer.writerows(failures)
        print(f'{len(failures)} replay(s) failed; details: {failure_path}')
    return saved


if __name__ == '__main__':
    main()
