"""Compute fitted-parameter envelopes from selected saved checkpoints.

The selectors combine: a checkpoint must satisfy every supplied criterion.
Model IDs and excluded checkpoint IDs are case-insensitive substrings.  Site
and origin are resolved from the current dataset, so they can be used with
saved checkpoints that only record sample and replica numbers.

For example, write the envelope for the best ordinary checkpoint in each run
plus the best checkpoint retained for every CMA-ES restart in those runs::

    python faster/parameter_envelope.py --model-id UQ5A --site Samoylov \\
        --origin core --best-restart --exclude 46WF \\
        --output-directory bounds/samoylov_core
"""

from __future__ import annotations

import argparse
import fnmatch
import json
import math
from collections import defaultdict
from pathlib import Path

import checkpoint_report
import data
import parameter_bounds
from USER_VARIABLES import RESULTS_DIRECTORY


def sample_metadata(dataset):
    """Return site/origin metadata keyed by saved sample and replica IDs."""
    return {
        (str(sample.sample_name), str(replica.replica_number)): {
            'site': sample.site,
            'origin': sample.origin,
        }
        for sample in dataset.samples
        for replica in sample.replicas
    }


def _matches(value, criteria):
    """Match a case-insensitive substring, or a shell-style glob when given."""
    if criteria is None:
        return True
    value = str(value).casefold()
    for criterion in criteria:
        criterion = str(criterion).casefold()
        if ('*' in criterion or '?' in criterion) and fnmatch.fnmatchcase(value, criterion):
            return True
        if criterion in value:
            return True
    return False


def checkpoint_replicas(checkpoint):
    """Yield every sample/replica represented by a checkpoint."""
    records = getattr(checkpoint, 'multi_replica_results', ())
    if records:
        for record in records:
            sample = record.get('sample')
            replica = record.get('replica_number', record.get('replica'))
            if sample is not None and replica is not None:
                yield str(sample), str(replica)
        return
    yield tuple(map(str, checkpoint.replica))


def select_checkpoints(checkpoints, *, model_ids=None, fit_modes=None, samples=None,
                       replicas=None, sites=None, origins=None, excluded_ids=None,
                       loss_upper=1e99, metadata=None):
    """Select checkpoints matching all requested criteria.

    For joint fits, a sample/replica/site/origin criterion is satisfied when
    one included replica satisfies the complete set of replica-level filters.
    """
    selected = []
    metadata = {} if metadata is None else metadata
    for checkpoint in checkpoints:
        if checkpoint.loss > loss_upper:
            continue
        if not _matches(checkpoint.model_id, model_ids):
            continue
        if not _matches(checkpoint.fit_mode, fit_modes):
            continue
        if excluded_ids and _matches(checkpoint.id, excluded_ids):
            continue
        matching_replicas = [
            (sample, replica) for sample, replica in checkpoint_replicas(checkpoint)
            if _matches(sample, samples) and _matches(replica, replicas)
            and _matches(metadata.get((sample, replica), {}).get('site'), sites)
            and _matches(metadata.get((sample, replica), {}).get('origin'), origins)
        ]
        if matching_replicas:
            selected.append(checkpoint)
    return selected


def restart_index(checkpoint):
    """Return a CMA restart index, or ``None`` for an ordinary checkpoint."""
    restart = (getattr(checkpoint, 'metadata', None) or {}).get('cma_restart')
    if not isinstance(restart, dict):
        return None
    index = restart.get('index')
    return None if index is None else str(index)


def best_with_restarts(checkpoints):
    """Keep the best checkpoint per run and the best one per restart per run."""
    best = {}
    for checkpoint in checkpoints:
        ordinary_key = (checkpoint.run_id, None)
        current = best.get(ordinary_key)
        if current is None or checkpoint.loss < current.loss:
            best[ordinary_key] = checkpoint
        index = restart_index(checkpoint)
        if index is not None:
            restart_key = (checkpoint.run_id, index)
            current = best.get(restart_key)
            if current is None or checkpoint.loss < current.loss:
                best[restart_key] = checkpoint
    return sorted({checkpoint.source: checkpoint for checkpoint in best.values()}.values(),
                  key=lambda checkpoint: (checkpoint.model_id, checkpoint.run_id,
                                          checkpoint.loss, checkpoint.id))


def envelope_rows(checkpoints):
    """Return one minimum/maximum record for every finite fitted parameter."""
    values = defaultdict(list)
    for checkpoint in checkpoints:
        for name, value in checkpoint.parameters.get_values().items():
            if isinstance(value, bool) or not isinstance(value, (int, float)):
                continue
            value = float(value)
            if math.isfinite(value):
                values[name].append((value, checkpoint))

    rows = []
    for name in sorted(values):
        entries = values[name]
        low, low_checkpoint = min(entries, key=lambda entry: entry[0])
        high, high_checkpoint = max(entries, key=lambda entry: entry[0])
        rows.append({
            'parameter': name,
            'low': low,
            'high': high,
            'value_count': len(entries),
            'low_checkpoint_id': low_checkpoint.id,
            'high_checkpoint_id': high_checkpoint.id,
        })
    return rows


def checkpoint_records(checkpoints):
    """Return concise provenance records for every checkpoint used in an envelope."""
    return [{
        'checkpoint_id': checkpoint.id,
        'model_id': checkpoint.model_id,
        'run_id': checkpoint.run_id,
        'source': checkpoint.source,
        'sample_replicas': [f'{sample}-{replica}'
                            for sample, replica in checkpoint_replicas(checkpoint)],
    } for checkpoint in checkpoints]


def write_envelope(path, bounds_id, rows, checkpoints, args):
    """Write bounds and an auditable list of contributing checkpoints."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    document = {
        'bounds_id': bounds_id,
        'selection': {
            'model_ids': args.model_ids,
            'fit_modes': args.fit_modes,
            'samples': args.samples,
            'replicas': args.replicas,
            'sites': args.sites,
            'origins': args.origins,
            'excluded_ids': args.excluded_ids,
            'best_restart': args.best_restart,
            'loss_upper': args.loss_upper,
        },
        'checkpoint_count': len(checkpoints),
        'checkpoint_ids': [checkpoint.id for checkpoint in checkpoints],
        'checkpoints': checkpoint_records(checkpoints),
        'parameters': rows,
    }
    with path.open('w') as handle:
        json.dump(document, handle, indent=2, sort_keys=True)


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results-directory', default=RESULTS_DIRECTORY,
                        help='Checkpoint root to scan (default: results/).')
    parser.add_argument('--model-id', dest='model_ids', nargs='+', metavar='MODEL_ID',
                        help='Model-ID substring(s) to include.')
    parser.add_argument('--fit-mode', dest='fit_modes', nargs='+',
                        choices=('single', 'split', 'multi'), metavar='FIT_MODE',
                        help='Fit mode(s) to include.')
    parser.add_argument('--sample', dest='samples', nargs='+', metavar='SAMPLE',
                        help='Sample number(s) to include.')
    parser.add_argument('--replica', dest='replicas', nargs='+', metavar='REPLICA',
                        help='Replica number(s) to include.')
    parser.add_argument('--site', dest='sites', nargs='+', metavar='SITE',
                        help='Site criterion or criteria to include.')
    parser.add_argument('--origin', dest='origins', nargs='+', metavar='ORIGIN',
                        help='Origin criterion or criteria to include.')
    parser.add_argument('--exclude', dest='excluded_ids', nargs='+', metavar='CP_ID',
                        help='Checkpoint-ID substring(s) to omit.')
    parser.add_argument('--best-restart', action='store_true',
                        help=('Keep only the lowest-loss checkpoint for each run, plus '
                              'the lowest-loss checkpoint for each available CMA-ES restart.'))
    parser.add_argument('--loss-upper', type=float, default=1e99, metavar='LOSS',
                        help='Discard checkpoints with loss above this value (default: 1e99).')
    parser.add_argument('--output-directory', type=Path,
                        default=parameter_bounds.BOUNDS_DIRECTORY,
                        help='Directory for ID-named bounds files (default: bounds/).')
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    if not math.isfinite(args.loss_upper):
        raise ValueError('--loss-upper must be finite.')
    checkpoints = checkpoint_report.collect_all(args.results_directory, verbose=False)
    metadata = sample_metadata(data.get_data_before_day()) if args.sites or args.origins else {}
    selected = select_checkpoints(
        checkpoints, model_ids=args.model_ids, fit_modes=args.fit_modes,
        samples=args.samples,
        replicas=args.replicas, sites=args.sites, origins=args.origins,
        excluded_ids=args.excluded_ids, loss_upper=args.loss_upper, metadata=metadata)
    if args.best_restart:
        selected = best_with_restarts(selected)
    rows = envelope_rows(selected)
    if not rows:
        raise ValueError('No finite parameter values found in the selected checkpoints.')
    output_directory = Path(args.output_directory)
    output_directory.mkdir(parents=True, exist_ok=True)
    bounds_id = parameter_bounds.new_bounds_id()
    output = parameter_bounds.bounds_path(bounds_id, output_directory)
    while output.exists():
        bounds_id = parameter_bounds.new_bounds_id()
        output = parameter_bounds.bounds_path(bounds_id, output_directory)
    write_envelope(output, bounds_id, rows, selected, args)
    print(f'Wrote {output}: {len(rows)} parameter envelope(s) from '
          f'{len(selected)} checkpoint(s).')
    print(f'Bounds ID: {bounds_id}')
    return rows


if __name__ == '__main__':
    main()
