"""Catalogue and characterize distinct local minima for one fitted sample.

The script reads saved checkpoints (including CMA-ES restart checkpoints),
replays them to obtain a comparable combined R², then identifies basins by
single-linkage clustering in normalized parameter space.  PCA is deliberately
used only as an explanation/visualisation of the resulting minima: it is not
used to decide whether two fits are different minima.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np

import checkpoint_report
from checkpoint_io import checkpoint_id
from fit_sample import run
from USER_VARIABLES import RESULTS_DIRECTORY


def rms_distance(left, right):
    """Return Euclidean distance per dimension for comparable vectors."""
    left = np.asarray(left, dtype=float)
    right = np.asarray(right, dtype=float)
    if left.shape != right.shape or left.size == 0:
        return float('inf')
    return float(np.sqrt(np.mean((left - right) ** 2)))


def connected_components(points, tolerance):
    """Label points using transitive single-linkage components."""
    labels = [0] * len(points)
    component = 0
    for origin in range(len(points)):
        if labels[origin]:
            continue
        component += 1
        labels[origin] = component
        pending = [origin]
        while pending:
            current = pending.pop()
            for candidate in range(len(points)):
                if labels[candidate] or rms_distance(points[current], points[candidate]) > tolerance:
                    continue
                labels[candidate] = component
                pending.append(candidate)
    return labels


def accepted_candidates(candidates, loss_relative, loss_absolute):
    """Keep complete fits within an explicit loss envelope of the best one."""
    complete = [candidate for candidate in candidates
                if candidate.get('status') == 'complete'
                and np.isfinite(candidate.get('total_loss', np.nan))]
    if not complete:
        return [], float('inf')
    best = min(candidate['total_loss'] for candidate in complete)
    threshold = best * (1 + loss_relative) + loss_absolute
    return [candidate for candidate in complete if candidate['total_loss'] <= threshold], threshold


def cluster_candidates(candidates, parameter_tolerance, loss_relative,
                       loss_absolute, prediction_tolerance):
    """Annotate candidates with acceptance and basin labels.

    Parameter proximity establishes a candidate basin.  When prediction
    signatures are available, they must also be close, preventing solutions
    with different trajectories from being merged merely because their fitted
    parameters happen to be nearby.
    """
    accepted, threshold = accepted_candidates(candidates, loss_relative, loss_absolute)
    accepted_ids = {id(candidate) for candidate in accepted}
    for candidate in candidates:
        candidate['accepted'] = id(candidate) in accepted_ids
        candidate['loss_threshold'] = threshold
        candidate['basin'] = None
    if not accepted:
        return [], threshold

    coordinates = [candidate['parameter_coordinates'] for candidate in accepted]
    labels = connected_components(coordinates, parameter_tolerance)
    # Split parameter-space components if model trajectories are materially
    # different.  This retains parameter-identifiability information without
    # claiming that PCA coordinates themselves define a minimum.
    next_label = max(labels, default=0)
    for label in sorted(set(labels)):
        positions = [index for index, value in enumerate(labels) if value == label]
        signatures = [accepted[index].get('prediction_signature') for index in positions]
        if any(signature is None for signature in signatures):
            continue
        sublabels = connected_components(signatures, prediction_tolerance)
        if max(sublabels, default=1) <= 1:
            continue
        for position, sublabel in zip(positions, sublabels):
            labels[position] = label if sublabel == 1 else next_label + sublabel - 1
        next_label += max(sublabels) - 1
    for candidate, label in zip(accepted, labels):
        candidate['basin'] = int(label)
    return accepted, threshold


def parameter_coordinates(checkpoint, names):
    values = checkpoint.parameters.get_values()
    ranges = checkpoint.run_config.get('range', {})
    coordinates = []
    for name in names:
        value = float(values[name])
        bounds = ranges.get(name, {}).get('range')
        if bounds is None or bounds[1] <= bounds[0]:
            coordinates.append(value)
        else:
            coordinates.append((value - bounds[0]) / (bounds[1] - bounds[0]))
    return coordinates


def combined_r2(checkpoint):
    """Replay one checkpoint to calculate the same combined fit R²."""
    _, objective, run_log = run(checkpoint.run_config, checkpoint.parameters,
                                log_prediction=True)
    return float(objective.R2(run_log))


def fitted_candidate(args, restart, seed):
    """Load one saved minimum as an analysis candidate.

    The name reflects a candidate produced by a previous CMA-ES restart.  The
    command never reruns an expensive optimization; it analyzes durable
    checkpoint results, making a later invocation reproducible.
    """
    checkpoint = args._checkpoints[restart]
    names = args._parameter_names
    try:
        fit_r2 = combined_r2(checkpoint)
    except Exception as error:
        return {'status': 'failed', 'restart': restart, 'seed': seed,
                'source': checkpoint.source, 'error': str(error)}
    return {
        'status': 'complete',
        'restart': restart,
        'seed': seed,
        'source': checkpoint.source,
        'checkpoint_id': checkpoint_id(checkpoint.source),
        'parameters': checkpoint.parameters.get_values(),
        'total_loss': float(checkpoint.loss),
        'fit_r2': fit_r2,
        'run_config': checkpoint.run_config,
        'parameter_names': names,
        'parameter_coordinates': parameter_coordinates(checkpoint, names),
        # A later extension can add trajectory features here.  Parameter
        # clustering remains valid without them.
        'prediction_signature': None,
    }


def _json_value(value):
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, dict):
        return {str(key): _json_value(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_value(item) for item in value]
    return value


def write_summary(destination, candidates, accepted, threshold):
    rows = []
    for candidate in candidates:
        rows.append({
            'restart': candidate.get('restart'),
            'checkpoint_id': candidate.get('checkpoint_id', ''),
            'basin': candidate.get('basin', ''),
            'accepted': candidate.get('accepted', False),
            'loss': candidate.get('total_loss', ''),
            'fit_r2': candidate.get('fit_r2', ''),
            'source': candidate.get('source', ''),
        })
    with (destination / 'candidates.csv').open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys() if rows else (
            'restart', 'checkpoint_id', 'basin', 'accepted', 'loss', 'fit_r2', 'source'))
        writer.writeheader()
        writer.writerows(rows)

    basin_rows = []
    for basin in sorted({candidate['basin'] for candidate in accepted}):
        members = [candidate for candidate in accepted if candidate['basin'] == basin]
        best = min(members, key=lambda candidate: candidate['total_loss'])
        basin_rows.append({'basin': basin, 'members': len(members),
                           'best_loss': best['total_loss'], 'best_r2': best['fit_r2'],
                           'representative_checkpoint': best.get('checkpoint_id', '')})
    with (destination / 'basin_summary.csv').open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=(
            'basin', 'members', 'best_loss', 'best_r2', 'representative_checkpoint'))
        writer.writeheader()
        writer.writerows(basin_rows)
    (destination / 'analysis.json').write_text(json.dumps(_json_value({
        'loss_threshold': threshold, 'candidate_count': len(candidates),
        'accepted_count': len(accepted), 'basin_count': len(basin_rows),
    }), indent=2, sort_keys=True) + '\n')


def write_pca(destination, accepted):
    if len(accepted) < 2:
        return
    matrix = np.asarray([candidate['parameter_coordinates'] for candidate in accepted], dtype=float)
    centered = matrix - matrix.mean(axis=0)
    _, singular_values, vectors = np.linalg.svd(centered, full_matrices=False)
    variance = singular_values ** 2
    explained = variance / variance.sum() if variance.sum() else variance
    with (destination / 'pca_axes.csv').open('w', newline='') as handle:
        writer = csv.writer(handle)
        writer.writerow(['parameter'] + [f'PC{index + 1}' for index in range(len(vectors))])
        for index, name in enumerate(accepted[0]['parameter_names']):
            writer.writerow([name] + vectors[:, index].tolist())
    with (destination / 'pca_scores.csv').open('w', newline='') as handle:
        writer = csv.writer(handle)
        writer.writerow(['restart', 'basin', 'loss', 'r2'] +
                        [f'PC{index + 1}' for index in range(len(vectors))])
        for candidate, score in zip(accepted, centered @ vectors.T):
            writer.writerow([candidate['restart'], candidate['basin'], candidate['total_loss'],
                             candidate['fit_r2']] + score.tolist())
    (destination / 'pca_explained_variance.json').write_text(json.dumps(
        {'explained_variance_ratio': explained.tolist()}, indent=2) + '\n')


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('sample')
    parser.add_argument('replica')
    parser.add_argument('--results-dir', default=RESULTS_DIRECTORY)
    parser.add_argument('--output-root', default='minima')
    parser.add_argument('--run-name', default='minima')
    parser.add_argument('--restarts', type=int,
                        help='Maximum number of saved minima to analyze (default: all).')
    parser.add_argument('--loss-relative', type=float, default=.05,
                        help='Keep fits within this relative loss envelope of the best (default: .05).')
    parser.add_argument('--loss-absolute', type=float, default=0.,
                        help='Additional absolute loss envelope (default: 0).')
    parser.add_argument('--parameter-tolerance', type=float, default=.05,
                        help='RMS distance in normalized parameter space defining one basin (default: .05).')
    parser.add_argument('--prediction-tolerance', type=float, default=.05,
                        help='RMS prediction-signature distance defining one basin (default: .05).')
    parser.add_argument('--seed', type=int, default=0)
    args = parser.parse_args(argv)
    if args.restarts is not None and args.restarts <= 0:
        parser.error('--restarts must be positive.')
    for name in ('loss_relative', 'loss_absolute', 'parameter_tolerance', 'prediction_tolerance'):
        if getattr(args, name) < 0:
            parser.error('--' + name.replace('_', '-') + ' must be non-negative.')

    checkpoints = [candidate for candidate in checkpoint_report.collect_all(
        args.results_dir, verbose=False)
        if candidate.replica == (str(args.sample), str(args.replica))]
    checkpoints.sort(key=lambda candidate: candidate.loss)
    if args.restarts is not None:
        checkpoints = checkpoints[:args.restarts]
    if not checkpoints:
        raise ValueError(f'No checkpoints found for sample {args.sample}, replica {args.replica}.')
    args._checkpoints = checkpoints
    args._parameter_names = sorted({name for candidate in checkpoints
                                    for name in candidate.parameters.get_values()})
    # Only compare solutions with a complete shared parameter set.
    args._parameter_names = [name for name in args._parameter_names if all(
        name in candidate.parameters.get_values() for candidate in checkpoints)]

    destination = Path(args.output_root) / args.run_name
    (destination / 'candidates').mkdir(parents=True, exist_ok=True)
    candidates = []
    for restart in range(len(checkpoints)):
        candidate = fitted_candidate(args, restart, args.seed + restart)
        candidates.append(candidate)
        (destination / 'candidates' / f'restart-{restart:03d}.json').write_text(
            json.dumps(_json_value(candidate), indent=2, sort_keys=True) + '\n')
    accepted, threshold = cluster_candidates(
        candidates, args.parameter_tolerance, args.loss_relative,
        args.loss_absolute, args.prediction_tolerance)
    # Rewrite the durable per-restart records after assigning the quality
    # envelope and basin label, so they remain useful without opening the
    # aggregate CSV files.
    for restart, candidate in enumerate(candidates):
        (destination / 'candidates' / f'restart-{restart:03d}.json').write_text(
            json.dumps(_json_value(candidate), indent=2, sort_keys=True) + '\n')
    write_summary(destination, candidates, accepted, threshold)
    write_pca(destination, accepted)
    return str(destination)


if __name__ == '__main__':
    main()
