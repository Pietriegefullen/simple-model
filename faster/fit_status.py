"""Summarise fitted sample/replica combinations and their best total R².

Checkpoints store the objective loss, but not its R².  For a fixed objective,
the total R² is ``1 - loss / observed_variance``.  This module calculates that
observed variance from the data and the saved objective configuration, so it
can report every checkpoint without rerunning the ODE model.

Run from the repository root::

    python faster/fit_status.py

Use ``--results`` to inspect another checkpoint directory and ``--output`` to
save the same terminal-friendly report to a file.
"""

import argparse
import json
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

import numpy as np

import data
import hashing
import parameters
from model_variants import model_variant_from_id
from USER_VARIABLES import RESULTS_DIRECTORY


MISSING = '--'
NOT_APPLICABLE = ' '


@dataclass(frozen=True)
class FitResult:
    """The best reported R² for one saved checkpoint."""

    sample: str
    replica: str
    fit_mode: str
    variant: str
    total_r2: float
    path: Path


def checkpoint_files(results_directory):
    """Yield parseable annotated checkpoint files below *results_directory*."""
    for path in Path(results_directory).rglob('*'):
        if not path.is_file() or path.name.startswith('.'):
            continue
        try:
            with path.open() as handle:
                checkpoint = json.load(handle)
        except (OSError, UnicodeDecodeError, json.JSONDecodeError):
            continue
        if isinstance(checkpoint, dict) and {
                'parameters', 'total_loss', 'run_config'}.issubset(checkpoint):
            yield path, checkpoint


def _transform_for(replica, pool, transform_names):
    """Build the same observed-data transform as ``build_objective_function``."""
    _, values = getattr(replica, pool)()
    values = np.asarray(values)
    transform = parameters.IdentityTransform()
    for name in transform_names:
        if name == 'log':
            transform = parameters.LogTransform(transform)
        elif name == 'normalize':
            transformed = transform(values)
            finite = transformed[np.isfinite(transformed)]
            if finite.size == 0:
                raise ValueError(f'{replica} {pool} has no finite values to normalize')
            transform = parameters.Normalization(
                transform, np.min(finite), np.max(finite))
        else:
            raise ValueError(f'Unsupported objective transform {name!r}')
    return transform


def observed_variance(dataset, chosen, objective_config):
    """Return the denominator used by the objective's aggregate R².

    ``Objective.weighted_mse`` adds the per-pool mean squared deviations from
    each fitted replica's transformed observed mean.  It is independent of
    model parameters, and therefore lets us calculate R² from saved loss.
    """
    sample = dataset[str(chosen['sample'])]
    fit_replicas = sample.get_split(
        chosen['validation_replica'], chosen['fit_mode'])['fit']
    total = 0.0
    for replica in fit_replicas:
        for pool, weight in objective_config['loss_weight'].items():
            times, values = getattr(replica, pool)()
            times = np.asarray(times)
            values = np.asarray(values)
            transform = _transform_for(
                replica, pool, objective_config['transform'][pool])
            transformed = transform(values)
            usable = np.isfinite(times) & np.isfinite(transformed) & (times != 0)
            if chosen.get('t_start') is not None:
                usable &= times >= chosen['t_start']
            if chosen.get('t_end') is not None:
                usable &= times <= chosen['t_end']
            observed = transformed[usable]
            if observed.size == 0:
                raise ValueError(f'{replica} {pool} has no usable observations')
            total += float(weight) * np.mean((observed - np.mean(observed)) ** 2)
    return total


def _variant(model_id):
    """Return the short, human-readable variant label for a model ID."""
    return model_variant_from_id(model_id) or f'unknown ({model_id})'


def best_results(results_directory, dataset):
    """Return the highest total R² for every mode/variant/sample/replica.

    The cache key includes the complete objective and split configuration: two
    old runs with different windows or objective weights get the correct R²
    before they are compared.
    """
    variances = {}
    best = {}
    errors = []
    for path, checkpoint in checkpoint_files(results_directory):
        try:
            run_config = checkpoint['run_config']
            chosen = run_config['chosen']
            objective_config = run_config['objective']
            model_id = hashing.build_model_id(run_config['model'])
            cache_key = json.dumps(
                {'chosen': chosen, 'objective': objective_config}, sort_keys=True,
                default=str)
            if cache_key not in variances:
                variances[cache_key] = observed_variance(
                    dataset, chosen, objective_config)
            denominator = variances[cache_key]
            if not np.isfinite(denominator) or denominator == 0:
                raise ValueError('total observed variance is zero or non-finite')
            r2 = 1 - float(checkpoint['total_loss']) / denominator
            key = (str(chosen['fit_mode']), _variant(model_id),
                   str(chosen['sample']), str(chosen['validation_replica']))
            result = FitResult(key[2], key[3], key[0], key[1], r2, path)
            if key not in best or result.total_r2 > best[key].total_r2:
                best[key] = result
        except (KeyError, TypeError, ValueError, FloatingPointError) as error:
            errors.append(f'{path}: {error}')
    return list(best.values()), errors


def expected_replicas(dataset):
    """Return every usable sample and its replica numbers in display order."""
    return [
        (str(sample.sample_name), tuple(str(replica.replica_number)
                                        for replica in sample.replicas))
        for sample in dataset.samples if sample.has_replicas()
    ]


def format_report(results, dataset, fit_modes=None, variants=None):
    """Format a matrix of best R² values, with ``--`` for fits still missing."""
    all_replicas = expected_replicas(dataset)
    values = {(result.fit_mode, result.variant, result.sample, result.replica):
              result.total_r2 for result in results}
    available_modes = sorted({result.fit_mode for result in results})
    available_variants = sorted({result.variant for result in results})
    modes = fit_modes if fit_modes is not None else available_modes
    labels = variants if variants is not None else available_variants
    if not modes or not labels:
        return 'No annotated checkpoints found.'

    lines = [
        'Best total R² by fitted validation replica',
        'Cell values are the highest total R² across saved checkpoints.',
        f'{MISSING} = not fitted; blank = replica is not present in that sample.',
    ]
    for mode in modes:
        lines.append(f'\nfit mode: {mode}')
        for variant in labels:
            fitted = 0
            expected = sum(len(replicas) for _, replicas in all_replicas)
            lines.append(f'  model variant: {variant}')
            lines.append('  sample      1      2      3      4      5      6')
            for sample, replicas in all_replicas:
                cells = []
                for replica in map(str, range(1, 7)):
                    if replica not in replicas:
                        cells.append(NOT_APPLICABLE)
                        continue
                    r2 = values.get((mode, variant, sample, replica))
                    if r2 is None:
                        cells.append(MISSING)
                    else:
                        fitted += 1
                        cells.append(f'{r2:.3f}')
                lines.append(f'  {sample:<6} ' + ' '.join(f'{cell:>6}' for cell in cells))
            lines.append(f'  coverage: {fitted}/{expected} fitted')
    return '\n'.join(lines)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results', default=RESULTS_DIRECTORY,
                        help='Checkpoint directory to inspect (default: results).')
    parser.add_argument('--output', type=Path,
                        help='Optional text file to receive the report.')
    parser.add_argument('--fit-mode', action='append', dest='fit_modes',
                        help='Only report this fit mode; may be repeated.')
    parser.add_argument('--variant', action='append', dest='variants',
                        help='Only report this model variant; may be repeated.')
    return parser.parse_args()


def main():
    args = parse_args()
    dataset = data.get_data_before_day()
    results, errors = best_results(args.results, dataset)
    report = format_report(results, dataset, args.fit_modes, args.variants)
    print(report)
    if errors:
        print(f'\nSkipped {len(errors)} checkpoint(s) that could not be interpreted.')
        for error in errors[:5]:
            print(f'  {error}')
        if len(errors) > 5:
            print(f'  ... and {len(errors) - 5} more')
    if args.output:
        args.output.write_text(report + '\n')


if __name__ == '__main__':
    main()
