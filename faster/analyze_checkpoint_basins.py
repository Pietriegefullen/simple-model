"""Diagnose parameter--fit relationships in retained fit checkpoints.

The optimiser normally retains a collection of its best evaluated points.  This
module turns those points into two complementary diagnostics:

* ``parameter_vs_loss.png`` plots every fitted parameter against total loss.
  It is useful for seeing a monotonic sensitivity as well as a broad, flat
  parameter direction.
* ``parameter_basin_envelope.png`` divides each parameter range into equal
  intervals (on the parameter's optimisation scale) and plots the *best* loss
  reached in each interval.  An interval with a high best-loss envelope did
  not produce a good fit in the retained search, and is therefore evidence of
  a poorer local basin.  It is evidence rather than proof: another parameter
  combination might still rescue that interval.

Only checkpoints with the same model, objective, parameter ranges, sample,
replica, and fit mode are compared.  That restriction is important: a raw
loss from two different objectives or observed data sets is not comparable.

For example::

    python faster/analyze_checkpoint_basins.py --sample 1366 --replica 6 \
        --fit-mode single --run-name 1366-6-existing

The command is read-only with respect to ``results/``.  Reports are written to
``checkpoint_basin_results/<run-name>/`` unless ``--output-root`` is given.
"""

import argparse
import csv
import json
import math
import re
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib

# The script is commonly run on a compute node without a display.
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import numpy as np
from scipy.stats import spearmanr

import parameters
from USER_VARIABLES import RESULTS_DIRECTORY


DEFAULT_BINS = 8
DEFAULT_MIN_SUPPORT = 5


def _loss(candidate):
    """Return the loss in either checkpoint or compact test-record form."""
    return float(candidate.get('total_loss', candidate.get('loss')))


def loss_envelope(candidates, relative=0.05, absolute=0.0):
    """Return candidates within a loss envelope and its inclusive threshold.

    The original candidate order is retained so callers can attach the result
    to existing rows without relying on array or dictionary comparisons.
    """
    if not candidates:
        return [], math.nan
    best = min(_loss(candidate) for candidate in candidates)
    threshold = best * (1 + relative) + absolute
    return ([candidate for candidate in candidates if _loss(candidate) <= threshold],
            threshold)


def rms_distance(left, right):
    """Root-mean-square distance in already scaled parameter coordinates."""
    left = np.asarray(left, dtype=float)
    right = np.asarray(right, dtype=float)
    return float(np.sqrt(np.mean((left - right) ** 2)))


def connected_components(coordinates, threshold):
    """Label transitive single-linkage components within *threshold* RMS.

    Labels start at one.  A small, dependency-free union-find is preferable to
    a hierarchical clustering dependency here and makes the linkage rule
    explicit.
    """
    coordinates = np.asarray(coordinates, dtype=float)
    count = len(coordinates)
    parent = list(range(count))

    def find(index):
        while parent[index] != index:
            parent[index] = parent[parent[index]]
            index = parent[index]
        return index

    def union(left, right):
        left, right = find(left), find(right)
        if left != right:
            parent[right] = left

    for left in range(count):
        for right in range(left + 1, count):
            if rms_distance(coordinates[left], coordinates[right]) <= threshold:
                union(left, right)

    labels = {}
    result = []
    for index in range(count):
        root = find(index)
        if root not in labels:
            labels[root] = len(labels) + 1
        result.append(labels[root])
    return result


def _range_config(candidate, name):
    """Find a parameter's saved range, falling back to the project defaults."""
    run_config = candidate.get('run_config', {})
    config = run_config.get('range', {}).get(name)
    if config and config.get('range') is not None:
        return config
    for default in parameters.default_model_parameters():
        if default.name == name and default.is_variable():
            return {'range': [default.low, default.high], 'scale': default.scale}
    return None


def _coordinate(value, config):
    """Map one physical value to [0, 1] on its optimisation scale."""
    low, high = map(float, config['range'])
    value = float(value)
    if config.get('scale') == 'log':
        if value <= 0 or low <= 0 or high <= 0:
            return math.nan
        value, low, high = np.log10([value, low, high])
    if high == low:
        return math.nan
    return float((value - low) / (high - low))


def scaled_coordinates(candidates):
    """Return sorted variable names and their normalised coordinates.

    Log-scaled variables use their logarithmic range, matching the optimiser's
    search space rather than making a several-orders-of-magnitude parameter
    appear concentrated at zero.
    """
    names = sorted({name for candidate in candidates
                    for name, value in candidate.get('parameters', {}).items()
                    if isinstance(value, (int, float, np.number))
                    and not isinstance(value, (bool, np.bool_))
                    and _range_config(candidate, name) is not None})
    coordinates = np.full((len(candidates), len(names)), np.nan, dtype=float)
    for row, candidate in enumerate(candidates):
        for column, name in enumerate(names):
            if name in candidate.get('parameters', {}):
                config = _range_config(candidate, name)
                coordinates[row, column] = _coordinate(candidate['parameters'][name], config)
    return names, coordinates


def checkpoint_candidates(results_directory):
    """Load annotated checkpoint files below *results_directory*."""
    candidates = []
    for path in Path(results_directory).rglob('*'):
        if not path.is_file() or path.name.startswith('.'):
            continue
        try:
            with path.open() as handle:
                candidate = json.load(handle)
        except (OSError, UnicodeDecodeError, json.JSONDecodeError):
            continue
        if not isinstance(candidate, dict) or not {
                'parameters', 'total_loss', 'run_config'}.issubset(candidate):
            continue
        candidate['_path'] = str(path)
        match = re.search(r'_run-([^_]+)_', path.name)
        candidate['_run'] = match.group(1) if match else '--'
        candidates.append(candidate)
    return candidates


def _compatibility_key(candidate):
    """Serialize every setting that must agree before losses are compared."""
    config = candidate['run_config']
    chosen = config['chosen']
    comparable = {
        'sample': chosen.get('sample'),
        'replica': chosen.get('validation_replica'),
        'fit_mode': chosen.get('fit_mode'),
        'pathways': chosen.get('pathways'),
        'parameter_override': chosen.get('parameter_override', {}),
        'model': config.get('model'),
        'objective': config.get('objective'),
        'range': config.get('range'),
    }
    return json.dumps(comparable, sort_keys=True, default=str)


def select_compatible_group(candidates, sample=None, replica=None, fit_mode=None):
    """Select the largest comparable checkpoint group after user filters."""
    selected = []
    for candidate in candidates:
        chosen = candidate['run_config'].get('chosen', {})
        if sample is not None and str(chosen.get('sample')) != str(sample):
            continue
        if replica is not None and str(chosen.get('validation_replica')) != str(replica):
            continue
        if fit_mode is not None and chosen.get('fit_mode') != fit_mode:
            continue
        selected.append(candidate)
    grouped = defaultdict(list)
    for candidate in selected:
        grouped[_compatibility_key(candidate)].append(candidate)
    if not grouped:
        return []
    # Deterministic tie breaking makes repeated reports reproducible.
    _, group = max(grouped.items(), key=lambda item: (len(item[1]), item[0]))
    return sorted(group, key=_loss)


def _display_value(value, config):
    """Format a physical parameter value without pretending false precision."""
    if config.get('scale') == 'log' or abs(float(value)) < .01 or abs(float(value)) >= 1e4:
        return f'{float(value):.2g}'
    return f'{float(value):.3g}'


def _physical_interval(config, lower, upper):
    """Convert a normalised interval back to physical parameter units."""
    low, high = map(float, config['range'])
    if config.get('scale') == 'log':
        low, high = np.log10([low, high])
        values = 10 ** np.array([low + lower * (high - low),
                                 low + upper * (high - low)])
    else:
        values = np.array([low + lower * (high - low),
                           low + upper * (high - low)])
    return _display_value(values[0], config), _display_value(values[1], config)


def parameter_statistics(candidates, names, coordinates, bins=DEFAULT_BINS,
                         min_support=DEFAULT_MIN_SUPPORT):
    """Calculate rank effects, good-fit ranges, and binned loss envelopes."""
    losses = np.asarray([_loss(candidate) for candidate in candidates], dtype=float)
    log_losses = np.log10(losses)
    good_limit = np.quantile(losses, .10)
    summaries = []
    envelopes = np.full((len(names), bins), np.nan)
    supports = np.zeros((len(names), bins), dtype=int)

    for index, name in enumerate(names):
        values = coordinates[:, index]
        usable = np.isfinite(values) & np.isfinite(log_losses)
        rho, p_value = spearmanr(values[usable], log_losses[usable])
        config = _range_config(candidates[0], name)
        good_values = values[usable & (losses <= good_limit)]
        poor_values = values[usable & (losses >= np.quantile(losses, .75))]
        for bin_index in range(bins):
            lower, upper = bin_index / bins, (bin_index + 1) / bins
            in_bin = usable & (values >= lower)
            in_bin &= values <= upper if bin_index == bins - 1 else values < upper
            supports[index, bin_index] = int(np.sum(in_bin))
            if supports[index, bin_index] >= min_support:
                envelopes[index, bin_index] = np.min(losses[in_bin])
        good_low, good_high = (np.quantile(good_values, [.1, .9])
                               if len(good_values) else (np.nan, np.nan))
        poor_low, poor_high = (np.quantile(poor_values, [.1, .9])
                               if len(poor_values) else (np.nan, np.nan))
        summaries.append({
            'parameter': name,
            'scale': config.get('scale', 'linear'),
            'spearman_rho_log_loss': float(rho),
            'spearman_p_value': float(p_value),
            'good_10_90_normalized': (float(good_low), float(good_high)),
            'poor_25_10_90_normalized': (float(poor_low), float(poor_high)),
        })
    return summaries, envelopes, supports, good_limit


def plot_parameter_vs_loss(candidates, names, coordinates, summaries, output_path):
    """Create all-parameter scatter panels on a stable logarithmic loss axis."""
    losses = np.asarray([_loss(candidate) for candidate in candidates])
    columns = 4
    rows = int(math.ceil(len(names) / columns))
    figure, axes = plt.subplots(rows, columns, figsize=(4.5 * columns, 2.8 * rows),
                                sharey=True, squeeze=False, constrained_layout=True)
    log_loss = np.log10(losses)
    color = plt.cm.viridis_r(Normalize(np.min(log_loss), np.quantile(log_loss, .95))(log_loss))
    y_low, y_high = np.quantile(log_loss, [.01, .99])
    padding = max(.1, (y_high - y_low) * .08)
    for index, (name, summary) in enumerate(zip(names, summaries)):
        axis = axes.flat[index]
        x = coordinates[:, index]
        usable = np.isfinite(x)
        axis.scatter(x[usable], log_loss[usable], c=color[usable], s=22,
                     edgecolor='none', alpha=.8)
        axis.axhline(np.min(log_loss), color='#222222', linewidth=.8, linestyle='--')
        axis.set_title(f"{name}\n$\\rho$={summary['spearman_rho_log_loss']:+.2f}", fontsize=10)
        axis.set_xlim(0, 1)
        axis.set_ylim(y_low - padding, y_high + padding)
        axis.grid(axis='y', color='.9', linewidth=.7)
        if index // columns == rows - 1:
            axis.set_xlabel('parameter position in fitted range')
        if index % columns == 0:
            axis.set_ylabel(r'$\log_{10}$(total loss; lower is better)')
    for axis in axes.flat[len(names):]:
        axis.remove()
    figure.suptitle(f'Parameter positions versus fit loss ({len(candidates)} comparable checkpoints)',
                    fontsize=14)
    figure.savefig(output_path, dpi=220, bbox_inches='tight')
    plt.close(figure)


def plot_basin_envelope(candidates, names, envelopes, supports, output_path):
    """Plot each interval's best loss relative to the global best loss."""
    losses = np.asarray([_loss(candidate) for candidate in candidates])
    relative = np.log10(envelopes / np.min(losses))
    # Rows with the largest span are most helpful at the top of the image.
    spans = np.nanmax(relative, axis=1) - np.nanmin(relative, axis=1)
    order = np.argsort(np.nan_to_num(spans, nan=-np.inf))[::-1]
    figure_height = max(5, .34 * len(names) + 2.4)
    figure, axis = plt.subplots(figsize=(11, figure_height), constrained_layout=True)
    image = axis.imshow(relative[order], aspect='auto', interpolation='nearest',
                        cmap='magma', vmin=0, vmax=1.0)
    axis.set_yticks(range(len(names)), [names[index] for index in order])
    axis.set_xticks(range(envelopes.shape[1]),
                    [f'{index / envelopes.shape[1]:.2g}–{(index + 1) / envelopes.shape[1]:.2g}'
                     for index in range(envelopes.shape[1])])
    axis.set_xlabel('parameter position in fitted range')
    axis.set_title('Best-loss envelope by parameter interval\n'
                   '(0 = global best; 1 = tenfold worse; grey = too few samples)')
    for row, source_index in enumerate(order):
        for column in range(envelopes.shape[1]):
            if supports[source_index, column] > 0:
                axis.text(column, row, str(supports[source_index, column]), ha='center', va='center',
                          fontsize=7, color='white' if relative[source_index, column] > .55 else 'black')
    colorbar = figure.colorbar(image, ax=axis, pad=.015)
    colorbar.set_label(r'$\log_{10}$(best interval loss / global best loss)')
    figure.savefig(output_path, dpi=220, bbox_inches='tight')
    plt.close(figure)


def write_summary(output_path, candidates, summaries, envelopes, supports, bins):
    """Write machine-readable ranges and the strongest marginal basin warning."""
    best_loss = min(_loss(candidate) for candidate in candidates)
    with Path(output_path).open('w', newline='') as handle:
        fields = ['parameter', 'scale', 'spearman_rho_log_loss', 'spearman_p_value',
                  'good_10_percent_range', 'worst_quartile_range',
                  'worst_interval', 'worst_interval_best_loss',
                  'worst_interval_loss_ratio', 'worst_interval_support']
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for index, summary in enumerate(summaries):
            config = _range_config(candidates[0], summary['parameter'])
            good = summary['good_10_90_normalized']
            poor = summary['poor_25_10_90_normalized']
            valid = np.flatnonzero(np.isfinite(envelopes[index]))
            worst = valid[np.argmax(envelopes[index, valid])] if len(valid) else None
            def interval(values):
                if not np.all(np.isfinite(values)):
                    return '--'
                low, high = _physical_interval(config, *values)
                return f'{low}–{high}'
            if worst is None:
                worst_interval, worst_loss, worst_ratio, support = '--', math.nan, math.nan, 0
            else:
                worst_interval = interval((worst / bins, (worst + 1) / bins))
                worst_loss = float(envelopes[index, worst])
                worst_ratio = worst_loss / best_loss
                support = int(supports[index, worst])
            writer.writerow({
                'parameter': summary['parameter'],
                'scale': summary['scale'],
                'spearman_rho_log_loss': f"{summary['spearman_rho_log_loss']:.5g}",
                'spearman_p_value': f"{summary['spearman_p_value']:.5g}",
                'good_10_percent_range': interval(good),
                'worst_quartile_range': interval(poor),
                'worst_interval': worst_interval,
                'worst_interval_best_loss': f'{worst_loss:.6g}',
                'worst_interval_loss_ratio': f'{worst_ratio:.4g}',
                'worst_interval_support': support,
            })


def _report_metadata(candidates, output_path, good_limit):
    chosen = candidates[0]['run_config']['chosen']
    losses = [_loss(candidate) for candidate in candidates]
    runs = Counter(candidate['_run'] for candidate in candidates)
    payload = {
        'sample': chosen.get('sample'),
        'replica': chosen.get('validation_replica'),
        'fit_mode': chosen.get('fit_mode'),
        'checkpoint_count': len(candidates),
        'run_count': len(runs),
        'checkpoints_per_run': dict(sorted(runs.items())),
        'best_loss': min(losses),
        'loss_10th_percentile': float(np.quantile(losses, .10)),
        'good_fit_loss_limit': float(good_limit),
        'interpretation': (
            'A high best-loss envelope is evidence that the retained search did not '
            'find a global-quality solution in that parameter interval. It does not '
            'prove that every joint parameter combination in the interval is a local basin.'),
    }
    with Path(output_path).open('w') as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)


def run_analysis(results_directory, output_directory, sample=None, replica=None,
                 fit_mode=None, bins=DEFAULT_BINS, min_support=DEFAULT_MIN_SUPPORT):
    """Generate reports and return the selected compatible checkpoints."""
    candidates = select_compatible_group(
        checkpoint_candidates(results_directory), sample, replica, fit_mode)
    if len(candidates) < 2:
        raise ValueError('Need at least two compatible annotated checkpoints.')
    names, coordinates = scaled_coordinates(candidates)
    if not names:
        raise ValueError('The selected checkpoints contain no ranged numeric parameters.')
    output_directory = Path(output_directory)
    output_directory.mkdir(parents=True, exist_ok=True)
    summaries, envelopes, supports, good_limit = parameter_statistics(
        candidates, names, coordinates, bins=bins, min_support=min_support)
    plot_parameter_vs_loss(candidates, names, coordinates, summaries,
                           output_directory / 'parameter_vs_loss.png')
    plot_basin_envelope(candidates, names, envelopes, supports,
                        output_directory / 'parameter_basin_envelope.png')
    write_summary(output_directory / 'parameter_fit_summary.csv', candidates,
                  summaries, envelopes, supports, bins)
    _report_metadata(candidates, output_directory / 'report.json', good_limit)
    return candidates


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results', default=RESULTS_DIRECTORY,
                        help='checkpoint tree to inspect (default: results/)')
    parser.add_argument('--sample', help='fitted sample number')
    parser.add_argument('--replica', help='validation replica number')
    parser.add_argument('--fit-mode', choices=('single', 'split'),
                        help='fit mode to inspect')
    parser.add_argument('--run-name', default='checkpoint-analysis',
                        help='subdirectory name below --output-root')
    parser.add_argument('--output-root', default='checkpoint_basin_results',
                        help='directory for reports; results/ is never modified')
    parser.add_argument('--bins', type=int, default=DEFAULT_BINS,
                        help=f'intervals per parameter (default: {DEFAULT_BINS})')
    parser.add_argument('--min-support', type=int, default=DEFAULT_MIN_SUPPORT,
                        help=f'checkpoints required in an interval (default: {DEFAULT_MIN_SUPPORT})')
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    if args.bins < 2:
        raise ValueError('--bins must be at least two')
    if args.min_support < 1:
        raise ValueError('--min-support must be positive')
    output = Path(args.output_root) / args.run_name
    candidates = run_analysis(args.results, output, args.sample, args.replica,
                              args.fit_mode, args.bins, args.min_support)
    chosen = candidates[0]['run_config']['chosen']
    print(f'Analysed {len(candidates)} compatible checkpoints for sample '
          f"{chosen['sample']}, replica {chosen['validation_replica']}, "
          f"{chosen['fit_mode']} fit.")
    print(f'Reports written to {output}')
    return output


if __name__ == '__main__':
    main()
