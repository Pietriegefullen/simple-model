"""Plot fitted-parameter distributions for successful single-replica fits.

The command first retains the lowest-loss checkpoint for every unique
``(loss definition, model ID, sample, replica)`` single fit.  It then replays
that checkpoint, keeps it only when its *combined fit* R² is strictly greater
than ``--r2-threshold`` (0.8 by default), and draws every retained parameter
grouped by model ID.  Different models may expose different parameters; a
parameter's plot includes only models for which that parameter was fitted.

Run from the repository root::

    python faster/plot_single_parameter_distributions.py

Use ``--r2-threshold`` to change the quality filter, or ``--model-id`` one or
more times to restrict the report to matching generated model IDs.
"""

from __future__ import annotations

import argparse
import contextlib
import csv
import io
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

import collect_best
from fit_sample import run
from USER_VARIABLES import PROJECT_DIRECTORY, RESULTS_DIRECTORY


DEFAULT_R2_THRESHOLD = 0.8
PLOTS_PER_PAGE = 12


def single_checkpoints(results_directory=RESULTS_DIRECTORY, model_ids=None):
    """Return the retained best checkpoint for every requested single fit.

    ``collect_best.collect_best`` deliberately keeps one checkpoint per model,
    loss definition, fit mode, and replica.  The additional model-ID filter is
    substring based, matching the convention used by ``collect_best``.
    """
    checkpoints = [
        checkpoint for checkpoint in collect_best.collect_best(results_directory)
        if checkpoint.fit_mode == 'single'
    ]
    if model_ids:
        checkpoints = collect_best.filter_checkpoints(
            checkpoints, model_ids=model_ids)
    return sorted(checkpoints, key=lambda checkpoint: (
        checkpoint.model_id, checkpoint.loss_id, checkpoint.replica))


def fit_r2(checkpoint):
    """Replay one checkpoint and return its combined R² on its fitted replica."""
    # Replaying a fit reports optimiser-style progress.  The report itself
    # needs only a concise summary, so retain no per-checkpoint chatter.
    with contextlib.redirect_stdout(io.StringIO()):
        _, objective, run_log = run(checkpoint.run_config, checkpoint.parameters)
        return float(objective.R2(run_log))


def scored_checkpoints(checkpoints, r2_threshold=DEFAULT_R2_THRESHOLD):
    """Return ``(checkpoint, r2)`` pairs that pass the strict R² filter."""
    if not np.isfinite(r2_threshold):
        raise ValueError('--r2-threshold must be finite.')

    accepted = []
    for checkpoint in checkpoints:
        r2 = fit_r2(checkpoint)
        if np.isfinite(r2) and r2 > r2_threshold:
            accepted.append((checkpoint, r2))
    return accepted


def parameter_rows(scored):
    """Flatten qualified checkpoints into stable, CSV-ready parameter rows."""
    rows = []
    for checkpoint, r2 in scored:
        sample, replica = checkpoint.replica
        for name, value in checkpoint.parameters.get_values().items():
            if value is None or not np.isfinite(value):
                continue
            rows.append({
                'model_id': checkpoint.model_id,
                'loss_id': checkpoint.loss_id,
                'sample': sample,
                'replica': replica,
                'checkpoint_id': checkpoint.id,
                'loss': checkpoint.loss,
                'r2': r2,
                'parameter': name,
                'value': float(value),
            })
    return sorted(rows, key=lambda row: (
        row['parameter'], row['model_id'], row['loss_id'],
        row['sample'], row['replica']))


def _use_log_scale(values):
    """Use log scale only when it makes a broad positive range readable."""
    values = np.asarray(values, dtype=float)
    positive = values[values > 0]
    return len(positive) == len(values) and np.max(positive) / np.min(positive) >= 100


def _draw_parameter(axis, parameter, by_model):
    model_ids = sorted(by_model)
    values = [by_model[model_id] for model_id in model_ids]
    checkpoint_count = sum(len(model_values) for model_values in values)
    positions = np.arange(1, len(model_ids) + 1)
    boxes = axis.boxplot(
        values, positions=positions, widths=.62, patch_artist=True,
        medianprops={'color': 'black', 'linewidth': 1.3},
        flierprops={'marker': 'o', 'markersize': 3, 'markerfacecolor': 'none'},
    )
    for box in boxes['boxes']:
        box.set(facecolor='tab:blue', alpha=.55)
    # Overlaying retained fits avoids hiding the sample count behind a box.
    for position, model_id, model_values in zip(positions, model_ids, values):
        offsets = np.linspace(-.08, .08, len(model_values))
        axis.scatter(position + offsets, model_values, s=14, color='tab:blue',
                     alpha=.72, zorder=3)
    all_values = [value for model_values in values for value in model_values]
    if _use_log_scale(all_values):
        axis.set_yscale('log')
    axis.set_title(f'{parameter} (n={checkpoint_count})')
    axis.set_xticks(
        positions,
        [f'{model_id}\n(n={len(model_values)})'
         for model_id, model_values in zip(model_ids, values)],
        rotation=45, ha='right', fontsize=8)
    axis.set_ylabel('Fitted value')
    axis.grid(axis='y', color='.88', linewidth=.7)
    axis.set_axisbelow(True)


def plot_distributions(rows, output, r2_threshold):
    """Write a paginated PDF of parameter distributions grouped by model ID."""
    by_parameter = defaultdict(lambda: defaultdict(list))
    for row in rows:
        by_parameter[row['parameter']][row['model_id']].append(row['value'])
    if not by_parameter:
        raise ValueError('No parameter values passed the R² filter.')

    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    parameter_names = sorted(by_parameter)
    with PdfPages(output) as pdf:
        for start in range(0, len(parameter_names), PLOTS_PER_PAGE):
            page_parameters = parameter_names[start:start + PLOTS_PER_PAGE]
            columns = 3
            rows_count = int(np.ceil(len(page_parameters) / columns))
            figure, axes = plt.subplots(rows_count, columns,
                                        figsize=(14, 3.8 * rows_count),
                                        squeeze=False)
            for axis, parameter in zip(axes.flat, page_parameters):
                _draw_parameter(axis, parameter, by_parameter[parameter])
            for axis in axes.flat[len(page_parameters):]:
                axis.set_visible(False)
            figure.suptitle(
                'Single-replica best checkpoints by model ID '
                f'(combined fit $R^2 > {r2_threshold:g}$)')
            figure.tight_layout(rect=(0, 0, 1, .96))
            pdf.savefig(figure)
            plt.close(figure)
    return len(parameter_names)


def write_rows(path, rows):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=(
            'model_id', 'loss_id', 'sample', 'replica', 'checkpoint_id',
            'loss', 'r2', 'parameter', 'value'))
        writer.writeheader()
        writer.writerows(rows)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results-directory', default=RESULTS_DIRECTORY,
                        help='Directory containing annotated checkpoints.')
    parser.add_argument('--r2-threshold', type=float, default=DEFAULT_R2_THRESHOLD,
                        help='Keep fits with combined R² strictly above this value (default: 0.8).')
    parser.add_argument('--model-id', dest='model_ids', action='append',
                        help='Generated model-ID substring to include; repeat as needed.')
    parser.add_argument('--output', type=Path,
                        default=Path(PROJECT_DIRECTORY) / 'single_parameter_distributions.pdf')
    parser.add_argument('--rows-csv', type=Path,
                        default=Path(PROJECT_DIRECTORY) / 'single_parameter_distribution_rows.csv')
    return parser.parse_args()


def main(args=None):
    args = parse_args() if args is None else args
    checkpoints = single_checkpoints(args.results_directory, args.model_ids)
    scored = scored_checkpoints(checkpoints, args.r2_threshold)
    rows = parameter_rows(scored)
    plot_count = plot_distributions(rows, args.output, args.r2_threshold)
    write_rows(args.rows_csv, rows)
    print(f'Wrote {args.output} ({plot_count} parameter distributions).')
    print(f'Wrote {args.rows_csv}.')
    print(f'Kept {len(scored)} of {len(checkpoints)} single-replica best checkpoints '
          f'with R² > {args.r2_threshold:g}.')
    return rows


if __name__ == '__main__':
    main()
