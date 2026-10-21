"""Plot held-out-replica R² distributions for all split-fit variant results.

The retained checkpoint in each fold is selected by calibration loss. This
report then replays that checkpoint and scores its *held-out* replica. Every
available split-fit fold for each requested model variant is included. Model
variants can therefore have different numbers of folds; these box plots are
descriptive distributions, not a paired A/B estimate.

Run from the repository root::

    python faster/plot_validation_r2_distribution.py
"""

from __future__ import annotations

import argparse
import csv
import contextlib
import io
from concurrent.futures import ProcessPoolExecutor
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import plot_best
import data
import optimizer
from fit_sample import run
from model_variants import model_variant_from_id
from USER_VARIABLES import PROJECT_DIRECTORY, RESULTS_DIRECTORY


METRICS = ('combined', 'CO2', 'CH4')


def variant_checkpoints(results_directory, variants, fit_mode='split'):
    """Select the lowest-calibration-loss checkpoint for every variant fold."""
    grouped = defaultdict(dict)
    for checkpoint in plot_best.collect_best(results_directory):
        if checkpoint.fit_mode != fit_mode:
            continue
        variant = model_variant_from_id(checkpoint.model_id)
        if variant not in variants:
            continue
        key = (checkpoint.loss_id, checkpoint.replica)
        previous = grouped[key].get(variant)
        if previous is None or checkpoint.loss < previous.loss:
            grouped[key][variant] = checkpoint

    selected = {variant: [] for variant in variants}
    for by_variant in grouped.values():
        for variant, checkpoint in by_variant.items():
            selected[variant].append(checkpoint)
    missing = [variant for variant, checkpoints in selected.items()
               if not checkpoints]
    if missing:
        raise ValueError('No split-fit folds found for: ' + ', '.join(missing))
    return selected


def validation_scores(checkpoint, dataset):
    """Replay a checkpoint and return its combined and per-gas held-out R²."""
    pathway_model, _, run_log = run(checkpoint.run_config, checkpoint.parameters)
    sample, replica = checkpoint.replica
    validation_replica = dataset[sample][replica]
    validation = optimizer.build_objective_function(
        pathway_model, validation_replica, checkpoint.run_config['objective'])
    result = {'combined': validation.R2(run_log)}
    for loss, _ in validation.loss_contributions():
        result[loss.pool] = loss.R2(run_log)
    return result


def _score_checkpoint(job):
    """Replay one checkpoint in a worker process and return its score row."""
    variant, checkpoint = job
    # Fit replay prints an optimiser-style progress line; it obscures the
    # distribution report, particularly when several workers run concurrently.
    with contextlib.redirect_stdout(io.StringIO()):
        dataset = data.get_data_before_day()
        scores = validation_scores(checkpoint, dataset)
    sample, replica = checkpoint.replica
    return {
        'loss_id': checkpoint.loss_id,
        'sample': sample,
        'replica': replica,
        'fit_mode': checkpoint.fit_mode,
        'variant': variant,
        **scores,
    }


def score_rows(checkpoints_by_variant, variants, workers=1):
    """Return one held-out score row per requested variant and available fold."""
    rows = []
    jobs = [
        (variant, checkpoint)
        for variant in variants
        for checkpoint in checkpoints_by_variant[variant]
    ]
    if workers == 1:
        rows = [_score_checkpoint(job) for job in jobs]
    else:
        with ProcessPoolExecutor(max_workers=workers) as executor:
            rows = list(executor.map(_score_checkpoint, jobs))
    return sorted(rows, key=lambda row: (row['variant'], row['sample'], row['replica']))


def summary_rows(rows, variants):
    """Compute conventional box-plot quantities for every metric and variant."""
    summaries = []
    for metric in METRICS:
        for variant in variants:
            values = np.asarray([
                row[metric] for row in rows if row['variant'] == variant], dtype=float)
            summaries.append({
                'metric': metric,
                'variant': variant,
                'n': len(values),
                'below_zero': np.count_nonzero(values < 0),
                'minimum': np.min(values),
                'q1': np.percentile(values, 25),
                'median': np.median(values),
                'q3': np.percentile(values, 75),
                'maximum': np.max(values),
                'mean': np.mean(values),
            })
    return summaries


def _summary_by_metric_variant(summaries):
    return {(row['metric'], row['variant']): row for row in summaries}


def plot_distribution(rows, summaries, variants, output):
    """Write three same-scale boxplot panels with Q1/median/Q3 annotations."""
    summary = _summary_by_metric_variant(summaries)
    colors = {'A': '#3975a4', 'B': '#d95f02'}

    figure, axes = plt.subplots(1, len(METRICS), figsize=(12, 5.7), sharey=True)
    for axis, metric in zip(axes, METRICS):
        values = [[row[metric] for row in rows if row['variant'] == variant]
                  for variant in variants]
        boxes = axis.boxplot(
            values, tick_labels=[f'Model {variant}' for variant in variants], widths=.55,
            patch_artist=True, showmeans=True,
            medianprops={'color': 'black', 'linewidth': 1.6},
            meanprops={'marker': 'D', 'markerfacecolor': 'white',
                       'markeredgecolor': 'black', 'markersize': 4.5},
        )
        for box, variant in zip(boxes['boxes'], variants):
            box.set(facecolor=colors.get(variant, '#777777'), alpha=.72)
        axis.axhline(0, color='.4', linewidth=.9, linestyle='--', zorder=0)
        axis.grid(axis='y', color='.88', linewidth=.8)
        axis.set_ylim(0, 1)
        axis.set_title({'combined': 'Combined objective',
                        'CO2': r'$\mathrm{CO_2}$',
                        'CH4': r'$\mathrm{CH_4}$'}[metric])
        axis.set_xlabel('Model variant')
        for position, variant in enumerate(variants, start=1):
            stats = summary[(metric, variant)]
            text = (f"n={stats['n']}\n"
                    f"below 0: {stats['below_zero']}\n"
                    f"Q1 {stats['q1']:.2f}\n"
                    f"med {stats['median']:.2f}\n"
                    f"Q3 {stats['q3']:.2f}")
            axis.text(position, -.28, text, transform=axis.get_xaxis_transform(),
                      ha='center', va='top', fontsize=8, clip_on=False,
                      bbox={'boxstyle': 'round,pad=.25', 'facecolor': 'white',
                            'edgecolor': '.7', 'alpha': .94})
    axes[0].set_ylabel(r'Held-out replicate $R^2$')
    figure.suptitle(
        'Split-fit validation: distribution across all available held-out replicas\n'
        r'Boxes span Q1--Q3; line = median; diamond = mean; whiskers use $1.5\times$IQR',
        y=.98)
    figure.subplots_adjust(left=.08, right=.99, bottom=.28, top=.80, wspace=.08)
    figure.savefig(output, dpi=300)
    plt.close(figure)


def write_csv(path, fieldnames, rows):
    with Path(path).open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results-directory', default=RESULTS_DIRECTORY)
    parser.add_argument('--variants', nargs='+', default=('A', 'B'),
                        help='Named variants to compare (default: A B).')
    parser.add_argument('--fit-mode', choices=('split', 'single'), default='split',
                        help='Fit mode to summarise (default: split).')
    parser.add_argument('--output', type=Path,
                        default=Path(PROJECT_DIRECTORY) / 'validation_r2_distribution.png')
    parser.add_argument('--scores-csv', type=Path,
                        default=Path(PROJECT_DIRECTORY) / 'validation_r2_scores.csv')
    parser.add_argument('--summary-csv', type=Path,
                        default=Path(PROJECT_DIRECTORY) / 'validation_r2_summary.csv')
    parser.add_argument('--workers', type=int, default=1,
                        help='Checkpoint replays to evaluate in parallel (default: 1).')
    return parser.parse_args()


def main(args=None):
    args = parse_args() if args is None else args
    variants = tuple(args.variants)
    if len(set(variants)) != len(variants):
        raise ValueError('Each requested variant must be listed once.')
    checkpoints_by_variant = variant_checkpoints(
        args.results_directory, variants, args.fit_mode)
    if args.workers < 1:
        raise ValueError('--workers must be at least 1.')
    rows = score_rows(checkpoints_by_variant, variants, args.workers)
    summaries = summary_rows(rows, variants)
    for path in (args.output, args.scores_csv, args.summary_csv):
        path.parent.mkdir(parents=True, exist_ok=True)
    plot_distribution(rows, summaries, variants, args.output)
    write_csv(args.scores_csv,
              ('loss_id', 'sample', 'replica', 'fit_mode', 'variant', *METRICS), rows)
    write_csv(args.summary_csv,
              ('metric', 'variant', 'n', 'below_zero', 'minimum', 'q1', 'median', 'q3',
               'maximum', 'mean'), summaries)
    print(f'Wrote {args.output}')
    print(f'Wrote {args.scores_csv}')
    print(f'Wrote {args.summary_csv}')
    print('Used ' + ', '.join(
        f'{len(checkpoints_by_variant[variant])} Model {variant} fold(s)'
        for variant in variants) + '.')
    return rows, summaries


if __name__ == '__main__':
    main()
