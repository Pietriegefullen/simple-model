"""Create the stand-alone, diagonal Model A/B validation-R² comparison.

By default it replays the best saved checkpoints for the requested fit mode.
Pass a fold-level CSV written by ``plot_validation_r2_distribution`` to reuse
scores already calculated. Only replica folds available for both Model A and
Model B can appear in a diagonal comparison. Scores outside [0, 1] are
deliberately cropped from the square plot and reported in its lower margin.

Run from the repository root::

    python faster/plot_validation_r2_scatter.py
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt

import compare_variants
import plot_validation_r2_distribution
from USER_VARIABLES import PROJECT_DIRECTORY, RESULTS_DIRECTORY


def paired_rows(score_rows, fit_mode):
    """Convert per-variant scores to paired A/B rows for one fit mode."""
    by_fold = {}
    for row in score_rows:
        if row.get('fit_mode', 'split') != fit_mode:
            continue
        key = (row['loss_id'], row['sample'], row['replica'])
        by_fold.setdefault(key, {})[row['variant']] = row

    rows = []
    for (_, sample, replica), scores in sorted(by_fold.items()):
        if {'A', 'B'} <= scores.keys():
            rows.append({
                'sample': sample,
                'replica': replica,
                'fit_mode': fit_mode,
                'r2_a': float(scores['A']['score']),
                'r2_b': float(scores['B']['score']),
            })
    if not rows:
        raise ValueError('No held-out replica folds contain both Model A and Model B.')
    return rows


def score_rows_from_csv(scores_csv, metric):
    """Read scores and return the Model A/B pairs sharing a held-out replica."""
    score_rows = []
    with Path(scores_csv).open(newline='') as handle:
        for row in csv.DictReader(handle):
            score_rows.append({**row, 'score': float(row[metric])})
    return score_rows


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--scores-csv', type=Path,
                        help='Optional existing fold-level score CSV to reuse.')
    parser.add_argument('--results-directory', default=RESULTS_DIRECTORY,
                        help='Checkpoint root when replaying scores.')
    parser.add_argument('--fit-mode', choices=('split', 'single'), default='split',
                        help='Fit mode to plot (default: split).')
    parser.add_argument('--metric', choices=('combined', 'CO2', 'CH4'),
                        default='combined')
    parser.add_argument('--workers', type=int, default=8,
                        help='Parallel checkpoint replays when no CSV is supplied (default: 8).')
    parser.add_argument('--output', type=Path,
                        default=Path(PROJECT_DIRECTORY) / 'validation_r2_scatter.png')
    return parser.parse_args()


def main(args=None):
    args = parse_args() if args is None else args
    if args.workers < 1:
        raise ValueError('--workers must be at least 1.')
    if args.scores_csv is not None:
        scores = score_rows_from_csv(args.scores_csv, args.metric)
    else:
        checkpoints = plot_validation_r2_distribution.variant_checkpoints(
            args.results_directory, ('A', 'B'), args.fit_mode)
        scores = plot_validation_r2_distribution.score_rows(
            checkpoints, ('A', 'B'), args.workers)
        for row in scores:
            row['score'] = row[args.metric]
    rows = paired_rows(scores, args.fit_mode)
    figure = compare_variants.plot_ab_r2_scatter(rows)
    figure.axes[0].set_title(
        f'A versus B {args.fit_mode}-fit validation ' +
        {'combined': r'combined $R^2$', 'CO2': r'$\mathrm{CO_2}\ R^2$',
         'CH4': r'$\mathrm{CH_4}\ R^2$'}[args.metric])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.output, dpi=300)
    plt.close(figure)
    print(f'Wrote {args.output} from {len(rows)} matched held-out replica fold(s).')
    return rows


if __name__ == '__main__':
    main()
