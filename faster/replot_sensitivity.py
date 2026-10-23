"""Recreate a sensitivity influence plot from previously saved results.

This command only reads ``sensitivity.csv`` and ``metadata.json``; it does not
replay the model.  Derivative plots use a signed logarithmic (``symlog``)
horizontal axis so positive and negative central differences remain visible.

Examples
--------
    python faster/replot_sensitivity.py results/sensitivity/1352-5/cp-ABCD_CO2_final
    python faster/replot_sensitivity.py results/sensitivity/1352-5/cp-ABCD_CO2_final \
        --plot elasticity --output parameter_elasticity.png
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from sensitivity import write_plot


NUMERIC_COLUMNS = frozenset({
    'baseline_parameter', 'minus_parameter', 'plus_parameter', 'baseline_target',
    'minus_target', 'plus_target', 'minus_percent_change', 'plus_percent_change',
    'elasticity', 'absolute_percent_effect', 'derivative', 'absolute_derivative',
    'absolute_elasticity', 'absolute_target_effect',
})


def load_rows(path):
    """Load CSV values back into the types consumed by ``write_plot``."""
    with path.open(newline='') as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError(f'Sensitivity table is empty: {path}')
    for row in rows:
        for name in NUMERIC_COLUMNS & set(row):
            if row[name] in ('', None):
                row[name] = None
            else:
                row[name] = float(row[name])
    return rows


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('result_directory', type=Path,
                        help='Directory containing sensitivity.csv and metadata.json.')
    parser.add_argument('--plot', choices=('derivative', 'elasticity'),
                        help='Metric to plot; defaults to the saved choice.')
    parser.add_argument('--output', type=Path,
                        help='Output PNG path; default: parameter_influence.png in the result directory.')
    args = parser.parse_args(argv)
    result_directory = args.result_directory.resolve()
    table = result_directory / 'sensitivity.csv'
    manifest = result_directory / 'metadata.json'
    if not table.is_file():
        parser.error(f'Missing sensitivity table: {table}')
    if not manifest.is_file():
        parser.error(f'Missing sensitivity metadata: {manifest}')
    with manifest.open() as handle:
        metadata = json.load(handle)
    metric = args.plot or metadata.get('plot_metric', 'derivative')
    rows = load_rows(table)
    output = (args.output or result_directory / 'parameter_influence.png').resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    expression = metadata.get('target_expression', metadata.get('quantity', 'quantity'))
    target = f"{metadata.get('statistic', 'final')}({expression})"
    baseline = float(metadata['baseline_target'])
    write_plot(rows, target, baseline, output, metric)
    print(f'Wrote {metric} sensitivity plot to {output}')


if __name__ == '__main__':
    main()
