"""One-at-a-time parameter sensitivity around a saved checkpoint.

The command replays a single replica with the fitted checkpoint values, then
replays it with every selected numeric parameter moved down and up by a small
relative amount.  It reports the centred finite-difference derivative
``dQ/dp ~= (Q(p+h)-Q(p-h))/(2h)``.  It never modifies the source checkpoint. A joint
multi-replica checkpoint is also supported: ``--replica`` selects which saved
local solution is combined with the checkpoint's shared solution.

Examples
--------
    python faster/sensitivity.py KRQY --quantity CO2 --stat final
    python faster/sensitivity.py JOINT --replica 1351 5 --scope shared
    python faster/sensitivity.py KRQY --quantity CH4 --stat integral \
        --parameters '*_v_max' '*_CUE' --fraction 0.05
"""

from __future__ import annotations

import argparse
import ast
import copy
import csv
import fnmatch
import json
import math
import re
from datetime import datetime, timezone
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import data
import model
import parameters
from checkpoint_io import checkpoint_id, replay_inputs
from USER_VARIABLES import PROJECT_DIRECTORY


REPLICA_PROVIDED = frozenset({'H2O', 'TOC', 'DOC', 'CO2', 'CH4', 'H2', 'Fe2'})


def _replica_key(record):
    return str(record.get('sample')), str(record.get('replica_number', record.get('replica')))[-1:]


def checkpoint_solution(source, contents, requested_replica=None):
    """Return config, values, roles, and the one replica to replay."""
    run_config = copy.deepcopy(contents['run_config'])
    chosen = run_config['chosen']
    multi = chosen.get('fit_mode') == 'multi'
    shared_names = set(run_config.get('multi_replica', {}).get('shared', ()))
    values = dict(contents['parameters'])

    if not multi:
        sample = str(chosen['sample'])
        replica = str(chosen['validation_replica'])
        if requested_replica and tuple(map(str, requested_replica)) != (sample, replica):
            raise ValueError('This is a single-replica checkpoint; --replica must match '
                             f'{sample} {replica}.')
        roles = {name: 'single' for name in values}
        return run_config, values, roles, (sample, replica), None

    if not requested_replica:
        raise ValueError('A joint multi-replica checkpoint needs --replica SAMPLE REPLICA.')
    wanted = tuple(map(str, requested_replica))
    record = next((item for item in contents.get('replicas', [])
                   if isinstance(item, dict) and _replica_key(item) == wanted), None)
    if record is None:
        available = ', '.join('/'.join(_replica_key(item)) for item in contents.get('replicas', [])
                              if isinstance(item, dict))
        raise ValueError(f'Replica {wanted[0]}/{wanted[1]} is not stored in this joint checkpoint '
                         f'(available: {available}).')
    values.update(record.get('parameters', {}))
    roles = {name: ('shared' if name in shared_names else 'local') for name in values}
    # ``model.predict`` needs only the structure and selected replica.  Turn
    # the stored joint placeholder into a normal, self-contained replay config.
    chosen.update({'sample': int(wanted[0]), 'validation_replica': int(wanted[1]),
                   'fit_mode': 'single'})
    chosen.setdefault('parameter_override', {}).update(
        {name: values[name] for name in shared_names if name in values})
    return run_config, values, roles, wanted, record


def selected_parameters(values, roles, patterns, scope):
    """Find perturbable, numeric model parameters in deterministic order."""
    selected = []
    for name, value in values.items():
        if name in REPLICA_PROVIDED or isinstance(value, bool):
            continue
        if not isinstance(value, (int, float)) or not math.isfinite(value):
            continue
        if scope != 'all' and roles.get(name) != scope:
            continue
        if patterns and not any(fnmatch.fnmatchcase(name, pattern) for pattern in patterns):
            continue
        selected.append(name)
    if not selected:
        raise ValueError('No numeric parameters matched the requested scope/patterns.')
    return sorted(selected)


def _combine_series(left, right, operator):
    """Evaluate two logged quantities on their shared time interval."""
    left_t, left_y = left
    right_t, right_y = right
    start, end = max(left_t[0], right_t[0]), min(left_t[-1], right_t[-1])
    if start > end:
        raise ValueError('Target quantities have no overlapping time interval.')
    times = np.union1d(left_t[(left_t >= start) & (left_t <= end)],
                       right_t[(right_t >= start) & (right_t <= end)])
    return times, operator(np.interp(times, left_t, left_y),
                           np.interp(times, right_t, right_y))


def target_series(run_log, expression):
    """Return a logged quantity or a safe arithmetic expression of quantities."""
    if expression in run_log.keys():
        times, values = run_log[expression]
        return np.asarray(times, dtype=float), np.asarray(values, dtype=float)
    try:
        tree = ast.parse(expression, mode='eval')
    except SyntaxError as error:
        raise ValueError(f'Invalid target expression {expression!r}: {error.msg}') from error

    def evaluate(node):
        if isinstance(node, ast.Name):
            if node.id not in run_log.keys():
                available = ', '.join(sorted(run_log.keys()))
                raise ValueError(f'Quantity {node.id!r} was not logged. Available quantities: {available}')
            times, values = run_log[node.id]
            return np.asarray(times, dtype=float), np.asarray(values, dtype=float)
        if isinstance(node, ast.BinOp) and isinstance(node.op, (ast.Add, ast.Sub, ast.Div)):
            operator = {
                ast.Add: np.add,
                ast.Sub: np.subtract,
                ast.Div: np.divide,
            }[type(node.op)]
            return _combine_series(evaluate(node.left), evaluate(node.right), operator)
        raise ValueError('Target expressions may contain logged quantity names joined with +, -, or /.')

    return evaluate(tree.body)


def replay(run_config, values, replica, target, stat):
    """Run one trajectory and reduce a target expression to one scalar."""
    pathway_model = model.build_model(run_config['chosen']['pathways'],
                                      run_config['chosen']['normalized_parameters'])
    pathway_model.parameters().set(parameters.ModelParameters(values))
    run_log = pathway_model.predict(replica)
    times, series = target_series(run_log, target)
    if stat == 'final':
        return float(series[-1])
    if stat == 'maximum':
        return float(np.max(series))
    if stat == 'integral':
        # ``trapz`` is available in the NumPy versions supported by this
        # project; ``trapezoid`` was added only in newer releases.
        return float(np.trapz(series, times))
    raise ValueError(f'Unknown statistic {stat!r}')


def perturbed(value, fraction, direction):
    """Relative OAT perturbation; zero-valued parameters cannot be varied this way."""
    if value == 0:
        return None
    return value * (1.0 + direction * fraction)


def write_plot(rows, target, baseline, output, metric):
    successful = [row for row in rows if row['minus_target'] is not None
                  and row['plus_target'] is not None]
    if not successful:
        figure, axis = plt.subplots(figsize=(8, 3))
        axis.text(.5, .5, 'No perturbation completed successfully.',
                  ha='center', va='center', transform=axis.transAxes)
        axis.set_axis_off()
        figure.tight_layout()
        figure.savefig(output, dpi=180)
        plt.close(figure)
        return
    relative = baseline != 0
    ranking = ({'derivative': 'absolute_derivative',
                'elasticity': 'absolute_elasticity'}[metric]
               if (metric == 'derivative' or relative) else 'absolute_target_effect')
    successful.sort(key=lambda row: row[ranking], reverse=True)
    labels = [f"{row['parameter']} ({row['role']})" for row in successful]
    if metric == 'derivative':
        plus = [row['derivative'] for row in successful]
        minus = plus
        x_label = f'Centred derivative d({target}) / d(parameter)'
    elif relative:
        plus = [row['plus_percent_change'] for row in successful]
        minus = [row['minus_percent_change'] for row in successful]
        x_label = f'Finite-difference elasticity d ln({target}) / d ln(parameter)'
    else:
        plus = [row['plus_target'] - baseline for row in successful]
        minus = [row['minus_target'] - baseline for row in successful]
        x_label = f'Absolute change in {target} from baseline'
    height = max(4, 0.32 * len(successful) + 1.8)
    figure, axis = plt.subplots(figsize=(10, height))
    y = np.arange(len(successful))
    if metric == 'derivative':
        axis.barh(y, plus, height=.60, label='centred difference')
    else:
        axis.barh(y - .19, minus, height=.36, label='− perturbation')
        axis.barh(y + .19, plus, height=.36, label='+ perturbation')
    axis.axvline(0, color='black', linewidth=.8)
    # A derivative can be positive or negative.  ``symlog`` is logarithmic
    # away from zero while retaining the sign, unlike an ordinary log axis.
    if metric == 'derivative':
        maximum = max(abs(value) for value in plus)
        if maximum > 0:
            axis.set_xscale('symlog', linthresh=maximum * 1e-3)
    axis.set_yticks(y, labels)
    axis.invert_yaxis()
    axis.set_xlabel(x_label)
    axis.set_title(f'One-at-a-time sensitivity: {target} (baseline {baseline:.6g})')
    axis.legend()
    figure.tight_layout()
    figure.savefig(output, dpi=180)
    plt.close(figure)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('checkpoint', help='Saved checkpoint path or unique checkpoint ID.')
    parser.add_argument('--replica', nargs=2, metavar=('SAMPLE', 'REPLICA'),
                        help='Replica to replay; required for a joint multi-replica checkpoint.')
    parser.add_argument('--target', '--quantity', dest='target', default='CO2',
                        help=('A logged quantity or + / - / expression of logged quantities, '
                              "e.g. 'CO2 + CH4' or 'CO2 / CH4' (default: CO2)."))
    parser.add_argument('--stat', choices=('final', 'maximum', 'integral'), default='final',
                        help='How to reduce the logged time series (default: final).')
    parser.add_argument('--fraction', type=float, default=1e-4,
                        help='Relative finite-difference step h/|p| (default: 1e-4).')
    parser.add_argument('--scope', choices=('all', 'shared', 'local'), default='all',
                        help='Parameter role to test; meaningful for joint fits (default: all).')
    parser.add_argument('--parameters', nargs='*', metavar='GLOB',
                        help='Optional parameter names or shell-style patterns, e.g. M_* Hydro_*.')
    parser.add_argument('--plot', choices=('derivative', 'elasticity'), default='derivative',
                        help=('Metric in parameter_influence.png (default: derivative). '
                              'Elasticity is dimensionless and useful for comparing unlike units.'))
    parser.add_argument('--output-dir', type=Path,
                        help='Override the default results/sensitivity/SAMPLE-REPLICA directory.')
    args = parser.parse_args(argv)
    if not 0 < args.fraction < 1:
        parser.error('--fraction must be greater than zero and less than one.')

    source, _, _, contents = replay_inputs(args.checkpoint)
    run_config, values, roles, replica_id, local_record = checkpoint_solution(
        source, contents, args.replica)
    dataset = data.get_data_before_day()
    replica = dataset[replica_id[0]][replica_id[1]]
    names = selected_parameters(values, roles, args.parameters, args.scope)
    baseline = replay(run_config, values, replica, args.target, args.stat)

    short_id = checkpoint_id(source) or source.name
    # This is deliberately derived from the project root, rather than the
    # process working directory: analyses always belong in
    # ``simple-model/results/``, never ``simple-model/faster/results/``.
    default_root = (Path(PROJECT_DIRECTORY) / 'results' / 'sensitivity' /
                    f'{replica_id[0]}-{replica_id[1]}')
    target_slug = re.sub(r'[^A-Za-z0-9_.-]+', '_', args.target).strip('_')
    output_dir = (args.output_dir or default_root /
                  f'cp-{short_id}_{target_slug}_{args.stat}').resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for name in names:
        print('perturbing', f'{name[:20]:20s}', end = '')
        value = float(values[name])
        low, high = perturbed(value, args.fraction, -1), perturbed(value, args.fraction, 1)
        row = {'parameter': name, 'role': roles.get(name, 'single'), 'baseline_parameter': value,
               'minus_parameter': low, 'plus_parameter': high, 'baseline_target': baseline,
               'minus_target': None, 'plus_target': None, 'minus_percent_change': None,
               'plus_percent_change': None, 'elasticity': None, 'absolute_percent_effect': None,
               'derivative': None, 'absolute_derivative': None,
               'absolute_elasticity': None, 'absolute_target_effect': None,
               'status': 'ok'}
        if low is None:
            row.update(status='skipped: zero-valued parameter')
            rows.append(row)
            continue
        try:
            lower_values, upper_values = dict(values), dict(values)
            lower_values[name], upper_values[name] = low, high
            row['minus_target'] = replay(run_config, lower_values, replica, args.target, args.stat)
            row['plus_target'] = replay(run_config, upper_values, replica, args.target, args.stat)
            row['absolute_target_effect'] = max(abs(row['minus_target'] - baseline),
                                                abs(row['plus_target'] - baseline))
            row['derivative'] = ((row['plus_target'] - row['minus_target'])
                                 / (high - low))
            print(f'{row["derivative"]:-.2e}')
            row['absolute_derivative'] = abs(row['derivative'])
            if baseline != 0:
                row['minus_percent_change'] = 100 * (row['minus_target'] - baseline) / baseline
                row['plus_percent_change'] = 100 * (row['plus_target'] - baseline) / baseline
                row['absolute_percent_effect'] = max(abs(row['minus_percent_change']),
                                                     abs(row['plus_percent_change']))
                row['elasticity'] = ((row['plus_target'] - row['minus_target']) / baseline
                                     / (2 * args.fraction))
                row['absolute_elasticity'] = abs(row['elasticity'])
            else:
                row.update(status='baseline target is zero; percent effects unavailable')
        except Exception as error:  # preserve results from the other parameters
            row.update(status=f'failed: {type(error).__name__}: {error}')
        rows.append(row)

    fields = list(rows[0])
    with (output_dir / 'sensitivity.csv').open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    write_plot(rows, f'{args.stat}({args.target})', baseline,
               output_dir / 'parameter_influence.png', args.plot)
    metadata = {
        'created_at_utc': datetime.now(timezone.utc).isoformat(),
        'source_checkpoint': str(source.resolve()), 'checkpoint_id': short_id,
        'checkpoint_total_loss': contents.get('total_loss'), 'replica': {
            'sample': replica_id[0], 'replica': replica_id[1]},
        'joint_checkpoint': run_config['chosen'].get('fit_mode') == 'multi' or local_record is not None,
        'target_expression': args.target, 'quantity': args.target,
        'statistic': args.stat, 'baseline_target': baseline,
        'method': 'centred finite difference',
        'derivative_definition': 'dQ/dp ~= (Q(p+h) - Q(p-h)) / (2h)',
        'relative_step_h_over_abs_p': args.fraction, 'plot_metric': args.plot,
        'derivative_plot_scale': 'symlog (signed logarithmic)',
        'scope': args.scope,
        'parameter_patterns': args.parameters, 'tested_parameters': names,
        'parameter_roles': {name: roles.get(name) for name in names},
        'selected_local_solution': local_record,
        'files': {'table': 'sensitivity.csv', 'plot': 'parameter_influence.png'},
    }
    with (output_dir / 'metadata.json').open('w') as handle:
        json.dump(metadata, handle, indent=2, default=str)
    print(f'Wrote sensitivity analysis for {len(names)} parameters to {output_dir}')


if __name__ == '__main__':
    main()
