"""Functional PCA of biomass-scaled pathway-activity signatures.

Each selected checkpoint is evaluated once on a common time grid.  Its
signature is the vector-valued function consisting of the six accepted-state
``ln(max(M * v, floor))`` reaction-rate curves in :data:`PATHWAYS` order.
The FPCA follows the documented scikit-fda discretized workflow: an
``FDataGrid`` is fitted with :class:`skfda.preprocessing.dim_reduction.FPCA`.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import checkpoint_report
import data
import parameter_envelope
import parameters
from model_factory import build_model
from plot import pathway_color
from USER_VARIABLES import PLOT_DIRECTORY, RESULTS_DIRECTORY


# This is deliberately independent of model construction order.  It defines
# the coordinates of every signature, output figure, and saved coefficient.
PATHWAYS = checkpoint_report.PATHWAY_PLOT_ORDER
DEFAULT_ACTIVITY_FLOOR = 1e-20


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results-directory', default=RESULTS_DIRECTORY)
    parser.add_argument('--model-id', dest='model_ids', nargs='+', metavar='MODEL_ID')
    parser.add_argument('--fit-mode', dest='fit_modes', nargs='+',
                        choices=('single', 'split', 'multi'), metavar='FIT_MODE')
    parser.add_argument('--sample', dest='samples', nargs='+', metavar='SAMPLE')
    parser.add_argument('--replica', dest='replicas', nargs='+', metavar='REPLICA')
    parser.add_argument('--site', dest='sites', nargs='+', metavar='SITE')
    parser.add_argument('--origin', dest='origins', nargs='+', metavar='ORIGIN')
    parser.add_argument('--exclude', dest='excluded_ids', nargs='+', metavar='CP_ID')
    parser.add_argument('--best-restart', action='store_true',
                        help=('Keep the best checkpoint per run and the best checkpoint '
                              'from each CMA-ES restart, as in parameter_envelope.py.'))
    parser.add_argument('--loss-upper', type=float, default=1e99, metavar='LOSS')
    parser.add_argument('--components', type=int, default=3,
                        help='Number of functional principal components (default: 3).')
    parser.add_argument('--component-factor', type=float, default=10.,
                        help=('Multiplier used for mean ± component mode plots '
                              '(default: 10).'))
    parser.add_argument('--time-end', type=float,
                        help=('End of the common evaluation grid in days. Defaults to the '
                              'shortest selected experimental time range.'))
    parser.add_argument('--time-points', type=int, default=256,
                        help='Number of equally spaced time points (default: 256).')
    parser.add_argument('--activity-floor', type=float, default=DEFAULT_ACTIVITY_FLOOR,
                        help=('Positive floor applied to M*v before its natural logarithm '
                              f'(default: {DEFAULT_ACTIVITY_FLOOR:g}).'))
    parser.add_argument('--output-directory', type=Path,
                        default=Path(PLOT_DIRECTORY) / 'fpca')
    parser.add_argument('--yes', action='store_true',
                        help='Skip the interactive checkpoint confirmation.')
    return parser.parse_args(argv)


def select_checkpoints(args):
    """Select checkpoints using the same criteria as parameter envelopes."""
    if not math.isfinite(args.loss_upper):
        raise ValueError('--loss-upper must be finite.')
    checkpoints = checkpoint_report.collect_all(args.results_directory, verbose=False)
    metadata = (parameter_envelope.sample_metadata(data.get_data_before_day())
                if args.sites or args.origins else {})
    selected = parameter_envelope.select_checkpoints(
        checkpoints, model_ids=args.model_ids, fit_modes=args.fit_modes,
        samples=args.samples, replicas=args.replicas, sites=args.sites,
        origins=args.origins, excluded_ids=args.excluded_ids,
        loss_upper=args.loss_upper, metadata=metadata)
    if args.best_restart:
        selected = parameter_envelope.best_with_restarts(selected)
    return selected


def replica_for_checkpoint(checkpoint, dataset):
    """Return the single experimental replica represented by a checkpoint.

    A joint checkpoint has multiple experimental initial states, hence it
    cannot have one unambiguous signature.  The checkpoint selectors can
    still include joint fits, but the caller must narrow them to a unique
    replica before this analysis is meaningful.
    """
    records = checkpoint.multi_replica_results
    if records:
        raise ValueError(
            f'{checkpoint.id} is a joint multi-replica checkpoint. FPCA needs one '
            'trajectory per checkpoint; select the corresponding saved inner '
            'checkpoints instead.')
    sample_id, replica_id = map(str, checkpoint.replica)
    try:
        sample = dataset[sample_id]
    except KeyError as error:
        raise ValueError(f'{checkpoint.id}: sample {sample_id!r} is unavailable.') from error
    for replica in sample.replicas:
        if str(replica.replica_number) == replica_id:
            return replica
    raise ValueError(f'{checkpoint.id}: replica {sample_id}/{replica_id} is unavailable.')


def common_time_grid(checkpoints, dataset, time_end, time_points):
    """Return a common grid without extrapolating by default."""
    if time_points < 2:
        raise ValueError('--time-points must be at least 2.')
    if time_end is None:
        ends = [float(np.max(replica_for_checkpoint(checkpoint, dataset)['days']))
                for checkpoint in checkpoints]
        time_end = min(ends)
    if not math.isfinite(time_end) or time_end <= 0:
        raise ValueError('--time-end must be finite and positive.')
    return np.linspace(0., time_end, time_points)


def pathway_model_for(checkpoint):
    """Build a checkpoint model without constructing an optimisation objective."""
    config = checkpoint.run_config
    chosen = config['chosen']
    model = build_model(str(config['model']['version']), chosen['pathways'],
                        normalized=chosen['normalized_parameters'])
    model.use_gas_limits = bool(chosen.get('gas_limits', False))
    model.parameters().set(parameters.ModelParameters(checkpoint.parameters))
    return model


def log_activity(signature, floor=DEFAULT_ACTIVITY_FLOOR):
    """Return ``ln(max(M*v, floor))`` with a physical zero-rate safeguard."""
    if not math.isfinite(floor) or floor <= 0:
        raise ValueError('--activity-floor must be finite and positive.')
    signature = np.asarray(signature, dtype=float)
    if not np.all(np.isfinite(signature)):
        raise ValueError('M*v signature contains non-finite values.')
    if np.any(signature < 0):
        raise ValueError('M*v signature contains negative reaction rates.')
    return np.log(np.maximum(signature, floor))


def signature_for_checkpoint(checkpoint, replica, time, *, activity_floor=DEFAULT_ACTIVITY_FLOOR):
    """Evaluate once and return its six-pathway log-activity signature."""
    model = pathway_model_for(checkpoint)
    # This single prediction also derives all accepted-state reaction-rate
    # curves; it never runs the saved optimiser or changes checkpoint files.
    run_log = model.predict_with_logging(
        replica, t_eval=time, quiet=True, reference='L', microbe_unit='micromol')
    curves = []
    for pathway in PATHWAYS:
        name = f'{pathway}_reaction_rate'
        if name not in run_log.keys():
            raise ValueError(f'{checkpoint.id}: required pathway {pathway!r} is absent.')
        logged_time, values = run_log[name]
        if not np.array_equal(logged_time, time):
            values = np.interp(time, logged_time, values)
        curves.append(np.asarray(values, dtype=float))
    signature = np.asarray(curves)
    try:
        return log_activity(signature, activity_floor)
    except ValueError as error:
        raise ValueError(f'{checkpoint.id}: {error}') from error


def trapezoid_weights(time):
    """Quadrature weights defining the L2 time inner product."""
    delta = np.diff(np.asarray(time, dtype=float))
    if np.any(delta <= 0):
        raise ValueError('FPCA time grid must be strictly increasing.')
    return np.r_[delta[0] / 2, (delta[:-1] + delta[1:]) / 2, delta[-1] / 2]


def functional_pca(signatures, time, components):
    """Fit scikit-fda discretized FPCA and return mean, components, and scores."""
    signatures = np.asarray(signatures, dtype=float)
    if signatures.ndim != 3 or signatures.shape[1] != len(PATHWAYS):
        raise ValueError('Signatures must have shape (checkpoint, pathway, time).')
    count = signatures.shape[0]
    if not 1 <= components <= min(count, signatures.shape[1] * signatures.shape[2]):
        raise ValueError('--components must be between 1 and the number of checkpoints.')
    try:
        from skfda import FDataGrid
        from skfda.preprocessing.dim_reduction import FPCA
    except ImportError as error:
        raise ImportError(
            'FPCA requires scikit-fda in the interpreter running this script '
            f'({sys.executable}). Install it with:\n'
            f'  {sys.executable} -m pip install -r requirements.txt\n'
            f'Underlying import error: {error}') from error

    # FPCA's grid implementation accepts scalar-valued FDataGrid objects.
    # Represent the vector-valued function on six disjoint, adjacent pathway
    # blocks and supply explicit trapezoidal weights for each block. Thus its
    # scalar inner product is exactly the sum of the six pathway L2 inner
    # products, while every score/component still belongs to one joint
    # six-dimensional signature (not six separate FPCAs).
    weights = np.tile(trapezoid_weights(time), len(PATHWAYS))
    fd = FDataGrid(
        data_matrix=signatures.reshape(count, -1, 1),
        grid_points=np.arange(signatures.shape[1] * signatures.shape[2], dtype=float),
    )
    fpca = FPCA(n_components=components, _weights=weights)
    scores = fpca.fit_transform(fd)
    mean = fpca.mean_.data_matrix[0, :, 0].reshape(len(PATHWAYS), -1)
    basis = fpca.components_.data_matrix[:, :, 0].reshape(
        components, len(PATHWAYS), -1)
    explained = np.asarray(fpca.explained_variance_ratio_, dtype=float)
    return mean, basis, scores, explained


def _signature_axis(axis, time, values, title, ylabel):
    for pathway, curve in zip(PATHWAYS, values):
        axis.plot(time, curve, color=pathway_color(pathway), label=pathway)
    axis.set_title(title)
    axis.set_xlabel('Time [d]')
    axis.set_ylabel(ylabel)
    # The values are already log-transformed; a log axis would apply the
    # transformation twice and cannot display legitimate negative log values.
    axis.set_yscale('linear')
    axis.spines[['top', 'right']].set_visible(False)
    axis.legend(loc='best', fontsize='small')


def shared_component_ylim(basis, margin=0.05):
    """Return one padded y-range that contains every component loading."""
    values = np.asarray(basis, dtype=float)
    low = float(np.min(values))
    high = float(np.max(values))
    if not math.isfinite(low) or not math.isfinite(high):
        raise ValueError('FPCA component loadings contain non-finite values.')
    span = high - low
    padding = (span * margin) if span else max(abs(low) * margin, 1.)
    return low - padding, high + padding


def _add_checkpoint_labels(axis, scores, checkpoint_ids, dimensions):
    """Label a manageable number of points without obscuring dense plots."""
    if len(checkpoint_ids) > 15:
        return
    for checkpoint_id, score in zip(checkpoint_ids, scores):
        label = f'cp-{checkpoint_id}'
        if dimensions == 2:
            axis.annotate(label, (score[0], score[1]), xytext=(4, 4),
                          textcoords='offset points', fontsize='x-small')
        else:
            axis.text(score[0], score[1], score[2], label, fontsize='x-small')


def write_checkpoint_scores(output_directory, scores, checkpoint_ids):
    """Plot FPCA checkpoint coefficients in the most legible coordinate system."""
    component_count = scores.shape[1]
    if component_count == 2:
        figure, axis = plt.subplots(figsize=(7, 5.5))
        axis.scatter(scores[:, 0], scores[:, 1], alpha=.8)
        axis.axhline(0, color='.7', linewidth=.8)
        axis.axvline(0, color='.7', linewidth=.8)
        axis.set_xlabel('PC 1 score')
        axis.set_ylabel('PC 2 score')
        axis.set_title('Checkpoint FPCA coefficients')
        axis.spines[['top', 'right']].set_visible(False)
        _add_checkpoint_labels(axis, scores, checkpoint_ids, dimensions=2)
        figure.tight_layout()
        figure.savefig(output_directory / 'checkpoint_coefficients_2d.png', dpi=200)
        plt.close(figure)
        return

    if component_count == 3:
        figure = plt.figure(figsize=(7, 6))
        axis = figure.add_subplot(projection='3d')
        axis.scatter(scores[:, 0], scores[:, 1], scores[:, 2], alpha=.8)
        axis.set_xlabel('PC 1 score')
        axis.set_ylabel('PC 2 score')
        axis.set_zlabel('PC 3 score')
        axis.set_title('Checkpoint FPCA coefficients')
        _add_checkpoint_labels(axis, scores, checkpoint_ids, dimensions=3)
        figure.tight_layout()
        figure.savefig(output_directory / 'checkpoint_coefficients_3d.png', dpi=200)
        plt.close(figure)
        return

    angles = np.linspace(0, 2 * np.pi, component_count, endpoint=False)
    angles = np.r_[angles, angles[0]]
    scales = np.std(scores, axis=0, ddof=1) if len(scores) > 1 else np.ones(component_count)
    scales[scales == 0] = 1
    figure, axis = plt.subplots(figsize=(7, 7), subplot_kw={'polar': True})
    for checkpoint_id, score in zip(checkpoint_ids, scores):
        radial = np.r_[score / scales, score[0] / scales[0]]
        axis.plot(angles, radial, linewidth=.9, alpha=.7,
                  label=f'cp-{checkpoint_id}')
    axis.set_xticks(angles[:-1], [f'PC {index}' for index in range(1, component_count + 1)])
    axis.set_title('Checkpoint FPCA coefficients (each PC scaled to 1 SD)')
    if len(checkpoint_ids) <= 15:
        axis.legend(loc='upper left', bbox_to_anchor=(1.05, 1.05), fontsize='x-small')
    figure.tight_layout()
    figure.savefig(output_directory / 'checkpoint_coefficients_radar.png', dpi=200)
    plt.close(figure)


def write_figures(output_directory, time, mean, basis, scores, explained, checkpoint_ids,
                  component_factor):
    output_directory.mkdir(parents=True, exist_ok=True)
    figure, axis = plt.subplots(figsize=(8, 4.5))
    _signature_axis(axis, time, mean, 'Mean log pathway-activity signature',
                    r'$\ln(M \times v)$ [reference: 1 reaction extent / L / d]')
    figure.tight_layout()
    figure.savefig(output_directory / 'mean_signature.png', dpi=200)
    plt.close(figure)

    component_ylim = shared_component_ylim(basis)
    mode_signatures = np.concatenate((
        mean[None, :, :] + component_factor * basis,
        mean[None, :, :] - component_factor * basis,
    ))
    mode_ylim = shared_component_ylim(mode_signatures, margin=0.05)
    for index, component in enumerate(basis, start=1):
        figure, axis = plt.subplots(figsize=(8, 4.5))
        _signature_axis(axis, time, component,
                        f'FPCA component {index} ({explained[index - 1]:.1%} variance)',
                        r'$\ln(M \times v)$ component loading')
        axis.set_ylim(component_ylim)
        figure.tight_layout()
        figure.savefig(output_directory / f'component_{index:02d}.png', dpi=200)
        plt.close(figure)

        upper = mean + component_factor * component
        lower = mean - component_factor * component
        figure, axis = plt.subplots(figsize=(8, 4.5))
        for pathway, mean_curve, plus_curve, minus_curve in zip(
                PATHWAYS, mean, upper, lower):
            color = pathway_color(pathway)
            axis.fill_between(time, np.minimum(plus_curve, minus_curve),
                              np.maximum(plus_curve, minus_curve),
                              color=color, alpha=.18)
            axis.plot(time, mean_curve, color=color, label=pathway)
        axis.set_title(
            f'FPCA mode {index}: mean ± {component_factor:g} × PC {index} '
            f'({explained[index - 1]:.1%} variance)')
        axis.set_xlabel('Time [d]')
        axis.set_ylabel(r'$\ln(M \times v)$')
        axis.set_yscale('linear')
        axis.set_ylim(mode_ylim)
        axis.spines[['top', 'right']].set_visible(False)
        axis.legend(loc='best', fontsize='small')
        figure.tight_layout()
        figure.savefig(output_directory / f'component_{index:02d}_mode.png', dpi=200)
        plt.close(figure)

    write_checkpoint_scores(output_directory, scores, checkpoint_ids)


def write_results(output_directory, args, checkpoints, time, scores, explained):
    with (output_directory / 'coefficients.csv').open('w', newline='') as handle:
        writer = csv.writer(handle)
        writer.writerow(['checkpoint_id', 'source', *[f'PC{index}' for index in range(1, scores.shape[1] + 1)]])
        for checkpoint, score in zip(checkpoints, scores):
            writer.writerow([checkpoint.id, checkpoint.source, *score])
    with (output_directory / 'manifest.json').open('w') as handle:
        json.dump({
            'pathway_order': PATHWAYS,
            'signature': 'ln(max(M * v, activity_floor))',
            'activity_floor': args.activity_floor,
            'component_factor': args.component_factor,
            'activity_reference_unit': '1 reaction extent / L / d',
            'time_grid_days': time.tolist(),
            'checkpoint_ids': [checkpoint.id for checkpoint in checkpoints],
            'explained_variance_ratio': explained.tolist(),
            'selection': {key: getattr(args, key) for key in (
                'model_ids', 'fit_modes', 'samples', 'replicas', 'sites', 'origins',
                'excluded_ids', 'best_restart', 'loss_upper', 'activity_floor',
                'component_factor')},
        }, handle, indent=2)


def confirm(count, *, input_function=input):
    try:
        response = input_function('Evaluate one log(M*v) signature for each included checkpoint? [y/N] ')
    except EOFError as error:
        raise RuntimeError('FPCA evaluation was not confirmed.') from error
    if response.strip().casefold() not in {'y', 'yes'}:
        raise RuntimeError('FPCA evaluation was not confirmed.')


def main(argv=None):
    args = parse_args(argv)
    checkpoints = select_checkpoints(args)
    if len(checkpoints) < 2:
        raise ValueError('FPCA requires at least two selected checkpoints.')
    if args.components < 1 or args.components > len(checkpoints):
        raise ValueError('--components must be between 1 and the selected checkpoint count.')
    if not math.isfinite(args.activity_floor) or args.activity_floor <= 0:
        raise ValueError('--activity-floor must be finite and positive.')
    if not math.isfinite(args.component_factor):
        raise ValueError('--component-factor must be finite.')
    print(f'Included checkpoints: {len(checkpoints)}')
    if not args.yes:
        confirm(len(checkpoints))
    dataset = data.get_data_before_day()
    time = common_time_grid(checkpoints, dataset, args.time_end, args.time_points)
    print(f'Common time grid: 0 to {time[-1]:g} d ({len(time)} points).')
    signatures = np.asarray([
        signature_for_checkpoint(
            checkpoint, replica_for_checkpoint(checkpoint, dataset), time,
            activity_floor=args.activity_floor)
        for checkpoint in checkpoints])
    mean, basis, scores, explained = functional_pca(signatures, time, args.components)
    write_figures(args.output_directory, time, mean, basis, scores, explained,
                  [checkpoint.id for checkpoint in checkpoints], args.component_factor)
    write_results(args.output_directory, args, checkpoints, time, scores, explained)
    print(f'Wrote FPCA figures and coefficients to {args.output_directory}')
    return mean, basis, scores


if __name__ == '__main__':
    main()
