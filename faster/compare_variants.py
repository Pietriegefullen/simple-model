"""Generate like-for-like visual comparisons of fitted model variants.

The command groups best checkpoints by loss definition, fit mode, sample, and
validation replica.  A fit comparison is therefore only made when the models
were calibrated against the same data split and objective.

Examples
--------
Generate every available A/B(/C) fit comparison and the focused A/C Homo
Gibbs-energy comparison::

    python compare_variants.py

Limit output to the one currently complete A/B/C comparison group::

    python compare_variants.py --sample 1351 --fit-mode split --replica 4
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

import collect_best
import data
import optimizer
import plot
from chemistry import GIBBS_MINIMUM
from fit_sample import run
from USER_VARIABLES import PROJECT_DIRECTORY, RESULTS_DIRECTORY


VARIANT_LINESTYLES = {
    'A': '-',       # full model
    'B': '--',      # no thermodynamic factor
    'C': ':',       # no Hydro pathway
}
VARIANT_COLORS = {
    'A': 'tab:blue',
    'B': 'tab:orange',
    'C': 'tab:green',
}
DEFAULT_VARIANTS = ('A', 'B', 'C')


def comparison_groups(checkpoints, variants=DEFAULT_VARIANTS):
    """Return comparable checkpoints, grouped by common data and objective.

    Some historical model IDs map to the same named variant.  In that case the
    lower-loss checkpoint is retained, rather than making duplicate curves for
    one conceptual variant.
    """
    wanted = set(variants)
    groups = defaultdict(dict)
    for checkpoint in checkpoints:
        variant = collect_best.model_variant(checkpoint)
        if variant not in wanted:
            continue
        key = (checkpoint.loss_id, checkpoint.fit_mode, checkpoint.replica)
        previous = groups[key].get(variant)
        if previous is None or checkpoint.loss < previous.loss:
            groups[key][variant] = checkpoint
    return dict(groups)


def comparison_suptitle(checkpoint, variants, suffix=''):
    """Return the fit-plot title with multiple model variants made explicit."""
    parts = [checkpoint.replica[0], checkpoint.fit_mode,
             'models ' + ', '.join(variants)]
    if suffix:
        parts.append(suffix)
    return ' | '.join(parts)


def _validate_comparison(checkpoints):
    """Reject model overlays that do not use the same fitting context."""
    reference = checkpoints[0]
    reference_key = (reference.loss_id, reference.fit_mode, reference.replica)
    for checkpoint in checkpoints[1:]:
        key = (checkpoint.loss_id, checkpoint.fit_mode, checkpoint.replica)
        if key != reference_key:
            raise ValueError('Variant comparisons require the same loss, fit mode, '
                             'sample, and validation replica.')


def _validation_objective(checkpoint, pathway_model, validation_replica):
    return optimizer.build_objective_function(
        pathway_model,
        validation_replica,
        checkpoint.run_config['objective'],
        t_start=0,
        t_end=None,
    )


def run_variants(checkpoints):
    """Replay each checkpoint once and index the resulting runs by variant."""
    results = {}
    for checkpoint in checkpoints:
        results[collect_best.model_variant(checkpoint)] = run(
            checkpoint.run_config, checkpoint.parameters)
    return results


def plot_fit_comparison(checkpoints, dataset, model_runs=None):
    """Overlay two or more fitted variants using the standard fit-plot style.

    Gas identity remains encoded by the CO2/CH4 panel colours; model identity
    is encoded by line style.  Measurements are drawn once, so identical data
    are not visually amplified for every compared model.
    """
    checkpoints = sorted(checkpoints, key=collect_best.model_variant)
    if len(checkpoints) < 2:
        raise ValueError('At least two model variants are required for a comparison.')
    _validate_comparison(checkpoints)
    if model_runs is None:
        model_runs = run_variants(checkpoints)

    reference = checkpoints[0]
    sample, validation_replica, fitted_replicas = collect_best.fit_replicas(
        reference, dataset)
    figure, axes = collect_best.create_fit_figure()
    handles = collect_best.plot_fit_measurements(reference, dataset, axes)

    variants = []
    for checkpoint in checkpoints:
        variant = collect_best.model_variant(checkpoint)
        pathway_model, objective, run_log = model_runs[variant]
        validation_objective = _validation_objective(
            checkpoint, pathway_model, validation_replica)
        linestyle = VARIANT_LINESTYLES.get(variant, '-')
        for pool, axis in axes.items():
            time, values = run_log[pool]
            axis.plot(time, values, color=plot.pool_color(pool),
                      linestyle=linestyle, linewidth=1.5)
        handles.append(Line2D(
            [], [], color='k', linestyle=linestyle, linewidth=1.5,
            label=collect_best.fit_score_label(
                checkpoint, objective, validation_objective, run_log),
        ))
        variants.append(variant)

    axes['CH4'].legend(
        handles=handles, fancybox=False, edgecolor='k', loc='lower right')
    collect_best.format_fit_axes(reference, sample, axes)
    collect_best.finish_fit_figure(
        figure, comparison_suptitle(reference, variants))
    return figure, axes


def plot_homo_delta_g_comparison(checkpoints, model_runs=None):
    """Compare Homo-pathway ΔG trajectories for variants A and C.

    This answers the mechanistic question directly: removing the Hydro pathway
    (model C) changes H2/CO2 availability, which shifts Homo's Gibbs energy.
    The horizontal reference is the model's Gibbs minimum.
    """
    checkpoints = sorted(checkpoints, key=collect_best.model_variant)
    _validate_comparison(checkpoints)
    variants = [collect_best.model_variant(checkpoint) for checkpoint in checkpoints]
    if set(variants) != {'A', 'C'}:
        raise ValueError('The Homo ΔG comparison requires exactly variants A and C.')
    if model_runs is None:
        model_runs = run_variants(checkpoints)

    reference = checkpoints[0]
    figure, axes = plot.create_figure(plot.FigureSpec(
        axis_names=('homo_delta_g',),
        figsize=(4, 4),
    ))
    axis = axes['homo_delta_g']
    for checkpoint in checkpoints:
        variant = collect_best.model_variant(checkpoint)
        _, _, run_log = model_runs[variant]
        try:
            time, values = run_log['Homo_deltaG_r']
        except KeyError as error:
            raise ValueError(
                f'Model {variant} did not log Homo_deltaG_r.') from error
        axis.plot(
            time, values,
            color=VARIANT_COLORS[variant],
            linestyle=VARIANT_LINESTYLES[variant],
            linewidth=1.5,
            label=rf'$\mathrm{{model\ {variant}}}$',
        )
    axis.axhline(
        GIBBS_MINIMUM, color='k', linestyle='--', linewidth=1,
        label='Gibbs minimum',
    )
    axis.set_title('Homo pathway')
    axis.set_xlabel(r'$\mathrm{time\ [d]}$')
    axis.set_ylabel(collect_best.quantity_ylabel('Homo_deltaG_r'))
    plot.format_ax(axis)
    axis.legend(fancybox=False, edgecolor='k', loc='best')
    figure.suptitle(
        comparison_suptitle(reference, variants, 'Homo $\\Delta G_\\mathrm{r}$'))
    # The single-panel Gibbs plot is deliberately square, so it needs more
    # room on the left for its vertical energy label than the two-panel fit.
    figure.subplots_adjust(left=0.25, right=0.95, bottom=0.16, top=0.78)
    return figure, axis


def comparison_target(output_directory, group):
    """Return a separate, deterministic target directory for one comparison."""
    loss_id, fit_mode, replica = group
    return Path(output_directory) / f'{loss_id}_{fit_mode}' / '-'.join(replica)


def save_comparisons(results_directory=RESULTS_DIRECTORY, output_directory=None,
                     variants=DEFAULT_VARIANTS, **filters):
    """Write variant fit and A/C Homo ΔG figures; return their paths."""
    output_directory = (Path(PROJECT_DIRECTORY) / 'variant_comparisons'
                        if output_directory is None else Path(output_directory))
    checkpoints = collect_best.filter_checkpoints(
        collect_best.collect_best(results_directory), **filters)
    groups = comparison_groups(checkpoints, variants)
    dataset = data.get_data_before_carex()
    paths = []
    for group, by_variant in sorted(groups.items()):
        selected = [by_variant[variant] for variant in variants
                    if variant in by_variant]
        model_runs = run_variants(selected) if len(selected) >= 2 else None
        if len(selected) >= 2:
            figure, _ = plot_fit_comparison(selected, dataset, model_runs=model_runs)
            selected_variants = [collect_best.model_variant(cp) for cp in selected]
            directory = comparison_target(output_directory, group)
            directory.mkdir(parents=True, exist_ok=True)
            path = directory / ('01_fit_' + '_vs_'.join(selected_variants) + '.png')
            figure.savefig(path, dpi=300)
            plt.close(figure)
            paths.append(path)

        if 'A' in by_variant and 'C' in by_variant:
            selected = [by_variant['A'], by_variant['C']]
            figure, _ = plot_homo_delta_g_comparison(selected, model_runs=model_runs)
            directory = comparison_target(output_directory, group)
            directory.mkdir(parents=True, exist_ok=True)
            path = directory / '02_homo_delta_g_A_vs_C.png'
            figure.savefig(path, dpi=300)
            plt.close(figure)
            paths.append(path)
    return paths


def parse_args():
    parser = argparse.ArgumentParser(
        description='Compare best fitted model variants on shared data splits.')
    parser.add_argument('--results-directory', default=RESULTS_DIRECTORY)
    parser.add_argument('--output-directory',
                        default=str(Path(PROJECT_DIRECTORY) / 'variant_comparisons'))
    parser.add_argument('--sample', dest='samples', nargs='+', metavar='SAMPLE')
    parser.add_argument('--fit-mode', dest='fit_modes', choices=('single', 'split'),
                        nargs='+', metavar='FIT_MODE')
    parser.add_argument('--replica', dest='replicas', nargs='+', metavar='REPLICA')
    parser.add_argument('--loss-id', dest='loss_ids', nargs='+', metavar='LOSS_ID')
    parser.add_argument('--variants', nargs='+', choices=DEFAULT_VARIANTS,
                        default=DEFAULT_VARIANTS)
    parser.add_argument('--dry-run', action='store_true')
    return parser.parse_args()


def main(args=None):
    args = parse_args() if args is None else args
    filters = {
        'samples': args.samples,
        'fit_modes': args.fit_modes,
        'replicas': args.replicas,
        'loss_ids': args.loss_ids,
        'model_ids': None,
    }
    if args.dry_run:
        checkpoints = collect_best.filter_checkpoints(
            collect_best.collect_best(args.results_directory), **filters)
        for group, by_variant in sorted(comparison_groups(
                checkpoints, args.variants).items()):
            available = [variant for variant in args.variants if variant in by_variant]
            if len(available) >= 2:
                print(group, 'models ' + ', '.join(available))
        return
    paths = save_comparisons(
        args.results_directory, args.output_directory, args.variants, **filters)
    for path in paths:
        print(path)
    return paths


if __name__ == '__main__':
    main()
