"""Generate like-for-like visual comparisons of fitted model variants.

The command groups best checkpoints by loss definition, fit mode, sample, and
validation replica.  A fit comparison is therefore only made when the models
were calibrated against the same data split and objective.

Examples
--------
Generate every available A/B plus A-derived C/D fit comparison and the focused A/C hydrogen
consumption comparison::

    python compare_variants.py

Limit output to the one currently complete A/B/C/D comparison group::

    python compare_variants.py --sample 1351 --fit-mode split --replica 4
"""

from __future__ import annotations

import argparse
import copy
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

import plot_best
import compare_checkpoints
import data
import optimizer
import parameters
import plot
import system
import variant_diagnostics
from fit_sample import run
from model_variants import model_variant_from_id
from USER_VARIABLES import PROJECT_DIRECTORY, RESULTS_DIRECTORY


DEFAULT_VARIANTS = ('A', 'B', 'C', 'D')


@dataclass(frozen=True)
class VariantRunSpec:
    """A model variant, its fitting metadata, and values used for its replay."""

    variant: str
    checkpoint: object
    run_config: dict
    parameter_values: parameters.ModelParameters


def variant_legend_label(variant):
    """Label A-derived counterfactuals as replays, rather than refits."""
    if variant == 'C':
        return r'$\mathrm{model\ {C}\ (A\ parameters)}$'
    if variant == 'D':
        return r'$\mathrm{model\ {D}\ (A\ parameters,\ no\ thermodynamics)}$'
    return rf'$\mathrm{{model\ {{{variant}}}}}$'


def comparison_groups(checkpoints, variants=DEFAULT_VARIANTS):
    """Return comparable checkpoints, grouped by common data and objective.

    Some historical model IDs map to the same named variant.  In that case the
    lower-loss checkpoint is retained, rather than making duplicate curves for
    one conceptual variant.
    """
    wanted = set(variants)
    groups = defaultdict(dict)
    for checkpoint in checkpoints:
        # Results may also contain exploratory models that do not have an
        # A/B/C label.  They are irrelevant to this report, not an error.
        variant = model_variant_from_id(checkpoint.model_id)
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


def counterfactual_c_spec(a_checkpoint):
    """Build no-Hydro model C from A's checkpoint and fitted values."""
    run_config = copy.deepcopy(a_checkpoint.run_config)
    chosen = run_config['chosen']
    chosen['pathways'] = [pathway for pathway in chosen['pathways']
                          if pathway != 'Hydro']
    # ``fit_sample.run`` rebuilds from ``chosen['pathways']``.  Keep the
    # stored structural metadata accurate too, for diagnostics and reuse.
    run_config['model'].pop('Hydro', None)
    run_config['range'] = {
        name: config for name, config in run_config['range'].items()
        if not name.startswith('Hydro_')
    }
    parameter_values = parameters.ModelParameters({
        name: config for name, config in a_checkpoint.parameters.get_config().items()
        if not name.startswith('Hydro_')
    })
    return VariantRunSpec('C', a_checkpoint, run_config, parameter_values)


def counterfactual_d_spec(a_checkpoint):
    """Build no-thermodynamics model D from A's fitted parameters.

    D retains every pathway and fitted biological parameter from model A. Its
    thermodynamics switches alone are fixed to ``False``, causing pathway
    rates to be evaluated without their thermodynamic multiplier.
    """
    run_config = copy.deepcopy(a_checkpoint.run_config)
    parameter_values = parameters.ModelParameters(
        a_checkpoint.parameters.get_config())
    thermodynamics = [
        name for name in parameter_values.get_values()
        if name.endswith('_thermodynamics')
    ]
    if not thermodynamics:
        raise ValueError('Model A has no thermodynamics switches to disable.')
    for name in thermodynamics:
        parameter_values[name].constant(False)
        pathway = name.removesuffix('_thermodynamics')
        if pathway in run_config['model']:
            run_config['model'][pathway]['use_thermodynamics'] = False
    run_config['chosen'].setdefault('parameter_override', {}).update({
        name: False for name in thermodynamics
    })
    return VariantRunSpec('D', a_checkpoint, run_config, parameter_values)


def variant_run_specs(by_variant, variants=DEFAULT_VARIANTS):
    """Return exactly the requested replay specifications.

    C and D are counterfactual replays rather than fitted checkpoint variants,
    so they are constructed from A only when their names were explicitly
    requested.  This makes ``--variants`` an exact selection rather than a
    suggestion.
    """
    specs = []
    for variant in variants:
        checkpoint = by_variant.get(variant)
        if variant not in ('C', 'D') and checkpoint is not None:
            specs.append(VariantRunSpec(
                variant, checkpoint, checkpoint.run_config, checkpoint.parameters))
    if 'A' in by_variant:
        if 'C' in variants:
            specs.append(counterfactual_c_spec(by_variant['A']))
        if ('D' in variants
                and any(name.endswith('_thermodynamics')
                        for name in by_variant['A'].parameters.get_values())):
            specs.append(counterfactual_d_spec(by_variant['A']))
    return specs


def _validate_comparison(specifications):
    """Reject model overlays that do not use the same fitting context."""
    reference = specifications[0].checkpoint
    reference_key = (reference.loss_id, reference.fit_mode, reference.replica)
    for specification in specifications[1:]:
        checkpoint = specification.checkpoint
        key = (checkpoint.loss_id, checkpoint.fit_mode, checkpoint.replica)
        if key != reference_key:
            raise ValueError('Variant comparisons require the same loss, fit mode, '
                             'sample, and validation replica.')


def _validation_objective(run_config, pathway_model, validation_replica):
    return optimizer.build_objective_function(
        pathway_model,
        validation_replica,
        run_config['objective'],
        t_start=0,
        t_end=None,
    )


def run_variants(specifications):
    """Replay variant specifications once each."""
    results = {}
    for specification in specifications:
        results[specification.variant] = run(
            specification.run_config, specification.parameter_values)
    return results


def plot_fit_comparison(specifications, dataset, model_runs=None):
    """Overlay two or more fitted variants using the standard fit-plot style.

    Gas identity remains encoded by the CO2/CH4 panel colours; model identity
    is encoded by line style.  Measurements are drawn once, so identical data
    are not visually amplified for every compared model.
    """
    specifications = sorted(specifications, key=lambda spec: spec.variant)
    if len(specifications) < 2:
        raise ValueError('At least two model variants are required for a comparison.')
    _validate_comparison(specifications)
    if model_runs is None:
        model_runs = run_variants(specifications)

    reference = specifications[0].checkpoint
    sample, validation_replica, fitted_replicas = plot_best.fit_replicas(
        reference, dataset)
    figure, axes = plot_best.create_fit_figure()
    handles = plot_best.plot_fit_measurements(reference, dataset, axes)

    variants = []
    for specification in specifications:
        checkpoint = specification.checkpoint
        variant = specification.variant
        pathway_model, objective, run_log = model_runs[variant]
        validation_objective = _validation_objective(
            specification.run_config, pathway_model, validation_replica)
        linestyle = plot.MODEL_VARIANT_LINESTYLES.get(variant, '-')
        for pool, axis in axes.items():
            time, values = run_log[pool]
            axis.plot(time, values, color=plot.pool_color(pool),
                      linestyle=linestyle, linewidth=1.5)
        handles.append(Line2D(
            [], [], color='k', linestyle=linestyle, linewidth=1.5,
            label=plot_best.fit_score_label(
                checkpoint, objective, validation_objective, run_log,
                model_label=variant_legend_label(variant)),
        ))
        variants.append(variant)

    axes['CH4'].legend(
        handles=handles, fancybox=False, edgecolor='k', loc='lower right')
    plot_best.format_fit_axes(reference, sample, axes)
    plot_best.finish_fit_figure(
        figure, comparison_suptitle(reference, variants))
    return figure, axes


def pathway_h2_consumption_rate(pathway_model, run_log, pathway_name):
    """Return the positive H2 uptake rate for one pathway in a replay.

    ``*_v`` is the pathway's biomass-specific rate. Multiplying it by the
    corresponding microbial biomass and the absolute H2 stoichiometric
    coefficient gives the amount of H2 consumed per dry-weight and day.
    """
    pathway = next((candidate for candidate in pathway_model.contributing_pathways
                    if candidate.__class__.__name__ == pathway_name), None)
    if pathway is None:
        raise ValueError(f'Model does not contain the {pathway_name} pathway.')

    h2_coefficient = -float(pathway.stoichiometry[system.index('H2')])
    if h2_coefficient <= 0:
        raise ValueError(f'{pathway_name} does not consume H2.')

    try:
        time, specific_rate = run_log[f'{pathway_name}_v']
        biomass_time, biomass = run_log[f'M_{pathway_name}']
    except KeyError as error:
        raise ValueError(
            f'Model did not log the values required for {pathway_name} H2 consumption.'
        ) from error

    time = np.asarray(time, dtype=float)
    specific_rate = np.asarray(specific_rate, dtype=float)
    biomass_time = np.asarray(biomass_time, dtype=float)
    biomass = np.asarray(biomass, dtype=float)
    order = np.argsort(biomass_time)
    biomass_on_rate_grid = np.interp(
        time, biomass_time[order], biomass[order])
    return time, h2_coefficient * biomass_on_rate_grid * specific_rate


def h2_consumption_label(variant, pathway):
    """Return a TeX legend label that identifies both model and pathway."""
    if variant == 'C':
        model = r'\mathrm{model\ C\ (A\ parameters)}'
    else:
        model = rf'\mathrm{{model\ {variant}}}'
    return rf'${model}:\ \mathrm{{{pathway}}}$'


def plot_h2_consumption_comparison(specifications, model_runs=None):
    """Compare H2 uptake by A's Homo/Hydro and C's Homo pathways.

    Model C is the no-Hydro replay using model A's fitted parameters. Showing
    all three pathway fluxes on one axes makes the pathway removed in C and
    Homo's response directly comparable.
    """
    specifications = sorted(specifications, key=lambda spec: spec.variant)
    _validate_comparison(specifications)
    variants = [specification.variant for specification in specifications]
    if set(variants) != {'A', 'C'}:
        raise ValueError('The H2-consumption comparison requires variants A and C.')
    if model_runs is None:
        model_runs = run_variants(specifications)

    reference = specifications[0].checkpoint
    figure, axes = plot.create_figure(plot.FigureSpec(
        axis_names=('h2_consumption',),
        figsize=(4, 4),
    ))
    axis = axes['h2_consumption']
    curves = (('A', 'Homo'), ('A', 'Hydro'), ('C', 'Homo'))
    for variant, pathway in curves:
        pathway_model, _, run_log = model_runs[variant]
        time, values = pathway_h2_consumption_rate(
            pathway_model, run_log, pathway)
        axis.plot(
            time, values,
            color=plot.pathway_color(pathway),
            linestyle=plot.MODEL_VARIANT_LINESTYLES[variant],
            linewidth=1.5,
            label=h2_consumption_label(variant, pathway),
        )
    axis.set_title(r'$\mathrm{H_2\ consumption\ rate}$')
    axis.set_xlabel(r'$\mathrm{time\ [d]}$')
    axis.set_ylabel(
        r'$\mathrm{H_2\ consumption\ [\mu mol/g\ dry\ weight/d]}$')
    plot.format_ax(axis)
    axis.legend(fancybox=False, edgecolor='k', loc='best')
    figure.suptitle(
        comparison_suptitle(reference, variants, r'$\mathrm{H_2\ consumption}$'))
    figure.subplots_adjust(left=0.25, right=0.95, bottom=0.16, top=0.78)
    return figure, axis


def _specification_by_variant(specifications, variant):
    """Return one named replay specification or raise a helpful error."""
    for specification in specifications:
        if specification.variant == variant:
            return specification
    raise ValueError(f'Model {variant} is not available for this comparison.')


def _replica_predictions(specifications, model_runs, replica, variants=('A', 'B')):
    """Replay selected variants on one replica at its measurement time points."""
    predictions = {}
    for variant in variants:
        specification = _specification_by_variant(specifications, variant)
        pathway_model, _, _ = model_runs[variant]
        predictions[variant] = (
            specification,
            pathway_model,
            pathway_model.predict(replica),
        )
    return predictions


def relative_percent_difference(reference, comparison):
    """Return the signed percentage change from *reference* to *comparison*."""
    reference = np.asarray(reference, dtype=float)
    comparison = np.asarray(comparison, dtype=float)
    result = np.full(reference.shape, np.nan)
    usable = (np.isfinite(reference) & np.isfinite(comparison)
              & (np.abs(reference) > np.finfo(float).tiny))
    result[usable] = 100 * (comparison[usable] - reference[usable]) / reference[usable]
    return result


def methane_log_relative_change(reference, comparison):
    """Return a symmetric, percent-like logarithmic change for positive CH4.

    The value is ``100 ln(B / A)``. It is approximately the ordinary percent
    change for small differences, while a doubling and a halving have equal
    magnitude with opposite signs.
    """
    reference = np.asarray(reference, dtype=float)
    comparison = np.asarray(comparison, dtype=float)
    result = np.full(reference.shape, np.nan)
    usable = (np.isfinite(reference) & np.isfinite(comparison)
              & (reference > 0) & (comparison > 0))
    result[usable] = 100 * np.log(comparison[usable] / reference[usable])
    return result


def ab_replica_detail_data(specifications, model_runs, replica):
    """Replay A/B on one replica and calculate their comparable differences."""
    _validate_comparison(specifications)
    predictions = _replica_predictions(specifications, model_runs, replica)
    _, _, a_log = predictions['A']
    _, _, b_log = predictions['B']
    differences = {}
    for pool in ('CO2', 'CH4'):
        a_time, a_values = a_log[pool]
        b_time, b_values = b_log[pool]
        b_on_a_grid = np.interp(a_time, b_time, b_values)
        if pool == 'CO2':
            values = relative_percent_difference(a_values, b_on_a_grid)
        else:
            values = methane_log_relative_change(a_values, b_on_a_grid)
        differences[pool] = (np.asarray(a_time, dtype=float), values)

    scores = {}
    for variant, (specification, pathway_model, run_log) in predictions.items():
        objective = _validation_objective(
            specification.run_config, pathway_model, replica)
        scores[variant] = objective.R2(run_log)
    return {
        'predictions': predictions,
        'differences': differences,
        'scores': scores,
    }


def _nice_symmetric_limit(value):
    """Round a positive data extent outward for an easily read zero-centred axis."""
    value = max(float(value), 1.0)
    scale = 10 ** np.floor(np.log10(value))
    return np.ceil(value / scale) * scale


def shared_difference_limits(detail_data):
    """Choose robust shared A/B-difference limits, one per measured pool.

    A representative within-replica 95-percent extent defines the common
    scale. A replica whose difference exceeds that scale receives its own
    expanded limit when it is plotted, rather than flattening every other
    replica's comparison.
    """
    limits = {}
    for pool in ('CO2', 'CH4'):
        extents = []
        for detail in detail_data:
            values = np.asarray(detail['differences'][pool][1], dtype=float)
            finite = np.abs(values[np.isfinite(values)])
            if finite.size:
                extents.append(np.percentile(finite, 95))
        extent = np.median(extents) if extents else 1.0
        limits[pool] = _nice_symmetric_limit(extent)
    return limits


def replica_difference_limits(detail, shared_limits):
    """Use shared limits unless one replica is a clear, readable outlier."""
    limits = {}
    for pool, shared_limit in shared_limits.items():
        values = np.asarray(detail['differences'][pool][1], dtype=float)
        finite = np.abs(values[np.isfinite(values)])
        extent = np.max(finite) if finite.size else 0.0
        limit = shared_limit
        if extent > 1.5 * shared_limit:
            limit = _nice_symmetric_limit(extent)
        limits[pool] = (-limit, limit)
    return limits


def plot_ab_replica_detail(specifications, model_runs, replica, *, detail=None,
                           difference_limits=None):
    """Plot A/B fits and the prediction difference for one replica.

    The lower CO2 panel is B's percent change relative to A. The lower CH4
    panel is ``100 ln(B / A)``, its symmetric log-scale equivalent. Both have
    a zero reference, so their sign directly says which model predicts more.
    """
    if detail is None:
        detail = ab_replica_detail_data(specifications, model_runs, replica)
    predictions = detail['predictions']
    a_specification, a_model, a_log = predictions['A']
    b_specification, b_model, b_log = predictions['B']
    sample_limits = plot_best.sample_y_limits(replica.sample)

    figure, grid = plt.subplots(2, 2, figsize=(10, 7), sharex='col')
    top_axes = dict(zip(('CO2', 'CH4'), grid[0]))
    difference_axes = dict(zip(('CO2', 'CH4'), grid[1]))
    for pool in ('CO2', 'CH4'):
        top_axis = top_axes[pool]
        observed_time, observed_values = getattr(replica, pool)()
        top_axis.plot(
            observed_time, observed_values, marker='o', markersize=3,
            linestyle='None', color=plot.pool_color(pool), clip_on=False)
        a_time, a_values = a_log[pool]
        b_time, b_values = b_log[pool]
        top_axis.plot(a_time, a_values, color=plot.pool_color(pool),
                      linestyle=plot.MODEL_VARIANT_LINESTYLES['A'], linewidth=1.5)
        top_axis.plot(b_time, b_values, color=plot.pool_color(pool),
                      linestyle=plot.MODEL_VARIANT_LINESTYLES['B'], linewidth=1.5)
        top_axis.set_title(rf'$\mathrm{{{pool[:-1]}_{pool[-1]}}}$')
        top_axis.set_ylabel(r'$\mathrm{substance\ [\mu mol/g\ dry\ weight]}$')
        plot.format_ax(
            top_axis,
            log_scale='log' in a_specification.run_config['objective']['transform'][pool],
        )
        if pool in sample_limits:
            top_axis.set_ylim(sample_limits[pool])

        difference_axis = difference_axes[pool]
        difference_time, difference = detail['differences'][pool]
        if pool == 'CO2':
            difference_axis.set_ylabel(r'$100(B-A)/A\ [\%]$')
        else:
            difference_axis.set_ylabel(r'$100\ln(B/A)$')
        difference_axis.plot(difference_time, difference, color=plot.pool_color(pool),
                             linewidth=1.5)
        difference_axis.axhline(0, color='k', linewidth=1, linestyle='--')
        difference_axis.set_xlabel(r'$\mathrm{time\ [d]}$')
        plot.format_ax(difference_axis)
        if difference_limits is not None:
            difference_axis.set_ylim(difference_limits[pool])

    top_axes['CH4'].legend(handles=[
        Line2D([], [], color='k', marker='o', linestyle='None', markersize=3,
               label='incubation data'),
        Line2D([], [], color='k', linestyle=plot.MODEL_VARIANT_LINESTYLES['A'],
               label='model A'),
        Line2D([], [], color='k', linestyle=plot.MODEL_VARIANT_LINESTYLES['B'],
               label='model B'),
    ], fancybox=False, edgecolor='k', loc='lower right')
    figure.suptitle(
        f'{replica.sample.sample_name} | replica {replica.replica_number} | '
        'models A and B')
    figure.subplots_adjust(left=.11, right=.97, bottom=.11, top=.84,
                           hspace=.35, wspace=.28)

    return figure, detail['scores']


def plot_ab_score_summary(rows):
    """Summarise A/B evaluation delta R-squared over every selected replica."""
    if not rows:
        return None
    rows = sorted(rows, key=lambda row: (row['sample'], row['replica'], row['fit_mode']))
    delta_r2 = np.asarray([row['r2_b'] - row['r2_a'] for row in rows])
    positions = np.arange(len(rows))
    width = max(8, .38 * len(rows) + 2)
    figure, axis = plt.subplots(figsize=(width, 5))
    colors = np.where(delta_r2 >= 0, 'tab:blue', 'tab:red')
    axis.vlines(positions, 0, delta_r2, colors=colors, linewidth=1.5)
    for position, value, color, row in zip(positions, delta_r2, colors, rows):
        marker = 'o' if row['fit_mode'] == 'split' else 's'
        axis.plot(position, value, marker=marker, color=color, markersize=5)
    axis.axhline(0, color='k', linewidth=1, linestyle='--')
    median = np.median(delta_r2)
    axis.axhline(median, color='tab:purple', linewidth=1.2, linestyle=':',
                 label=rf'median $\Delta R^2 = {median:.2f}$')
    axis.set_xticks(positions, [
        f"{row['sample']}-{row['replica']}\n{row['fit_mode']}" for row in rows],
        rotation=90)
    axis.set_ylabel(r'$\Delta R^2 = R^2_B - R^2_A$')
    axis.set_title('A vs B evaluation difference by replica')
    axis.legend(handles=[
        Line2D([], [], marker='o', color='k', linestyle='None', label='split validation'),
        Line2D([], [], marker='s', color='k', linestyle='None', label='single fit'),
        Line2D([], [], color='tab:purple', linestyle=':',
               label=rf'median $\Delta R^2 = {median:.2f}$'),
    ], fancybox=False, edgecolor='k', loc='best')
    plot.format_ax(axis)
    figure.subplots_adjust(left=.10, right=.98, bottom=.27, top=.90)
    return figure


def plot_ab_r2_scatter(rows):
    """Compare per-replica A/B evaluation scores against the equality line."""
    if not rows:
        return None
    scores_a = np.asarray([row['r2_a'] for row in rows], dtype=float)
    scores_b = np.asarray([row['r2_b'] for row in rows], dtype=float)
    usable = np.isfinite(scores_a) & np.isfinite(scores_b)
    if not np.any(usable):
        return None

    figure, axis = plt.subplots(figsize=(5, 5))
    markers = {'split': 'o', 'single': 's'}
    in_display_range = ((scores_a >= 0) & (scores_a <= 1)
                        & (scores_b >= 0) & (scores_b <= 1))
    for fit_mode, marker in markers.items():
        chosen = np.asarray([
            row['fit_mode'] == fit_mode for row in rows
        ]) & usable & in_display_range
        if np.any(chosen):
            axis.scatter(scores_a[chosen], scores_b[chosen], marker=marker,
                         color='k', s=16,
                         label=('split validation' if fit_mode == 'split'
                                else 'single fit'))

    limits = (0, 1)
    axis.plot(limits, limits, color='k', linestyle='--', linewidth=1,
              label=r'$R^2_\mathrm{A} = R^2_\mathrm{B}$')
    axis.set_xlim(limits)
    axis.set_ylim(limits)
    axis.set_aspect('equal', adjustable='box')
    axis.set_xlabel(r'$R^2_\mathrm{A}$')
    axis.set_ylabel(r'$R^2_\mathrm{B}$')
    cropped = int(np.count_nonzero(usable & ~in_display_range))
    nonfinite = int(np.count_nonzero(~usable))
    crop_label = f'{cropped} paired point' + ('' if cropped == 1 else 's')
    crop_label += r' outside $[0, 1]$ cropped'
    if nonfinite:
        crop_label += f'; {nonfinite} non-finite point(s) omitted'
    axis.text(.5, .02, crop_label, transform=axis.transAxes,
              ha='center', va='bottom', fontsize=8)
    axis.set_title('A versus B validation score by replica')
    axis.legend(fancybox=False, edgecolor='k', loc='best')
    plot.format_ax(axis)
    figure.subplots_adjust(left=.18, right=.96, bottom=.15, top=.90)
    return figure


def comparison_target(output_directory, group):
    """Return a separate, deterministic target directory for one comparison."""
    loss_id, fit_mode, replica = group
    return Path(output_directory) / f'{loss_id}_{fit_mode}' / '-'.join(replica)


def save_comparisons(results_directory=RESULTS_DIRECTORY, output_directory=None,
                     variants=DEFAULT_VARIANTS, **filters):
    """Write variant fit and A/C H2-consumption figures; return their paths."""
    output_directory = (Path(PROJECT_DIRECTORY) / 'variant_comparisons'
                        if output_directory is None else Path(output_directory))
    checkpoints = plot_best.filter_checkpoints(
        plot_best.collect_best(results_directory), **filters)
    groups = comparison_groups(checkpoints, variants)
    dataset = data.get_data_before_carex()
    paths = []
    ab_detail_jobs = []
    for group, by_variant in sorted(groups.items()):
        specifications = variant_run_specs(by_variant, variants)
        model_runs = run_variants(specifications) if len(specifications) >= 2 else None
        directory = comparison_target(output_directory, group)
        if 'A' in by_variant and 'B' in by_variant:
            a_checkpoint = by_variant['A']
            b_checkpoint = by_variant['B']
            if a_checkpoint.source is None or b_checkpoint.source is None:
                raise ValueError('Parameter comparison requires checkpoint source paths.')
            directory.mkdir(parents=True, exist_ok=True)
            path = directory / '00_parameters_A_reference.png'
            compare_checkpoints.plot_comparison(
                [a_checkpoint.source, b_checkpoint.source], path,
                reference_index=0, sort='group', show_constant=True,
                show_r2=False)
            paths.append(path)
        if len(specifications) >= 2:
            figure, _ = plot_fit_comparison(
                specifications, dataset, model_runs=model_runs)
            selected_variants = [specification.variant for specification in specifications]
            directory.mkdir(parents=True, exist_ok=True)
            path = directory / ('01_fit_' + '_vs_'.join(selected_variants) + '.png')
            figure.savefig(path, dpi=300)
            plt.close(figure)
            paths.append(path)

        by_variant_specification = {
            specification.variant: specification for specification in specifications}
        if {'A', 'C'}.issubset(by_variant_specification):
            selected = [by_variant_specification['A'], by_variant_specification['C']]
            figure, _ = plot_h2_consumption_comparison(selected, model_runs=model_runs)
            directory.mkdir(parents=True, exist_ok=True)
            path = directory / '02_h2_consumption_A_vs_C.png'
            figure.savefig(path, dpi=300)
            plt.close(figure)
            paths.append(path)

        if len(specifications) >= 2:
            diagnostic_title = comparison_suptitle(
                specifications[0].checkpoint,
                [specification.variant for specification in specifications],
                r'$\mathrm{mechanistic\ pathway\ diagnostics}$')
            for name, figure in variant_diagnostics.diagnostic_figures(
                    specifications, model_runs, diagnostic_title):
                path = directory / f'{name}.png'
                figure.savefig(path, dpi=300)
                plt.close(figure)
                paths.append(path)

        if {'A', 'B'}.issubset(by_variant_specification):
            reference = by_variant_specification['A'].checkpoint
            _, validation_replica, _ = plot_best.fit_replicas(reference, dataset)
            detail = ab_replica_detail_data(
                [by_variant_specification['A'], by_variant_specification['B']],
                model_runs, validation_replica)
            ab_detail_jobs.append({
                'directory': directory,
                'reference': reference,
                'specifications': [by_variant_specification['A'],
                                   by_variant_specification['B']],
                'model_runs': model_runs,
                'replica': validation_replica,
                'detail': detail,
            })

    shared_limits = shared_difference_limits(
        [job['detail'] for job in ab_detail_jobs])
    ab_score_rows = []
    for job in ab_detail_jobs:
        figure, scores = plot_ab_replica_detail(
            job['specifications'], job['model_runs'], job['replica'],
            detail=job['detail'],
            difference_limits=replica_difference_limits(
                job['detail'], shared_limits))
        path = job['directory'] / '03_ab_replica_detail.png'
        figure.savefig(path, dpi=300)
        plt.close(figure)
        paths.append(path)
        reference = job['reference']
        ab_score_rows.append({
            'sample': reference.replica[0],
            'replica': reference.replica[1],
            'fit_mode': reference.fit_mode,
            'r2_a': scores['A'],
            'r2_b': scores['B'],
        })
    summary = plot_ab_score_summary(ab_score_rows)
    if summary is not None:
        output_directory.mkdir(parents=True, exist_ok=True)
        path = output_directory / '04_ab_dataset_delta_r2.png'
        summary.savefig(path, dpi=300)
        plt.close(summary)
        paths.append(path)
    scatter = plot_ab_r2_scatter(ab_score_rows)
    if scatter is not None:
        output_directory.mkdir(parents=True, exist_ok=True)
        path = output_directory / '05_ab_dataset_r2_scatter.png'
        scatter.savefig(path, dpi=300)
        plt.close(scatter)
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
        checkpoints = plot_best.filter_checkpoints(
            plot_best.collect_best(args.results_directory), **filters)
        for group, by_variant in sorted(comparison_groups(
                checkpoints, args.variants).items()):
            available = [specification.variant for specification
                         in variant_run_specs(by_variant, args.variants)]
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
