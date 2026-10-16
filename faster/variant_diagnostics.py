"""Mechanistic diagnostics for replayed model-variant comparisons.

The functions here deliberately work from the replay logs and pathway
stoichiometry rather than from fitting scores.  They therefore expose how a
variant changes pathway interactions after it has been fitted.
"""

from __future__ import annotations

import math

import matplotlib.pyplot as plt
import numpy as np

import plot
import system
from pathways import pathway_color


H2_PATHWAYS = ('Homo', 'Hydro')
ACETATE_PATHWAYS = ('Fe3', 'Aceto')
LIMITATION_PATHWAYS = ('Homo', 'Hydro', 'Aceto', 'Fe3')
MISSING_DIAGNOSTIC_COLOR = '#d9d9d9'


def _specification(specifications, variant):
    return next(specification for specification in specifications
                if specification.variant == variant)


def _pathway(pathway_model, pathway_name):
    result = next((candidate for candidate in pathway_model.contributing_pathways
                   if candidate.__class__.__name__ == pathway_name), None)
    if result is None:
        raise ValueError(f'Model does not contain the {pathway_name} pathway.')
    return result


def _microbe_pool(pathway_name):
    return {'Homo': 'M_Homo', 'Hydro': 'M_Hydro', 'Fe3': 'M_Fe3',
            'Aceto': 'M_Ac', 'Fermentation': 'M_Ferm',
            'Hydrolysis': 'M_Ferm'}[pathway_name]


def available_pathways(pathway_model, names):
    """Return requested pathway names that are active in a model."""
    present = {candidate.__class__.__name__
               for candidate in pathway_model.contributing_pathways}
    return tuple(name for name in names if name in present)


def pathway_activity(pathway_model, run_log, pathway_name):
    """Return ``M v``: pathway flux before applying stoichiometry."""
    _pathway(pathway_model, pathway_name)
    time, specific_rate = run_log[f'{pathway_name}_v']
    biomass_time, biomass = run_log[_microbe_pool(pathway_name)]
    time = np.asarray(time, dtype=float)
    biomass_time = np.asarray(biomass_time, dtype=float)
    order = np.argsort(biomass_time)
    biomass_on_rate_grid = np.interp(time, biomass_time[order],
                                     np.asarray(biomass, dtype=float)[order])
    return time, biomass_on_rate_grid * np.asarray(specific_rate, dtype=float)


def pathway_pool_rate(pathway_model, run_log, pathway_name, pool):
    """Return one pathway's signed contribution to a pool's change rate."""
    pathway = _pathway(pathway_model, pathway_name)
    time, activity = pathway_activity(pathway_model, run_log, pathway_name)
    coefficient = float(pathway.pathway_vector[system.index(pool)])
    return time, activity * coefficient


def pathway_consumption_rate(pathway_model, run_log, pathway_name, pool):
    """Return positive net consumption of *pool* by a pathway."""
    time, rate = pathway_pool_rate(pathway_model, run_log, pathway_name, pool)
    #if np.all(rate >= 0):
    #    raise ValueError(f'{pathway_name} does not consume {pool}.')
    return time, -rate


def _curve_label(variant, pathway):
    return rf'$\mathrm{{model\ {variant}:\ {pathway}}}$'


def _finish(figure, axes, suptitle, *, legend=True):
    for axis in np.ravel(axes):
        plot.format_ax(axis)
        axis.set_xlabel(r'$\mathrm{time\ [d]}$')
        if legend:
            handles, _ = axis.get_legend_handles_labels()
            if handles:
                axis.legend(fancybox=False, edgecolor='k', loc='best', fontsize=8)
    figure.suptitle(suptitle)


def plot_flux_partition(specifications, model_runs, suptitle):
    """Plot H2 and acetate consumption partitioned by pathway and variant."""
    figure, axes = plt.subplots(1, 2, figsize=(9, 4), sharex=False)
    for axis, pool, pathways in zip(axes, ('H2', 'Acetate'),
                                    (H2_PATHWAYS, ACETATE_PATHWAYS)):
        for specification in specifications:
            variant = specification.variant
            pathway_model, _, run_log = model_runs[variant]
            for pathway in available_pathways(pathway_model, pathways):
                try:
                    time, values = pathway_consumption_rate(
                        pathway_model, run_log, pathway, pool)
                except ValueError:
                    continue
                axis.plot(time, values, color=pathway_color(pathway),
                          linestyle=plot.MODEL_VARIANT_LINESTYLES[variant], linewidth=1.4,
                          label=_curve_label(variant, pathway))
        axis.set_title(rf'$\mathrm{{{pool}\ consumption}}$')
        axis.set_ylabel(
            rf'$\mathrm{{{pool}\ consumption\ [\mu mol/g\ dry\ weight/d]}}$')
    _finish(figure, axes, suptitle)
    figure.subplots_adjust(left=.10, right=.98, bottom=.16, top=.78, wspace=.32)
    return figure


def _thermodynamic_factor(run_log, pathway_name):
    """Return a logged thermodynamic factor, or ``None`` when it is absent."""
    key = f'{pathway_name}_thermodynamic_factor'
    if key in run_log.keys():
        return run_log[key]
    return None


def _logged_series(run_log, key):
    """Return one logged series, distinguishing missing values from zero."""
    return run_log[key] if key in run_log.keys() else None


def _missing_diagnostic_note(axis, variant):
    """Annotate a line diagnostic whose series was not written to the log."""
    count = sum(text.get_gid() == 'missing-diagnostic'
                for text in axis.texts)
    axis.text(.02, .96 - .12 * count,
              f'model {variant}: not logged',
              color='0.45', transform=axis.transAxes,
              ha='left', va='top', fontsize=8, gid='missing-diagnostic')


def plot_ab_thermodynamic_effect(specifications, model_runs, suptitle):
    """Separate thermodynamic, kinetic, and biomass effects by selection."""
    variants = [specification.variant for specification in specifications]
    if len(variants) < 2:
        return None
    reference_variant = 'A' if 'A' in model_runs else variants[0]
    reference_model, _, _ = model_runs[reference_variant]
    pathways = available_pathways(reference_model, LIMITATION_PATHWAYS)
    if not pathways:
        return None
    figure, axes = plt.subplots(3, len(pathways), figsize=(3 * len(pathways), 8),
                                sharex='col', squeeze=False)
    row_titles = (r'$f_\mathrm{thermo}$', r'$v$', r'$M v$')
    for column, pathway in enumerate(pathways):
        for variant in variants:
            pathway_model, _, run_log = model_runs[variant]
            if pathway not in available_pathways(pathway_model, (pathway,)):
                continue
            factor_series = _thermodynamic_factor(run_log, pathway)
            v_time, v_values = run_log[f'{pathway}_v']
            activity_time, activity_values = pathway_activity(
                pathway_model, run_log, pathway)
            series = (factor_series, (v_time, v_values),
                      (activity_time, activity_values))
            for axis, values in zip(axes[:, column], series):
                if values is None:
                    _missing_diagnostic_note(axis, variant)
                    continue
                time, values = values
                axis.plot(time, values, color=pathway_color(pathway),
                          linestyle=plot.MODEL_VARIANT_LINESTYLES[variant], linewidth=1.3,
                          label=rf'$\mathrm{{model\ {variant}}}$')
        axes[0, column].set_title(rf'$\mathrm{{{pathway}\ pathway}}$')
        axes[0, column].set_ylim(0, 1)
    for row, label in enumerate(row_titles):
        axes[row, 0].set_ylabel(label)
    _finish(figure, axes, suptitle)
    figure.subplots_adjust(left=.08, right=.99, bottom=.09, top=.88,
                            hspace=.28, wspace=.25)
    return figure


def _relative_percent(reference, comparison):
    reference = np.asarray(reference, dtype=float)
    comparison = np.asarray(comparison, dtype=float)
    result = np.full(reference.shape, np.nan)
    usable = (np.isfinite(reference) & np.isfinite(comparison)
              & (np.abs(reference) > np.finfo(float).tiny))
    result[usable] = 100 * (comparison[usable] - reference[usable]) / reference[usable]
    return result


def plot_ac_pathway_removal_response(specifications, model_runs, suptitle):
    """Show the selected A-derived counterfactual responses to model A."""
    variants = [specification.variant for specification in specifications]
    counterfactuals = [variant for variant in ('C', 'D') if variant in variants]
    if 'A' not in variants or not counterfactuals:
        return None
    figure, axes = plt.subplots(2, 2, figsize=(8, 7), sharex=True)
    a_model, _, a_log = model_runs['A']
    for axis, pool in zip(axes.flat[:3], ('H2', 'CO2', 'Acetate')):
        time, a_values = a_log[pool]
        for variant in counterfactuals:
            _, _, counterfactual_log = model_runs[variant]
            comparison_time, comparison_values = counterfactual_log[pool]
            comparison_on_a_grid = np.interp(
                time, comparison_time, comparison_values)
            axis.plot(
                time, _relative_percent(a_values, comparison_on_a_grid),
                color='k', linestyle=plot.MODEL_VARIANT_LINESTYLES[variant],
                linewidth=1.5, label=rf'$\mathrm{{model\ {variant}}}$')
        axis.axhline(0, color='k', linestyle='--', linewidth=1)
        axis.set_title(rf'$\mathrm{{{pool}:\ counterfactual\ relative\ to\ A}}$')
        axis.set_ylabel(r'$100(X-A)/A\ [\%]$')
    axis = axes.flat[3]
    for variant in ('A', *counterfactuals):
        pathway_model, _, run_log = model_runs[variant]
        time, values = pathway_consumption_rate(pathway_model, run_log, 'Homo', 'H2')
        axis.plot(time, values, color=pathway_color('Homo'),
                  linestyle=plot.MODEL_VARIANT_LINESTYLES[variant], linewidth=1.5,
                  label=rf'$\mathrm{{model\ {variant}}}$')
    axis.set_title(r'$\mathrm{Homo\ H_2\ consumption}$')
    axis.set_ylabel(r'$\mathrm{H_2\ consumption\ [\mu mol/g\ dry\ weight/d]}$')
    _finish(figure, axes, suptitle)
    figure.subplots_adjust(left=.12, right=.98, bottom=.11, top=.84,
                            hspace=.32, wspace=.30)
    return figure


def plot_substrate_flux_phase(specifications, model_runs, suptitle):
    """Relate H2/acetate availability directly to each consuming pathway flux."""
    figure, axes = plt.subplots(1, 2, figsize=(9, 4))
    for axis, pool, pathways in zip(axes, ('H2', 'Acetate'),
                                    (H2_PATHWAYS, ACETATE_PATHWAYS)):
        for specification in specifications:
            variant = specification.variant
            pathway_model, _, run_log = model_runs[variant]
            state_time, state_values = run_log[pool]
            for pathway in available_pathways(pathway_model, pathways):
                try:
                    time, rate = pathway_consumption_rate(
                        pathway_model, run_log, pathway, pool)
                except ValueError:
                    continue
                state = np.interp(time, state_time, state_values)
                axis.plot(state, rate, color=pathway_color(pathway),
                          linestyle=plot.MODEL_VARIANT_LINESTYLES[variant], linewidth=1.4,
                          label=_curve_label(variant, pathway))
        axis.set_xlabel(rf'$\mathrm{{{pool}\ [\mu mol/g\ dry\ weight]}}$')
        axis.set_ylabel(
            rf'$\mathrm{{{pool}\ consumption\ [\mu mol/g\ dry\ weight/d]}}$')
        axis.set_title(rf'$\mathrm{{{pool}\ availability\ and\ pathway\ flux}}$')
        plot.format_ax(axis)
        axis.legend(fancybox=False, edgecolor='k', loc='best', fontsize=8)
    figure.suptitle(suptitle)
    figure.subplots_adjust(left=.10, right=.98, bottom=.18, top=.78, wspace=.34)
    return figure


def _cumulative(time, rate):
    time = np.asarray(time, dtype=float)
    rate = np.asarray(rate, dtype=float)
    order = np.argsort(time)
    time, rate = time[order], rate[order]
    time, unique_indices = np.unique(time, return_index=True)
    rate = rate[unique_indices]
    if time.size < 2:
        return 0.0
    return float(np.sum(np.diff(time) * (rate[:-1] + rate[1:]) / 2))


def plot_mass_balance_attribution(specifications, model_runs, suptitle):
    """Attribute final pool turnover to pathway-level cumulative fluxes."""
    pools = ('H2', 'Acetate', 'CO2', 'CH4')
    variants = [specification.variant for specification in specifications]
    pathway_names = sorted({candidate.__class__.__name__
                            for variant in variants
                            for candidate in model_runs[variant][0].contributing_pathways})
    figure, axes = plt.subplots(2, 2, figsize=(9, 7))
    positions = np.arange(len(variants))
    for axis, pool in zip(axes.flat, pools):
        positive = np.zeros(len(variants))
        negative = np.zeros(len(variants))
        for pathway in pathway_names:
            values = []
            for variant in variants:
                pathway_model, _, run_log = model_runs[variant]
                if pathway not in available_pathways(pathway_model, (pathway,)):
                    values.append(0.0)
                    continue
                time, rate = pathway_pool_rate(pathway_model, run_log, pathway, pool)
                values.append(_cumulative(time, rate))
            values = np.asarray(values)
            axis.bar(positions, np.maximum(values, 0), bottom=positive,
                     color=pathway_color(pathway), label=pathway)
            axis.bar(positions, np.minimum(values, 0), bottom=negative,
                     color=pathway_color(pathway))
            positive += np.maximum(values, 0)
            negative += np.minimum(values, 0)
        axis.axhline(0, color='k', linewidth=.8)
        axis.set_xticks(positions, [rf'$\mathrm{{model\ {variant}}}$'
                                    for variant in variants])
        axis.set_title(rf'$\mathrm{{{pool}\ turnover}}$')
        axis.set_ylabel(r'$\mathrm{cumulative\ change\ [\mu mol/g\ dry\ weight]}$')
        plot.format_ax(axis)
    handles, labels = axes.flat[0].get_legend_handles_labels()
    figure.legend(handles, labels, loc='lower center', ncol=min(4, len(handles)),
                  fancybox=False, edgecolor='k')
    figure.suptitle(suptitle)
    figure.subplots_adjust(left=.11, right=.98, bottom=.18, top=.84,
                            hspace=.36, wspace=.32)
    return figure


def _window_means(time, values, edges):
    time = np.asarray(time, dtype=float)
    values = np.asarray(values, dtype=float)
    centres = (edges[:-1] + edges[1:]) / 2
    return np.interp(centres, time, values)


def plot_ab_limitation_fingerprint(specifications, model_runs, suptitle):
    """Compare selected models' normalized pathway limitation signatures."""
    variants = [specification.variant for specification in specifications]
    if len(variants) < 2:
        return None
    reference_variant = 'A' if 'A' in model_runs else variants[0]
    reference_model, _, reference_log = model_runs[reference_variant]
    pathways = available_pathways(reference_model, LIMITATION_PATHWAYS)
    if not pathways:
        return None
    max_time = max(float(np.max(reference_log[pathway + '_v'][0]))
                   for pathway in pathways)
    edges = np.linspace(0, max_time, 7)
    diagnostics = []
    labels = []
    for pathway in pathways:
        for diagnostic in ('thermodynamic factor', 'MM factor', 'biomass', 'activity'):
            diagnostics.append((pathway, diagnostic))
            labels.append(rf'$\mathrm{{{pathway}:\ {diagnostic}}}$')

    values = {variant: [] for variant in variants}
    for pathway, diagnostic in diagnostics:
        for variant in variants:
            pathway_model, _, run_log = model_runs[variant]
            if pathway not in available_pathways(pathway_model, (pathway,)):
                logged_series = None
            elif diagnostic == 'thermodynamic factor':
                logged_series = _thermodynamic_factor(run_log, pathway)
            elif diagnostic == 'MM factor':
                logged_series = _logged_series(run_log, f'{pathway}_MM')
            elif diagnostic == 'biomass':
                logged_series = _logged_series(run_log, _microbe_pool(pathway))
            else:
                try:
                    logged_series = pathway_activity(
                        pathway_model, run_log, pathway)
                except KeyError:
                    logged_series = None
            if logged_series is None:
                values[variant].append(np.full(len(edges) - 1, np.nan))
            else:
                time, series = logged_series
                values[variant].append(_window_means(time, series, edges))

    normalized = {}
    for variant in variants:
        normalized[variant] = np.asarray(values[variant], dtype=float)
    combined = np.stack([normalized[variant] for variant in variants])
    finite = np.isfinite(combined)
    valid_diagnostics = np.any(finite, axis=(0, 2))
    lower = np.min(np.where(finite, combined, np.inf), axis=(0, 2))
    upper = np.max(np.where(finite, combined, -np.inf), axis=(0, 2))
    lower[~valid_diagnostics] = np.nan
    upper[~valid_diagnostics] = np.nan
    # Thermodynamic and Michaelis--Menten factors are already normalized to
    # the physical [0, 1] interval.  Preserve that scale so a constant factor
    # of one is displayed as one, rather than as zero after min--max scaling.
    factor_rows = np.asarray([
        diagnostic in ('thermodynamic factor', 'MM factor')
        for _, diagnostic in diagnostics
    ])
    lower[factor_rows & valid_diagnostics] = 0.
    upper[factor_rows & valid_diagnostics] = 1.
    scale = np.where(upper > lower, upper - lower, 1)[:, None]
    for variant in variants:
        normalized[variant] = (normalized[variant] - lower[:, None]) / scale

    figure, axes = plt.subplots(
        1, len(variants),
        figsize=(5 * len(variants), max(5, .36 * len(labels) + 2)),
        sharey=True, squeeze=False)
    axes = axes[0]
    diagnostic_map = plt.get_cmap('viridis').copy()
    diagnostic_map.set_bad(MISSING_DIAGNOSTIC_COLOR)
    for axis, variant in zip(axes, variants):
        image = axis.imshow(normalized[variant], aspect='auto', vmin=0, vmax=1,
                            cmap=diagnostic_map, interpolation='nearest')
        axis.set_title(rf'$\mathrm{{model\ {variant}}}$')
        axis.set_xticks(np.arange(len(edges) - 1),
                         [f'{left:.0f}-{right:.0f}'
                          for left, right in zip(edges[:-1], edges[1:])],
                         rotation=45, ha='right')
        axis.set_xlabel(r'$\mathrm{time\ window\ [d]}$')
        plot.format_ax(axis)
    axes[0].set_yticks(np.arange(len(labels)), labels, fontsize=8)
    colorbar = figure.colorbar(image, ax=axes.tolist(), fraction=.04, pad=.03)
    colorbar.set_label(
        r'$\mathrm{factor\ value\ or\ within\ diagnostic\ normalized\ value;\ grey:\ not\ logged}$')
    figure.suptitle(suptitle)
    figure.subplots_adjust(left=.30, right=.94, bottom=.20, top=.88, wspace=.08)
    return figure


def diagnostic_figures(specifications, model_runs, suptitle):
    """Return every supported mechanistic figure for one comparison group."""
    figures = [('05_flux_partition', plot_flux_partition(
        specifications, model_runs, suptitle))]
    optional = (
        ('06_ab_thermodynamic_effect', plot_ab_thermodynamic_effect),
        ('07_ac_pathway_removal_response', plot_ac_pathway_removal_response),
        ('08_substrate_flux_phase', plot_substrate_flux_phase),
        ('09_mass_balance_attribution', plot_mass_balance_attribution),
        ('10_ab_limitation_fingerprint', plot_ab_limitation_fingerprint),
    )
    for name, builder in optional:
        figure = builder(specifications, model_runs, suptitle)
        if figure is not None:
            figures.append((name, figure))
    return figures
