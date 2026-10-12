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
    key = f'{pathway_name}_thermodynamic_factor'
    if key in run_log.keys():
        return run_log[key]
    time, values = run_log[f'{pathway_name}_v']
    return time, np.ones_like(np.asarray(values, dtype=float))


def plot_ab_thermodynamic_effect(specifications, model_runs, suptitle):
    """Separate thermodynamic, kinetic, and biomass effects in A versus B."""
    variants = {specification.variant for specification in specifications}
    if not {'A', 'B'}.issubset(variants):
        return None
    a_model, _, _ = model_runs['A']
    pathways = available_pathways(a_model, LIMITATION_PATHWAYS)
    if not pathways:
        return None
    figure, axes = plt.subplots(3, len(pathways), figsize=(3 * len(pathways), 8),
                                sharex='col', squeeze=False)
    row_titles = (r'$f_\mathrm{thermo}$', r'$v$', r'$M v$')
    for column, pathway in enumerate(pathways):
        for variant in ('A', 'B'):
            pathway_model, _, run_log = model_runs[variant]
            if pathway not in available_pathways(pathway_model, (pathway,)):
                continue
            factor_time, factor_values = _thermodynamic_factor(run_log, pathway)
            v_time, v_values = run_log[f'{pathway}_v']
            activity_time, activity_values = pathway_activity(
                pathway_model, run_log, pathway)
            for axis, time, values in zip(axes[:, column],
                                          (factor_time, v_time, activity_time),
                                          (factor_values, v_values, activity_values)):
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
    """Show how removal of Hydro redistributes A's pools and Homo activity."""
    variants = {specification.variant for specification in specifications}
    if not {'A', 'C'}.issubset(variants):
        return None
    figure, axes = plt.subplots(2, 2, figsize=(8, 7), sharex=True)
    a_model, _, a_log = model_runs['A']
    c_model, _, c_log = model_runs['C']
    for axis, pool in zip(axes.flat[:3], ('H2', 'CO2', 'Acetate')):
        time, a_values = a_log[pool]
        c_time, c_values = c_log[pool]
        c_on_a_grid = np.interp(time, c_time, c_values)
        axis.plot(time, _relative_percent(a_values, c_on_a_grid),
                  color='k', linestyle=plot.MODEL_VARIANT_LINESTYLES['C'], linewidth=1.5)
        axis.axhline(0, color='k', linestyle='--', linewidth=1)
        axis.set_title(rf'$\mathrm{{{pool}:\ model\ C\ relative\ to\ A}}$')
        axis.set_ylabel(r'$100(C-A)/A\ [\%]$')
    axis = axes.flat[3]
    for variant, pathway_model, run_log in (
            ('A', a_model, a_log), ('C', c_model, c_log)):
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
    """Compare normalized pathway limitation signatures for A and B."""
    variants = {specification.variant for specification in specifications}
    if not {'A', 'B'}.issubset(variants):
        return None
    a_model, _, a_log = model_runs['A']
    pathways = available_pathways(a_model, LIMITATION_PATHWAYS)
    if not pathways:
        return None
    max_time = max(float(np.max(a_log[pathway + '_v'][0])) for pathway in pathways)
    edges = np.linspace(0, max_time, 7)
    diagnostics = []
    labels = []
    for pathway in pathways:
        for diagnostic in ('thermodynamic factor', 'MM factor', 'biomass', 'activity'):
            diagnostics.append((pathway, diagnostic))
            labels.append(rf'$\mathrm{{{pathway}:\ {diagnostic}}}$')

    values = {variant: [] for variant in ('A', 'B')}
    for pathway, diagnostic in diagnostics:
        for variant in ('A', 'B'):
            pathway_model, _, run_log = model_runs[variant]
            if diagnostic == 'thermodynamic factor':
                time, series = _thermodynamic_factor(run_log, pathway)
            elif diagnostic == 'MM factor':
                time, series = run_log[f'{pathway}_MM']
            elif diagnostic == 'biomass':
                time, series = run_log[_microbe_pool(pathway)]
            else:
                time, series = pathway_activity(pathway_model, run_log, pathway)
            values[variant].append(_window_means(time, series, edges))

    normalized = {}
    for variant in ('A', 'B'):
        normalized[variant] = np.asarray(values[variant], dtype=float)
    combined = np.stack((normalized['A'], normalized['B']))
    lower = np.nanmin(combined, axis=(0, 2))[:, None]
    upper = np.nanmax(combined, axis=(0, 2))[:, None]
    scale = np.where(upper > lower, upper - lower, 1)
    for variant in ('A', 'B'):
        normalized[variant] = (normalized[variant] - lower) / scale

    figure, axes = plt.subplots(1, 2, figsize=(10, max(5, .36 * len(labels) + 2)),
                                sharey=True)
    for axis, variant in zip(axes, ('A', 'B')):
        image = axis.imshow(normalized[variant], aspect='auto', vmin=0, vmax=1,
                            cmap='viridis', interpolation='nearest')
        axis.set_title(rf'$\mathrm{{model\ {variant}}}$')
        axis.set_xticks(np.arange(len(edges) - 1),
                         [f'{left:.0f}-{right:.0f}'
                          for left, right in zip(edges[:-1], edges[1:])],
                         rotation=45, ha='right')
        axis.set_xlabel(r'$\mathrm{time\ window\ [d]}$')
        plot.format_ax(axis)
    axes[0].set_yticks(np.arange(len(labels)), labels, fontsize=8)
    colorbar = figure.colorbar(image, ax=axes, fraction=.04, pad=.03)
    colorbar.set_label(r'$\mathrm{within\ diagnostic\ normalized\ value}$')
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
