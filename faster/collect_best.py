"""Collect the best checkpoint for each model, replica, and fit mode.

Running this module writes a data/model comparison plot for every selected
checkpoint below ``best_<loss-id>_<model-id>_<fit-mode>``.
"""

import argparse
import os
from dataclasses import dataclass
from functools import cached_property
from itertools import cycle
import traceback

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

import data
import hashing
import optimizer
import parameters
import plot
import system
from chemistry import GIBBS_MINIMUM
from fit_sample import run
from pathways import pathway_color
from model_variants import model_variant_from_id
from USER_VARIABLES import PROJECT_DIRECTORY, RESULTS_DIRECTORY

GIBBS_ENERGY_SCALE = 1e-3  # Internal reaction energies are J/mol; plots use kJ/mol.
GIBBS_ENERGY_UNIT = r'\mathrm{kJ}\,\mathrm{mol}^{-1}'

# ``00`` is reserved for the data/model fit (CO2 and CH4).  The remaining
# system pools are grouped with the pathway diagnostics below, rather than
# duplicating the gases or a microbe pool in a separate group.
POOL_PLOT_ORDER = ('TOC', 'DOC', 'Acetate', 'Fe3', 'H2', 'H2O', 'Fe2')
PATHWAY_PLOT_ORDER = (
    'Hydrolysis', 'Fermentation', 'Fe3', 'Aceto', 'Hydro', 'Homo',
)
PATHWAY_MICROBE_POOL = {
    'Hydrolysis': 'M_Ferm',
    'Fermentation': 'M_Ferm',
    'Fe3': 'M_Fe3',
    'Aceto': 'M_Ac',
    'Hydro': 'M_Hydro',
    'Homo': 'M_Homo',
}
PATHWAY_V_MAX_PARAMETER = {
    'Hydrolysis': 'Hydrolysis_v_max',
    'Fermentation': 'Ferm_v_max',
    'Fe3': 'Fe3_v_max',
    'Aceto': 'Ac_v_max',
    'Hydro': 'Hydro_v_max',
    'Homo': 'Homo_v_max',
}
PATHWAY_QUANTITY_ORDER = (
    'thermodynamic_factor', 'deltaG_r', 'MM', 'microbe_pool',
)

@dataclass
class Checkpoint:
    """A loaded checkpoint and the metadata used to group it."""

    parameters: parameters.ModelParameters
    loss: float
    run_config: dict
    id: str
    source: str = None

    @cached_property
    def model_id(self):
        return hashing.build_model_id(self.run_config['model'])

    @cached_property
    def loss_id(self):
        return hashing.build_loss_id(self.run_config['objective'])

    @cached_property
    def replica(self):
        chosen = self.run_config['chosen']
        return str(chosen['sample']), str(chosen['validation_replica'])

    @cached_property
    def replica_name(self):
        return '-'.join(self.replica)

    @cached_property
    def fit_mode(self):
        return self.run_config['chosen']['fit_mode']

    @cached_property
    def key(self):
        return '_'.join((self.loss_id, self.model_id))
    
    def __str__(self):
        return self.key

def collect_all(results_directory=RESULTS_DIRECTORY):
    """Load every annotated checkpoint in ``results_directory``."""
    collected = []
    for root, _, files in os.walk(results_directory):
        for filename in files:
            if filename.startswith('.') or filename.endswith('.log') or filename.endswith('.json'):
                continue

            file_path = os.path.join(root, filename)
            try:
                loaded_parameters, loss, run_config = parameters.load_parameter_file(file_path)
            except Exception as ex:
                print(traceback.format_exc())
                input()
            checkpoint_ids = [part for part in filename.split('_')
                              if part.startswith('cp-')]
            if not checkpoint_ids:
                raise ValueError(f'Checkpoint ID missing from {file_path}')
            if not run_config:
                raise ValueError(f'Run configuration missing from {file_path}')

            collected.append(Checkpoint(loaded_parameters, loss, run_config,
                                        checkpoint_ids[0], file_path))
    return collected


def collect_best(results_directory=RESULTS_DIRECTORY):
    """Return the lowest-loss checkpoint per model/loss/replica/fit-mode."""
    best_by_group = {}
    for checkpoint in collect_all(results_directory):
        group = (checkpoint.key, checkpoint.fit_mode, checkpoint.replica)
        best = best_by_group.get(group)
        if best is None or checkpoint.loss < best.loss:
            best_by_group[group] = checkpoint
    return list(best_by_group.values())


def filter_checkpoints(checkpoints, samples=None, model_ids=None, fit_modes=None,
                       replicas=None, loss_ids=None):
    """Keep checkpoints matching every supplied command-line filter.

    Model-ID filters are substrings, so a short model hash can select all of
    its checkpoints without requiring the full generated identifier.
    """
    filters = {
        'sample': None if samples is None else set(map(str, samples)),
        'model_id': None if model_ids is None else set(model_ids),
        'fit_mode': None if fit_modes is None else set(fit_modes),
        'replica': None if replicas is None else set(map(str, replicas)),
        'loss_id': None if loss_ids is None else set(loss_ids),
    }
    selected = []
    for checkpoint in checkpoints:
        sample, replica = checkpoint.replica
        if filters['sample'] is not None and sample not in filters['sample']:
            continue
        if (filters['model_id'] is not None and not any(
                model_id in checkpoint.model_id
                for model_id in filters['model_id'])):
            continue
        if (filters['fit_mode'] is not None and
                checkpoint.fit_mode not in filters['fit_mode']):
            continue
        if filters['replica'] is not None and replica not in filters['replica']:
            continue
        if filters['loss_id'] is not None and checkpoint.loss_id not in filters['loss_id']:
            continue
        selected.append(checkpoint)
    return selected


def plot_target(checkpoint, number=0, quantity=None):
    """Return the directory and basename for a numbered checkpoint plot."""
    directory = os.path.join(
        PROJECT_DIRECTORY,
        '_'.join(('best', checkpoint.key, checkpoint.fit_mode)),
        checkpoint.replica_name,
    )
    filename_parts = [
        f'{number:02d}',
    ]
    if quantity is not None:
        filename_parts.append(quantity)
    filename_parts.extend((
        checkpoint.replica_name, checkpoint.model_id,
        checkpoint.fit_mode, checkpoint.id,
    ))
    return directory, '_'.join(filename_parts)


def _run_log_keys(run_log):
    """Return run-log keys without depending on its private storage."""
    return set(run_log.keys())


def pathway_plot_quantities(run_log):
    """Yield ``(name, time, values)`` in the required plot-file order.

    Quantities added to a pathway log in the future are appended after the
    specified diagnostics for that pathway.  ``v_max`` is shown as a dashed
    reference line on the corresponding ``v`` plot, not as its own plot.
    """
    keys = _run_log_keys(run_log)

    for pool in POOL_PLOT_ORDER:
        if pool in keys:
            yield pool, *run_log[pool]

    for pathway in PATHWAY_PLOT_ORDER:
        prefix = f'{pathway}_'
        pathway_keys = {key for key in keys if key.startswith(prefix)}
        if not pathway_keys:
            continue

        for quantity in PATHWAY_QUANTITY_ORDER:
            if quantity == 'microbe_pool':
                microbe_pool = PATHWAY_MICROBE_POOL[pathway]
                if microbe_pool in keys:
                    yield microbe_pool, *run_log[microbe_pool]
            else:
                name = prefix + quantity
                if name in pathway_keys:
                    yield name, *run_log[name]

        handled = {
            prefix + quantity
            for quantity in PATHWAY_QUANTITY_ORDER
            if quantity != 'microbe_pool'
        }
        # Add every remaining pathway diagnostic (for example ``v``,
        # inhibition, and methane production) in a stable order.
        for name in sorted(pathway_keys - handled):
            yield name, *run_log[name]


def pathway_v_max_values(run_log):
    """Return the configured maximum rate for each active pathway."""
    parameters = getattr(run_log, '_parameters', None)
    if parameters is None:
        return {}

    keys = _run_log_keys(run_log)
    return {
        pathway: float(parameters[PATHWAY_V_MAX_PARAMETER[pathway]])
        for pathway in PATHWAY_PLOT_ORDER
        if any(key.startswith(f'{pathway}_') for key in keys)
    }


def pathway_reaction_rate_curves(run_log, v_max_values):
    """Yield actual and maximum biomass-scaled rates for each pathway.

    ``<pathway>_v`` is a specific rate, whereas the reaction rate affecting
    the model pools is the specific rate multiplied by the pathway's microbial
    biomass. Pathway and pool logs can have different time grids, so biomass
    is interpolated onto the specific-rate grid before both curves are made.
    """
    keys = _run_log_keys(run_log)
    for pathway in PATHWAY_PLOT_ORDER:
        specific_rate_key = f'{pathway}_v'
        microbe_pool = PATHWAY_MICROBE_POOL[pathway]
        if (specific_rate_key not in keys or microbe_pool not in keys or
                pathway not in v_max_values):
            continue

        time, specific_rate = run_log[specific_rate_key]
        biomass_time, biomass = run_log[microbe_pool]
        time = np.asarray(time, dtype=float)
        specific_rate = np.asarray(specific_rate, dtype=float)
        biomass_time = np.asarray(biomass_time, dtype=float)
        biomass = np.asarray(biomass, dtype=float)

        valid_biomass = np.isfinite(biomass_time) & np.isfinite(biomass)
        if not np.any(valid_biomass):
            continue
        biomass_time = biomass_time[valid_biomass]
        biomass = biomass[valid_biomass]
        order = np.argsort(biomass_time)
        biomass_time = biomass_time[order]
        biomass = biomass[order]
        biomass_time, unique_indices = np.unique(biomass_time,
                                                  return_index=True)
        biomass = biomass[unique_indices]

        biomass_on_rate_grid = np.interp(time, biomass_time, biomass)
        actual_rate = biomass_on_rate_grid * specific_rate
        maximum_rate = biomass_on_rate_grid * v_max_values[pathway]
        yield pathway, time, actual_rate, maximum_rate


def quantity_ylabel(name):
    """Return a compact, quantity-appropriate y-axis label."""
    if name.endswith('_deltaG_r'):
        return rf'$\Delta G_\mathrm{{r}}\ [{GIBBS_ENERGY_UNIT}]$'
    if name.endswith(('_thermodynamic_factor', '_MM', '_inhib')):
        return r'$\mathrm{factor}$'
    if name.endswith('_v') or name.endswith('_v_max'):
        return r'$\mathrm{rate}$'
    return r'$\mathrm{substance\ [\mu mol/g\ dry\ weight]}$'


def figure_suptitle(checkpoint, fit_replicas=(), curve_replica=None):
    """Return shared metadata, including the data behind a model curve."""
    fit_replica_numbers = ', '.join(
        str(replica.replica_number) for replica in fit_replicas)
    curve_replica_number = getattr(curve_replica, 'replica_number', None)
    first_line = ' | '.join((
        checkpoint.replica[0],
        checkpoint.fit_mode,
        f'model {model_variant(checkpoint)}',
    ))
    return first_line


def plot_quantity(name, time, values, checkpoint, v_max=None, suptitle=None):
    """Plot one modeled pool or pathway quantity in the fit-plot style."""
    figure, axis = plt.subplots(figsize=(4, 4))
    if name.endswith('_deltaG_r'):
        values = np.asarray(values, dtype=float) * GIBBS_ENERGY_SCALE
    axis.plot(time, values, '-', color='k', label=name)
    if v_max is not None:
        axis.axhline(v_max, color='k', linestyle='--', label='v_max')
        if name.endswith('_v'):
            axis.set_ylim(0, v_max * 1.05)
        axis.legend(fancybox=False, edgecolor='k', loc='best')
    if name.endswith('_deltaG_r'):
        axis.axhline(GIBBS_MINIMUM * GIBBS_ENERGY_SCALE,
                     color='k', linestyle='--', linewidth=1,
                     label=r'$\mathrm{Gibbs\ minimum}$')
        axis.legend(fancybox=False, edgecolor='k', loc='best')
    axis.set_title(name)
    axis.set_xlabel(r'$\mathrm{time\ [d]}$')
    axis.set_ylabel(quantity_ylabel(name))
    if name.endswith(('_thermodynamic_factor', '_MM', '_inhib')):
        axis.set_ylim(0, 1)
    plot.format_ax(axis)
    figure.suptitle(figure_suptitle(checkpoint) if suptitle is None else suptitle)
    figure.tight_layout(rect=(0, 0, 1, 0.85))
    return figure, axis


def plot_reaction_rate(pathway, time, actual_rate, maximum_rate, checkpoint,
                       *, suptitle=None):
    """Plot a pathway's biomass-scaled rate and its biomass-scaled maximum."""
    figure, axis = plt.subplots(figsize=(4, 4))
    axis.plot(time, actual_rate, '-', color='k', label=r'$M(t)\,v(t)$')
    axis.plot(time, maximum_rate, '--', color='k',
              label=r'$M(t)\,v_{\max}$')
    axis.set_title(f'{pathway} reaction rate')
    axis.set_xlabel(r'$\mathrm{time\ [d]}$')
    axis.set_ylabel(r'$\mathrm{reaction\ rate\ [\mu mol/g\ dry\ weight/d]}$')
    axis.legend(fancybox=False, edgecolor='k', loc='best')
    plot.format_ax(axis)
    figure.suptitle(figure_suptitle(checkpoint) if suptitle is None else suptitle)
    figure.tight_layout(rect=(0, 0, 1, 0.85))
    return figure, axis


def _factor_product(factor_series):
    """Return the product of factor series on the first series' time grid."""
    time = np.asarray(factor_series[0][0], dtype=float)
    product = np.ones_like(time)
    for factor_time, factor_values in factor_series:
        factor_time = np.asarray(factor_time, dtype=float)
        factor_values = np.asarray(factor_values, dtype=float)
        valid = np.isfinite(factor_time) & np.isfinite(factor_values)
        if not np.any(valid):
            return time, np.full_like(time, np.nan)
        factor_time = factor_time[valid]
        factor_values = factor_values[valid]
        order = np.argsort(factor_time)
        factor_time = factor_time[order]
        factor_values = factor_values[order]
        factor_time, unique_indices = np.unique(factor_time,
                                                 return_index=True)
        factor_values = factor_values[unique_indices]
        product *= np.interp(time, factor_time, factor_values)
    return time, product


def plot_pathway_limitation_factors(pathway, mm_series, thermodynamic_series,
                                    inhibition_series, checkpoint, *, suptitle=None):
    """Plot all active limitation factors and their product for one pathway."""
    figure, axis = plt.subplots(figsize=(4, 4))
    if mm_series is None:
        axis.text(.02, .96, 'Michaelis-Menten factor not logged', color='0.45',
                  transform=axis.transAxes, ha='left', va='top', fontsize=8)
    else:
        time, values = mm_series
        axis.plot(time, values, color='tab:blue', linewidth=1.5,
                  label=r'$\mathrm{Michaelis\! -\! Menten}$')
    if thermodynamic_series is None:
        axis.text(.02, .84, 'Thermodynamic factor not logged', color='0.45',
                  transform=axis.transAxes, ha='left', va='top', fontsize=8)
    else:
        time, values = thermodynamic_series
        axis.plot(time, values, color='tab:orange', linestyle='--', linewidth=1.5,
                  label=r'$\mathrm{thermodynamic}$')
    if inhibition_series is not None:
        time, values = inhibition_series
        axis.plot(time, values, color='tab:purple', linestyle='-.', linewidth=1.5,
                  label=r'$\mathrm{inhibition}$')
    if mm_series is None:
        axis.text(.02, .72, 'Total factor not available', color='0.45',
                  transform=axis.transAxes, ha='left', va='top', fontsize=8)
    else:
        # Missing inhibition or thermodynamic logs mean that mechanism was not
        # active in the model, i.e. its factor is one. The product is therefore
        # v / v_max, equivalently M(t)v(t) / [M(t)v_max].
        active_factors = [mm_series]
        if thermodynamic_series is not None:
            active_factors.append(thermodynamic_series)
        if inhibition_series is not None:
            active_factors.append(inhibition_series)
        time, total_factor = _factor_product(active_factors)
        axis.plot(time, total_factor, color='k', linewidth=2,
                  label=r'$\mathrm{all\ factors}\ [M(t)v(t)]/[M(t)v_{\max}]$')
    axis.axhline(1, color='0.65', linestyle=':', linewidth=1,
                 label=r'$\mathrm{no\ limitation}$')
    axis.set_title(f'{pathway} limitation factors')
    axis.set_xlabel(r'$\mathrm{time\ [d]}$')
    axis.set_ylabel(r'$\mathrm{factor}$')
    axis.set_ylim(0, 1.05)
    plot.format_ax(axis)
    axis.legend(fancybox=False, edgecolor='k', loc='best')
    figure.suptitle(figure_suptitle(checkpoint) if suptitle is None else suptitle)
    figure.tight_layout(rect=(0, 0, 1, 0.85))
    return figure, axis


def pathway_limitation_factor_figures(checkpoint, run_log, *, suptitle=None):
    """Return one all-factor limitation plot for each active pathway."""
    keys = _run_log_keys(run_log)
    figures = []
    for pathway in PATHWAY_PLOT_ORDER:
        if f'{pathway}_v' not in keys:
            continue
        mm_key = f'{pathway}_MM'
        thermodynamic_key = f'{pathway}_thermodynamic_factor'
        inhibition_key = f'{pathway}_inhib'
        figures.append((f'{pathway}_limitation_factors',
                        plot_pathway_limitation_factors(
                            pathway,
                            run_log[mm_key] if mm_key in keys else None,
                            (run_log[thermodynamic_key]
                             if thermodynamic_key in keys else None),
                            run_log[inhibition_key] if inhibition_key in keys else None,
                            checkpoint, suptitle=suptitle)))
    return figures


def _pathway_from_model(pathway_model, pathway_name):
    """Return a named pathway object from a built model."""
    pathway = next(
        (candidate for candidate in pathway_model.contributing_pathways
         if candidate.__class__.__name__ == pathway_name),
        None,
    )
    if pathway is None:
        raise ValueError(f'Model does not contain the {pathway_name} pathway.')
    return pathway


def pathway_substrate_consumption_rate(pathway_model, run_log, pathway_name,
                                       substrate):
    """Return a pathway's positive net consumption rate of one substrate.

    ``pathway_vector`` includes both reaction stoichiometry and the carbon
    assimilated into biomass, so acetate uptake is reported as the total model
    consumption attributed to that pathway.
    """
    pathway = _pathway_from_model(pathway_model, pathway_name)
    coefficient = -float(pathway.pathway_vector[system.index(substrate)])
    if coefficient <= 0:
        raise ValueError(f'{pathway_name} does not consume {substrate}.')

    microbe_pool = PATHWAY_MICROBE_POOL[pathway_name]
    try:
        time, specific_rate = run_log[f'{pathway_name}_v']
        biomass_time, biomass = run_log[microbe_pool]
    except KeyError as error:
        raise ValueError(
            f'Model did not log the values required for {pathway_name} '
            f'{substrate} consumption.') from error

    time = np.asarray(time, dtype=float)
    specific_rate = np.asarray(specific_rate, dtype=float)
    biomass_time = np.asarray(biomass_time, dtype=float)
    biomass = np.asarray(biomass, dtype=float)
    order = np.argsort(biomass_time)
    biomass_on_rate_grid = np.interp(time, biomass_time[order], biomass[order])
    return time, coefficient * biomass_on_rate_grid * specific_rate


def cumulative_pathway_production(time, rate):
    """Integrate a pathway production rate after consolidating solver times."""
    time = np.asarray(time, dtype=float)
    rate = np.asarray(rate, dtype=float)
    valid = np.isfinite(time) & np.isfinite(rate)
    time, rate = time[valid], rate[valid]
    order = np.argsort(time)
    time, rate = time[order], rate[order]
    time, inverse = np.unique(time, return_inverse=True)
    rate_sum = np.bincount(inverse, weights=rate)
    rate_count = np.bincount(inverse)
    rate = np.maximum(rate_sum / rate_count, 0.)
    if len(time) < 2:
        return time, np.zeros_like(time)
    increments = np.diff(time) * (rate[:-1] + rate[1:]) / 2
    return time, np.concatenate(([0.], np.cumsum(increments)))


def ch4_pathway_production_figure(checkpoint, run_log, *, suptitle=None):
    """Stack cumulative methane production by pathway against total methane.

    Each area is the cumulative methane produced by an active CH4-producing
    pathway. The legend reports each pathway's final cumulative production,
    while the dashed line is the modeled total CH4 pool increase.
    """
    keys = _run_log_keys(run_log)
    curves = []
    for pathway in PATHWAY_PLOT_ORDER:
        key = f'{pathway}_CH4 from {PATHWAY_MICROBE_POOL[pathway]}'
        if key in keys:
            time, cumulative = cumulative_pathway_production(*run_log[key])
            if len(time):
                curves.append((pathway, time, cumulative))
    if not curves or 'CH4' not in keys:
        return None

    time = np.unique(np.concatenate([curve_time for _, curve_time, _ in curves]))
    contributions = [
        np.interp(time, curve_time, cumulative, left=0., right=cumulative[-1])
        for _, curve_time, cumulative in curves
    ]
    figure, axes = plot.create_figure(plot.FigureSpec(
        axis_names=('diagnostic',), figsize=(5, 4)))
    axis = axes['diagnostic']
    labels = [
        rf'$\mathrm{{{pathway} pathway}}$' + f': {contribution[-1]:.3g}'
        for (pathway, _, _), contribution in zip(curves, contributions)
    ]
    axis.stackplot(
        time, *contributions,
        labels=labels,
        colors=[pathway_color(pathway) for pathway, _, _ in curves],
        alpha=.75,
    )
    total_time, total_ch4 = run_log['CH4']
    total_ch4 = np.asarray(total_ch4, dtype=float)
    total_production = total_ch4 - total_ch4[0]
    axis.plot(total_time, total_production, 'k--', linewidth=1.5,
              label=(r'$\mathrm{total CH_4 produced}$'
                     f': {total_production[-1]:.3g}'))
    axis.set_title(r'$\mathrm{cumulative CH_4 production}$')
    axis.set_xlabel(r'$\mathrm{time [d]}$')
    axis.set_ylabel(r'$\mathrm{cumulative CH_4 [\mu mol/g dry weight]}$')
    plot.format_ax(axis)
    axis.legend(
        title=r'$\mathrm{final cumulative CH_4 [\mu mol/g]}$',
        fancybox=False, edgecolor='k', loc='best')
    figure.suptitle(figure_suptitle(checkpoint) if suptitle is None else suptitle)
    figure.subplots_adjust(left=.19, right=.95, bottom=.16, top=.80)
    return figure, axis


def plot_pathway_pair(checkpoint, curves, title, ylabel, *, suptitle=None,
                      gibbs_minimum=False, factor_limits=False,
                      secondary_curve=None, secondary_ylabel=None):
    """Plot two pathway diagnostics with shared styling and TeX labels."""
    figure, axes = plot.create_figure(plot.FigureSpec(
        axis_names=('diagnostic',), figsize=(4, 4)))
    axis = axes['diagnostic']
    for pathway, time, values in curves:
        if gibbs_minimum:
            values = np.asarray(values, dtype=float) * GIBBS_ENERGY_SCALE
        axis.plot(time, values, color=pathway_color(pathway), linewidth=1.5,
                  label=rf'$\mathrm{{{pathway}\ pathway}}$')
    if gibbs_minimum:
        axis.axhline(GIBBS_MINIMUM * GIBBS_ENERGY_SCALE,
                     color='k', linestyle='--', linewidth=1,
                     label=r'$\mathrm{Gibbs\ minimum}$')
    axis.set_title(title)
    axis.set_xlabel(r'$\mathrm{time\ [d]}$')
    axis.set_ylabel(ylabel)
    if factor_limits:
        axis.set_ylim(0, 1)
    plot.format_ax(axis)

    handles, labels = axis.get_legend_handles_labels()
    if secondary_curve is not None:
        time, values, label = secondary_curve
        secondary_axis = axis.twinx()
        secondary_axis.plot(time, values, color='k', linestyle='--',
                            linewidth=1.5, label=label)
        secondary_axis.set_ylabel(secondary_ylabel)
        secondary_axis.spines['top'].set_visible(False)
        extra_handles, extra_labels = secondary_axis.get_legend_handles_labels()
        handles.extend(extra_handles)
        labels.extend(extra_labels)
    axis.legend(handles, labels, fancybox=False, edgecolor='k', loc='best')
    figure.suptitle(figure_suptitle(checkpoint) if suptitle is None else suptitle)
    figure.subplots_adjust(
        left=.24, right=.76 if secondary_curve is not None else .95,
        bottom=.16, top=.80)
    return figure, axis


def fe3_aceto_diagnostic_figures(checkpoint, model_run, *, suptitle=None):
    """Return shared-axis Fe3/Aceto pathway diagnostics when both are active."""
    pathway_model, _, run_log = model_run
    keys = _run_log_keys(run_log)
    pathways = ('Fe3', 'Aceto')
    figures = []

    microbe_keys = tuple(PATHWAY_MICROBE_POOL[pathway] for pathway in pathways)
    if all(name in keys for name in microbe_keys):
        figures.append(('Fe3_Aceto_microbes', plot_pathway_pair(
            checkpoint,
            [(pathway, *run_log[PATHWAY_MICROBE_POOL[pathway]])
             for pathway in pathways],
            r'$\mathrm{living\ microbes}$',
            r'$\mathrm{microbial\ biomass\ [\mu mol/g\ dry\ weight]}$',
            suptitle=suptitle)))

    delta_g_keys = tuple(f'{pathway}_deltaG_r' for pathway in pathways)
    if all(name in keys for name in delta_g_keys):
        figures.append(('Fe3_Aceto_deltaG_r', plot_pathway_pair(
            checkpoint,
            [(pathway, *run_log[f'{pathway}_deltaG_r']) for pathway in pathways],
            r'$\Delta G_\mathrm{r}$',
            rf'$\Delta G_\mathrm{{r}}\ [{GIBBS_ENERGY_UNIT}]$',
            suptitle=suptitle,
            gibbs_minimum=True)))

    try:
        acetate_curves = [
            (pathway, *pathway_substrate_consumption_rate(
                pathway_model, run_log, pathway, 'Acetate'))
            for pathway in pathways
        ]
    except ValueError:
        acetate_curves = []
    if acetate_curves:
        acetate_pool = None
        if 'Acetate' in keys:
            acetate_pool = (*run_log['Acetate'],
                            r'$\mathrm{total\ Acetate\ pool}$')
        figures.append(('Fe3_Aceto_Acetate_consumption', plot_pathway_pair(
            checkpoint, acetate_curves,
            r'$\mathrm{Acetate\ consumption}$',
            r'$\mathrm{Acetate\ consumption\ [\mu mol/g\ dry\ weight/d]}$',
            suptitle=suptitle,
            secondary_curve=acetate_pool,
            secondary_ylabel=r'$\mathrm{total\ Acetate\ [\mu mol/g\ dry\ weight]}$')))

    mm_keys = tuple(f'{pathway}_MM' for pathway in pathways)
    if all(name in keys for name in mm_keys):
        figures.append(('Fe3_Aceto_MM', plot_pathway_pair(
            checkpoint,
            [(pathway, *run_log[f'{pathway}_MM']) for pathway in pathways],
            r'$\mathrm{Michaelis\! -\! Menten\ factor}$',
            r'$\mathrm{factor}$',
            suptitle=suptitle,
            factor_limits=True)))
    return figures

def model_variant(cp):
    variant = model_variant_from_id(cp.model_id)
    if variant is not None:
        return variant
    raise Exception('Model variant not identified: ' + cp.model_id)


def sample_y_limits(sample, relative_padding=0.05):
    """Return fixed CO2 and CH4 limits from all measured sample replicas."""
    limits = {}
    for pool in ('CO2', 'CH4'):
        values = []
        for replica in sample.replicas:
            _, replica_values = getattr(replica, pool)()
            replica_values = np.asarray(replica_values, dtype=float)
            usable = np.isfinite(replica_values)
            if pool == 'CH4':
                usable &= replica_values > 0
            values.extend(replica_values[usable])

        if not values:
            continue
        values = np.asarray(values)
        if pool == 'CO2':
            low, high = np.min(values), np.max(values)
            padding = max((high - low) * relative_padding,
                          max(abs(low), abs(high), 1.) * relative_padding)
            limits[pool] = (low - padding, high + padding)
        else:
            log_values = np.log10(values)
            low, high = np.min(log_values), np.max(log_values)
            padding = max((high - low) * relative_padding, relative_padding)
            limits[pool] = (10**(low - padding), 10**(high + padding))
    return limits


def fit_replicas(checkpoint, dataset):
    """Return the sample, validation replica, and fitted replicas for a plot."""
    sample_number, validation_number = checkpoint.replica
    sample = dataset[sample_number]
    validation_replica = sample[validation_number]
    replicas = sample.get_split(validation_number, checkpoint.fit_mode)['fit']
    return sample, validation_replica, replicas


def create_fit_figure():
    """Create the standard two-pool figure used for fitted-model plots."""
    return plot.create_figure(plot.FigureSpec(
        ncols=2,
        axis_names=('CO2', 'CH4'),
        figsize=(10, 5),
        sharex=True,
    ))


def plot_fit_measurements(checkpoint, dataset, axes):
    """Draw the fit/validation measurements and return their legend handles."""
    _, validation_replica, fitted_replicas = fit_replicas(checkpoint, dataset)
    markersize = 3
    markers = cycle(('^', 's', 'D', 'v', 'P', 'X'))
    handles = []
    for fitted_replica in fitted_replicas:
        marker = next(markers)
        plot.plot_replica_pools(
            fitted_replica, axes, marker=marker, ms=markersize, mfc='none')
        handles.append(Line2D(
            [], [], color='k', marker=marker, linestyle='None',
            markerfacecolor='none', markersize=markersize,
            label=f'fit replica {fitted_replica.replica_number}',
        ))

    if checkpoint.fit_mode == 'split':
        plot.plot_replica_pools(
            validation_replica, axes, marker='o', ms=markersize)
        handles.append(Line2D(
            [], [], color='k', marker='o', linestyle='None', markersize=markersize,
            label=f'validation replica {validation_replica.replica_number}',
        ))
    return handles


def fit_score_label(checkpoint, objective, validation_objective, run_log,
                    model_label=None):
    """Return the TeX R-squared label used for a model line in a fit plot."""
    r2_values = [(objective.R2(run_log), 'fit')]
    if checkpoint.fit_mode != 'single':
        r2_values.append((validation_objective.R2(run_log), 'val'))
    r2_text = ', '.join(
        rf'$R^2_{{\mathrm{{{name}}}}} = {value:.2f}$'
        for value, name in r2_values
    )
    if model_label is None:
        model_label = rf'$\mathrm{{model\ {{{model_variant(checkpoint)}}}}}$'
    return ' '.join((model_label, r2_text))


def format_fit_axes(checkpoint, sample, axes):
    """Apply the standard fit-plot scales, limits, and axes style."""
    y_limits = sample_y_limits(sample)
    for pool, axis in axes.items():
        plot.format_ax(
            axis,
            log_scale='log' in checkpoint.run_config['objective']['transform'][pool],
        )
        if pool in y_limits:
            axis.set_ylim(y_limits[pool])


def finish_fit_figure(figure, suptitle):
    """Apply the fixed report layout shared by fit comparison figures."""
    figure.suptitle(suptitle)
    figure.subplots_adjust(
        left=0.12, right=0.97, bottom=0.14, top=0.78, wspace=0.30)



def plot_fit(checkpoint, dataset, model_run=None):
    """Plot fitted replicas, validation data, and a model run for a checkpoint."""
    if model_run is None:
        model_run = run(checkpoint.run_config, checkpoint.parameters)
    pathway_model, objective, run_log = model_run
    sample, validation_replica, fitted_replicas = fit_replicas(checkpoint, dataset)
    validation_objective = optimizer.build_objective_function(
        pathway_model,
        validation_replica,
        checkpoint.run_config['objective'],
        t_start=0,
        t_end=None,
    )
    figure, axes = create_fit_figure()
    handles = plot_fit_measurements(checkpoint, dataset, axes)
    plot.plot_run_pools(run_log, axes)

    handles.append(Line2D(
        [], [], color='k', linestyle='-',
        label=fit_score_label(checkpoint, objective, validation_objective, run_log),
    ))

    axes['CH4'].legend(
        handles=handles, fancybox=False, edgecolor='k', loc='lower right')
    format_fit_axes(checkpoint, sample, axes)
    finish_fit_figure(figure, figure_suptitle(
        checkpoint, fitted_replicas, getattr(run_log, 'replica', None)))

    return figure, axes


def plot_all(results_directory=RESULTS_DIRECTORY, **filters):
    """Generate plots for all best checkpoints and return their paths."""
    dataset = data.get_data_before_carex()
    paths = []
    checkpoints = filter_checkpoints(collect_best(results_directory), **filters)
    for checkpoint in checkpoints:
        model_run = run(checkpoint.run_config, checkpoint.parameters)
        figure, _ = plot_fit(checkpoint, dataset, model_run=model_run)
        directory, filename = plot_target(checkpoint)

        os.makedirs(directory, exist_ok=True)
        path = os.path.join(directory, filename)
        #print('saving to  ', os.path.split(directory)[-1], filename)
        figure.savefig(path, dpi=300)
        plt.close(figure)
        paths.append(path)

        _, _, run_log = model_run
        sample = dataset[checkpoint.replica[0]]
        fit_replicas = sample.get_split(
            checkpoint.replica[1], checkpoint.fit_mode)['fit']
        suptitle = figure_suptitle(
            checkpoint, fit_replicas, getattr(run_log, 'replica', None))
        v_max_values = pathway_v_max_values(run_log)
        last_number = 0
        for number, (quantity, time, values) in enumerate(
                pathway_plot_quantities(run_log), start=1):
            pathway = quantity[:-2] if quantity.endswith('_v') else None
            figure, _ = plot_quantity(
                quantity, time, values, checkpoint,
                v_max=v_max_values.get(pathway), suptitle=suptitle)
            _, filename = plot_target(checkpoint, number, quantity)
            path = os.path.join(directory, filename)
            # print('saving to  ', os.path.split(directory)[-1], filename)
            figure.savefig(path, dpi=300)
            plt.close(figure)
            paths.append(path)
            last_number = number

        for number, (pathway, time, actual_rate, maximum_rate) in enumerate(
                pathway_reaction_rate_curves(run_log, v_max_values),
                start=last_number + 1):
            quantity = f'{pathway}_reaction_rate'
            figure, _ = plot_reaction_rate(
                pathway, time, actual_rate, maximum_rate, checkpoint,
                suptitle=suptitle)
            _, filename = plot_target(checkpoint, number, quantity)
            path = os.path.join(directory, filename)
            figure.savefig(path, dpi=300)
            plt.close(figure)
            paths.append(path)

        for number, (quantity, (figure, _)) in enumerate(
                pathway_limitation_factor_figures(
                    checkpoint, run_log, suptitle=suptitle), start=80):
            _, filename = plot_target(checkpoint, number, quantity)
            path = os.path.join(directory, filename)
            figure.savefig(path, dpi=300)
            plt.close(figure)
            paths.append(path)

        for number, (quantity, (figure, _)) in enumerate(
                fe3_aceto_diagnostic_figures(
                    checkpoint, model_run, suptitle=suptitle), start=90):
            _, filename = plot_target(checkpoint, number, quantity)
            path = os.path.join(directory, filename)
            figure.savefig(path, dpi=300)
            plt.close(figure)
            paths.append(path)

        ch4_production = ch4_pathway_production_figure(
            checkpoint, run_log, suptitle=suptitle)
        if ch4_production is not None:
            quantity = 'CH4_pathway_production'
            figure, _ = ch4_production
            _, filename = plot_target(checkpoint, 94, quantity)
            path = os.path.join(directory, filename)
            figure.savefig(path, dpi=300)
            plt.close(figure)
            paths.append(path)
    return paths


def parse_args():
    parser = argparse.ArgumentParser(
        description='Plot the best checkpoints, optionally limited to selected fits.')
    parser.add_argument('--results-directory', default=RESULTS_DIRECTORY,
                        help='Checkpoint root to scan.')
    parser.add_argument('--sample', dest='samples', nargs='+', metavar='SAMPLE',
                        help='Only plot these sample IDs.')
    parser.add_argument('--model-id', dest='model_ids', nargs='+', metavar='MODEL_ID',
                        help='Only plot model IDs containing these values.')
    parser.add_argument('--fit-mode', dest='fit_modes', choices=('single', 'split'),
                        nargs='+', metavar='FIT_MODE',
                        help='Only plot these fit modes.')
    parser.add_argument('--replica', dest='replicas', nargs='+', metavar='REPLICA',
                        help='Only plot these validation replica numbers.')
    parser.add_argument('--loss-id', dest='loss_ids', nargs='+', metavar='LOSS_ID',
                        help='Only plot these exact loss IDs.')
    parser.add_argument('--dry-run', action='store_true',
                        help='Print selected checkpoints without creating plots.')
    return parser.parse_args()


def main(args=None):
    args = parse_args() if args is None else args
    filters = {
        'samples': args.samples,
        'model_ids': args.model_ids,
        'fit_modes': args.fit_modes,
        'replicas': args.replicas,
        'loss_ids': args.loss_ids,
    }
    if args.dry_run:
        checkpoints = filter_checkpoints(
            collect_best(args.results_directory), **filters)
        for checkpoint in checkpoints:
            print(checkpoint.model_id, checkpoint.loss_id,
                  checkpoint.replica_name, checkpoint.fit_mode)
        print(f'{len(checkpoints)} checkpoint(s) selected')
        return checkpoints
    return plot_all(args.results_directory, **filters)


if __name__ == '__main__':
    main()
