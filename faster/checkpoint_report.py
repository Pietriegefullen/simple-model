"""Checkpoint catalogue and reusable standard report figures.

The command interface lives in :mod:`plot_checkpoint`; this module only
loads saved checkpoints and draws their standard report figures.
"""

import json
import os
from dataclasses import dataclass
from functools import cached_property
from itertools import cycle
from pathlib import Path

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
import fit_replicas as joint_fit
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
    'thermodynamic_factor', 'deltaG_r', 'MM', 'microbe_pool', 'v',
)

@dataclass
class Checkpoint:
    """A loaded checkpoint and the metadata used to group it."""

    parameters: parameters.ModelParameters
    loss: float
    run_config: dict
    id: str
    source: str = None
    metadata: dict = None

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

    @property
    def is_multi_replica(self):
        """Whether this checkpoint is an outer joint-replica fit."""
        return self.fit_mode == 'multi'

    @property
    def multi_replica_results(self):
        """Return the saved per-replica results for a joint fit.

        These are written by :class:`fit_replicas.SharedReplicaObjective` at
        the same outer evaluation as the checkpoint.  In particular, their
        local parameter values belong to the outer shared parameter values in
        this checkpoint.
        """
        if not self.is_multi_replica:
            return ()
        records = (self.metadata or {}).get('replicas', ())
        if not isinstance(records, list):
            return ()
        return tuple(record for record in records if isinstance(record, dict))

    @cached_property
    def run_id(self):
        """The full configuration identity, including the multi-fit setup."""
        return hashing.build_run_id(self.run_config)
    
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
            except Exception as error:
                raise ValueError(f'Could not load checkpoint {file_path}') from error
            checkpoint_ids = [part for part in filename.split('_')
                              if part.startswith('cp-')]
            if not checkpoint_ids:
                raise ValueError(f'Checkpoint ID missing from {file_path}')
            if not run_config:
                raise ValueError(f'Run configuration missing from {file_path}')

            # ``load_parameter_file`` deliberately exposes the common
            # parameter/loss/config triple.  Joint fits additionally persist
            # a per-replica solution in the annotated JSON, which is needed
            # to reproduce every curve in the combined fit.
            contents = parameters.load_file(file_path)
            metadata = {'replicas': contents.get('replicas', [])} \
                if isinstance(contents, dict) else {}
            collected.append(Checkpoint(loaded_parameters, loss, run_config,
                                        checkpoint_ids[0], file_path, metadata))
    return collected


def choose_best(checkpoints):
    """Keep the lowest-loss checkpoint for each model, fit mode, and replica."""
    best_by_group = {}
    for checkpoint in checkpoints:
        # A multi fit's placeholder ``0/0`` replica is not an identity.  Its
        # complete run configuration is: it includes selected replicas,
        # shared parameters, parameter ranges, objective, and algorithms.
        group = (checkpoint.model_id, checkpoint.fit_mode,
                 checkpoint.run_id if checkpoint.is_multi_replica
                 else checkpoint.replica)
        best = best_by_group.get(group)
        if best is None or checkpoint.loss < best.loss:
            best_by_group[group] = checkpoint
    return list(best_by_group.values())


def collect_best(results_directory=RESULTS_DIRECTORY):
    """Return the best checkpoint for each model, fit mode, and replica."""
    return choose_best(collect_all(results_directory))


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
        # Keep this public filtering helper usable with the lightweight
        # checkpoint doubles used by callers and tests.
        multi_records = getattr(checkpoint, 'multi_replica_results', ())
        multi_samples = {str(record.get('sample')) for record in multi_records}
        multi_replicas = {
            str(value) for record in multi_records
            for value in (record.get('replica_number'), record.get('replica'))
            if value is not None
        }
        if filters['sample'] is not None and (
                (not multi_records and sample not in filters['sample'])
                or (multi_records and not multi_samples & filters['sample'])):
            continue
        if (filters['model_id'] is not None and not any(
                model_id in checkpoint.model_id
                for model_id in filters['model_id'])):
            continue
        if (filters['fit_mode'] is not None and
                checkpoint.fit_mode not in filters['fit_mode']):
            continue
        if filters['replica'] is not None and (
                (not multi_records and replica not in filters['replica'])
                or (multi_records and not multi_replicas & filters['replica'])):
            continue
        if filters['loss_id'] is not None and checkpoint.loss_id not in filters['loss_id']:
            continue
        selected.append(checkpoint)
    return selected


class ReplicaPlotCheckpoint:
    """A single-replica view of a joint checkpoint for shared plot helpers."""

    def __init__(self, checkpoint, replica):
        self.checkpoint = checkpoint
        self.replica = (str(replica.sample.sample_name),
                        str(replica.replica_number))

    @property
    def parameters(self):
        return self.checkpoint.parameters

    @property
    def run_config(self):
        return self.checkpoint.run_config

    @property
    def model_id(self):
        return self.checkpoint.model_id

    @property
    def loss_id(self):
        return self.checkpoint.loss_id

    @property
    def id(self):
        return self.checkpoint.id

    @property
    def fit_mode(self):
        # Each curve is fitted directly to its corresponding replica.
        return 'single'


def multi_plot_directory(checkpoint):
    """Return an unambiguous result directory for one joint checkpoint."""
    return os.path.join(
        PROJECT_DIRECTORY,
        '_'.join(('best', checkpoint.key, 'multi', checkpoint.run_id,
                  checkpoint.id)),
    )


def multi_plot_target(checkpoint, replica, number=0, quantity=None):
    """Return a target below the joint-fit directory for one replica plot."""
    replica_name = f'{replica.sample.sample_name}-{replica.replica_number}'
    directory = os.path.join(multi_plot_directory(checkpoint), replica_name)
    parts = [f'{number:02d}']
    if quantity is not None:
        parts.append(quantity)
    parts.extend((replica_name, checkpoint.model_id, 'multi',
                  checkpoint.run_id, checkpoint.id))
    return directory, '_'.join(parts)


def write_multi_fit_manifest(directory, checkpoint):
    """Write human-readable provenance beside a joint fit's plots."""
    records = checkpoint.multi_replica_results
    lines = [
        'JOINT MULTI-REPLICA FIT',
        '=' * 24,
        f'Outer checkpoint: {checkpoint.source}',
        f'Checkpoint ID: {checkpoint.id}',
        f'Run ID: {checkpoint.run_id}',
        f'Model ID: {checkpoint.model_id}',
        f'Loss ID: {checkpoint.loss_id}',
        f'Total loss: {checkpoint.loss:.12g}',
        '',
        f'Replicas in fit ({len(records)}):',
    ]
    for record in records:
        replica = record.get('replica',
                             f"{record.get('sample', '?')}{record.get('replica_number', '?')}")
        loss = record.get('loss')
        r2 = record.get('r2')
        lines.append(
            f'  - {replica}: loss={loss!r}, R2={r2!r}; '
            f"inner checkpoint: {record.get('checkpoint', 'not recorded')}")
    lines.extend((
        '',
        'FULL FIT CONFIGURATION',
        '-' * 22,
        json.dumps(checkpoint.run_config, indent=2, sort_keys=True),
        '',
        'PER-REPLICA SAVED SOLUTIONS',
        '-' * 27,
        json.dumps(records, indent=2, sort_keys=True),
        '',
    ))
    os.makedirs(directory, exist_ok=True)
    path = os.path.join(directory, 'fit_metadata.txt')
    with open(path, 'w') as handle:
        handle.write('\n'.join(lines))
    return path


def _run_log_keys(run_log):
    """Return run-log keys without depending on its private storage."""
    return set(run_log.keys())


def dissolved_h2_series(pathway_model, run_log):
    """Return the dissolved-H2 series evaluated on the total-H2 time grid.

    The system-state log stores gas pools as totals.  The pathway model uses
    ``compute_dissolved`` internally to partition them between the gas and
    water phases, so use that same calculation for the diagnostic instead of
    duplicating its Henry-law implementation here.
    """
    if pathway_model is None or not hasattr(pathway_model, 'compute_dissolved'):
        return None
    keys = _run_log_keys(run_log)
    if not {'H2', 'H2O'} <= keys:
        return None

    time, total_h2 = run_log['H2']
    water_time, water = run_log['H2O']
    time = np.asarray(time, dtype=float)
    total_h2 = np.asarray(total_h2, dtype=float)
    water_time = np.asarray(water_time, dtype=float)
    water = np.asarray(water, dtype=float)
    valid_water = np.isfinite(water_time) & np.isfinite(water)
    if not np.any(valid_water):
        return None
    water_time = water_time[valid_water]
    water = water[valid_water]
    order = np.argsort(water_time)
    water_time = water_time[order]
    water = water[order]
    water_time, unique_indices = np.unique(water_time, return_index=True)
    water = water[unique_indices]

    # Only H2 and H2O affect H2's dissolved fraction.  Supplying the full
    # state vector keeps this calculation identical to the model's own use.
    states = np.zeros((len(system.SYSTEM), len(time)))
    states[system.index('H2')] = total_h2
    states[system.index('H2O')] = np.interp(time, water_time, water)
    # ``compute_dissolved`` accepts one state vector at a time: its Henry-law
    # constants are one value per pool rather than per time point.
    dissolved_states = np.column_stack([
        pathway_model.compute_dissolved(state) for state in states.T
    ])
    dissolved_h2 = dissolved_states[system.index('H2')]
    return time, dissolved_h2


def pathway_plot_quantities(run_log, include_pathway=False, pathway_model=None):
    """Yield modeled quantities in the required plot-file order.

    Quantities added to a pathway log in the future are appended after the
    specified diagnostics for that pathway.  ``v_max`` is shown as a dashed
    reference line on the corresponding ``v`` plot, not as its own plot.
    With ``include_pathway=True``, each result also carries the pathway that
    owns it.  This distinguishes the shared ``M_Ferm`` pool when it represents
    a capacity for Hydrolysis versus Fermentation.
    """
    def emitted_quantity(name, time, values, pathway=None):
        result = (name, time, values)
        return result + (pathway,) if include_pathway else result

    keys = _run_log_keys(run_log)

    for pool in POOL_PLOT_ORDER:
        if pool in keys:
            yield emitted_quantity(pool, *run_log[pool])
            if pool == 'H2':
                dissolved_h2 = dissolved_h2_series(pathway_model, run_log)
                if dissolved_h2 is not None:
                    yield emitted_quantity('H2_dissolved', *dissolved_h2)

    for pathway in PATHWAY_PLOT_ORDER:
        prefix = f'{pathway}_'
        pathway_keys = {key for key in keys if key.startswith(prefix)}
        if not pathway_keys:
            continue

        for quantity in PATHWAY_QUANTITY_ORDER:
            if quantity == 'microbe_pool':
                microbe_pool = PATHWAY_MICROBE_POOL[pathway]
                if microbe_pool in keys:
                    yield emitted_quantity(
                        microbe_pool, *run_log[microbe_pool], pathway)
            else:
                name = prefix + quantity
                if name in pathway_keys:
                    yield emitted_quantity(name, *run_log[name], pathway)

        handled = {
            prefix + quantity
            for quantity in PATHWAY_QUANTITY_ORDER
            if quantity != 'microbe_pool'
        }
        # Add every remaining pathway diagnostic (for example ``v``,
        # inhibition, and methane production) in a stable order.
        for name in sorted(pathway_keys - handled):
            yield emitted_quantity(name, *run_log[name], pathway)


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


def pathway_uses_capacity_parameter(run_log, pathway):
    """Whether a pathway's ``M`` value is a capacity rather than biomass.

    When ``v_max`` has been removed from the fit and fixed at one, multiplying
    it by ``M`` leaves ``M`` as the pathway capacity for that pathway's
    configured reactant/product stoichiometry.  A fitted value of one does not
    have that interpretation, so both the value and constant status matter.
    """
    model_parameters = getattr(run_log, '_parameters', None)
    parameter_name = PATHWAY_V_MAX_PARAMETER[pathway]
    if model_parameters is None or parameter_name not in model_parameters:
        return False

    v_max = model_parameters[parameter_name]
    return v_max.is_constant() and np.isclose(float(v_max), 1.)


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


def quantity_ylabel(name, capacity_pathway=None):
    """Return a compact, quantity-appropriate y-axis label."""
    if capacity_pathway is not None:
        return (r'$\mathrm{pathway\ capacity\ '
                r'[\mu mol/g\ dry\ weight/d]}$')
    if name.endswith('_deltaG_r'):
        return rf'$\Delta G_\mathrm{{r}}\ [{GIBBS_ENERGY_UNIT}]$'
    if name.endswith(('_thermodynamic_factor', '_MM', '_inhib')):
        return r'$\mathrm{factor}$'
    if name.endswith('_v') or name.endswith('_v_max'):
        return r'$\mathrm{rate}$'
    return r'$\mathrm{substance\ [\mu mol/g\ dry\ weight]}$'


def quantity_title(name, capacity_pathway=None):
    """Return the plot title for one modeled quantity."""
    if capacity_pathway is not None:
        return f'{capacity_pathway} pathway capacity ({name})'
    if name == 'H2_dissolved':
        return 'Dissolved H2'
    return name


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


def plot_quantity(name, time, values, checkpoint, v_max=None, suptitle=None,
                  capacity_pathway=None, 
                  log_scale = False):
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
    if capacity_pathway is not None:
        axis.set_title(quantity_title(name, capacity_pathway))
        axis.text(
            .02, .96,
            '$v_{\\max} = 1$ (constant)\n'
            'M represents capacity for this\n'
            'pathway stoichiometry',
            transform=axis.transAxes, ha='left', va='top', fontsize=8,
            bbox={'facecolor': 'white', 'edgecolor': '0.7', 'alpha': .9},
        )
    else:
        axis.set_title(quantity_title(name))
    axis.set_xlabel(r'$\mathrm{time\ [d]}$')
    axis.set_ylabel(quantity_ylabel(name, capacity_pathway))
    if name.endswith(('_thermodynamic_factor', '_MM', '_inhib')):
        axis.set_ylim(0, 1)
    plot.format_ax(axis, log_scale)
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


def pathway_substance_production_rate(pathway_model, run_log, pathway_name,
                                      substance):
    """Return a pathway's positive net production rate of one substance.

    When complete system-state logs are available, rates are recalculated on
    accepted solver states.  This matches the total-pool curve and avoids
    integrating diagnostic rates from rejected solver trial states.
    """
    pathway = _pathway_from_model(pathway_model, pathway_name)
    coefficient = float(pathway.pathway_vector[system.index(substance)])
    if coefficient <= 0:
        raise ValueError(f'{pathway_name} does not produce {substance}.')

    keys = _run_log_keys(run_log)
    if (hasattr(pathway_model, 'compute_dissolved') and
            all(pool in keys for pool in system.SYSTEM)):
        time = np.asarray(run_log[substance][0], dtype=float)
        states = np.empty((len(system.SYSTEM), len(time)))
        for index, pool in enumerate(system.SYSTEM):
            pool_time, values = run_log[pool]
            pool_time = np.asarray(pool_time, dtype=float)
            values = np.asarray(values, dtype=float)
            valid = np.isfinite(pool_time) & np.isfinite(values)
            if not np.any(valid):
                break
            order = np.argsort(pool_time[valid])
            pool_time = pool_time[valid][order]
            values = values[valid][order]
            pool_time, unique_indices = np.unique(pool_time,
                                                  return_index=True)
            states[index] = np.interp(time, pool_time, values[unique_indices])
        else:
            state_logger = getattr(pathway, 'state_logger', None)
            try:
                pathway.state_logger = None
                rates = np.asarray([
                    pathway(t, pathway_model.compute_dissolved(state))[
                        system.index(substance)]
                    for t, state in zip(time, states.T)
                ])
            finally:
                pathway.state_logger = state_logger
            return time, rates

    # Lightweight logs used by older outputs and tests omit complete system
    # states.  Retain the rate-log reconstruction as a compatibility fallback.
    microbe_pool = PATHWAY_MICROBE_POOL[pathway_name]
    try:
        time, specific_rate = run_log[f'{pathway_name}_v']
        biomass_time, biomass = run_log[microbe_pool]
    except KeyError as error:
        raise ValueError(
            f'Model did not log the values required for {pathway_name} '
            f'{substance} production.') from error

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


def ch4_pathway_production_figure(checkpoint, run_log, *, suptitle=None,
                                  log_scale=False):
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
    plot.format_ax(axis, log_scale=log_scale)
    axis.legend(
        title=r'$\mathrm{final cumulative CH_4 [\mu mol/g]}$',
        fancybox=False, edgecolor='k', loc='best')
    figure.suptitle(figure_suptitle(checkpoint) if suptitle is None else suptitle)
    figure.subplots_adjust(left=.19, right=.95, bottom=.16, top=.80)
    return figure, axis


def co2_pathway_production_figure(checkpoint, pathway_model, run_log, *,
                                  suptitle=None, log_scale=False):
    """Stack cumulative CO2 production by every CO2-producing pathway.

    The dashed curve is the net change in the total CO2 pool.  It can differ
    from the stack because other pathways may consume CO2.
    """
    keys = _run_log_keys(run_log)
    curves = []
    for pathway in PATHWAY_PLOT_ORDER:
        try:
            rate_time, rate = pathway_substance_production_rate(
                pathway_model, run_log, pathway, 'CO2')
        except ValueError:
            continue
        time, cumulative = cumulative_pathway_production(rate_time, rate)
        if len(time):
            curves.append((pathway, time, cumulative))
    if not curves or 'CO2' not in keys:
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
    total_time, total_co2 = run_log['CO2']
    total_co2 = np.asarray(total_co2, dtype=float)
    total_change = total_co2 - total_co2[0]
    axis.plot(total_time, total_change, 'k--', linewidth=1.5,
              label=(r'$\mathrm{net\ total\ CO_2\ change}$'
                     f': {total_change[-1]:.3g}'))
    axis.set_title(r'$\mathrm{cumulative\ CO_2\ production}$')
    axis.set_xlabel(r'$\mathrm{time [d]}$')
    axis.set_ylabel(r'$\mathrm{cumulative\ CO_2\ [\mu mol/g dry weight]}$')
    plot.format_ax(axis, log_scale=log_scale)
    axis.legend(
        title=r'$\mathrm{final cumulative\ CO_2\ [\mu mol/g]}$',
        fancybox=False, edgecolor='k', loc='best')
    figure.suptitle(figure_suptitle(checkpoint) if suptitle is None else suptitle)
    figure.subplots_adjust(left=.19, right=.95, bottom=.16, top=.80)
    return figure, axis


def plot_pathway_pair(checkpoint, curves, title, ylabel, *, suptitle=None,
                      gibbs_minimum=False, factor_limits=False,
                      secondary_curve=None, secondary_ylabel=None, log_scale = False):
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
    plot.format_ax(axis, log_scale=log_scale)

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


def fe3_aceto_diagnostic_figures(checkpoint, model_run, *, suptitle=None, log_scale = False):
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
            suptitle=suptitle,
            log_scale = log_scale)))

    delta_g_keys = tuple(f'{pathway}_deltaG_r' for pathway in pathways)
    if all(name in keys for name in delta_g_keys):
        figures.append(('Fe3_Aceto_deltaG_r', plot_pathway_pair(
            checkpoint,
            [(pathway, *run_log[f'{pathway}_deltaG_r']) for pathway in pathways],
            r'$\Delta G_\mathrm{r}$',
            rf'$\Delta G_\mathrm{{r}}\ [{GIBBS_ENERGY_UNIT}]$',
            suptitle=suptitle,
            gibbs_minimum=True,
            log_scale = log_scale)))

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
            secondary_ylabel=r'$\mathrm{total\ Acetate\ [\mu mol/g\ dry\ weight]}$',
            log_scale = log_scale)))

    mm_keys = tuple(f'{pathway}_MM' for pathway in pathways)
    if all(name in keys for name in mm_keys):
        figures.append(('Fe3_Aceto_MM', plot_pathway_pair(
            checkpoint,
            [(pathway, *run_log[f'{pathway}_MM']) for pathway in pathways],
            r'$\mathrm{Michaelis\! -\! Menten\ factor}$',
            r'$\mathrm{factor}$',
            suptitle=suptitle,
            factor_limits=True,
            log_scale = log_scale)))
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
        log_applied = plot.format_ax(
            axis,
            log_scale='log' in checkpoint.run_config['objective']['transform'][pool],
        )
        if pool in y_limits and not log_applied:
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


def _save_diagnostic_plots(checkpoint, model_run, directory, target, suptitle,
                           *, overwrite=True):
    """Save the non-fit plots shared by single, split, and multi fits.

    ``target(number, quantity)`` supplies the filename while this function
    retains the established ``plot_best`` diagnostic ordering and styling.
    """
    paths = []

    def save(figure, path):
        if overwrite or not Path(path).exists():
            figure.savefig(path, dpi=300)
            paths.append(path)
        plt.close(figure)

    pathway_model, _, run_log = model_run
    v_max_values = pathway_v_max_values(run_log)
    last_number = 0
    for number, (quantity, time, values, pathway_context) in enumerate(
            pathway_plot_quantities(
                run_log, include_pathway=True, pathway_model=pathway_model),
            start=1):
        pathway = quantity[:-2] if quantity.endswith('_v') else None
        capacity_pathway = (
            pathway_context if quantity.startswith('M_')
            and pathway_context is not None
            and pathway_uses_capacity_parameter(run_log, pathway_context)
            else None
        )
        for log_scale in (False, True):
            figure, axis = plot_quantity(
                quantity, time, values, checkpoint,
                v_max=v_max_values.get(pathway), suptitle=suptitle,
                capacity_pathway=capacity_pathway, 
                log_scale = log_scale)
            if log_scale and axis.get_yscale() != 'log':
                plt.close(figure)
                continue
            filename_quantity = (
                f'{capacity_pathway}_{quantity}'
                if capacity_pathway is not None else quantity)
            path = target(number, filename_quantity + ('_log' if log_scale else ''))
            save(figure, path)
            last_number = number

    for number, (pathway, time, actual_rate, maximum_rate) in enumerate(
            pathway_reaction_rate_curves(run_log, v_max_values),
            start=last_number + 1):
        quantity = f'{pathway}_reaction_rate'
        figure, _ = plot_reaction_rate(
            pathway, time, actual_rate, maximum_rate, checkpoint,
            suptitle=suptitle)
        path = target(number, quantity)
        save(figure, path)

    for number, (quantity, (figure, _)) in enumerate(
            pathway_limitation_factor_figures(
                checkpoint, run_log, suptitle=suptitle), start=80):
        path = target(number, quantity)
        save(figure, path)

    for log_scale in (False, True):
        for number, (quantity, (figure, _)) in enumerate(
                fe3_aceto_diagnostic_figures(
                    checkpoint, model_run, suptitle=suptitle,
                    log_scale=log_scale), start=90):
            if log_scale and figure.axes[0].get_yscale() != 'log':
                plt.close(figure)
                continue
            path = target(number, quantity + ('_log' if log_scale else ''))
            save(figure, path)

    for log_scale in (False, True):
        ch4_production = ch4_pathway_production_figure(
            checkpoint, run_log, suptitle=suptitle, log_scale=log_scale)
        if ch4_production is None:
            continue
        figure, axis = ch4_production
        if log_scale and axis.get_yscale() != 'log':
            plt.close(figure)
            continue
        quantity = 'CH4_pathway_production' + ('_log' if log_scale else '')
        path = target(94, quantity)
        save(figure, path)

    for log_scale in (False, True):
        co2_production = co2_pathway_production_figure(
            checkpoint, pathway_model, run_log, suptitle=suptitle,
            log_scale=log_scale)
        if co2_production is None:
            continue
        figure, axis = co2_production
        if log_scale and axis.get_yscale() != 'log':
            plt.close(figure)
            continue
        quantity = 'CO2_pathway_production' + ('_log' if log_scale else '')
        path = target(95, quantity)
        save(figure, path)
    return paths


def multi_replica_model_run(checkpoint, record, dataset):
    """Recreate one saved local solution from a joint checkpoint."""
    sample = dataset[str(record['sample'])]
    replica = sample[str(record['replica_number'])]
    chosen = checkpoint.run_config['chosen']
    pathway_model = joint_fit.build_parameter_template(
        chosen['pathways'], chosen['normalized_parameters'])
    # The outer checkpoint contains shared values.  The record stores all
    # remaining values from the same outer evaluation, including initial pools
    # and per-replica local parameters.
    pathway_model.parameters().set(checkpoint.parameters)
    pathway_model.parameters().set(record['parameters'])
    objective = optimizer.build_objective_function(
        pathway_model, [replica], checkpoint.run_config['objective'])
    objective._call(pathway_model.parameters(), transformed=False)
    return replica, (pathway_model, objective, pathway_model.system_state_log)


def plot_multi_checkpoint(checkpoint, dataset):
    """Plot every replica in one joint fit beneath a shared provenance folder."""
    if not checkpoint.multi_replica_results:
        raise ValueError(
            f'Joint checkpoint {checkpoint.source} has no saved replica solutions.')
    root = multi_plot_directory(checkpoint)
    paths = [write_multi_fit_manifest(root, checkpoint)]
    replica_names = ', '.join(
        str(record.get('replica', '?')) for record in checkpoint.multi_replica_results)
    for record in checkpoint.multi_replica_results:
        replica, model_run = multi_replica_model_run(checkpoint, record, dataset)
        plot_checkpoint = ReplicaPlotCheckpoint(checkpoint, replica)
        figure, _ = plot_fit(plot_checkpoint, dataset, model_run=model_run)
        suptitle = ' | '.join((
            f'joint multi fit ({replica_names})',
            f'replica {replica.sample.sample_name}-{replica.replica_number}',
            f'model {model_variant(plot_checkpoint)}',
        ))
        finish_fit_figure(figure, suptitle)
        directory, filename = multi_plot_target(checkpoint, replica)
        os.makedirs(directory, exist_ok=True)
        path = os.path.join(directory, filename)
        figure.savefig(path, dpi=300)
        plt.close(figure)
        paths.append(path)

        def target(number, quantity):
            _, diagnostic_filename = multi_plot_target(
                checkpoint, replica, number, quantity)
            return os.path.join(directory, diagnostic_filename)

        paths.extend(_save_diagnostic_plots(
            plot_checkpoint, model_run, directory, target, suptitle))
    return paths
