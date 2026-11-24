"""Checkpoint catalogue and reusable standard report figures.

The command interface lives in :mod:`plot_checkpoint`; this module only
loads saved checkpoints and draws their standard report figures.
"""

import json
import os
from dataclasses import dataclass
from fnmatch import fnmatchcase
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
from CONSTANTS import GIBBS_MINIMUM, MOLAR_MASS_BIOMASS
from fit_sample import run
import fit_replicas as joint_fit
from plot import pathway_color
from model_variants import model_variant_from_id
from USER_VARIABLES import PROJECT_DIRECTORY, RESULTS_DIRECTORY

GIBBS_ENERGY_SCALE = 1e-3  # Internal reaction energies are J/mol; plots use kJ/mol.
GIBBS_ENERGY_UNIT = r'\mathrm{kJ}/\mathrm{mol}'
GIBBS_Y_LIMITS = (-1000., 500.)

# ``00`` is reserved for the data/model fit (CO2 and CH4).  The remaining
# system pools are grouped with the pathway diagnostics below, rather than
# duplicating the gases or a microbe pool in a separate group.
MODELLED_GAS_POOLS = ('CO2', 'CH4')
POOL_PLOT_ORDER = (*MODELLED_GAS_POOLS,
                   'TOC', 'DOC', 'Acetate', 'Fe3', 'H2', 'H2O', 'Fe2')
PATHWAY_PLOT_ORDER = (
    'Hydrolysis', 'Fermentation', 'Fe3', 'Aceto', 'Hydro', 'Homo',
)

# The complete standard plot suite. Edit these flags to suppress an entire
# plot family without changing the report code. ``PATHWAY_PLOT_FLAGS`` is a
# second, deliberately short layer: setting a pathway to ``False`` excludes
# it from every pathway-specific family (including joint signature plots).
PLOT_FLAGS = {
    'fit': True,
    'quantities': True,
    'gas_pools': True,
    'pathway_capacity': True,
    'reaction_rate_signature': True,
    'limitation_factors': True,
    'pathway_delta_g': True,
    'fe3_aceto_diagnostics': True,
    'ch4_pathway_production': True,
    'co2_pathway_production': True,
}
PATHWAY_PLOT_FLAGS = dict.fromkeys(PATHWAY_PLOT_ORDER, True)

def plot_enabled(name, pathway=None):
    """Whether one standard plot family is enabled for an optional pathway."""
    return (PLOT_FLAGS.get(name, False) and
            (pathway is None or PATHWAY_PLOT_FLAGS.get(pathway, True)))


def enabled_pathways():
    """Return active plot pathways in the shared display order."""
    return tuple(pathway for pathway in PATHWAY_PLOT_ORDER
                 if PATHWAY_PLOT_FLAGS.get(pathway, True))


PATHWAY_DISPLAY_NAMES = {
    'Hydrolysis': 'Hydrolysis',
    'Fermentation': 'Fermentation',
    'Ferm': 'Fermentation',
    'Fe3': 'Iron\ reduction',
    'Aceto': 'Acetoclastic\ methanogenesis',
    'Hydro': 'Hydrogenotrophic\ methanogenesis',
    'Homo': 'Homoacetogenesis',
}
PATHWAY_MICROBE_POOL = {
    'Hydrolysis': 'M_Ferm',
    'Fermentation': 'M_Ferm',
    'Fe3': 'M_Fe3',
    'Aceto': 'M_Ac',
    'Hydro': 'M_Hydro',
    'Homo': 'M_Homo',
}
MICROBIAL_POOLS = frozenset(PATHWAY_MICROBE_POOL.values())
MICROBIAL_BIOMASS_UNIT = r'\mu\mathrm{mol}/\mathrm{L}'
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


def checkpoint_filename_metadata(file_path):
    """Extract catalogue fields embedded in a canonical checkpoint filename.

    Checkpoint filenames are generated from the same configuration that is
    stored in their JSON payload.  Reading these fields avoids normalising and
    hashing that payload merely to populate or filter the catalogue.  Older
    or manually renamed files simply return no value for the missing fields
    and fall back to the stored configuration.
    """
    parts = os.path.basename(file_path).split('_')
    metadata = {}
    for part in parts:
        if part.startswith('model-'):
            metadata['model_id'] = part
        elif part.startswith('run-'):
            metadata['run_id'] = part
        elif part.startswith('cp-'):
            metadata['checkpoint_id'] = part
    for position, part in enumerate(parts):
        if part in {'single', 'split', 'multi'} and position >= 2:
            metadata.update({
                'sample': parts[position - 2],
                'replica': parts[position - 1],
                'fit_mode': part,
            })
            break
    return metadata

@dataclass
class Checkpoint:
    """A loaded checkpoint and the metadata used to group it."""

    parameters: parameters.ModelParameters
    loss: float
    run_config: dict
    id: str
    source: str = None
    metadata: dict = None
    filename_metadata: dict = None

    @cached_property
    def model_id(self):
        if self.filename_metadata and self.filename_metadata.get('model_id'):
            return self.filename_metadata['model_id']
        return hashing.build_model_id(self.run_config['model'])

    @cached_property
    def loss_id(self):
        return hashing.build_loss_id(self.run_config['objective'])

    @cached_property
    def replica(self):
        if (self.filename_metadata and self.filename_metadata.get('sample')
                and self.filename_metadata.get('replica')):
            return (str(self.filename_metadata['sample']),
                    str(self.filename_metadata['replica']))
        chosen = self.run_config['chosen']
        return str(chosen['sample']), str(chosen['validation_replica'])

    @cached_property
    def replica_name(self):
        return '-'.join(self.replica)

    @cached_property
    def fit_mode(self):
        if self.filename_metadata and self.filename_metadata.get('fit_mode'):
            return self.filename_metadata['fit_mode']
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
        if self.filename_metadata and self.filename_metadata.get('run_id'):
            return self.filename_metadata['run_id']
        return hashing.build_run_id(self.run_config)
    
    def __str__(self):
        return self.key


@dataclass(frozen=True)
class CheckpointLoadFailure:
    """A checkpoint-shaped file that could not be deserialised."""

    source: str
    checkpoint_id: str
    error: str
    filename_metadata: dict


def collect_all_with_failures(results_directory=RESULTS_DIRECTORY, verbose=True):
    """Load valid checkpoints and retain errors for each unreadable file."""
    collected = []
    failures = []
    for root, _, files in os.walk(results_directory):
        for filename in files:
            legacy_restart = (
                Path(root).name == 'cma_restart_minima'
                and fnmatchcase(filename, 'restart-[0-9][0-9][0-9]-*.json'))
            # Persisted checkpoints are deliberately extensionless.  Report
            # artefacts can contain checkpoint-like IDs.  The exception is
            # the older JSON layout for retained CMA restart minima.
            if filename.startswith('.') or (Path(filename).suffix and not legacy_restart):
                continue

            if verbose:
                print(f'{len(collected):6d}', filename)

            # Results directories also contain report artefacts (for example
            # sensitivity.csv and parameter_influence.png).  A persisted
            # checkpoint has its generated ``cp-...`` identifier as a filename
            # component, except for legacy restart minima whose ID is derived
            # from their contents below.
            checkpoint_ids = [part for part in filename.split('_')
                              if part.startswith('cp-')]
            if legacy_restart:
                checkpoint_ids = [Path(filename).stem]
            if not checkpoint_ids:
                continue

            file_path = os.path.join(root, filename)
            # Read JSON exactly once.  The former implementation first read
            # this file through ``load_parameter_file`` and then opened it a
            # second time just to obtain joint-fit metadata.
            try:
                contents = parameters.load_file(file_path)
                loaded_parameters = parameters.ModelParameters(contents['parameters'])
                loss = float(contents['total_loss'])
                run_config = contents['run_config']
                if legacy_restart:
                    checkpoint_ids = [hashing.build_checkpoint_id(contents)]
            except Exception as error:
                failures.append(CheckpointLoadFailure(
                    file_path, checkpoint_ids[0], str(error),
                    checkpoint_filename_metadata(file_path)))
                continue
            if not run_config:
                failures.append(CheckpointLoadFailure(
                    file_path, checkpoint_ids[0], 'Run configuration missing',
                    checkpoint_filename_metadata(file_path)))
                continue

            # ``load_parameter_file`` deliberately exposes the common
            # parameter/loss/config triple.  Joint fits additionally persist
            # a per-replica solution in the annotated JSON, which is needed
            # to reproduce every curve in the combined fit.
            metadata = ({'replicas': contents.get('replicas', []),
                         'cma_restart': contents.get('cma_restart')}
                        if isinstance(contents, dict) else {})
            collected.append(Checkpoint(loaded_parameters, loss, run_config,
                                        checkpoint_ids[0], file_path, metadata,
                                        checkpoint_filename_metadata(file_path)))
    if verbose:
        print('done')
    return collected, failures


def collect_all(results_directory=RESULTS_DIRECTORY, verbose=True):
    """Load every valid annotated checkpoint in ``results_directory``.

    Consumers that need to present malformed files can use
    :func:`collect_all_with_failures` alongside the usable checkpoints.
    """
    collected, failures = collect_all_with_failures(results_directory, verbose)
    if failures:
        failure = failures[0]
        raise ValueError(f'Could not load checkpoint {failure.source}: {failure.error}')
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
    its checkpoints without requiring the full generated identifier. Sample
    filters also accept shell-style patterns, such as ``137*``.
    """
    filters = {
        'sample': None if samples is None else set(map(str, samples)),
        'model_id': None if model_ids is None else set(model_ids),
        'fit_mode': None if fit_modes is None else set(fit_modes),
        'replica': None if replicas is None else set(map(str, replicas)),
        'loss_id': None if loss_ids is None else set(loss_ids),
    }
    selected = []
    total = len(checkpoints)
    count = 0
    for checkpoint in checkpoints:
        count += 1
        print(f'\r{count:6d}/{total:d}, selected {len(selected):5d}', end = '')
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
        matches_sample = lambda value: any(
            fnmatchcase(value, pattern) for pattern in filters['sample'])
        if filters['sample'] is not None and (
                (not multi_records and not matches_sample(sample))
                or (multi_records and not any(matches_sample(value)
                                               for value in multi_samples))):
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


def litres_per_gdw(run_log):
    """Return the exact volume reference used by this run's replica."""
    replica = getattr(run_log, 'replica', None)
    if replica is None:
        return None
    return (float(replica.water_content) * 1e-3) / float(replica.dry_weight)


def litre_pool_values(run_log, pool, values):
    """Convert a gdw/mass display pool back to the model's µmol/L basis."""
    scale = litres_per_gdw(run_log)
    values = np.asarray(values, dtype=float)
    reference = getattr(run_log, '_reference', 'gdw')
    microbe_unit = getattr(run_log, '_microbe_unit', 'mg')
    if reference == 'gdw':
        if scale is None:
            return values
        values = values / scale
    elif reference != 'L':
        raise ValueError(f'Unsupported run-log reference unit: {reference!r}')
    if pool.startswith('M_') and microbe_unit == 'mg':
        values = values / MOLAR_MASS_BIOMASS
    elif pool.startswith('M_') and microbe_unit != 'micromol':
        raise ValueError(f'Unsupported run-log microbial unit: {microbe_unit!r}')
    return values


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
    total_h2 = litre_pool_values(run_log, 'H2', total_h2)
    water_time = np.asarray(water_time, dtype=float)
    water = litre_pool_values(run_log, 'H2O', water)
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

    :data:`PATHWAY_QUANTITY_ORDER` is an explicit allowlist, rather than an
    invitation to draw every pathway entry in ``run_log``. ``v_max`` is shown
    as a dashed reference line on the corresponding ``v`` plot, not as its
    own plot. With ``include_pathway=True``, each result also carries the
    pathway that owns it. This distinguishes the shared ``M_Ferm`` pool when
    it represents a capacity for Hydrolysis versus Fermentation.
    """
    def emitted_quantity(name, time, values, pathway=None):
        if name in system.SYSTEM:
            values = litre_pool_values(run_log, name, values)
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
        biomass = litre_pool_values(run_log, microbe_pool, biomass)

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
        reaction_rate_key = f'{pathway}_reaction_rate'
        actual_rate = (np.asarray(run_log[reaction_rate_key][1], dtype=float)
                       if reaction_rate_key in keys else
                       biomass_on_rate_grid * specific_rate)
        if reaction_rate_key in keys:
            reaction_time = np.asarray(run_log[reaction_rate_key][0], dtype=float)
            if not np.array_equal(reaction_time, time):
                actual_rate = np.interp(time, reaction_time, actual_rate)
        maximum_rate = biomass_on_rate_grid * v_max_values[pathway]
        yield pathway, time, actual_rate, maximum_rate


def quantity_ylabel(name, capacity_pathway=None):
    """Return a compact, quantity-appropriate y-axis label."""
    if capacity_pathway is not None:
        return rf'$\mathrm{{pathway\ capacity\ [{plot.REACTION_RATE_UNIT}]}}$'
    if name in MICROBIAL_POOLS:
        return rf'$\mathrm{{microbial\ biomass\ [{MICROBIAL_BIOMASS_UNIT}]}}$'
    if name.endswith('_deltaG_r'):
        return rf'$\Delta G_\mathrm{{r}}\ [{GIBBS_ENERGY_UNIT}]$'
    if name.endswith('_MM'):
        return r'$f_{\mathrm{MM}}\ [-]$'
    if name.endswith('_inhib'):
        return r'$f_{\mathrm{I}}\ [-]$'
    if name.endswith('_thermodynamic_factor'):
        return r'$f_{\mathrm{\Delta G}}\ [-]$'
    if name.endswith('_v') or name.endswith('_v_max'):
        return r'$v\ [\mathrm{d}^{-1}]$'
    return rf'$\mathrm{{substance\ [{plot.SUBSTANCE_UNIT}]}}$'


def pathway_display_name(pathway):
    """Return the descriptive name used in report-facing pathway labels."""
    return PATHWAY_DISPLAY_NAMES.get(pathway, pathway)


def _microbe_pathway(name):
    """Return the default pathway name for a microbial-pool diagnostic."""
    return {
        'M_Ferm': 'Fermentation', 'M_Fe3': 'Fe3', 'M_Ac': 'Aceto',
        'M_Hydro': 'Hydro', 'M_Homo': 'Homo',
    }.get(name)


def quantity_title(name, capacity_pathway=None, pathway=None):
    """Return the plot title for one modeled quantity."""
    if capacity_pathway is not None:
        display_name = pathway_display_name(capacity_pathway)
        return rf'$\mathrm{{{display_name}}}\ \mathrm{{capacity\ rate}}$'
    if name in MICROBIAL_POOLS:
        display_name = pathway_display_name(pathway or _microbe_pathway(name))
        return rf'$\mathrm{{biomass}}\ \mathrm{{{display_name}}}$'
    if name == 'H2_dissolved':
        return r'$\mathrm{dissolved\ H_2}$'
    if name.endswith('_MM'):
        pathway = name.removesuffix('_MM')
        display_name = pathway_display_name(pathway)
        return rf'$\mathrm{{{display_name}}}$'
    if name.endswith('_inhib'):
        pathway = name.removesuffix('_inhib')
        display_name = pathway_display_name(pathway)
        return rf'$\mathrm{{{display_name}}}\ \text{{inhibition}}$'
    if name.endswith('_thermodynamic_factor'):
        pathway = name.removesuffix('_thermodynamic_factor')
        display_name = pathway_display_name(pathway)
        return rf'$\mathrm{{{display_name}}}\ \text{{thermodynamic factor}}$'
    if name.endswith('_v'):
        pathway = name.removesuffix('_v')
        display_name = pathway_display_name(pathway)
        return rf'$\mathrm{{{display_name}}}\ \mathrm{{turnover\ rate}}$'
    if ' from ' in name:
        pathway, produced = name.split('_', 1)
        gas = produced.split(' from ', 1)[0]
        display_name = PATHWAY_DISPLAY_NAMES.get(pathway, pathway)
        return (rf'${plot.chemical_formula(gas)}\ \mathrm{{from}}\ '
                rf'\mathrm{{{display_name}}}$')
    return plot.pool_title(name)


def delta_g_legend_label(name):
    """Return a readable legend label for one pathway's Gibbs energy."""
    pathway = name.removesuffix('_deltaG_r')
    displayed = (plot.chemical_formula(pathway)
                 if pathway in {'Fe3', 'Fe2'}
                 else rf'\mathrm{{{PATHWAY_DISPLAY_NAMES.get(pathway, pathway)}}}')
    return rf'$\Delta G_\mathrm{{r}}\ ({displayed})$'


def gibbs_lower_limit(time, values, *, start_day=50., relative_padding=.05):
    """Set the lower bound from the stable, post-transient Gibbs trajectory."""
    time = np.asarray(time, dtype=float)
    values = np.asarray(values, dtype=float)
    usable = np.isfinite(time) & np.isfinite(values) & (time >= start_day)
    if not np.any(usable):
        usable = np.isfinite(values)
    if not np.any(usable):
        return None
    lower, upper = np.min(values[usable]), np.max(values[usable])
    padding = max((upper - lower) * relative_padding,
                  max(abs(lower), abs(upper), 1.) * relative_padding)
    return lower - padding


def constrain_gibbs_axis(axis):
    """Keep Gibbs diagnostic axes within the standard readable range."""
    lower, upper = axis.get_ylim()
    lower_limit, upper_limit = GIBBS_Y_LIMITS
    lower, upper = max(lower, lower_limit), min(upper, upper_limit)
    if lower >= upper:
        lower, upper = GIBBS_Y_LIMITS
    axis.set_ylim(lower, upper)


def figure_suptitle(checkpoint, fit_replicas=(), curve_replica=None):
    """Return shared metadata, including the data behind a model curve."""
    fit_replica_numbers = ', '.join(
        str(replica.replica_number) for replica in fit_replicas)
    curve_replica_number = getattr(curve_replica, 'replica_number', None)
    first_line = ' | '.join((
        f'sample {checkpoint.replica[0]}',
        f'replica {checkpoint.replica[1]}',
        checkpoint.fit_mode,
        f'model {model_variant(checkpoint)}',
    ))
    return first_line


def plot_quantity(name, time, values, checkpoint, v_max=None, suptitle=None,
                  capacity_pathway=None, pathway=None,
                  log_scale = False):
    """Plot one modeled pool or pathway quantity in the fit-plot style."""
    figure, axis = plt.subplots(figsize=(4, 4))
    if name.endswith('_deltaG_r'):
        values = np.asarray(values, dtype=float) * GIBBS_ENERGY_SCALE
    if capacity_pathway is not None:
        if v_max is None:
            raise ValueError('A pathway-capacity plot requires v_max.')
        values = np.asarray(values, dtype=float) * v_max
        axis.plot(time, values, '-', color='k', label=r'$M(t)\,v_{\max}$')
    else:
        axis.plot(time, values, '-', color='k',
                  label=(delta_g_legend_label(name) if name.endswith('_deltaG_r')
                         else r'$v(t)$' if name.endswith('_v') else name))
    if v_max is not None and capacity_pathway is None:
        axis.axhline(v_max, color='k', linestyle='--',
                     label=r'$v_{\mathrm{max}}$')
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
        axis.legend(handles=[
            Line2D([], [], color='k', linestyle='-', label=r'$M(t)\,v_{\max}$'),
            Line2D([], [], color='none', label=r'$v_{\max}=1\ \text{(constant)}$'),
        ], fancybox=False, edgecolor='k', loc='best')
    else:
        axis.set_title(quantity_title(name, pathway=pathway))
    axis.set_xlabel(r'$\mathrm{time\ [d]}$')
    axis.set_ylabel(quantity_ylabel(name, capacity_pathway))
    if name.endswith(('_thermodynamic_factor', '_MM', '_inhib')):
        axis.set_ylim(0, 1)
    plot.format_ax(axis, log_scale)
    if name.endswith('_deltaG_r'):
        lower = gibbs_lower_limit(time, values)
        if lower is not None:
            axis.set_ylim(lower, axis.get_ylim()[1])
        constrain_gibbs_axis(axis)
    if name.endswith('_MM') and axis.get_yscale() == 'log':
        axis.set_ylim(1e-6, 1.)
    figure.suptitle(figure_suptitle(checkpoint) if suptitle is None else suptitle)
    figure.tight_layout(rect=(0, 0, 1, 0.85))
    return figure, axis


def reaction_rate_signature_figure(checkpoint, run_log, *, suptitle=None,
                                   log_scale=False, pathways=None):
    """Plot the FPCA-style, biomass-scaled pathway-activity signature.

    A signature contains one actual ``M*v`` curve for every active pathway.
    It deliberately omits the separate ``M*v_max`` reference curves: those
    are capacities rather than elements of the reaction-rate signature.
    """
    figure, axis = plt.subplots(figsize=(8, 4.5))
    v_max_values = pathway_v_max_values(run_log)
    pathways = set(enabled_pathways() if pathways is None else pathways)
    for pathway, time, actual_rate, _ in pathway_reaction_rate_curves(
            run_log, v_max_values):
        if pathway not in pathways:
            continue
        axis.plot(time, actual_rate, color=pathway_color(pathway),
                  label=pathway_display_name(pathway))
    axis.set_title(r'$M(t)\,v(t)$ pathway-activity signature')
    axis.set_xlabel(r'$\mathrm{time\ [d]}$')
    axis.set_ylabel(rf'$\mathrm{{reaction\ rate\ [{plot.REACTION_RATE_UNIT}]}}$')
    axis.legend(fancybox=False, edgecolor='k', loc='best')
    plot.format_ax(axis, log_scale=log_scale)
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


def format_michaelis_menten_km(value):
    """Format a checkpoint's half-saturation constant for a legend."""
    value = float(value)
    return f'{value:.0f}' if abs(value) > 100 else f'{value:.2f}'


def plot_pathway_limitation_factors(pathway, mm_series, thermodynamic_series,
                                    inhibition_series, checkpoint, *,
                                    mm_component_series=(), suptitle=None):
    """Plot all active limitation factors and their product for one pathway."""
    figure, axis = plt.subplots(figsize=(4, 4))
    if mm_series is None:
        axis.text(.02, .96, 'Michaelis-Menten factor not logged', color='0.45',
                  transform=axis.transAxes, ha='left', va='top', fontsize=8)
    else:
        time, values = mm_series
        for (substance, (component_time, component_values), km), color in zip(
                mm_component_series,
                plt.rcParams['axes.prop_cycle'].by_key()['color']):
            km_label = (rf',\ K_{{M}} = '
                        rf'{format_michaelis_menten_km(km)}'
                        if km is not None else '')
            axis.plot(component_time, component_values, color=color, linewidth=1.25,
                      linestyle='--',
                      label=(rf'$f_{{\mathrm{{MM}}}}('
                             rf'{plot.chemical_formula(substance)}){km_label}$'))
    if thermodynamic_series is None:
        axis.text(.02, .84, 'Thermodynamic factor not logged', color='0.45',
                  transform=axis.transAxes, ha='left', va='top', fontsize=8)
    else:
        time, values = thermodynamic_series
        axis.plot(time, values, color='tab:orange', linestyle='-.', linewidth=1.5,
                  label=r'$f_{\mathrm{G}}$')
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
        axis.plot(time, total_factor, color='k', linestyle='-', linewidth=2,
                  label=r'$\mathrm{all\ factors}\ [M(t)v(t)]/[M(t)v_{\max}]$')
    axis.axhline(1, color='0.65', linestyle=':', linewidth=1,
                 label=r'$\mathrm{no\ limitation}$')
    axis.set_title(rf'${plot.pathway_name(pathway)}\ \mathrm{{limitation\ factors}}$')
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
        component_prefix = f'{mm_key}_'
        mm_components = []
        for key in sorted(keys):
            if not key.startswith(component_prefix):
                continue
            substance = key.removeprefix(component_prefix)
            km_key = f'{pathway}_Km_{substance}'
            km = None
            if km_key in keys:
                _, km_values = run_log[km_key]
                finite_values = np.asarray(km_values, dtype=float)
                finite_values = finite_values[np.isfinite(finite_values)]
                if finite_values.size:
                    km = finite_values[0]
            mm_components.append((substance, run_log[key], km))
        figures.append((f'{pathway}_limitation_factors',
                        plot_pathway_limitation_factors(
                            pathway,
                            run_log[mm_key] if mm_key in keys else None,
                            (run_log[thermodynamic_key]
                             if thermodynamic_key in keys else None),
                            run_log[inhibition_key] if inhibition_key in keys else None,
                            checkpoint, mm_component_series=tuple(mm_components),
                            suptitle=suptitle)))
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


def _activity_stoichiometry(pathway):
    """Return the pool change per unit pathway activity.

    Model v0.4 represents the reaction and growth-associated anabolism as
    separate vectors.  Older checkpoint models retain their combined vector.
    """
    if hasattr(pathway, 'stoichiometry'):
        return (np.asarray(pathway.stoichiometry, dtype=float)
                + np.asarray(getattr(pathway, 'anabolism', 0.), dtype=float))
    return np.asarray(pathway.pathway_vector, dtype=float)


def pathway_substrate_consumption_rate(pathway_model, run_log, pathway_name,
                                       substrate):
    """Return a pathway's positive net consumption rate of one substrate.

    Growth-associated carbon assimilation is included, so acetate uptake is
    reported as the total model consumption attributed to that pathway.
    """
    pathway = _pathway_from_model(pathway_model, pathway_name)
    coefficient = -float(_activity_stoichiometry(pathway)[system.index(substrate)])
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


def pathway_substance_rate(pathway_model, run_log, pathway_name, substance):
    """Return a pathway's signed net rate of change for one substance.

    When complete system-state logs are available, rates are recalculated on
    accepted solver states.  This matches the total-pool curve and avoids
    integrating diagnostic rates from rejected solver trial states.
    """
    pathway = _pathway_from_model(pathway_model, pathway_name)
    coefficient = float(_activity_stoichiometry(pathway)[system.index(substance)])

    keys = _run_log_keys(run_log)
    if (hasattr(pathway_model, 'compute_dissolved') and
            all(pool in keys for pool in system.SYSTEM)):
        time = np.asarray(run_log[substance][0], dtype=float)
        states = np.empty((len(system.SYSTEM), len(time)))
        for index, pool in enumerate(system.SYSTEM):
            pool_time, values = run_log[pool]
            pool_time = np.asarray(pool_time, dtype=float)
            values = np.asarray(values, dtype=float)
            values = litre_pool_values(run_log, pool, values)
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
    biomass = litre_pool_values(run_log, microbe_pool, biomass)
    order = np.argsort(biomass_time)
    biomass_on_rate_grid = np.interp(time, biomass_time[order], biomass[order])
    return time, coefficient * biomass_on_rate_grid * specific_rate


def pathway_substance_production_rate(pathway_model, run_log, pathway_name,
                                      substance):
    """Return a pathway's positive net production rate of one substance."""
    pathway = _pathway_from_model(pathway_model, pathway_name)
    if float(_activity_stoichiometry(pathway)[system.index(substance)]) <= 0:
        raise ValueError(f'{pathway_name} does not produce {substance}.')
    return pathway_substance_rate(pathway_model, run_log, pathway_name,
                                  substance)


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


def cumulative_pathway_consumption(time, rate):
    """Integrate the negative portion of a signed pathway rate."""
    return cumulative_pathway_production(time, -np.asarray(rate, dtype=float))


def cumulative_pathway_contributions(curves, total_time, total_pool):
    """Integrate signed pathway rates and close them to a modeled pool.

    The solver advances the pool through internal time steps, whereas pathway
    diagnostics are evaluated only at the report time grid.  Independently
    applying the trapezoidal rule to each sparse diagnostic curve can therefore
    fail to reproduce the solved pool.  Here every pathway is integrated on a
    shared grid and the small interval-wise quadrature residual is distributed
    in proportion to the pathways' absolute fluxes.  Consequently, the signed
    cumulative contributions sum *exactly* to ``total_pool - total_pool[0]``
    at every plotted time while preserving the rate-derived pathway shares.
    """
    total_time = np.asarray(total_time, dtype=float)
    total_pool = np.asarray(total_pool, dtype=float)
    valid_total = np.isfinite(total_time) & np.isfinite(total_pool)
    total_time, total_pool = total_time[valid_total], total_pool[valid_total]
    order = np.argsort(total_time)
    total_time, total_pool = total_time[order], total_pool[order]
    total_time, indices = np.unique(total_time, return_index=True)
    total_pool = total_pool[indices]

    time = np.unique(np.concatenate(
        [total_time] + [np.asarray(curve_time, dtype=float)
                        for _, curve_time, _ in curves]))
    pool = np.interp(time, total_time, total_pool)
    rates = np.asarray([
        np.interp(time, np.asarray(curve_time, dtype=float),
                  np.asarray(rate, dtype=float), left=0., right=0.)
        for _, curve_time, rate in curves
    ])
    increments = np.diff(time) * (rates[:, :-1] + rates[:, 1:]) / 2
    residual = np.diff(pool) - np.sum(increments, axis=0)
    weights = np.abs(increments)
    weight_sum = np.sum(weights, axis=0)
    nonzero = weight_sum > 0
    if np.any(nonzero):
        increments[:, nonzero] += (
            residual[nonzero] * weights[:, nonzero] / weight_sum[nonzero])
    if np.any(~nonzero & (np.abs(residual) > 0)):
        # This should be unreachable for a pool governed by these pathways,
        # but retaining the residual makes an incomplete legacy log explicit.
        increments[0, ~nonzero] += residual[~nonzero]
    cumulative = np.column_stack((np.zeros(len(curves)), np.cumsum(increments, axis=1)))
    return time, cumulative, pool


def gas_pathway_production_figure(checkpoint, pathway_model, run_log, gas, *,
                                  suptitle=None, log_scale=False, pathways=None):
    """Show an absolute gas pool and its cumulative pathway contributions.

    Product pathways stack above the initial pool.  Each consuming pathway is
    a same-colour band back down from that production stack, so its lower edge
    reconstructs the total pool shown by the dashed line.
    """
    keys = _run_log_keys(run_log)
    if gas not in keys:
        return None

    producers, consumers = [], []
    pathways = enabled_pathways() if pathways is None else pathways
    for pathway in pathways:
        try:
            rate_time, rate = pathway_substance_rate(
                pathway_model, run_log, pathway, gas)
        except ValueError:
            continue
        volume_scale = litres_per_gdw(run_log)
        if volume_scale is not None:
            rate = np.asarray(rate, dtype=float) * volume_scale
        coefficient = float(_activity_stoichiometry(_pathway_from_model(
            pathway_model, pathway))[system.index(gas)])
        if coefficient > 0:
            producers.append((pathway, rate_time, rate))
        elif coefficient < 0:
            consumers.append((pathway, rate_time, rate))

    if not producers:
        return None

    total_time, total_pool = run_log[gas]
    total_time = np.asarray(total_time, dtype=float)
    total_pool = np.asarray(total_pool, dtype=float)
    curves = producers + consumers
    time, cumulative, pool = cumulative_pathway_contributions(
        curves, total_time, total_pool)
    contributions = [curve for curve in cumulative[:len(producers)]]
    consumer_contributions = [
        -curve for curve in cumulative[len(producers):]
    ]
    initial_pool = float(total_pool[0])
    production_top = initial_pool + np.sum(contributions, axis=0)

    figure, axes = plot.create_figure(plot.FigureSpec(
        axis_names=('diagnostic',), figsize=(5, 4)))
    axis = axes['diagnostic']
    labels = [
        rf'$\mathrm{{{pathway} pathway}}$' + f': {contribution[-1]:.3g}'
        for (pathway, _, _), contribution in zip(producers, contributions)
    ]
    production_boundary = initial_pool * np.ones_like(time)
    for (pathway, _, _), contribution, label in zip(
            producers, contributions, labels):
        next_boundary = production_boundary + contribution
        axis.fill_between(
            time, production_boundary, next_boundary,
            color=pathway_color(pathway), alpha=.75, linewidth=0, label=label)
        production_boundary = next_boundary
    consumption_boundary = production_top.copy()
    for (pathway, _, _), contribution in zip(consumers, consumer_contributions):
        next_boundary = consumption_boundary - contribution
        axis.fill_between(
            time, next_boundary, consumption_boundary,
            color=pathway_color(pathway), alpha=.75, linewidth=0,
            label=(rf'$\mathrm{{{pathway}\ pathway\ ({gas}\ consumed)}}$'
                   f': {-contribution[-1]:.3g}'))
        consumption_boundary = next_boundary
    axis.plot(time, pool, 'k--', linewidth=1.5,
              label=(rf'$\mathrm{{total\ {plot.chemical_formula(gas)}}}$'
                     f': {total_pool[-1]:.3g}'))
    axis.set_title(rf'$\mathrm{{{plot.chemical_formula(gas)}\ pool\ and\ pathway\ contributions}}$')
    axis.set_xlabel(r'$\mathrm{time [d]}$')
    axis.set_ylabel(rf'$\mathrm{{{plot.chemical_formula(gas)}\ [{plot.GAS_SUBSTANCE_UNIT}]}}$')
    plot.format_ax(axis, log_scale=log_scale)
    axis.legend(
        title=rf'$\mathrm{{final\ {plot.chemical_formula(gas)}\ [{plot.GAS_SUBSTANCE_UNIT}]}}$',
        fancybox=False, edgecolor='k', loc='best')
    figure.suptitle(figure_suptitle(checkpoint) if suptitle is None else suptitle)
    figure.subplots_adjust(left=.19, right=.95, bottom=.16, top=.80)
    return figure, axis


def ch4_pathway_production_figure(checkpoint, run_log, *, suptitle=None,
                                  log_scale=False, pathways=None):
    """Show the absolute CH4 pool with cumulative pathway contributions."""
    keys = _run_log_keys(run_log)
    curves = []
    pathways = enabled_pathways() if pathways is None else pathways
    for pathway in pathways:
        key = f'{pathway}_CH4 from {PATHWAY_MICROBE_POOL[pathway]}'
        if key in keys:
            rate_time, rate = run_log[key]
            volume_scale = litres_per_gdw(run_log)
            if volume_scale is not None:
                rate = np.asarray(rate, dtype=float) * volume_scale
            if len(rate_time):
                curves.append((pathway, rate_time, rate))
    if not curves or 'CH4' not in keys:
        return None

    total_time, total_ch4 = run_log['CH4']
    total_ch4 = np.asarray(total_ch4, dtype=float)
    time, cumulative, pool = cumulative_pathway_contributions(
        curves, total_time, total_ch4)
    contributions = [curve for curve in cumulative]
    figure, axes = plot.create_figure(plot.FigureSpec(
        axis_names=('diagnostic',), figsize=(5, 4)))
    axis = axes['diagnostic']
    labels = [
        rf'$\mathrm{{{pathway} pathway}}$' + f': {contribution[-1]:.3g}'
        for (pathway, _, _), contribution in zip(curves, contributions)
    ]
    production_boundary = float(total_ch4[0]) * np.ones_like(time)
    for (pathway, _, _), contribution, label in zip(curves, contributions, labels):
        next_boundary = production_boundary + contribution
        axis.fill_between(
            time, production_boundary, next_boundary,
            color=pathway_color(pathway), alpha=.75, linewidth=0, label=label)
        production_boundary = next_boundary
    axis.plot(time, pool, 'k--', linewidth=1.5,
              label=(r'$\mathrm{total\ CH_4}$' f': {total_ch4[-1]:.3g}'))
    axis.set_title(r'$\mathrm{CH_4\ pool\ and\ pathway\ contributions}$')
    axis.set_xlabel(r'$\mathrm{time\ [d]}$')
    axis.set_ylabel(rf'$\mathrm{{CH_4\ [{plot.GAS_SUBSTANCE_UNIT}]}}$')
    plot.format_ax(axis, log_scale=log_scale)
    axis.legend(title=rf'$\mathrm{{final\ CH_4\ [{plot.GAS_SUBSTANCE_UNIT}]}}$',
                fancybox=False, edgecolor='k', loc='best')
    figure.suptitle(figure_suptitle(checkpoint) if suptitle is None else suptitle)
    figure.subplots_adjust(left=.19, right=.95, bottom=.16, top=.80)
    return figure, axis


def co2_pathway_production_figure(checkpoint, pathway_model, run_log, *,
                                  suptitle=None, log_scale=False, pathways=None):
    """Show the absolute CO2 pool with producing and consuming pathways."""
    return gas_pathway_production_figure(
        checkpoint, pathway_model, run_log, 'CO2', suptitle=suptitle,
        log_scale=log_scale, pathways=pathways)


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
    if factor_limits and log_scale and axis.get_yscale() == 'log':
        axis.set_ylim(1e-6, 1.)
    if gibbs_minimum:
        lower_limits = [gibbs_lower_limit(time, values)
                        for _, time, values in curves]
        lower_limits = [lower for lower in lower_limits if lower is not None]
        if lower_limits:
            axis.set_ylim(min(lower_limits), axis.get_ylim()[1])
        constrain_gibbs_axis(axis)

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


def pathway_delta_g_figures(checkpoint, run_log, *, suptitle=None,
                            log_scale=False, pathways=None):
    """Return combined ΔGᵣ plots for all active pathways and Hydro/Homo.

    The Fe3/Aceto diagnostic remains intentionally separate: it is useful for
    examining their competing acetate use.  Avoid adding a duplicate "all"
    plot to models that only logged that same two-pathway pair.
    """
    keys = _run_log_keys(run_log)

    def figure(pathways, name, title):
        return (name, plot_pathway_pair(
            checkpoint,
            [(pathway, *run_log[f'{pathway}_deltaG_r']) for pathway in pathways],
            title,
            rf'$\Delta G_\mathrm{{r}}\ [{GIBBS_ENERGY_UNIT}]$',
            suptitle=suptitle,
            gibbs_minimum=True,
            log_scale=log_scale))

    pathways = set(enabled_pathways() if pathways is None else pathways)
    active = tuple(pathway for pathway in PATHWAY_PLOT_ORDER
                   if pathway in pathways
                   if f'{pathway}_deltaG_r' in keys)
    figures = []
    if len(active) > 2:
        figures.append(figure(
            active, 'All_pathways_deltaG_r',
            r'$\Delta G_\mathrm{r}$ (all pathways)'))

    hydro_homo = ('Hydro', 'Homo')
    if all(pathway in pathways and f'{pathway}_deltaG_r' in keys
           for pathway in hydro_homo):
        figures.append(figure(
            hydro_homo, 'Hydro_Homo_deltaG_r',
            r'$\Delta G_\mathrm{r}$ (Hydro and Homo pathways)'))
    return figures


def fe3_aceto_diagnostic_figures(checkpoint, model_run, *, suptitle=None, log_scale = False):
    """Return shared-axis Fe3/Aceto pathway diagnostics when both are active."""
    pathway_model, _, run_log = model_run
    keys = _run_log_keys(run_log)
    pathways = ('Fe3', 'Aceto')
    if not all(pathway in enabled_pathways() for pathway in pathways):
        return []
    figures = []

    microbe_keys = tuple(PATHWAY_MICROBE_POOL[pathway] for pathway in pathways)
    if all(name in keys for name in microbe_keys):
        figures.append(('Fe3_Aceto_microbes', plot_pathway_pair(
            checkpoint,
            [(pathway, run_log[PATHWAY_MICROBE_POOL[pathway]][0],
              litre_pool_values(run_log, PATHWAY_MICROBE_POOL[pathway],
                                run_log[PATHWAY_MICROBE_POOL[pathway]][1]))
             for pathway in pathways],
            r'$\mathrm{living\ microbes}$',
            rf'$\mathrm{{microbial\ biomass\ [{MICROBIAL_BIOMASS_UNIT}]}}$',
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
            acetate_pool = (run_log['Acetate'][0],
                            litre_pool_values(run_log, 'Acetate', run_log['Acetate'][1]),
                            r'$\mathrm{total\ Acetate\ pool}$')
        figures.append(('Fe3_Aceto_Acetate_consumption', plot_pathway_pair(
            checkpoint, acetate_curves,
            r'$\mathrm{Acetate\ consumption}$',
            rf'$\mathrm{{Acetate\ consumption\ [{plot.REACTION_RATE_UNIT}]}}$',
            suptitle=suptitle,
            secondary_curve=acetate_pool,
            secondary_ylabel=rf'$\mathrm{{total\ Acetate\ [{plot.SUBSTANCE_UNIT}]}}$',
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


def sample_time_limits(sample):
    """Return one padded time range spanning every displayed measurement."""
    times = []
    for replica in sample.replicas:
        time, _ = replica.CO2()
        times.append(time)
    return plot.time_limits(times)


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


def plot_fit_measurements(checkpoint, dataset, axes, pools=('CO2', 'CH4')):
    """Draw selected fit/validation measurements and return legend handles.

    ``pools`` lets an interactive viewer render either gas on its own while
    preserving the marker and legend conventions of the standard fit figure.
    """
    _, validation_replica, fitted_replicas = fit_replicas(checkpoint, dataset)
    markersize = 3
    markers = cycle(('^', 's', 'D', 'v', 'P', 'X'))
    handles = []
    for fitted_replica in fitted_replicas:
        marker = next(markers)
        for pool in pools:
            plot.plot_pool_data(
                fitted_replica, pool, axes[pool], marker=marker, ms=markersize,
                mfc='none')
        handles.append(Line2D(
            [], [], color='k', marker=marker, linestyle='None',
            markerfacecolor='none', markersize=markersize,
            label=f'fit replica {fitted_replica.replica_number}',
        ))

    if checkpoint.fit_mode == 'split':
        for pool in pools:
            plot.plot_pool_data(
                validation_replica, pool, axes[pool], marker='o', ms=markersize)
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


def format_fit_axes(checkpoint, sample, axes, *, y_limits=None, x_limits=None):
    """Apply data-defined scales and limits to a CO2/CH4 fit plot.

    The fit display is independent of the objective's loss transform: CO2 is
    always linear and CH4 is always logarithmic.  Limits deliberately come
    only from measured replicas, so an implausible model trajectory cannot
    expand or compress the measurement scale.
    """
    del checkpoint  # Kept in the signature for existing callers.
    if y_limits is None:
        y_limits = sample_y_limits(sample)
    if x_limits is None:
        x_limits = sample_time_limits(sample)
    for pool, axis in axes.items():
        plot.format_ax(axis)
        if pool == 'CH4':
            axis.set_yscale('log')
        else:
            axis.set_yscale('linear')
        if pool in y_limits:
            lower, upper = y_limits[pool]
            axis.set_ylim(lower if pool == 'CH4' else min(0., lower), upper)
        plot.format_time_axis(axis, x_limits)


def finish_fit_figure(figure, suptitle):
    """Apply the fixed report layout shared by fit comparison figures."""
    figure.suptitle(suptitle)
    figure.subplots_adjust(
        left=0.12, right=0.97, bottom=0.14, top=0.78, wspace=0.30)



def plot_fit(checkpoint, dataset, model_run=None):
    """Plot fitted replicas, validation data, and a model run for a checkpoint."""
    if model_run is None:
        model_run = run(checkpoint.run_config, checkpoint.parameters,
                        log_prediction=True)
    pathway_model, objective, run_log = model_run
    sample, validation_replica, fitted_replicas = fit_replicas(checkpoint, dataset)
    print('plot_fit', sample, validation_replica)
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
                           *, x_limits=None, overwrite=True):
    """Save the non-fit plots shared by single, split, and multi fits.

    ``target(number, quantity)`` supplies the filename while this function
    retains the established ``plot_best`` diagnostic ordering and styling.
    """
    paths = []

    def save(figure, path):
        for axis in figure.axes:
            plot.format_time_axis(axis, x_limits)
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
        if not plot_enabled('quantities', pathway_context):
            continue
        if quantity in MODELLED_GAS_POOLS and not plot_enabled('gas_pools'):
            continue
        pathway = (pathway_context if quantity.startswith('M_')
                   else quantity[:-2] if quantity.endswith('_v') else None)
        capacity_pathway = (
            pathway_context if quantity.startswith('M_')
            and pathway_context is not None
            and pathway_uses_capacity_parameter(run_log, pathway_context)
            else None
        )
        if (capacity_pathway is not None and
                not plot_enabled('pathway_capacity', capacity_pathway)):
            continue
        for log_scale in (False, True):
            figure, axis = plot_quantity(
                quantity, time, values, checkpoint,
                v_max=(v_max_values.get(pathway)
                       if quantity.endswith('_v') or capacity_pathway is not None
                       else None),
                suptitle=suptitle, capacity_pathway=capacity_pathway,
                pathway=pathway_context,
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

    rate_pathways = enabled_pathways()
    if (plot_enabled('reaction_rate_signature') and rate_pathways and
            any(pathway in rate_pathways
                for pathway, *_ in pathway_reaction_rate_curves(run_log, v_max_values))):
        # The joint pathway-activity signature replaces the former six separate
        # reaction-rate panels.  Keep both scales: the linear view exposes the
        # dominant pathway, while the log view retains the low-activity curves.
        for log_scale in (False, True):
            figure, axis = reaction_rate_signature_figure(
                checkpoint, run_log, suptitle=suptitle, log_scale=log_scale,
                pathways=rate_pathways)
            if log_scale and axis.get_yscale() != 'log':
                plt.close(figure)
                continue
            quantity = 'reaction_rate_signature' + ('_log' if log_scale else '')
            save(figure, target(last_number + 1, quantity))

    if plot_enabled('limitation_factors'):
        for number, (quantity, (figure, _)) in enumerate(
                pathway_limitation_factor_figures(
                    checkpoint, run_log, suptitle=suptitle), start=80):
            pathway = quantity.removesuffix('_limitation_factors')
            if not plot_enabled('limitation_factors', pathway):
                plt.close(figure)
                continue
            path = target(number, quantity)
            save(figure, path)

    if plot_enabled('pathway_delta_g'):
        for log_scale in (False, True):
            for number, (quantity, (figure, _)) in enumerate(
                    pathway_delta_g_figures(
                        checkpoint, run_log, suptitle=suptitle,
                        log_scale=log_scale), start=88):
                if log_scale and figure.axes[0].get_yscale() != 'log':
                    plt.close(figure)
                    continue
                path = target(number, quantity + ('_log' if log_scale else ''))
                save(figure, path)

    if plot_enabled('fe3_aceto_diagnostics'):
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
        if not plot_enabled('ch4_pathway_production'):
            break
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
        if not plot_enabled('co2_pathway_production'):
            break
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
        checkpoint.run_config['model']['version'], chosen['pathways'],
        chosen['normalized_parameters'])
    pathway_model.use_gas_limits = bool(chosen.get('gas_limits', False))
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
        suptitle = ' | '.join((
            f'joint multi fit ({replica_names})',
            f'sample {replica.sample.sample_name}',
            f'replica {replica.replica_number}',
            f'model {model_variant(plot_checkpoint)}',
        ))
        finish_fit_figure(figure, suptitle)
        directory, filename = multi_plot_target(checkpoint, replica)
        os.makedirs(directory, exist_ok=True)
        if plot_enabled('fit'):
            figure, _ = plot_fit(plot_checkpoint, dataset, model_run=model_run)
            path = os.path.join(directory, filename)
            figure.savefig(path, dpi=300)
            plt.close(figure)
            paths.append(path)

        def target(number, quantity):
            _, diagnostic_filename = multi_plot_target(
                checkpoint, replica, number, quantity)
            return os.path.join(directory, diagnostic_filename)

        paths.extend(_save_diagnostic_plots(
            plot_checkpoint, model_run, directory, target, suptitle,
            x_limits=sample_time_limits(replica.sample)))
    return paths
