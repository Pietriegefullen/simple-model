import os
import cProfile
import multiprocessing
import queue
import uuid
from pathlib import Path

import numpy as np
import scipy.integrate
import matplotlib.pyplot as plt
import json

import system
import optimizer
from optimizer import r2
import parameters
from parameters import Parameter
from model_run import ModelRun

import USER_VARIABLES
import CONSTANTS
from CONSTANTS import EPS

from . import pathways
from . import chemistry

# A candidate that cannot be integrated within this budget is not useful to
# the optimiser.  ``solve_ivp`` runs in a child process so this is a hard
# deadline, including for LSODA's Fortran internals.
INTEGRATION_TIMEOUT_SECONDS = 10.

# When set, every integration watchdog process writes one cProfile data file
# here.  It is deliberately environment-controlled so normal fits pay no
# profiling or file-system cost.
INTEGRATION_PROFILE_DIRECTORY_ENV = 'SOIL_INTEGRATION_PROFILE_DIRECTORY'

# Stop only trajectories that are wildly outside the scale of the particular
# replica.  The floor avoids a zero limit for all-zero gas measurements.
GAS_RUNAWAY_MULTIPLIER = 100.0
GAS_RUNAWAY_FLOOR = 1.0


version_major = '0'
version_minor = '4.2' # minor 2: Henry's law adjusted
                      #       3: allow zero pool values.
                      #       4: use micromoles per litre of solvent internally.
                      #     4.1: retain evaluated parameter values across rebuilds.
                      #     4.2: keep replica pools out of the optimiser state.
__version__ = '.'.join([version_major, version_minor])

# exposed:
# build_model
# integrate
# pathway_by_name
# Model

def build_model(pathways, normalized):
    pathway_model = Model(pathways, normalized)
    return pathway_model

def pathway_by_name(name):
    pwys = {'Homo': pathways.Homo,
            'Aceto': pathways.Aceto,
            'Hydrolysis': pathways.Hydrolysis,
            'Fe3': pathways.Fe3,
            'Hydro': pathways.Hydro,
            'Fermentation': pathways.Fermentation}
    return pwys[name]


def gas_upper_limits(replica, initial_state):
    """Return conservative, replica-specific CO2 and CH4 abort limits."""
    limits = {}
    L_per_gdw = (replica.water_content*1e-3)/replica.dry_weight # L/gdw
    for pool in ('CO2', 'CH4'):
        _, measured = getattr(replica, pool)()
        measured = np.asarray(measured, dtype=float)
        finite_measured = measured[np.isfinite(measured)]
        observed_peak = (np.max(np.abs(finite_measured))
                         if finite_measured.size else 0.0)/L_per_gdw
        initial_value = initial_state[system.index(pool)]
        scale = max(GAS_RUNAWAY_FLOOR, observed_peak, abs(initial_value))
        limits[pool] = GAS_RUNAWAY_MULTIPLIER * scale
    return limits


def gas_runaway_events(limits):
    """Build terminal solve_ivp events for excessive gas concentrations."""
    events = []
    for pool in ('CO2', 'CH4'):
        pool_index = system.index(pool)
        upper_limit = limits[pool]

        def event(_time, state, index=pool_index, limit=upper_limit):
            return limit - state[index]

        event.terminal = True
        event.direction = -1
        events.append(event)
    return events


def _solve_segment(f, t_span, initial_state, t_eval, events):
    """Solve one continuous segment and translate every abnormal exit."""
    try:
        result = scipy.integrate.solve_ivp(
            f, t_span, initial_state,
            t_eval=t_eval,
            events=events,
            method='LSODA',
            rtol = 1e-6,
            atol = 1e-9,
            max_step=1.0,
            #first_step=1e-6,
        )
    except optimizer.IntegrationFailure:
        raise
    except Exception as error:
        raise optimizer.IntegrationFailure(
            f'Integration raised {type(error).__name__}: {error}') from error

    if not result.success:
        raise optimizer.IntegrationFailure(f'Integration failed: {result.message}')
    if any(event_times.size for event_times in result.t_events):
        triggered = 'CO2' if result.t_events[0].size else 'CH4'
        raise optimizer.IntegrationFailure(
            f'Integration stopped because {triggered} exceeded its runaway limit')
    return result


def integrate(f, t, S0, solver_result, reset_Fe3, gas_limits):
    #print('solving IVP')
    events = [] if gas_limits is None else gas_runaway_events(gas_limits)

    t_after = None
    if isinstance(reset_Fe3, int):
        if reset_Fe3 >= max(t):
            raise Exception('Trying to reset iron after run end.')
        t = np.sort(np.unique(np.concatenate([t, (reset_Fe3,)])))
        t_before = t[t <= reset_Fe3]
        t_after = t[t >= reset_Fe3]
        if not reset_Fe3 in t:
            t_before = np.concatenate([t_before, [reset_Fe3]])
            t_after = np.concatenate([[reset_Fe3], t_after], axis = 0)
    else:
        t_before = t

    result = _solve_segment(f, (0, max(t_before)), S0, t_before, events)

    solver_result.append(result)

    if not t_after is None:
        S1 = result.y[:,-1] # system state on last day before resetting
        fe3_index = system.index('Fe3')
        print('resetting Fe3 on day', min(t_after), 'from', S1[fe3_index], 'to', S0[fe3_index])
        S1[fe3_index] = S0[fe3_index]
        print('solving IVP after reset')
        t_eval = np.sort(np.unique(np.concatenate([t_after,
                                                   np.arange(min(t_after), max(t_after), 10) ])))
        result_after = _solve_segment(
            f, (min(t_after), max(t_after)), S1, t_eval, events)

        solver_result[0].y = np.concatenate([solver_result[0].y,
                                             result_after.y], axis = -1)

        solver_result[0].t = np.concatenate([solver_result[0].t,
                                             result_after.t])

    return solver_result


def _integration_worker(f, t, S0, result_queue, reset_Fe3, gas_limits):
    """Return integration failures to the parent instead of printing a traceback."""
    profiler = None
    profile_directory = os.environ.get(INTEGRATION_PROFILE_DIRECTORY_ENV)
    if profile_directory:
        profiler = cProfile.Profile()
    try:
        solver_result = []
        if profiler is None:
            integrate(f, t, S0, solver_result, reset_Fe3, gas_limits)
        else:
            profiler.runcall(integrate, f, t, S0, solver_result, reset_Fe3, gas_limits)
        result_queue.put(solver_result[0])
    except optimizer.IntegrationFailure as error:
        result_queue.put(error)
    finally:
        if profiler is not None:
            # A unique filename matters because the outer optimiser can start
            # several integration workers concurrently.  If the worker is
            # killed by the hard timeout it cannot flush a profile; successful
            # workers still provide representative RHS/solver measurements.
            directory = Path(profile_directory)
            directory.mkdir(parents=True, exist_ok=True)
            profiler.dump_stats(str(directory / (
                f'integration-{os.getpid()}-{uuid.uuid4().hex}.pstats')))

def get_pathways(model_type):
    basic = ['Hydrolysis',
             'Fermentation',
             'Hydro',
             'Aceto']
    if model_type == 'complex':
        return basic + ['Homo',
                        'Fe3']
    elif model_type == 'simple':
        #return basic
    #else:
        print('model type:', model_type)
        raise NotImplementedError()

class Model():
    @classmethod
    def hypothesis_model_parameters(cls, normalize):
        # Ferm_Km: choose half pool size of minimum TOC content
        #           => ensure f_MM = 1 at t = 0 for ferm 

        # Hydrolysis_Kmb: if M_0 == 1 (i.e. 100%), Hydrolysis_Kmb is relative,
        #                 i.e. 200% for 2.0, 1000% for 10 etc.
        #                 Hence, the actual value depends on the actual amount
        #                 of microbes in each replica, which will remain unknown.


        # Acetate range inspired by bhz data
        # 

        HYDROLYSIS_KMB = None #[1_113, 59_222_344]
        FERM_KM =  None #[0.47, 315_000]

        HYDRO_KM_CO2 = [0.2,2_650_000]
        HOMO_KM_CO2 = [1.88, 2_230_000]

         # sehr affin zu H2 -> eher kleiner Km-Wert?
        HYDRO_KM_H2 =  [13, 3_210_000]
        HOMO_KM_H2 = [2, 1_970_000]

        ACETO_KM_AC = None #[170, 1_010_000]

        FE3_KM_AC = None #[0.315, 3_720_000]
        FE3_KM_FE3 = None #[0.285, 3_730_000]

        p = [
         Parameter('Hydrolysis_Kmb', 10, HYDROLYSIS_KMB, 'log',  normalize = normalize), 
         Parameter('Ferm_Km',        9495/2, FERM_KM, scale= 'log', normalize = normalize),

         Parameter('Fe3_Km_Fe3',     100, FE3_KM_FE3, scale = 'log', normalize = normalize),
         Parameter('Fe3_Km_Ac',      100, FE3_KM_AC, scale = 'log',  normalize = normalize),
         Parameter('Aceto_Km_Ac',    100, ACETO_KM_AC, 'log', normalize = normalize),

         Parameter('Hydro_Km_CO2',   10_000, HYDRO_KM_CO2, scale = 'log', normalize = normalize),
         Parameter('Hydro_Km_H2',    10_000, HYDRO_KM_H2, scale = 'log', normalize = normalize),
         Parameter('Homo_Km_H2',     10_000, HOMO_KM_H2, scale = 'log', normalize = normalize),
         Parameter('Homo_Km_CO2',    10_000, HOMO_KM_CO2, scale = 'log', normalize = normalize),
         Parameter('Ferm_inhibition', 1., [0.698, 754_000], 'log', normalize = normalize), # no source


         Parameter('Hydrolysis_v_max', 0.01, [0.01, 0.1], 'log', normalize = normalize), 

         Parameter('death_rate',  8.3e-5, scale = 'log', normalize = normalize),

         Parameter('Ferm_v_max',     0.01, [0.01, 0.1], 'log', normalize = normalize),
         Parameter('Ferm_CUE',        .0812, [2e-9, 0.15], scale = 'linear', normalize = normalize),
         
         
         Parameter('Hydro_v_max',    0.01, [0.01, 0.1], scale = 'log', normalize = normalize),
         Parameter('Hydro_CUE',       .08, [0.01, 0.5], scale = 'linear', normalize = normalize), 
         
         Parameter('Homo_v_max',     0.01, [0.01, 0.1], scale = 'log', normalize = normalize), # no source
         Parameter('Homo_CUE',        .2, [0.03, 0.5],  scale = 'linear', normalize = normalize),
         
         Parameter('Ac_v_max',       0.01, [0.01, 0.1], scale = 'log', normalize = normalize),
         Parameter('Ac_CUE',          0.00876, [0.00069, 0.5], scale = 'linear', normalize = normalize), 
         
         Parameter('Fe3_v_max',      0.01, [0.01, 0.1], scale = 'log', normalize = normalize), 
         Parameter('Fe3_CUE',        0.05, scale ='linear', normalize = normalize),
         
         Parameter('Acetate',         200, [100, 2_000], scale = 'log', normalize = normalize),
         Parameter('Fe3',            150, [100, 100_000], scale = 'log', normalize = normalize),
        
        # if biomass is variable, it must be in micromol per liter
         Parameter('M_Ferm',        1., scale = 'log',  normalize = normalize),
         Parameter('M_Hydro',       1., scale = 'log', normalize = normalize),
         Parameter('M_Fe3',         1., scale = 'log', normalize = normalize),
         Parameter('M_Homo',        1., scale = 'log', normalize = normalize),
         Parameter('M_Ac',          1., scale = 'log', normalize = normalize),
         
         Parameter('Hydro_thermodynamics', True, {True}, scale = 'linear', normalize = False),
         Parameter('Homo_thermodynamics', True, {True}, scale = 'linear', normalize = False),
         Parameter('Aceto_thermodynamics', True, {True}, scale = 'linear', normalize = False),
         Parameter('Fe3_thermodynamics', True, {True}, scale = 'linear', normalize = False),
            ]
            
        par = parameters.ModelParameters(p)
        return par


    @classmethod
    def envelope_converted_model_parameters(cls, normalize):


        # (nearly) unidentifiable:
        # M_Ac * Ac_v_max * f_MM(Ac; Aceto_Km_Ac)
        # M_Hydro * Hydro_v_max * f_MM(CO2; Hydro_Km_CO2) * f_MM(H2; Hydro_Km_H2)
        # M_Homo * Homo_v_max * f_MM(CO2; Homo_Km_CO2) * f_MM(H2; Homo_Km_H2)
        # M_Fe3 * Fe3_v_max * f_MM(Fe3; Fe3_Km_Fe3) * f_MM(Ac; Fe3_Km_Ac)

        # M_Ferm * Ferm_v_max * f_MM(DOC; Ferm_Km)
        # M_Ferm * Hydrolysis_v_max

        # => fixed Km's: eliminates 8 unidentifiable parameters, 
        #                11 remain. Can they be sufficiently restricted?


        # => merge M and v_max: eliminates 6 and 
        #                   => identifiable for fixed Km.


        p = [
         Parameter('death_rate',  8.3e-5, scale = 'log', normalize = normalize),

         Parameter('Hydrolysis_v_max', 0.01, [0.01, 0.1], 'log', normalize = normalize), 
         Parameter('Hydrolysis_Kmb', 10_000, [1_113, 59_222_344], 'log',  normalize = normalize), # no source


         Parameter('Ferm_v_max',     0.01, [0.01, 0.1], 'log', normalize = normalize),
         Parameter('Ferm_Km',        1_000, [0.47, 315_000], scale= 'log', normalize = normalize),
         Parameter('Ferm_inhibition', 1., [0.698, 754_000], 'log', normalize = normalize), # no source
         Parameter('Ferm_CUE',        .0812, [2e-9, 0.5], scale = 'linear', normalize = normalize),
         
         
         Parameter('Hydro_v_max',    0.01, [0.01, 0.1], scale = 'log', normalize = normalize),
         Parameter('Hydro_CUE',       .08, [0.01, 0.5], scale = 'linear', normalize = normalize), 
         Parameter('Hydro_Km_CO2',   100, [0.2,2_650_000], scale = 'log', normalize = normalize),
         Parameter('Hydro_Km_H2',    50, [13, 3_210_000], scale = 'log', normalize = normalize),
         
         Parameter('Homo_Km_H2',     10, [2, 1_970_000], scale = 'log', normalize = normalize),
         Parameter('Homo_Km_CO2',    10, [1.88, 2_230_000], scale = 'log', normalize = normalize),
         Parameter('Homo_v_max',     0.01, [0.01, 0.1], scale = 'log', normalize = normalize), # no source
         Parameter('Homo_CUE',        .2, [0.03, 0.5],  scale = 'linear', normalize = normalize),
         
         Parameter('Aceto_Km_Ac',    10, [170, 1_010_000], 'log', normalize = normalize),
         Parameter('Ac_v_max',       0.01, [0.01, 0.1], scale = 'log', normalize = normalize),
         Parameter('Ac_CUE',          0.00876, [0.00069, 0.5], scale = 'linear', normalize = normalize), 
         
         Parameter('Fe3_Km_Ac',      11.99, [0.315, 3_720_000], scale = 'log',  normalize = normalize),
         Parameter('Fe3_Km_Fe3',     230, [0.285, 3_730_000], scale = 'log', normalize = normalize),
         Parameter('Fe3_v_max',      0.01, [0.01, 0.1], scale = 'log', normalize = normalize), 
         Parameter('Fe3_CUE',        0.05, scale ='linear', normalize = normalize),
         
         Parameter('Acetate',         20, [0.01, 162_000], scale = 'linear', normalize = normalize), # no source
         Parameter('Fe3',            150, [1_110, 659_000], scale = 'linear', normalize = normalize), # no source
        
        # if biomass is variable, it must be in micromol per liter
         Parameter('M_Ferm',        1., [5e-4, 6e5], scale = 'log',  normalize = normalize),
         Parameter('M_Hydro',       5., [2e-4, 9e5], scale = 'log', normalize = normalize),
         Parameter('M_Fe3',         1., [1e-2, 9e5], scale = 'log', normalize = normalize),
         Parameter('M_Homo',        1., [4e-4, 2e6], scale = 'log', normalize = normalize),
         Parameter('M_Ac',          4., [3e-4, 5e3], scale = 'log', normalize = normalize),
         
         Parameter('Hydro_thermodynamics', True, {True}, scale = 'linear', normalize = False),
         Parameter('Homo_thermodynamics', True, {True}, scale = 'linear', normalize = False),
         Parameter('Aceto_thermodynamics', True, {True}, scale = 'linear', normalize = False),
         Parameter('Fe3_thermodynamics', True, {True}, scale = 'linear', normalize = False),
            ]
            
        par = parameters.ModelParameters(p)
        return par

    @classmethod
    def default_model_parameters(cls, normalize):
        # no source for:
        # =============
        # Ferm_inhibition
        # death_rate
        # Hydrolysis_Kmb
        # Hydrolysis_v_max
        # Ferm_CUE
        # Homo_CUE
        # Homo_v_max
        # Fe3
        # Acetate
        # 
        # Ferm_inhibition has a large influence on inhibition!
        # for large values, there is almost full inhibition.


        p = [
         Parameter('death_rate',  1.29e-4, [1.29e-4,6e-2], scale = 'log', normalize = normalize),

         Parameter('Hydrolysis_v_max', .0083, [1e-8, 0.1], 'log', normalize = normalize), # no source
         Parameter('Hydrolysis_Kmb', 20, [1e-10, 200], 'log',  normalize = normalize), # no source


         Parameter('Ferm_v_max',     0.5533, [0.5533, 97.7322], 'log', normalize = normalize),
         Parameter('Ferm_Km',        999, [10, 15_000], scale= 'log', normalize = normalize),
         Parameter('Ferm_inhibition', 1e-3, [1e-6, 1e-1], 'log', normalize = normalize), # no source
         Parameter('Ferm_CUE',        .0812, scale = 'linear', normalize = normalize),
         
         
         Parameter('Hydro_v_max',    1_728, scale = 'log', normalize = normalize),
         Parameter('Hydro_CUE',       .08, scale = 'linear', normalize = normalize), 
         Parameter('Hydro_Km_CO2',   50, [29,100], scale = 'log', normalize = normalize),
         Parameter('Hydro_Km_H2',    50, [.1, 100], scale = 'log', normalize = normalize),
         
         Parameter('Homo_Km_H2',     10, [10, 100], scale = 'log', normalize = normalize),
         Parameter('Homo_Km_CO2',    10, [10, 100], scale = 'log', normalize = normalize),
         Parameter('Homo_v_max',     10, [1e-3, 1e3], scale = 'log', normalize = normalize), # no source
         Parameter('Homo_CUE',        .2, scale = 'linear', normalize = normalize),
         
         Parameter('Aceto_Km_Ac',    10, [10, 40], scale = 'log', normalize = normalize),
         Parameter('Ac_v_max',       0.977322, [0.977322, 97.7322], scale = 'log', normalize = normalize),
         Parameter('Ac_CUE',          0.00876, [0.00876, 0.3], scale = 'linear', normalize = normalize), 
         
         Parameter('Fe3_Km_Ac',      11.99, [11.99, 230], 'log',  normalize = normalize),
         Parameter('Fe3_Km_Fe3',     230, [230, 5_400_000], 'log', normalize = normalize), # Sulman2022: 1e-7
         Parameter('Fe3_v_max',      0.1296, [0.1296, 47.6], 'log', normalize = normalize), 
         Parameter('Fe3_CUE',        0.05, scale ='linear', normalize = normalize),
         
         Parameter('Acetate',         20, [0, 100], 'linear', normalize = normalize), # no source
         Parameter('Fe3',            150, [0, 300], 'linear', normalize = normalize), # no source
        
        # if biomass is variable, it must be in micromol per liter
         Parameter('M_Ferm',        1., [0.3, 300], scale = 'log',  normalize = normalize),
         Parameter('M_Hydro',       5., [5, 500], scale = 'log', normalize = normalize),
         Parameter('M_Fe3',         1., [0.0002, 360], scale = 'log', normalize = normalize),
         Parameter('M_Homo',        1., [.0002, 500], scale = 'log', normalize = normalize),
         Parameter('M_Ac',          4., [4, 1_000], scale = 'log', normalize = normalize),
         
         Parameter('Hydro_thermodynamics', True, {True}, scale = 'linear', normalize = False),
         Parameter('Homo_thermodynamics', True, {True}, scale = 'linear', normalize = False),
         Parameter('Aceto_thermodynamics', True, {True}, scale = 'linear', normalize = False),
         Parameter('Fe3_thermodynamics', True, {True}, scale = 'linear', normalize = False),
            ]
            
        par = parameters.ModelParameters(p)
        return par

    def __init__(self, pwys, normalized_parameters):
        self.system_state_log = ModelRun()
        self.system_change_log = ModelRun()
        self.system_change_log.rate = True
        self.model_parameters = parameters.ModelParameters()

        pathway_classes = [pathway_by_name(p) if isinstance(p, str) else p
                           for p in pwys]
        self._unbuilt_contributing_pathways = pathway_classes
        self.contributing_pathways = None

        self.t = None
        self.replica = None
        self.use_gas_limits = False

        self.normalize_parameters = normalized_parameters

        self.HCC_s = chemistry.dimensionless_henry_solubility(system.SYSTEM)

        self.build(quiet = True)

    def get_config(self, only_structure = False):
        cfg = {}
        for pathway in self.contributing_pathways:
            cfg[pathway.__class__.__name__] = pathway.get_config(only_structure)
        cfg.update({'version':__version__})
        return cfg

    def build(self, quiet = False):
        self.contributing_pathways = [p(self.model_parameters)
                                      for p in self._unbuilt_contributing_pathways]
        for p in self.contributing_pathways:
            p.inject_logger(self.system_state_log)

        # to initialize model parameters used in initial state
        _ = system.initial_state(None, self.model_parameters)

        normalize = self.normalize_parameters
        default_values = Model.default_model_parameters(normalize)
        for used_parameter in self.model_parameters:
            # ``predict`` rebuilds pathways so that numerical arrays such as
            # Km and anabolism reflect the current candidate.  Rebuilding
            # must never reset a candidate back to its catalogue default:
            # optimisation and checkpointing both depend on the parameter
            # collection retaining the evaluated values.
            if used_parameter.is_unset():
                used_parameter.set(default_values[used_parameter.name])

        for pool in ('CO2', 'H2', 'CH4', 'Fe2'):
            if self.model_parameters[pool].is_unset():
                self.model_parameters[pool].constant(0)

    def compute_dissolved(self, S):
        dissolved_S = S.copy()

        gas_volume = 1e-3 * self.replica.total_gas_volume  # L
        water_volume = 1e-3*self.replica.water_content # L
        lam =  self.HCC_s * water_volume/gas_volume
        dissolved_fractions = lam/(1+lam)
        for pool in ['CO2', 'CH4', 'H2']:
            i = system.index(pool)
            dissolved_S[i] = dissolved_fractions[i] * S[i]

        return dissolved_S

    def __call__(self, t, S):
        S = np.maximum(S, 0)
        dissolved_S = self.compute_dissolved(S) # micromol dissolved per Liter

        # compute rates
        dSj_dt = np.stack([pathway_j(t, dissolved_S)
                            for pathway_j in self.contributing_pathways], axis = -1)
        dS_dt = np.sum(dSj_dt, axis = -1)
        dS_dt = np.where(S == 0, np.maximum(0, dS_dt), dS_dt)

        #print()
        #y = system.SYSTEM
        ##y = ['CH4']
        #for i, p in enumerate(y):
        #    print(f'{t:8.3}', f'{p:20s}',  f'{S[system.index(p)]:8.3}', f'{dS_dt[system.index(p)]:8.3}')
        return dS_dt

    def add_t(self, t):
        if self.t is None:
            self.t = np.squeeze(t)
        else:
            self.t = np.concatenate([np.squeeze(t), self.t])
        self._clean_t()

    def _clean_t(self):
        self.t = np.sort(list(set(self.t)))

    def log_accepted_kinetics(self, t, states):
        # ``solve_ivp`` stores state as ``(pool, time)``.  Convert each time
        # column independently, then rebuild the same layout for ``deltaG_r``.
        dissolved_states = np.column_stack([
            self.compute_dissolved(state) for state in states.T
        ])
        for pathway in self.contributing_pathways:
            pathway_name = pathway.__class__.__name__
            effective_km = pathway.kinetic_half_saturation_constants(dissolved_states)
            mm_factors = pathway.kinetic_factors(dissolved_states)
            mm_factor = np.prod(mm_factors, axis=0)
            self.system_state_log.log(pathway.__class__.__name__ + '_MM',
                                      t, mm_factor)
            km_total = np.zeros_like(mm_factor, dtype=float)
            for pool_name in pathway.kinetic_factor_names():
                pool_index = system.index(pool_name)
                self.system_state_log.log(
                    f'{pathway_name}_MM_{pool_name}', t,
                    mm_factors[pool_index])
                self.system_state_log.log(
                    f'{pathway_name}_Km_{pool_name}', t,
                    effective_km[pool_index])
                km_total += effective_km[pool_index]
            # A convenient total for inspection; individual Km curves above
            # are the values actually used for each substrate.
            self.system_state_log.log(f'{pathway_name}_Km_total', t, km_total)
            if pathway.use_inhib:
                self.system_state_log.log(pathway_name + '_inhib', t,
                                          pathway.kinetic_inhibition(dissolved_states))
            if float(pathway.use_thermodynamics) == 1:
                self.system_state_log.log(pathway_name + '_thermodynamic_factor', t,
                                          pathway.thermodynamics(t, dissolved_states))
            specific_rate = pathway.compute_v(t, dissolved_states)
            self.system_state_log.log(pathway.__class__.__name__ + '_v', t, specific_rate)
            biomass = dissolved_states[pathway.microbe_index]
            self.system_state_log.log(
                pathway_name + '_reaction_rate', t, biomass * specific_rate)
            self.system_state_log.log(
                pathway_name + '_capacity_rate', t, biomass * float(pathway.v_max))
            for pool_index, pool_name in enumerate(system.SYSTEM):
                coefficient = (pathway.stoichiometry + pathway.anabolism)[pool_index]
                if coefficient != 0:
                    self.system_state_log.log(
                        f'{pathway_name}_rate_{pool_name}', t,
                        coefficient * biomass * specific_rate)
            if 'CH4' in pathway.products:
                rates = biomass * (
                    specific_rate * (pathway.stoichiometry + pathway.anabolism)[:, np.newaxis]
                    - pathway.death_rate[:, np.newaxis])
                self.system_state_log.log(
                    pathway_name + '_CH4 from ' + pathway.microbe.name, t,
                    rates[system.index('CH4')])

    def log_accepted_states(self, t, states, *, record_pools=True):
        """Record all accepted pools and their dissolved gas equivalents."""
        dissolved_states = np.column_stack([
            self.compute_dissolved(state) for state in states.T
        ])
        for pool_index, pool_name in enumerate(system.SYSTEM):
            if record_pools:
                self.system_state_log.log(pool_name, t, states[pool_index])
            self.system_state_log.log(
                f'{pool_name}_dissolved', t, dissolved_states[pool_index])


    def log_accepted_thermodynamics(self, t, states):
        """Log reaction Gibbs energies on accepted solver states only.

        Pathway calls occur at LSODA trial states as well as accepted states.
        Those trial states must not appear in diagnostic time series because
        they can yield nonphysical reaction quotients.
        """
        # ``solve_ivp`` stores state as ``(pool, time)``.  Convert each time
        # column independently, then rebuild the same layout for ``deltaG_r``.
        dissolved_states = np.column_stack([
            self.compute_dissolved(state) for state in states.T
        ])
        for pathway in self.contributing_pathways:
            if float(pathway.use_thermodynamics) == 1:
                self.system_state_log.log(
                    pathway.__class__.__name__ + '_deltaG_r',
                    t,
                    pathway.deltaG_r(dissolved_states))

    def predict(self, replica, t = None, t_eval = None, quiet = False, parallel = False,
                reset_Fe3 = None, days_beyond_reset = 1000,
                initial_mean_days = 0, reference = 'gdw', microbe_unit = 'mg',
                gas_limits = None):
        if t is not None and t_eval is not None:
            raise ValueError('Specify either t or t_eval, not both.')

        if t_eval is None:
            t_eval = replica['days']
            if t is not None:
                t_eval = np.concatenate([t_eval, np.atleast_1d(t)])
        t_eval = np.sort(np.unique(np.atleast_1d(t_eval)))

        if not reset_Fe3 is None:
            last_day = reset_Fe3 + days_beyond_reset
            t_eval = np.sort(np.unique(np.append(t_eval, last_day)))

        # prepare for solving
        self.build(quiet = quiet)
        # The v0.4 RHS operates on micromol/L states.  Experimental starting
        # pools are recorded per gdw, so convert them before integration.
        S0 = system.initial_state(
            replica, self.parameters(), initial_mean_days, reference='L')
        self.replica = replica
        self.parameters().check(ignore=('TOC', 'DOC', 'H2O'))
        self.system_state_log.reset(replica)
        self.system_change_log.reset(replica)
        use_gas_limits = self.use_gas_limits if gas_limits is None else bool(gas_limits)
        runaway_limits = gas_upper_limits(replica, S0) if use_gas_limits else None

        # A process watchdog is necessary because LSODA offers neither a
        # wall-clock nor a function-evaluation limit.  The ``parallel``
        # argument is retained for API compatibility; every prediction now
        # receives the same five-second protection.
        result_queue = multiprocessing.Queue(maxsize=1)
        try:
            process = multiprocessing.Process(
                target=_integration_worker,
                args=(self, t_eval, S0, result_queue, reset_Fe3, runaway_limits))
            process.daemon = True
            process.start()
            try:
                solver_result = result_queue.get(timeout=INTEGRATION_TIMEOUT_SECONDS)
            except queue.Empty:
                process.terminate()
                process.join()
                raise optimizer.IntegrationFailure(
                    f'Integration exceeded {INTEGRATION_TIMEOUT_SECONDS:g} seconds')
            # Receipt from the queue means the solver has returned; joining
            # now only reaps the completed worker.
            process.join()
        finally:
            result_queue.close()
            result_queue.join_thread()

        if isinstance(solver_result, optimizer.IntegrationFailure):
            raise solver_result

        # prevent nonpositive pool values
        solver_result.y = np.maximum(solver_result.y, 1e-24)

        # add pool values to log
        L_per_gdw = (replica.water_content*1e-3)/replica.dry_weight # L/gdw
        biomass_molar_mass = CONSTANTS.MOLAR_MASS_BIOMASS # mg/micromol
        for Si, pool_name in zip(solver_result.y, system.SYSTEM):
            if reference == 'gdw':
                Si *= L_per_gdw # transform pools from micromol L-1 to micromol gdw-1
                
            elif reference == 'L':
                pass

            else:
                raise NotImplementedError()

            if pool_name.startswith('M_'):
                if microbe_unit == 'mg':
                    Si *= biomass_molar_mass

                elif microbe_unit == 'micromol':
                    pass

                else:
                    raise NotImplementedError()
                
            self.system_state_log.log(pool_name, solver_result.t, Si)

        self.system_state_log._parameters = self.parameters()
        self.system_state_log._reference = reference
        self.system_state_log._microbe_unit = microbe_unit

        return self.system_state_log

    def predict_with_logging(self, replica, t = None, t_eval = None, quiet = False,
                             parallel = False, reset_Fe3 = None,
                             days_beyond_reset = 1000, initial_mean_days = 0, 
                             reference = 'gdw', microbe_unit = 'mg'):
        """Predict a trajectory and add diagnostics used by replay plots.

        ``predict`` deliberately records only pool trajectories, keeping
        optimisation evaluations free of diagnostic work.  This companion
        method derives pathway kinetics and thermodynamics once from the
        accepted solver states, so every logged series shares the pool time
        grid and is safe to use in replay reports.
        """
        run_log = self.predict(
            replica, t=t, t_eval=t_eval, quiet=quiet, parallel=parallel,
            reset_Fe3=reset_Fe3, days_beyond_reset=days_beyond_reset,
            initial_mean_days=initial_mean_days, reference=reference,
            microbe_unit=microbe_unit)
        missing_pools = [pool for pool in system.SYSTEM if pool not in run_log.keys()]
        if missing_pools:
            raise RuntimeError(
                'Diagnostic prediction is missing accepted pool trajectories: '
                + ', '.join(missing_pools))
        pool_time, _ = run_log[system.SYSTEM[0]]
        states = np.vstack([run_log[pool][1] for pool in system.SYSTEM])
        # Diagnostics evaluate the RHS in its native micromol/L basis, while
        # the returned pool trajectories honour the caller's display units.
        if microbe_unit == 'mg':
            biomass_molar_mass = CONSTANTS.MOLAR_MASS_BIOMASS
            for pool_index, pool_name in enumerate(system.SYSTEM):
                if pool_name.startswith('M_'):
                    states[pool_index] /= biomass_molar_mass
        elif microbe_unit != 'micromol':
            raise NotImplementedError()

        if reference == 'gdw':
            L_per_gdw = (replica.water_content * 1e-3) / replica.dry_weight
            states /= L_per_gdw
        elif reference != 'L':
            raise NotImplementedError()

        # The lean prediction has already logged pools in the requested output
        # units.  Add dissolved concentrations without overwriting them with
        # the internal solver-state representation.
        self.log_accepted_states(pool_time, states, record_pools=False)
        self.log_accepted_thermodynamics(pool_time, states)
        self.log_accepted_kinetics(pool_time, states)

        #run_log = self.system_state_log

        return run_log

    def parameters(self):
        return self.model_parameters

    def __str__(self):
        model_string = f'Model with {len(self.contributing_pathways)} Pathways:\n'
        model_string += len(model_string)*'=' + '\n'
        model_string += '\n'.join([str(p) for p in self.contributing_pathways])
        model_string += '\n'*2
        model_string += str(self.model_parameters)
        model_string += '\n'
        var_cnt = sum([int(p.is_variable()) for p in self.model_parameters])
        nan_cnt = sum([int(p.is_unset()) for p in self.model_parameters])
        model_string += '\n'
        model_string += f'{var_cnt:d} variables\n'
        model_string += f'{nan_cnt:d} unset'
        return model_string

    def save(self, target_directory, file_name):
        cfg = {'pathways': [p.__class__.__name__
                            for p in self.contributing_pathways],
               'parameters': self.parameters().get_config()}
        if not os.path.isdir(target_directory):
            os.makedirs(target_directory)
        with open(os.path.join(target_directory, file_name + '.json'), 'w') as df:
            json.dump(cfg, df, indent = 4)

    def model_type(self):
        simple = sorted(get_pathways('simple'))
        complex = sorted(get_pathways('complex'))
        pwys = sorted([p.__class__.__name__
                      for p in self.contributing_pathways])
        if len(simple) == len(pwys):
            for s,p in zip(simple, pwys):
                if not s == p:
                    raise Exception('Unknown model type')
            return 'simple'
        elif len(complex) == len(pwys):
            for c, p in zip(complex, pwys):
                if not c == p:
                    raise Exception('Unknown model type')
            return 'complex'
        raise Exception('Unknown model type')

    def load(self, file):
        with open(file, 'r') as df:
            cfg = json.load(df)

        self.__init__(cfg['pathways'])
        self.model_parameters.set(cfg['parameters'])
