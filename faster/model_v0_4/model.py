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
version_minor = '4' # minor 2: Henry's law adjusted
                    #       3: allow zero pool values.
                    #       4: convert Km values to micromoles per liter water per dry weight
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
    for pool in ('CO2', 'CH4'):
        _, measured = getattr(replica, pool)()
        measured = np.asarray(measured, dtype=float)
        finite_measured = measured[np.isfinite(measured)]
        observed_peak = (np.max(np.abs(finite_measured))
                         if finite_measured.size else 0.0)
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
    events = gas_runaway_events(gas_limits)

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
    def default_model_parameters(cls, normalize_parameters):
        p = [
             Parameter('Hydrolysis_v_max', .0083, [1e-8, 0.1], 'log', normalize = normalize_parameters),
             Parameter('Hydrolysis_Kmb', 20, [1e-10, 200], 'log',  normalize = normalize_parameters),
             
             Parameter('Ferm_v_max',     .525, [0.0001, 5], 'log', normalize = normalize_parameters),
             #Variable('Ferm_Kmb',       890, [0.0005, 2000]),
             Parameter('Ferm_Km',        83, [0.0001, 100], 'log', normalize = normalize_parameters), #833
             Parameter('Ferm_inhibition', 16, [0.001, 200], 'log', normalize = normalize_parameters),
             Parameter('Ferm_CUE',        .5, [0, 1], 'linear', normalize = normalize_parameters),
             
             Parameter('death_rate',  8.3e-5, scale = 'log', normalize = normalize_parameters),
             
             Parameter('Hydro_Km_CO2',   500, [.0005, 1000], 'log', normalize = normalize_parameters),
             Parameter('Hydro_v_max',    .17, [0.003, 1.], 'log', normalize = normalize_parameters),
             Parameter('Hydro_CUE',       .5, [0, 1], 'linear', normalize = normalize_parameters),
             Parameter('Hydro_Km_H2',    500, [.0005, 1000], 'log', normalize = normalize_parameters),
             
             Parameter('Homo_Km_H2',     500, [0.0005, 1000], 'log', normalize = normalize_parameters),
             Parameter('Homo_Km_CO2',    500, [0.0005, 1000], 'log', normalize = normalize_parameters),
             Parameter('Homo_v_max',      .5, [0.005, 1.], 'log', normalize = normalize_parameters),
             Parameter('Homo_CUE',        .5, [0, 1], 'linear', normalize = normalize_parameters),
             
             Parameter('Aceto_Km_Ac',    25, [0.5, 500], 'log', normalize = normalize_parameters), #166
             Parameter('Ac_v_max',       .0083, [0.005, 1.], 'log', normalize = normalize_parameters),
             Parameter('Ac_CUE',          .5, [0, 1], 'linear', normalize = normalize_parameters), 
             
             Parameter('Fe3_Km_Ac',      50, [0.0005, 1000], 'log',  normalize = normalize_parameters),#500
             Parameter('Fe3_Km_Fe3',     500, [0.0005, 1000], 'log', normalize = normalize_parameters),
             Parameter('Fe3_v_max',      1.5, [0.2, 3.], 'log', normalize = normalize_parameters), 
             Parameter('Fe3_CUE',        0.5, [0, 1], 'linear', normalize = normalize_parameters),
             
             Parameter('Acetate',         20, [0, 100], 'linear', normalize = normalize_parameters),
             Parameter('Fe3',            150, [0, 300], 'linear', normalize = normalize_parameters),
            
            # biomass in units milligrams per gram dry weight of soil!!!
             Parameter('M_Ferm',         .42, [1e-8, 5],'log',  normalize = normalize_parameters),
             Parameter('M_Hydro',       .083, [1e-8, 5], 'log', normalize = normalize_parameters),
             Parameter('M_Fe3',          .25, [1e-8, 5], 'log',  normalize = normalize_parameters),
             Parameter('M_Homo',         .25, [1e-8, 5], 'log', normalize = normalize_parameters),
             Parameter('M_Ac',         .0033, [1e-8, 5], 'log', normalize = normalize_parameters),
             
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
            used_parameter.set(default_values[used_parameter.name])

        self.model_parameters['CO2'].constant(0)
        self.model_parameters['H2'].constant(0)
        self.model_parameters['CH4'].constant(0)
        self.model_parameters['Fe2'].constant(0)

    def compute_dissolved(self, S):
        dissolved_S = S.copy()
        rho_H2O = CONSTANTS.DENSITY_WATER_4_DEGC # g/L
        m_H2O = CONSTANTS.MOLAR_MASS_H2O # g/mol
        total_gas_volume = self.replica.total_gas_volume # ml
        sample_weight = self.replica.dry_weight # g
        gas_volume = 1e-3 * total_gas_volume/sample_weight # L/gdw

        water_volume = 10**-6 * S[system.index('H2O')]*m_H2O/rho_H2O #  L/gdw
        lam =  self.HCC_s * water_volume/gas_volume
        dissolved_fractions = lam/(1+lam)
        for pool in ['CO2', 'CH4', 'H2']:
            i = system.index(pool)
            dissolved_S[i] = dissolved_fractions[i] * S[i]
        return dissolved_S

    def __call__(self, t, S):
        S = np.maximum(S, 0) 
        
        dissolved_S = self.compute_dissolved(S)
        dSj_dt = np.stack([pathway_j(t, dissolved_S)
                            for pathway_j in self.contributing_pathways], axis = -1)

        dS_dt = np.sum(dSj_dt, axis = -1)
        dS_dt = np.where(S == 0, np.maximum(0, dS_dt), dS_dt)

        #for dSi_dt, pool_name in zip(dS_dt, system.SYSTEM):
        #    self.system_change_log.log_snap(pool_name, t, dSi_dt)
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
                coefficient = pathway.pathway_vector[pool_index]
                if coefficient != 0:
                    self.system_state_log.log(
                        f'{pathway_name}_rate_{pool_name}', t,
                        coefficient * biomass * specific_rate)
            if 'CH4' in pathway.products:
                rates = biomass * (
                    specific_rate * pathway.pathway_vector[:, np.newaxis]
                    - pathway.death_rate[:, np.newaxis])
                self.system_state_log.log(
                    pathway_name + '_CH4 from ' + pathway.microbe.name, t,
                    rates[system.index('CH4')])

    def log_accepted_states(self, t, states):
        """Record all accepted pools and their dissolved gas equivalents."""
        dissolved_states = np.column_stack([
            self.compute_dissolved(state) for state in states.T
        ])
        for pool_index, pool_name in enumerate(system.SYSTEM):
            # Re-log raw pools in the diagnostic path, making the complete-log
            # contract explicit even if a future lean predictor changes.
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
                initial_mean_days = 0):
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
        S0 = system.initial_state(replica, self.parameters(), initial_mean_days)
        self.replica = replica
        self.parameters().check()
        self.system_state_log.reset(replica)
        self.system_change_log.reset(replica)
        gas_limits = gas_upper_limits(replica, S0)

        # A process watchdog is necessary because LSODA offers neither a
        # wall-clock nor a function-evaluation limit.  The ``parallel``
        # argument is retained for API compatibility; every prediction now
        # receives the same five-second protection.
        result_queue = multiprocessing.Queue(maxsize=1)
        try:
            process = multiprocessing.Process(
                target=_integration_worker,
                args=(self, t_eval, S0, result_queue, reset_Fe3, gas_limits))
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
        for Si, pool_name in zip(solver_result.y, system.SYSTEM):
            self.system_state_log.log(pool_name, solver_result.t, Si)

        #self.log_accepted_thermodynamics(solver_result.t, solver_result.y)
        #self.log_accepted_kinetics(solver_result.t, solver_result.y)

        #self.system_state_log._log['CO2_on_measured'] = t_eval, predicted_CO2_on_measured
        #self.system_state_log._log['CH4_on_measured'] = t_eval, predicted_CH4_on_measured
        #self.system_state_log._log['R2'] = {'CO2': co2_r2,
        #                                    'CH4': ch4_r2}

        #add_to_log = []
        #for k, val in self.system_state_log._log.items():
        #    if 'CH4 from' in k:
        #        t, v = val
        #        delta_t = np.diff(t)
        #        integral = delta_t*(v[:-1]+v[1:])*0.5 # trapezoidal rule
        #        n = k + ' (integrated)'
        #        add_to_log.append((n, t[1:], np.cumsum(integral)))

        #for name, ts, vs in add_to_log:
        #    self.system_state_log.log(name, ts, vs)

        self.system_state_log._parameters = self.parameters()

        return self.system_state_log

    def predict_with_logging(self, replica, t = None, t_eval = None, quiet = False,
                             parallel = False, reset_Fe3 = None,
                             days_beyond_reset = 1000, initial_mean_days = 0):
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
            initial_mean_days=initial_mean_days)
        missing_pools = [pool for pool in system.SYSTEM if pool not in run_log.keys()]
        if missing_pools:
            raise RuntimeError(
                'Diagnostic prediction is missing accepted pool trajectories: '
                + ', '.join(missing_pools))
        pool_time, _ = run_log[system.SYSTEM[0]]
        states = np.vstack([run_log[pool][1] for pool in system.SYSTEM])
        self.log_accepted_states(pool_time, states)
        self.log_accepted_thermodynamics(pool_time, states)
        self.log_accepted_kinetics(pool_time, states)
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

