import os
import json
import traceback
import multiprocessing
from concurrent.futures import ProcessPoolExecutor, TimeoutError as FuturesTimeout
import warnings
import numpy as np
import scipy.optimize
import matplotlib.pyplot as plt
import cma
from datetime import datetime
from abc import ABC, abstractmethod

from data import compute_rate
import parameters
import USER_VARIABLES
import checkpoint


# Optimisers require finite objective values.  This penalty is large enough
# that a valid candidate always wins, without risking overflow when losses are
# combined across replicas.
INVALID_OBJECTIVE_PENALTY = 1e100


class NonFiniteLoss(ValueError):
    """A single candidate produced no usable finite loss value."""


class IntegrationFailure(NonFiniteLoss):
    """The model trajectory terminated before it could produce a fit."""


class _IndexedObjective:
    """Return a candidate index with its loss from a pool worker."""

    def __init__(self, objective):
        self.objective = objective

    def __call__(self, indexed_candidate):
        index, candidate = indexed_candidate
        return index, self.objective(candidate)


def r2(predicted, measured, log = False):
    predicted = np.squeeze(predicted)
    measured = np.squeeze(measured)
    if log:
        predicted = np.log(predicted)
        measured = np.log(measured)
    usable = np.logical_and(np.isfinite(predicted), np.isfinite(measured))
    measured_mean = np.nanmean(measured[usable])
    SS_res = np.nansum((predicted[usable] - measured[usable])**2)
    SS_total = np.nansum((measured[usable] - measured_mean)**2)

    r2_value = 1 - SS_res/SS_total
    
    return r2_value


def get(algorithm_name):
    algo_classes = {'differential_evolution': DifferentialEvolution,
            'cma_es': CMAEvolutionStrategy,
            'Powell': Powell,
            'initial_sampling': InitialSampling}
    chosen_class = algo_classes[algorithm_name]
    algo_instance = chosen_class()
    return algo_instance

def mse(true, pred):
    with np.errstate(over='ignore', invalid='ignore'):
        return np.mean((true - pred)**2)

def build_objective_function(pathway_model, replicas, objective_config, t_start = 0, t_end = None):
    import parameters
    if not isinstance(replicas, (list, tuple)):
        replicas = [replicas]
    replica_objectives = []
    for replica in replicas:
        replica_objective = Objective(pathway_model, replica)
        for pool in ['CO2', 'CH4']:
            tf = parameters.IdentityTransform()
            for t in objective_config['transform'][pool]:
                if t == 'normalize':
                    if pool == 'CO2':
                        _,pool_values = replica.CO2()
                    elif pool == 'CH4':
                        _,pool_values = replica.CH4()
                    else:
                        raise NotImplementedError()
                    values = tf.transform(pool_values)
                    finite_values = values[np.isfinite(values)]
                    replica_low = np.min(finite_values)
                    replica_high = np.max(finite_values)
                    tf = parameters.Normalization(tf, replica_low, replica_high)
                elif t == 'log':
                    tf = parameters.LogTransform(tf)
            
            pool_loss = get_loss_function(pool, 
                                        objective_config['reduction'][pool], 
                                        tf,
                                        t_start = t_start,
                                        t_end = t_end)
            replica_objective.add_loss(pool_loss, objective_config['loss_weight'][pool])
        replica_objectives.append(replica_objective)
    total_objective = sum(replica_objectives)

    return total_objective

class Loss():
    def __init__(self, pool, reduction_function, 
                 transform = None, t_start = None, t_end = None):
        self.pool = pool
        self.reduction_function = reduction_function
        self.transform = transform
        self.t_start = t_start
        self.t_end = t_end
        
        self._model = None
        self._replica = None
    
    def set_model(self, model, replica):
        self._model = model
        self._replica = replica
        
        self._model.add_t(self._replica.incubation['days'])
        
    def squared_errors(self, replica, run_log):
        predicted, measured = self.get_values(replica, run_log)
        usable = np.logical_and(np.isfinite(predicted), np.isfinite(measured))
        measured_mean = np.nanmean(measured[usable])
        residual = (predicted[usable] - measured[usable])**2
        total = (measured[usable] - measured_mean)**2
        return residual, total
        
    def R2(self, run_log):
        S_res, S_total = self.squared_errors(self._replica, run_log)
        r2_value = 1 - np.nansum(S_res)/np.nansum(S_total)
        return r2_value

    def t_eval(self, replica):
        if self.pool == 'CO2':
            t, values = replica.CO2()
        elif self.pool == 'CH4':
            t, values = replica.CH4()
        else:
            raise NotImplementedError()

        #if not self.t_end is None:
        #    t = np.unique(np.concatenate([t, np.arange(np.floor(np.max(t)), self.t_end)], axis = 0))
        #if not self.t_start is None:
        #    t = np.unique(np.concatenate([np.arange(0, self.t_start), t], axis = 0))
        #t = np.sort(t)

        usable = np.isfinite(t)
        if callable(self.transform):
            usable &= np.isfinite(self.transform(values))

        usable &= t != 0
        if self.t_start is not None:
            usable &= t >= self.t_start
        if self.t_end is not None:
            usable &= t <= self.t_end

        return t[usable]
    
    def get_values(self, replica, run_log):
        t_pred, pool_pred = run_log[self.pool]
        pool_true = replica[self.pool]
        t_true = replica.incubation['days']
  
        if self.pool == 'CO2':
            t_true, pool_true = replica.CO2()
        elif self.pool == 'CH4':
            t_true, pool_true = replica.CH4()
        else:
            raise NotImplementedError()
            
        assert pool_true.size == t_true.size
        
        # Keep singleton selections one-dimensional.  A narrow replay window
        # can legitimately contain one observation; ``squeeze`` would turn
        # its index into a scalar and make the subsequent time filtering try
        # to iterate over a numpy float.
        ind = np.asarray([
            np.nonzero(np.isclose(t_pred, t))[0].item()
            for t in t_true if np.any(np.isclose(t, t_pred))
        ], dtype=int)
        t_pred = t_pred[ind]
        pool_pred = pool_pred[ind]

        ind = np.asarray([
            np.nonzero(np.isclose(t_true, t))[0].item()
            for t in t_pred if np.any(np.isclose(t, t_true))
        ], dtype=int)

        t_true = t_true[ind]
        pool_true = pool_true[ind]
        
        assert np.all(np.isclose(t_pred, t_true)), str(t_pred) +'\n' + str(t_true)


        if 0 in t_pred:
            pool_pred = pool_pred[t_pred != 0]
            pool_true = pool_true[t_pred != 0]
            t_true = t_true[t_pred != 0]
            t_pred = t_pred[t_pred != 0]
        
        if not self.t_start is None or not self.t_end is None:
            idx = np.arange(len(t_pred))
            if not self.t_start is None:
                idx = np.intersect1d(idx, np.nonzero(t_pred >= self.t_start)[0])
                
            if not self.t_end is None:
                idx = np.intersect1d(idx, np.nonzero(t_pred <= self.t_end)[0])
            
            pool_pred = pool_pred[idx]
            pool_true = pool_true[idx]

        before_pred = np.array(pool_pred)        
        if callable(self.transform):
            pool_pred = self.transform(pool_pred)
            pool_true = self.transform(pool_true)
        
        usable = np.logical_and(np.isfinite(pool_pred), np.isfinite(pool_true))
        
        sz = pool_pred.size
        
        unusable_count = np.count_nonzero(np.logical_not(usable))
        if unusable_count > .5*sz:
            pass
            #print()
            #print('WARNING: less than 50% usable values!', replica, self.pool)
            
        pool_pred = pool_pred[usable]
        pool_true = pool_true[usable]
        
        return pool_pred, pool_true
    
    def __str__(self):
        tf = self.pool
        if callable(self.transform):
            l, r = self.transform.operator()
            tf = l + tf + r
        f = str(self.reduction_function.__name__) + '(' + tf + ')'
        return f
    
    def __call__(self, replica, run_log):
        pool_pred, pool_true = self.get_values(replica, run_log)
        if pool_pred.size == 0 or pool_true.size == 0:
            raise NonFiniteLoss(f'{self} has no usable values')
        if not np.all(np.isfinite(pool_pred)):
            raise NonFiniteLoss(f'{self} has non-finite predicted values')
        if not np.all(np.isfinite(pool_true)):
            raise NonFiniteLoss(f'{self} has non-finite measured values')
        loss_value = self.reduction_function(pool_true, pool_pred)
        if not np.isfinite(loss_value):
            raise NonFiniteLoss(f'{self} returned {loss_value}')
        return loss_value

loss_functions = {'mse': mse}

def get_loss_function(pool, reduction = 'mse', transform = None, 
                  t_start = None, t_end = None):
    return Loss(pool, loss_functions[reduction], 
                transform = transform, 
                t_start = t_start, t_end = t_end)

class Addable():
    def __init__(self, left = None, right = None):
        self._left = left
        self._right = right
        
        self._callbacks = []
        self._call_log = []
        self._global_call_counter = None
        self._global_best_loss = None
        self._global_best_r2 = None
        self._global_call_lock = None
        self._final_call_count = None
        self._final_best_loss = None
        self._final_best_r2 = None
        self._best_local_loss = None
        self._best_local_r2 = None
        self._generation = None

    def _record_evaluation_result(self, value):
        """Allow leaf objectives to account for completed evaluations."""

    def loss_contributions(self):
        return self._left.loss_contributions() + self._right.loss_contributions()
        
    def model(self):
        if hasattr(self, '_model'):
            return self._model
        return self._left.model()
    
    def __radd__(self, other):
        if other == 0:
            return self
        return self.__add__(other)
    
    def __add__(self, other):
        assert other._model is self.model()
        obj = Addable(self, other)
        return obj

    def weighted_mse(self, run_log):
        mse = [objective.weighted_mse(run_log)
               for objective in (self._left, self._right)]
        if any(contribution is None for contribution, _ in mse):
            return None, None
        return tuple(sum(values) for values in zip(*mse))

    def R2(self, run_log):
        mse_res, mse_tot = self.weighted_mse(run_log)
        return None if mse_res is None else 1 - mse_res/mse_tot

    def last_weighted_mse(self):
        """Return residual and total variance from the most recent evaluation."""
        if self._left is None or self._right is None:
            return None, None
        mse = [objective.last_weighted_mse()
               for objective in (self._left, self._right)]
        if any(contribution is None or any(value is None for value in contribution)
               for contribution in mse):
            return None, None
        return tuple(sum(values) for values in zip(*mse))

    def last_R2(self):
        """Return R² for the exact model evaluations that produced last loss."""
        mse_res, mse_tot = self.last_weighted_mse()
        return None if mse_res is None or mse_tot is None else 1 - mse_res/mse_tot

    def fit_replicas(self):
        return list({self._left.fit_replicas(), self._right.fit_replicas()})

    def _call(self, args, **kwargs):
        return self._left(args, **kwargs) + self._right(args, **kwargs)

    def set_global_call_counter(self, counter, lock, best_loss = None, best_r2 = None):
        self._global_call_counter = counter
        self._global_best_loss = best_loss
        self._global_best_r2 = best_r2
        self._global_call_lock = lock
        self._final_call_count = None
        self._final_best_loss = None
        self._final_best_r2 = None

    def set_generation(self, generation):
        """Expose the optimiser generation to callbacks.

        ``generation`` may be a multiprocessing manager value, allowing
        callbacks running in worker processes to see the current generation.
        """
        self._generation = generation

    def generation(self):
        if self._generation is None:
            return None
        return getattr(self._generation, 'value', self._generation)

    def finalize_global_call_counter(self):
        with self._global_call_lock:
            self._final_call_count = self._global_call_counter.value
            if self._global_best_loss is not None:
                self._final_best_loss = self._global_best_loss.value
            if self._global_best_r2 is not None:
                self._final_best_r2 = self._global_best_r2.value
        self._global_call_counter = None
        self._global_best_loss = None
        self._global_best_r2 = None
        self._global_call_lock = None
    
    def __call__(self, args, **kwargs):
        if self._global_call_counter is not None:
            with self._global_call_lock:
                self._global_call_counter.value += 1
        value = self._call(args, **kwargs)
        self._record_evaluation_result(value)
        self._call_log.append((args, kwargs, value))
        last_r2 = self.last_R2()
        if self._best_local_loss is None or value < self._best_local_loss:
            self._best_local_loss = value
            self._best_local_r2 = last_r2

        # Each worker has its own call log.  Publish the minimum before
        # callbacks run so they all report the same best loss/R² pair.
        if self._global_best_loss is not None:
            with self._global_call_lock:
                if value < self._global_best_loss.value:
                    self._global_best_loss.value = value
                    if self._global_best_r2 is not None:
                        self._global_best_r2.value = (
                            float('nan') if last_r2 is None else float(last_r2))

        for callback in self._callbacks:
            callback()
            
        return value
    
    def call_count(self):
        if self._global_call_counter is not None:
            with self._global_call_lock:
                return self._global_call_counter.value
        if self._final_call_count is not None:
            return self._final_call_count
        return len(self._call_log)
    
    def last_call(self):
        return self._call_log[-1]
    
    def best_call(self):
        return sorted(self._call_log, key = lambda k: k[-1])[0]

    def best_loss(self):
        if self._global_best_loss is not None:
            with self._global_call_lock:
                return self._global_best_loss.value
        if self._final_best_loss is not None:
            return self._final_best_loss
        return self.best_call()[-1]

    def best_R2(self):
        if self._global_best_r2 is not None:
            with self._global_call_lock:
                value = self._global_best_r2.value
            return None if np.isnan(value) else value
        if self._final_best_r2 is not None:
            return None if np.isnan(self._final_best_r2) else self._final_best_r2
        return self._best_local_r2
    
    def add_callback(self, callback):
        assert callable(callback)
        callback.set_objective(self)
        self._callbacks.append(callback)

    def set_parameters(self, parameters, values, transformed = False):
        """Delegate parameter updates through a combined objective tree.

        A split fit is represented by ``Addable`` rather than an ``Objective``.
        InitialSampling selects its final candidate through this method, so the
        aggregate needs the same parameter-setting interface as a leaf.
        """
        if self._left is None:
            raise AttributeError('objective has no parameter receiver')
        self._left.set_parameters(parameters, values, transformed)

    def __str__(self):
        return ' + '.join([str(s) for s in [self._left, self._right]
                          if not s is None])
    
class Objective(Addable):
    def __init__(self, model, replica, replica_weight = 1.0):
        super().__init__()
        self._model = model

        self._replica = replica
        self._replica_weight = replica_weight
        
        self._loss_contributions = []
        self._last_weighted_mse = None
        self._penalty_count = 0
        self._penalty_counter = None
        self._last_evaluation_was_penalized = False
        self._suppress_penalty_warnings = False

    def set_penalty_warning_output(self, enabled):
        """Enable or suppress one warning per invalid objective evaluation."""
        self._suppress_penalty_warnings = not enabled

    def penalty_count(self):
        """Return invalid evaluations in this objective's current fit."""
        return (self._penalty_count if self._penalty_counter is None
                else self._penalty_counter.value())

    def evaluation_count(self):
        """Return completed evaluations, including those in worker copies."""
        return (self.call_count() if self._penalty_counter is None
                else self._penalty_counter.call_count())

    def set_penalty_counter(self, counter):
        """Use a process-safe counter shared by parallel worker copies."""
        self._penalty_counter = counter

    def _record_evaluation_result(self, value):
        if self._last_evaluation_was_penalized:
            self._penalty_count += 1
        if self._penalty_counter is not None:
            self._penalty_counter.record(self._last_evaluation_was_penalized)
            
    def loss_contributions(self):
        l = self._loss_contributions
        for _l, _ in l:
            _l._replica = self._replica 
        return l

    def fit_replicas(self):
        return self._replica

    def variables(self):
        return self.model().parameters().variables()

    def add_loss(self, loss_function, loss_weight):
        assert callable(loss_function)
        self.loss_contributions().append((loss_function, loss_weight))

    def t_eval(self):
        return np.unique(np.concatenate([
            loss_function.t_eval(self._replica)
            for loss_function, _ in self.loss_contributions()
        ]))

    def weighted_mse(self, run_log):
        weighted_model_mse = 0
        weighted_total_mse = 0
        for loss_function, loss_weight in self.loss_contributions():
            s_residual, s_total = loss_function.squared_errors(self._replica, run_log)
            weighted_model_mse += loss_weight*np.nanmean(s_residual)
            weighted_total_mse += loss_weight*np.nanmean(s_total)
        return weighted_model_mse, weighted_total_mse

    def set_parameters(self, parameters, values, transformed = False):
        if transformed:
            values = [var.inverse_transform(p) 
                            for var, p in zip(parameters, values)]

        _ = [v.set(p) for v, p in zip(parameters, np.atleast_1d(values))]
        
        
    def _call(self, transformed_parameters, transformed = True, t_eval = None):
        self._last_evaluation_was_penalized = False
        if isinstance(transformed_parameters, parameters.ModelParameters):
            p_dict = transformed_parameters.as_dict()
            p_values = list(p_dict.values())
            pars = self.model().parameters()
            assert all([p in pars for p in transformed_parameters])
            self.set_parameters([pars[n] for n in list(p_dict.keys())], 
                                [p.value for p in p_values],
                                transformed = transformed)
        else:
            self.set_parameters(self.variables(), transformed_parameters, 
                                transformed = transformed)

        try:
            compute_loss = t_eval is None
            t_eval = self.t_eval() if t_eval is None else t_eval
            run_log = self._model.predict(self._replica, t_eval = t_eval)

            replica_loss = 0
            if compute_loss:
                for loss_function, loss_weight in self.loss_contributions():
                    loss_value = loss_function(self._replica, run_log)
                    replica_loss += loss_weight*loss_value
        except NonFiniteLoss as error:
            self._last_evaluation_was_penalized = True
            if not self._suppress_penalty_warnings:
                print(f'\nWARNING: Replica {self._replica}: {error}; '
                      f'using penalty {INVALID_OBJECTIVE_PENALTY:g}.')
            self._last_weighted_mse = None, None
            return self._replica_weight * INVALID_OBJECTIVE_PENALTY
        except Exception as ex:
            ex.args = ('Replica ' + str(self._replica) + ': ' + ex.args[0], ) + ex.args[1:]
            raise

        # ``predict`` reuses the model's mutable run log.  Retain the derived
        # metrics now, while they still belong to this replica evaluation.
        self._last_weighted_mse = self.weighted_mse(run_log)
        
        objective_value = self._replica_weight * replica_loss
        return objective_value

    def last_weighted_mse(self):
        return self._last_weighted_mse
   

    def __str__(self):
        def format_weight(w):
            return '' if w == 1 else f'{w:.2g} * '
        
        def format_subscript(s):
            s = str(s)
            subs = {'0': '₀',
                    '1': '₁',
                    '2': '₂',
                    '3': '₃',
                    '4': '₄',
                    '5': '₅', 
                    '6': '₆', 
                    '7': '₇', 
                    '8': '₈',
                    '9': '₉', 
                    '-': '₋'}
            return ''.join([subs[i] if i in subs else i for i in s ])
            
            
        s_weight = format_weight(self._replica_weight)
        s_loss = ' + '.join([format_weight(w) + str(l) for l,w in self.loss_contributions()])
        replica = str(self._replica.sample.sample_name) + '-' + str(self._replica.replica_number)
        s_replica = format_subscript(replica)
        
        return s_weight + '[ '+ s_loss + ' ]' + s_replica

class Algorithm(ABC):
    def __init__(self, **default_kwargs):
        self._kwargs = dict(default_kwargs)
        self._valid_kwargs = frozenset(default_kwargs)

    def configure(self, kwargs):
        """Apply supported algorithm settings and warn about ignored keys.

        Callers commonly retain descriptive metadata alongside an algorithm
        configuration (for example, a staged-fit ``name``).  Passing that
        metadata through to SciPy produces an unhelpful ``unexpected keyword``
        error.  Each algorithm declares its supported settings through the
        defaults supplied to this base class, so only those keys are accepted.
        """
        if not isinstance(kwargs, dict):
            raise TypeError('Algorithm configuration must be a dictionary.')

        unused = sorted(set(kwargs) - self._valid_kwargs)
        if unused:
            warnings.warn(
                f'{self.__class__.__name__} ignored unsupported '
                f'configuration key(s): {", ".join(unused)}',
                RuntimeWarning, stacklevel=2)
        self._kwargs.update({
            key: value for key, value in kwargs.items()
            if key in self._valid_kwargs
        })
    
    def get_bounds(self, variables):
        lower_bounds = np.reshape([v.transform(v.lower()) for v in variables], (-1,))
        upper_bounds = np.reshape([v.transform(v.upper()) for v in variables], (-1,))
        bounds = list(zip(lower_bounds, upper_bounds))
        return bounds

    def minimize(self, objective, initial_parameters, verbose = True):
        objective.model().parameters().set(initial_parameters)
        if verbose:
            print('Fitting with ' + self.__class__.__name__)
            print(objective.model())
            print()
            print('search space average fraction:',
                  objective.model().parameters().search_space()[-1])
            print()
            print('Objective = ' + str(objective))
            print()
            print('Callbacks:')
            for cb in objective._callbacks:
                pth = cb.target_directory()
                if not pth is None:
                    print('\n'.join(set([cb.__class__.__name__ + ': ' + os.path.join(*pth.split(os.sep)[-2:])
                                    for cb in objective._callbacks])))
                else:
                    print(cb.__class__.__name__)

            print()
        if len(objective.model().parameters().variables()) == 0:
            raise Exception('Model has no variables. Check initial parameters.')
        self._minimize(objective, initial_parameters)
    
    @abstractmethod
    def _minimize(self, objective, initial_parameters):
        raise NotImplementedError()
    
    def get_config(self):
        return self._kwargs

class DifferentialEvolution(Algorithm):    
    def __init__(self):
        defaults = {'strategy': 'rand1bin',
                        'updating': 'deferred',#'immediate',
                        'popsize': 10,
                        'workers': -1,
                        'tol': 1e-4,
                        'init': 'sobol',
                        'polish': False,
                        'recombination': .7, # CR
                        'mutation': (.5,1.), # F
                        'seed': None,
                        # Optional transformed incumbent supplied by a
                        # checkpoint-resume caller.
                        'x0': None,
                        # SciPy has no per-evaluation budget.  The generation
                        # callback below honours this cap after each completed
                        # generation, which keeps staged restarts bounded.
                        'max_evaluations': None,
                        }
        super().__init__(**defaults)
        self.generation = 0
        self._process_manager = None
        self._worker_map = None

    def set_process_manager(self, manager):
        """Reuse a caller-owned manager for repeated parallel DE runs.

        Nested optimisation may invoke many short DE runs.  Starting a fresh
        manager for every one leaks file descriptors on some multiprocessing
        backends before the operating system can reclaim them.
        """
        self._process_manager = manager

    def set_worker_map(self, worker_map):
        """Use a caller-owned pool's ``map`` method for SciPy evaluation."""
        self._worker_map = worker_map
        
    def _minimize(self, objective, initial_parameters):
        self.generation = 1
        self._objective = objective
        self._shared_generation = None
        self._max_evaluations = self._kwargs.get('max_evaluations')
        if self._max_evaluations is not None:
            if (not isinstance(self._max_evaluations, (int, np.integer))
                    or self._max_evaluations <= 0):
                raise ValueError('max_evaluations must be a positive integer or None')
        variables = objective.model().parameters().variables()
        manager = None
        owns_manager = False
        # A caller-owned pool serialises the objective anew for every map
        # request.  Manager lock proxies retain a per-process connection and
        # can become stale across those repeated requests on Python 3.12.
        # The counter is only progress reporting, so omit it for that stable
        # reusable-pool path.  Checkpoints and stop callbacks remain local to
        # each worker and do not need the manager.
        if self._kwargs['workers'] != 1 and self._worker_map is None:
            manager = self._process_manager
            if manager is None:
                manager = multiprocessing.Manager()
                owns_manager = True
            objective.set_global_call_counter(
                manager.Value('i', 0), manager.Lock(),
                manager.Value('d', float('inf')),
                manager.Value('d', float('nan')))
            self._shared_generation = manager.Value('i', self.generation)
            objective.set_generation(self._shared_generation)
        else:
            objective.set_generation(self.generation)

        try:
            # ``max_evaluations`` is a repository-level control, not a SciPy
            # differential_evolution keyword.
            scipy_kwargs = dict(self._kwargs)
            del scipy_kwargs['max_evaluations']
            if self._worker_map is not None:
                scipy_kwargs['workers'] = self._worker_map
            _ = scipy.optimize.differential_evolution(objective,
                                                      bounds = self.get_bounds(variables),
                                                      callback = self.generation_counter,
                                                      **scipy_kwargs
                                                      )
        #except KeyboardInterrupt:
        #    return

        finally:
            if manager is not None:
                objective.finalize_global_call_counter()
                if owns_manager:
                    manager.shutdown()
        
    def generation_counter(self, *args, **kwargs):
        self.generation += 1
        if self._shared_generation is None:
            self._objective.set_generation(self.generation)
        else:
            self._shared_generation.value = self.generation
        return (self._max_evaluations is not None
                and self._objective.call_count() >= self._max_evaluations)

class CMAEvolutionStrategy(Algorithm):
    """Bounded CMA-ES backed by the maintained :mod:`pycma` package."""
    def __init__(self):
        super().__init__(popsize=None, sigma0=0.30, restarts=0,
                         restart_multiplier=2, seed=None, workers=1,
                         max_evaluations=None, x0=None, tolx=1e-8)
        self._process_manager = None
        self._worker_map = None
        self._state_checkpoint = None

    def set_process_manager(self, manager):
        self._process_manager = manager

    def set_worker_map(self, worker_map):
        self._worker_map = worker_map

    def set_state_checkpoint(self, state_checkpoint):
        """Use ``state_checkpoint`` to resume complete CMA-ES generations.

        The object intentionally has a tiny protocol (``load`` and ``save``)
        so the optimiser remains independent of the on-disk checkpoint
        implementation.  ``None`` disables state persistence.
        """
        self._state_checkpoint = state_checkpoint

    def _minimize(self, objective, initial_parameters):
        variables = objective.model().parameters().variables()
        dimension = len(variables)
        bounds = np.asarray(self.get_bounds(variables), dtype=float)
        widths = bounds[:, 1] - bounds[:, 0]
        if np.any(~np.isfinite(widths)) or np.any(widths <= 0):
            raise ValueError('CMA-ES requires finite, non-zero parameter bounds')

        max_evaluations = self._kwargs['max_evaluations']
        if max_evaluations is not None:
            if (not isinstance(max_evaluations, (int, np.integer))
                    or max_evaluations <= 0):
                raise ValueError('max_evaluations must be a positive integer')
        restarts = self._kwargs['restarts']
        if not isinstance(restarts, (int, np.integer)) or restarts < 0:
            raise ValueError('restarts must be a non-negative integer')
        restart_multiplier = self._kwargs['restart_multiplier']
        if restart_multiplier < 1:
            raise ValueError('restart_multiplier must be at least 1')
        sigma0 = float(self._kwargs['sigma0'])
        if not np.isfinite(sigma0) or sigma0 <= 0:
            raise ValueError('sigma0 must be a finite positive number')
        tolx = float(self._kwargs['tolx'])
        if not np.isfinite(tolx) or tolx <= 0:
            raise ValueError('tolx must be a finite positive number')
        workers = self._kwargs['workers']
        if not isinstance(workers, (int, np.integer)) or workers == 0 or workers < -1:
            raise ValueError('CMA-ES workers must be 1, -1, or a positive integer')

        configured_popsize = self._kwargs['popsize']
        if configured_popsize is not None and (
                not isinstance(configured_popsize, (int, np.integer))
                or configured_popsize <= 0):
            raise ValueError('popsize must be a positive integer')
        base_popsize = (int(configured_popsize) if configured_popsize is not None
                        else 4 + int(3 * np.log(dimension)))

        transformed_initial = np.asarray(
            [variable.transform(variable.value) for variable in variables], dtype=float)
        x0 = self._kwargs['x0']
        if x0 is not None:
            x0 = np.asarray(x0, dtype=float).reshape(-1)
            if x0.size != dimension or not np.all(np.isfinite(x0)):
                raise ValueError('x0 must contain one finite transformed value per variable')
            transformed_initial = x0
        mean0 = np.clip((transformed_initial - bounds[:, 0]) / widths, 0.0, 1.0)

        manager = None
        owns_manager = False
        # See DifferentialEvolution: do not send manager proxies through a
        # reusable worker pool on every inner optimisation.
        if workers != 1 and self._worker_map is None:
            manager = self._process_manager
            if manager is None:
                manager = multiprocessing.Manager()
                owns_manager = True
            objective.set_global_call_counter(
                manager.Value('i', 0), manager.Lock(),
                manager.Value('d', float('inf')),
                manager.Value('d', float('nan')))

        best_loss = float('inf')
        best_x = None
        evaluations = 0
        generation = 0
        start_restart = 0
        resumed_strategy = None
        pending_population = None
        pending_losses = None
        pending_restart = None
        if self._state_checkpoint is not None:
            saved_state = self._state_checkpoint.load()
            if isinstance(saved_state, dict):
                expected_bounds = bounds.tolist()
                if (saved_state.get('version') == 1
                        and saved_state.get('dimension') == dimension
                        and saved_state.get('bounds') == expected_bounds
                        and saved_state.get('base_popsize') == base_popsize
                        and saved_state.get('restart_multiplier') == restart_multiplier):
                    resumed_strategy = saved_state.get('strategy')
                    start_restart = saved_state.get('restart', 0)
                    evaluations = saved_state.get('evaluations', 0)
                    generation = saved_state.get('generation', 0)
                    best_loss = saved_state.get('best_loss', best_loss)
                    best_x = saved_state.get('best_x')
                    pending_population = saved_state.get('pending_population')
                    pending_losses = saved_state.get('pending_losses')
                    pending_restart = saved_state.get('pending_restart')
                    # pycma's default sampler uses NumPy's legacy global RNG.
                    # Restoring it makes the next ``ask`` identical to an
                    # uninterrupted run, rather than merely warm-started.
                    random_state = saved_state.get('numpy_random_state')
                    if random_state is not None:
                        np.random.set_state(random_state)
                    print('Resuming CMA-ES state at generation '
                          + str(generation) + ', evaluation ' + str(evaluations))

        def save_state(strategy, restart, pending=None, losses=None):
            if self._state_checkpoint is None:
                return
            self._state_checkpoint.save({
                'version': 1,
                'dimension': dimension,
                'bounds': bounds.tolist(),
                'base_popsize': base_popsize,
                'restart_multiplier': restart_multiplier,
                'restart': restart,
                'evaluations': evaluations,
                'generation': generation,
                'best_loss': best_loss,
                'best_x': best_x,
                'numpy_random_state': np.random.get_state(),
                'strategy': strategy,
                # ``strategy`` is saved after ``ask`` but before ``tell``
                # when a population is in flight.  Its RNG has already
                # advanced, so retaining the sampled points is essential:
                # asking again would produce a different generation.
                'pending_population': pending,
                'pending_losses': losses,
                'pending_restart': (restart if pending is not None else None),
            })

        try:
            for restart in range(start_restart, restarts + 1):
                if max_evaluations is not None and evaluations >= max_evaluations:
                    break
                population_size = int(base_popsize * (restart_multiplier ** restart))
                options = {
                    # pycma performs its own boundary transformation instead
                    # of piling out-of-bounds samples onto box corners.
                    'bounds': [0.0, 1.0],
                    'popsize': population_size,
                    'tolx': tolx,
                    'verbose': -9,
                    'verb_disp': 0,
                    'verb_log': 0,
                }
                if self._kwargs['seed'] is not None and resumed_strategy is None:
                    options['seed'] = int(self._kwargs['seed']) + restart
                mean = best_x if best_x is not None else mean0
                if resumed_strategy is not None:
                    strategy = resumed_strategy
                    resumed_strategy = None
                else:
                    strategy = cma.CMAEvolutionStrategy(mean, sigma0, options)

                while not strategy.stop():
                    continuing_pending = (
                        pending_population is not None and pending_restart == restart)
                    if (not continuing_pending and max_evaluations is not None
                            and evaluations >= max_evaluations):
                        break
                    if not continuing_pending:
                        generation += 1
                    if hasattr(objective, 'set_generation'):
                        objective.set_generation(generation)
                    if continuing_pending:
                        population = np.asarray(pending_population)
                        losses = list(pending_losses)
                        pending_population = None
                        pending_losses = None
                        pending_restart = None
                    else:
                        remaining = (max_evaluations - evaluations
                                     if max_evaluations is not None else population_size)
                        batch_size = min(population_size, remaining)
                        # ``ask`` returns bounded points in the unit box.  A
                        # partial final generation is evaluated but not fed
                        # back: pycma correctly adapts only complete
                        # populations.
                        population = np.asarray(strategy.ask())[:batch_size]
                        losses = [None] * batch_size
                        # This makes even an interruption during the first
                        # population resumable.  In a parallel map the whole
                        # unfinished population is retried; serial evaluation
                        # also retains each completed candidate below.
                        if batch_size == population_size:
                            save_state(strategy, restart, population, losses)
                    batch_size = len(population)
                    candidates = bounds[:, 0] + population * widths
                    if workers == 1:
                        for index, candidate in enumerate(candidates):
                            if losses[index] is not None:
                                continue
                            loss = objective(candidate)
                            losses[index] = loss
                            evaluations += 1
                            if np.isfinite(loss) and loss < best_loss:
                                best_loss = float(loss)
                                best_x = population[index].copy()
                            if batch_size == population_size:
                                save_state(strategy, restart, population, losses)
                    elif self._worker_map is not None:
                        missing = [index for index, loss in enumerate(losses)
                                   if loss is None]
                        mapped_losses = self._worker_map(
                            objective, [candidates[index] for index in missing])
                        for index, loss in zip(missing, mapped_losses):
                            losses[index] = loss
                        evaluations += len(missing)
                        if batch_size == population_size:
                            save_state(strategy, restart, population, losses)
                    else:
                        # ProcessPoolExecutor workers are non-daemonic, so a
                        # model evaluation can run its short-lived integration
                        # watchdog process.
                        with ProcessPoolExecutor(
                                max_workers=None if workers == -1 else workers) as pool:
                            missing = [index for index, loss in enumerate(losses)
                                       if loss is None]
                            mapped_losses = pool.map(
                                objective, [candidates[index] for index in missing])
                        for index, loss in zip(missing, mapped_losses):
                            losses[index] = loss
                        evaluations += len(missing)
                        if batch_size == population_size:
                            save_state(strategy, restart, population, losses)
                    losses = np.where(np.isfinite(losses), losses, INVALID_OBJECTIVE_PENALTY)
                    chosen = int(np.argmin(losses))
                    if losses[chosen] < best_loss:
                        best_loss = float(losses[chosen])
                        best_x = population[chosen].copy()
                    if batch_size == population_size:
                        strategy.tell(population, losses.tolist())
                        # Objective callbacks only persist model results.  A
                        # CMA state is persisted here, after ``tell`` updates
                        # its distribution for the whole population.
                        save_state(strategy, restart)
                    else:
                        # A partial final generation is intentionally not
                        # told to pycma.  Keep the preceding complete state
                        # checkpoint: a distribution after ``ask`` but before
                        # ``tell`` cannot safely be resumed without a full
                        # population of losses.
                        break
        finally:
            if manager is not None:
                objective.finalize_global_call_counter()
                if owns_manager:
                    manager.shutdown()

        if best_x is not None:
            best_transformed = bounds[:, 0] + best_x * widths
            objective.set_parameters(variables, best_transformed, transformed=True)


class InitialSampling(Algorithm):
    """Evaluate a space-filling initial population without evolving it.

    This is deliberately separate from :class:`DifferentialEvolution`: it
    produces candidate points using the same families of initial designs, but
    never calls SciPy's differential-evolution solver or creates a generation
    of trial vectors.  It is useful when the first population already gives
    the desired fit, and for staged space-narrowing runs based on checkpoints.
    """

    _SAMPLERS = ('sobol', 'latinhypercube', 'halton', 'random')

    def __init__(self):
        super().__init__(init='sobol', samples=128, seed=None, workers=1,
                         stall_timeout=None)
        self.best_x = None
        self.best_fun = None

    def _unit_population(self, dimension, samples):
        method = self._kwargs['init']
        seed = self._kwargs.get('seed')

        if method == 'sobol':
            # ``random_base2`` preserves Sobol's balance property.  Trimming
            # the final block lets callers request an exact evaluation budget.
            sampler = scipy.stats.qmc.Sobol(d=dimension, scramble=True, seed=seed)
            return sampler.random_base2(int(np.ceil(np.log2(samples))))[:samples]
        if method == 'latinhypercube':
            return scipy.stats.qmc.LatinHypercube(d=dimension, seed=seed).random(samples)
        if method == 'halton':
            return scipy.stats.qmc.Halton(d=dimension, scramble=True, seed=seed).random(samples)
        if method == 'random':
            return np.random.default_rng(seed).random((samples, dimension))
        raise ValueError(
            "initial_sampling init must be one of "
            + ', '.join(self._SAMPLERS) + f'; got {method!r}')

    def _minimize(self, objective, initial_parameters):
        samples = self._kwargs['samples']
        if not isinstance(samples, (int, np.integer)) or samples <= 0:
            raise ValueError('initial_sampling samples must be a positive integer')

        variables = objective.model().parameters().variables()
        if hasattr(objective, 'set_generation'):
            objective.set_generation(1)
        bounds = np.asarray(self.get_bounds(variables), dtype=float)
        unit_population = self._unit_population(len(variables), int(samples))
        population = bounds[:, 0] + unit_population * (bounds[:, 1] - bounds[:, 0])

        # Evaluating candidates directly (rather than using DE with maxiter=0)
        # means the call count is exactly ``samples`` and no mutation or
        # crossover is performed.
        workers = self._kwargs['workers']
        stall_timeout = self._kwargs['stall_timeout']
        if stall_timeout is not None and stall_timeout <= 0:
            raise ValueError('initial_sampling stall_timeout must be positive or None')
        if workers == 1 and stall_timeout is None:
            losses = np.asarray([objective(candidate) for candidate in population])
        else:
            if not isinstance(workers, (int, np.integer)) or workers == 0 or workers < -1:
                raise ValueError('initial_sampling workers must be 1, -1, or a positive integer')
            manager = multiprocessing.Manager()
            objective.set_global_call_counter(
                manager.Value('i', 0), manager.Lock(),
                manager.Value('d', float('inf')),
                manager.Value('d', float('nan')))
            # Unlike multiprocessing.Pool, these workers may start the
            # per-integration watchdog process used by Model.predict.
            pool = ProcessPoolExecutor(
                max_workers=None if workers == -1 else workers)
            losses = np.empty(len(population), dtype=float)
            try:
                indexed_objective = _IndexedObjective(objective)
                futures = [pool.submit(indexed_objective, indexed_candidate)
                           for indexed_candidate in enumerate(population)]
                for future in futures:
                    try:
                        if stall_timeout is None:
                            index, loss = future.result()
                        else:
                            index, loss = future.result(timeout=stall_timeout)
                    except FuturesTimeout as error:
                        raise checkpoint.StalledEvaluationInterrupt(
                            stall_timeout) from error
                    losses[index] = loss
            except BaseException:
                pool.shutdown(wait=False, cancel_futures=True)
                raise
            else:
                pool.shutdown()
            finally:
                objective.finalize_global_call_counter()
                manager.shutdown()
        best_index = int(np.argmin(losses))
        self.best_x = population[best_index].copy()
        self.best_fun = float(losses[best_index])

        # Leave the model in the selected state without another objective call:
        # callbacks have already checkpointed every evaluated candidate.
        objective.set_parameters(variables, self.best_x, transformed=True)


class Powell(Algorithm):
    def _minimize(self, objective, initial_parameters):
        variables = objective.model().parameters().variables()
        x0 = np.reshape([v.transform(v.value) for v in variables], (-1,))
        _ = scipy.optimize.minimize(objective, x0, method = 'Powell', 
                                    bounds = self.get_bounds(variables),
                                    **self._kwargs)
