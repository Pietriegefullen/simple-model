import os
import json
import traceback
import multiprocessing
import numpy as np
import scipy.optimize
import matplotlib.pyplot as plt
from datetime import datetime
from abc import ABC, abstractmethod

from data import compute_rate
import parameters
import USER_VARIABLES


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
            'Powell': Powell}
    chosen_class = algo_classes[algorithm_name]
    algo_instance = chosen_class()
    return algo_instance

def mse(true, pred):
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
        
        ind = np.array([np.nonzero(np.isclose(t_pred,t))[0].item() for t in t_true
                        if np.any(np.isclose(t,t_pred))]).squeeze()
        t_pred = t_pred[ind]
        pool_pred = pool_pred[ind]

        ind = np.array([np.nonzero(np.isclose(t_true,t))[0].item()
                        for t in t_pred if np.any(np.isclose(t,t_true))]).squeeze()

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
            print()
            print('WARNING: less than 50% usable values!', replica, self.pool)
            
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
        if not np.all(np.isfinite(pool_pred)):
            raise Exception('Non-finite pred values')
        if not np.all(np.isfinite(pool_true)):
            raise Exception('Non-finite true values')
        loss_value = self.reduction_function(pool_true, pool_pred)
        if not np.isfinite(loss_value):
            raise Exception(str(self) + ' returned ' + str(loss_value))
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
        if any(contribution is None for contribution, _ in mse):
            return None, None
        return tuple(sum(values) for values in zip(*mse))

    def last_R2(self):
        """Return R² for the exact model evaluations that produced last loss."""
        mse_res, mse_tot = self.last_weighted_mse()
        return None if mse_res is None else 1 - mse_res/mse_tot

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
        
        
    def _call(self, transformed_parameters, transformed = True):
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

        run_log = self._model.predict(self._replica, t_eval = self.t_eval())
        
        replica_loss = 0
        for loss_function, loss_weight in self.loss_contributions():
            try:
                loss_value = loss_function(self._replica, run_log)
            except Exception as ex:
                ex.args = ('Replica ' + str(self._replica) + ': ' + ex.args[0], ) + ex.args[1:]
                raise 
            replica_loss += loss_weight*loss_value

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
        self._kwargs = default_kwargs

    def configure(self, kwargs):
        self._kwargs.update(kwargs)
    
    def get_bounds(self, variables):
        lower_bounds = np.reshape([v.transform(v.lower()) for v in variables], (-1,))
        upper_bounds = np.reshape([v.transform(v.upper()) for v in variables], (-1,))
        bounds = list(zip(lower_bounds, upper_bounds))
        return bounds

    def minimize(self, objective, initial_parameters):
        objective.model().parameters().set(initial_parameters)
        print('Fitting with ' + self.__class__.__name__)
        print(objective.model())
        print()
        print('search space average fraction:', 
              objective.model().parameters().search_space()[-1])
        print()
        print('Objective = ' + str(objective))
        print()
        for cb in objective._callbacks:
            print()
        print('\n'.join(set([os.path.join(*cb.target_directory().split(os.sep)[-2:])
                        for cb in objective._callbacks])))
        
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
                        'mutation': (.5,1.)  # F
                        }
        super().__init__(**defaults)
        self.generation = 0
        
    def _minimize(self, objective, initial_parameters):
        self.generation = 1
        variables = objective.model().parameters().variables()
        manager = None
        if self._kwargs['workers'] != 1:
            manager = multiprocessing.Manager()
            objective.set_global_call_counter(
                manager.Value('i', 0), manager.Lock(),
                manager.Value('d', float('inf')),
                manager.Value('d', float('nan')))

        try:
            _ = scipy.optimize.differential_evolution(objective,
                                                      bounds = self.get_bounds(variables),
                                                      callback = self.generation_counter,
                                                      **self._kwargs
                                                      )
        except KeyboardInterrupt:
            return
        finally:
            if manager is not None:
                objective.finalize_global_call_counter()
                manager.shutdown()
        
    def generation_counter(self, *args, **kwargs):
        self.generation += 1
        return False


class Powell(Algorithm):
    def _minimize(self, objective, initial_parameters):
        variables = objective.model().parameters().variables()
        x0 = np.reshape([v.transform(v.value) for v in variables], (-1,))
        _ = scipy.optimize.minimize(objective, x0, method = 'Powell', 
                                    bounds = self.get_bounds(variables),
                                    **self._kwargs)
