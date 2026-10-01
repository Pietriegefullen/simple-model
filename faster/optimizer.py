import os
import json
import traceback
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
    return np.sqrt(np.sum((true - pred)**2))

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
        
    def R2(self, replica, run_log):
        predicted, measured = self.get_values()
        usable = np.logical_and(np.isfinite(predicted), np.isfinite(measured))
        measured_mean = np.nanmean(measured[usable])
        SS_res = np.nansum((predicted[usable] - measured[usable])**2)
        SS_total = np.nansum((measured[usable] - measured_mean)**2)
    
        r2_value = 1 - SS_res/SS_total
        
        return r2_value
    
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
        
        ind = np.squeeze([np.nonzero(t_pred==t)[0] for t in t_true])
        t_pred = t_pred[ind]
        pool_pred = pool_pred[ind]
        
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
    
    def _call(self, args, **kwargs):
        return self._left(args, **kwargs) + self._right(args, **kwargs)
    
    def __call__(self, args, **kwargs):
        value = self._call(args, **kwargs)
        self._call_log.append((args, kwargs, value))
        
        for callback in self._callbacks:
            callback()
            
        return value
    
    def call_count(self):
        return len(self._call_log)
    
    def last_call(self):
        return self._call_log[-1]
    
    def best_call(self):
        return sorted(self._call_log, key = lambda k: k[-1])[0]
    
    def add_callback(self, callback):
        assert callable(callback)
        callback.set_objective(self)
        self._callbacks.append(callback)

    def __str__(self):
        return ' + '.join([str(s) for s in [self._left, self._right]
                          if not s is None])
    
class Objective(Addable):
    def __init__(self, model, replica, replica_weight = 1.0, val_replica = None):
        super().__init__()
        self._model = model
        self._variables = model.parameters().variables()

        self._replica = replica
        self._replica_weight = replica_weight
        self._val_replica = val_replica
        
        self._loss_contributions = []
            
    def add_loss(self, loss_function, loss_weight):
        assert callable(loss_function)
        #loss_function.set_model(self._model, self._replica)
        self._model.add_t(self._replica.incubation['days'])
        self._loss_contributions.append((loss_function, loss_weight))

    def set_parameters(self, parameters, values, transformed = False):
        if transformed:
            values = [var.inverse_transform(p) 
                            for var, p in zip(parameters, values)]
        _ = [v.set(p) for v, p in zip(parameters, np.squeeze(values))]
        
    def _call(self, transformed_parameters, transformed = True):
        # TODO: variables are empty, if setting specific parameters!
        # => allow calling with ModelParameters()?
        # => but for minimize, must be plain array of numbers!
        if isinstance(transformed_parameters, parameters.ModelParameters):
            p_dict = transformed_parameters.as_dict()
            p_values = list(p_dict.values())
            pars = self.model().parameters()
            self.set_parameters([pars[n] for n in list(p_dict.keys())], 
                                [p.value for p in p_values], transformed = transformed)
        else:
            self.set_parameters(self._variables, transformed_parameters, 
                                transformed = transformed)
        
        run_log = self._model.predict(self._replica)
        
        replica_loss = 0
        for loss_function, loss_weight in self._loss_contributions:
            try:
                loss_value = loss_function(self._replica, run_log)
            except Exception as ex:
                ex.args = ('Replica ' + str(self._replica) + ': ' + ex.args[0], ) + ex.args[1:]
                raise 
            replica_loss += loss_weight*loss_value
        
        objective_value = self._replica_weight * replica_loss
        return objective_value

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
        s_loss = ' + '.join([format_weight(w) + str(l) for l,w in self._loss_contributions])
        replica = str(self._replica.sample.sample_name) + '-' + str(self._replica.replica_number)
        s_replica = format_subscript(replica)
        
        return s_weight + '[ '+ s_loss + ' ]' + s_replica
        

class Algorithm(ABC):
    def __init__(self, **default_kwargs):
        self._kwargs = default_kwargs
        self._configure(**default_kwargs)

    def _configure(self, **kwargs):
        for name, value in kwargs.items():
            setattr(self, name, value)
    
    def configure(self, kwargs):
        self._configure(**kwargs)
    
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

        try:
            _ = scipy.optimize.differential_evolution(objective,
                                                      bounds = self.get_bounds(variables),
                                                      callback = self.generation_counter,
                                                      strategy = self.strategy,
                                                      )
            """
                                                      updating = self.updating,
                                                      popsize = self.popsize,
                                                      workers = self.workers,
                                                      tol = self.tol,
                                                      init = self.init,
                                                      polish = self.polish,
                                                      recombination = self.recombination,
                                                      mutation = self.mutation
                                                      )
    """
            #,
                                                      #**self._kwargs)
        except KeyboardInterrupt:
            return
        
    def generation_counter(self, args, **kwargs):
        self.generation += 1
        return False


class Powell(Algorithm):
    def _minimize(self, objective, initial_parameters):
        variables = objective.model().parameters().variables()
        x0 = np.reshape([v.transform(v.value) for v in variables], (-1,))
        _ = scipy.optimize.minimize(objective, x0, method = 'Powell', 
                                    bounds = self.get_bounds(variables))


