import gc
import os
import json
import traceback
import numpy as np
import scipy.optimize
import matplotlib.pyplot as plt
from datetime import datetime

from data import compute_rate

import USER_VARIABLES

import linecache
import os


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
        
        if callable(self.transform):
            pool_pred = self.transform(pool_pred)
            pool_true = self.transform(pool_true)
            
        return pool_pred, pool_true
    
    def __str__(self):
        s = f'{self.pool}'
        return s
    
    def __call__(self, replica, run_log):
        pool_pred, pool_true = self.get_values(replica, run_log)
        return self.reduction_function(pool_true, pool_pred)

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
        return '\n'.join([str(s) for s in [self._left, self._right]
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

    def _call(self, transformed_parameters):
        parameter_values = [var.inverse_transform(p) 
                            for var, p in zip(self._variables, transformed_parameters)]
        _ = [v.set(p) for v, p in zip(self._variables, np.squeeze(parameter_values))]
        
        run_log = self._model.predict(self._replica)
        
        replica_loss = 0
        for loss_function, loss_weight in self._loss_contributions:
            loss_value = loss_function(self._replica, run_log)
            replica_loss += loss_weight*loss_value
        
        objective_value = self._replica_weight * replica_loss
        return objective_value

    def __str__(self):
        r_name = '-'.join([self._replica.sample.sample_name,self._replica.replica_number])
        lstr = '  '.join([str(loss_function)
                              for loss_function, _ in self._loss_contributions])
        return f'{r_name}: {lstr}'


class Algorithm():
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
        print('Fitting with ' + self.__class__.__name__)
        print(objective.model())
        print(objective)
        self._minimize(objective, initial_parameters)
    
    def _minimize(self, objective, initial_parameters):
        raise NotImplementedError()

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
                                                      callback = self.generation_counter)
            #,
                                                      #**self._kwargs)
        except KeyboardInterrupt:
            return
        
    def generation_counter(self, args, **kwargs):
        self.generation += 1
        return False

class Powell(Algorithm):
    def _minimize(self, objective, initial_parameters):
        vairables = objective.model().parameters().variables()
        x0 = np.reshape([v.transform(v.value) for v in variables], (-1,))
        _ = scipy.optimize.minimize(objective, x0, method = 'Powell', 
                                    bounds = self.get_bounds(variables))



def target_directory_path(replica_objectives, model_type, suffix = ''):
    str_model_type = '_' + model_type + '_'
    timestamp = datetime.now().strftime('%Y-%m-%d--%H-%M-%S')
    name = 'fit_' + '_'.join([str(r.replica)
                              for r in replica_objectives]) + str_model_type + timestamp
    cp_path = os.path.join(USER_VARIABLES.LOG_DIRECTORY + suffix, name) 
    return cp_path

class Objective2():
    def __init__(self, replica_objectives, variables, model, suffix = '', keep_only_best = False):
        self.model = model
        self.replica_objectives = replica_objectives
        self.variables = variables
        self.generation = None
        self._call_count = 0
        self._best_call = None
        model_type = model.model_type()
        self.cp_path = target_directory_path(replica_objectives, model_type, suffix)
        
        print('checkpoint path', self.cp_path)
        
        if not os.path.isdir(self.cp_path):
            os.makedirs(self.cp_path)
        self._keep_only_best = keep_only_best

    def transform(self, parameter_values):
        transformed_parameter_values = [var.transform(p)
                                        for var, p in zip(self.variables, parameter_values)]
        return np.squeeze(transformed_parameter_values)
    
    def inverse_transform(self, transformed_parameter_values):
        try:
            parameter_values = [var.inverse_transform(p) 
                            for var, p in zip(self.variables, transformed_parameter_values)]
        except Exception as ex:
            raise Exception(str(transformed_parameter_values))
        return np.squeeze(parameter_values)
    
    def __call__(self, transformed_parameter_values):
        self._call_count += 1
        parameter_values = self.inverse_transform(transformed_parameter_values)
        _ = [v.set(p) for v, p in zip(self.variables, np.squeeze(parameter_values))]
        total_loss = sum([obj() for obj in self.replica_objectives])   
        if np.isnan(total_loss):
            return 999.

        if self._best_call is None or total_loss < self.best_call()[0]:
            print()
            print('calls', f'{self._call_count:6d}', 'best total loss', total_loss)
            parameter_dict = self.model.parameters().get_config()
            #parameter_dict = {var.name:p
            #                for var, p in zip(self.variables, parameter_values)}
            replica_objectives = self.replica_objectives
            if not isinstance(replica_objectives, list):
                replica_objectives = [replica_objectives]
            checkpoint_file = os.path.join(self.cp_path, self.file_name(total_loss))
            with open(checkpoint_file, 'w') as cf:
                json.dump(parameter_dict, cf, indent = 4)
                
            if self._keep_only_best:
                for file in os.listdir(self.cp_path):
                    path = os.path.join(self.cp_path, file)
                    if path == checkpoint_file: continue
                    os.remove(path)
                
            self._best_call = (total_loss, parameter_dict)
            
            for ro in self.replica_objectives:
                predicted_CO2_on_measured, predicted_CH4_on_measured = ro.last_call
                
                used_days = ro.days[ro.used_indices]
                replica_days, replica_CO2 = ro.replica.CO2()
                _, replica_CH4 = ro.replica.CH4()
                
               # used_indices = np.squeeze([np.nonzero(replica_days == t)[0] for t in used_days])
                used_CO2 = replica_CO2[ro.used_indices]
                used_CH4 = replica_CH4[ro.used_indices]
                
                co2_r2 = r2(predicted_CO2_on_measured, used_CO2, log = ro.log_co2)
                ch4_r2 = r2(predicted_CH4_on_measured, used_CH4, log = ro.log_ch4)
                print('replica ', ro.replica, 'R2:', 'CO2', f'{co2_r2:.3f}', 'CH4', f'{ch4_r2:.3f}')
        best, _ = self._best_call
        genstr = 'None (local)'
        if not self.generation is None:
            genstr = str(self.generation)
        print('generation', genstr, 'calls', f'{self._call_count:6d}','total loss', f'{total_loss:7.2f}', 'current best', f'{best:7.2f}')

        gc.collect()
        return total_loss
  
    def file_name(self, total_loss):
        return f'call_{self._call_count:03d}_loss_{total_loss:7.4f}'

    def get_callback(self):        
        def callback(*args, **kwargs):
            self.generation += 1
            print('callback')
            return False
        return callback

    def best_call(self):
        return self._best_call
    
    def __str__(self):
        return 'Objective function: sum of loss from\n' + '\n'.join([str(s) 
                                                         for s in self.replica_objectives])

class Loss2():
    def __init__(self, predicted, measured):
        self.predicted = predicted
        self.measured = measured
    
    def MSE(self, weights = None):
        residuals = self.predicted - self.measured
        if not weights is None:
            return np.mean(weights*residuals**2)
        return np.mean(residuals**2)
    
    def RMSE(self):
        return np.sqrt(self.MSE())
    
    def R2(self):
        raise NotImplementedError()
        
    def NRMSE(self):
        # normalizing RMSE by value range
        rmse = self.RMSE()
        rng = np.max(self.measured) - np.min(self.measured)
        return rmse/rng
    
    def compute(self, name, **kwargs):
        if name.lower() == 'mse':
            return self.MSE(**kwargs)
        elif name.lower() =='rmse':
            return self.RMSE()
        else:
            raise NotImplementedError()

class ReplicaObjective():
    def __init__(self, replica, model, log_co2 = False, log_ch4 = False, 
                 fit_from = 0, fit_to = None,
                 loss_weight_CO2 = 1, loss_weight_CH4 = 1,
                 weighted_measurements = False,
                 loss_function_co2 = 'mse', 
                 loss_function_ch4 = 'mse',
                 rate_penalty = 0,
                 normalized = False,
                 initial_mean_days = 0):
        self.model = model
        self.replica = replica
        self.last_call = None
        self.log_co2 = log_co2
        self.log_ch4 = log_ch4
        self.loss_function_co2 = loss_function_co2
        self.loss_function_ch4 = loss_function_ch4
        self.w_CO2 = loss_weight_CO2
        self.w_CH4 = loss_weight_CH4
        self.weighted_measurements = weighted_measurements
        self._weights = None
        self._normalized = normalized
        self.rate_penalty = rate_penalty
        
        self.initial_mean_days = initial_mean_days
        
        if not fit_to is None and fit_from >= fit_to:
            raise ValueError('"from" value >= "to" value')
        self._fit_from = fit_from
        
        if not fit_to is None and fit_to <= 0:
            raise Exception('"to" value must be > 0')
        self._fit_to = fit_to
        
        days = self.replica.incubation['days']
        if not self.replica.last_day is None:
            days = days[days <= self.replica.last_day]
        self.days = days
        
        self.used_indices = np.nonzero(self.days >= self._fit_from)[0]
        if not self._fit_to is None:
            self.used_indices = np.intersect1d(self.used_indices, np.nonzero(self.days <= self._fit_to)[0])
        
        
    def __call__(self):        
        try:
            results = self.model.predict(self.replica, 
                                         self.days[self.used_indices], 
                                         quiet = True,
                                         initial_mean_days = self.initial_mean_days)
        except KeyboardInterrupt:
            while True:
                inp = input('continue ? [Y/n]> ')
                if inp == '' or inp == 'y':
                    return np.nan # continue
                elif inp == 'n':
                    raise Exception()
        except Exception as ex:
            print(traceback.format_exc())
            input()
            # This also catches timeout Exception or KeyboardInterrupt, should LSODA get stuck
            return np.nan
        
        
        _, predicted_CO2 = results['CO2']
        t, predicted_CH4 = results['CH4']
        self.last_call = predicted_CO2, predicted_CH4
    
        if self._weights is None:
            if not self.weighted_measurements:
                self._weights = np.ones((predicted_CO2.size,))
            delta_t = np.diff(t)
            weights = np.concatenate([[delta_t[0]/2],
                                            (delta_t[0:-1]+delta_t[1:])/2,
                                            [delta_t[-1]/2]], axis = 0)
            self._weights = weights/np.max(weights)
        w = np.array(self._weights)

        measured_CO2 = self.replica.incubation['CO2'][self.used_indices]
        measured_CH4 = self.replica.incubation['CH4'][self.used_indices]
        if 0 in t:
            measured_CO2 = measured_CO2[t != 0]
            predicted_CO2 = predicted_CO2[t != 0]
            measured_CH4 = measured_CH4[t != 0]
            predicted_CH4 = predicted_CH4[t != 0]
            w = w[t != 0]
            t = t[t!=0]

        if self.log_co2:
            measured_CO2 = np.log(measured_CO2)
            predicted_CO2 = np.log(predicted_CO2)
            measured_CO2 = np.where(np.isfinite(measured_CO2), measured_CO2, -12)
            predicted_CO2 = np.where(np.isfinite(predicted_CO2), predicted_CO2, -12)
        CO2_loss = Loss(predicted_CO2, measured_CO2)
        CO2_loss_value = CO2_loss.compute(self.loss_function_co2)

        if self.log_ch4:
            measured_CH4 = np.log(measured_CH4)
            predicted_CH4 = np.log(predicted_CH4)
            measured_CH4 = np.where(np.isfinite(measured_CH4), measured_CH4, -12)
            predicted_CH4 = np.where(np.isfinite(predicted_CH4), predicted_CH4, -12)
        CH4_loss = Loss(predicted_CH4, measured_CH4)
        CH4_loss_value = CH4_loss.compute(self.loss_function_ch4)

        w_CO2 = self.w_CO2
        w_CH4 = self.w_CH4
        if self._normalized:
            w_CO2 *= 1/(np.max(measured_CO2) - np.min(measured_CO2))
            w_CH4 *= 1/(np.max(measured_CH4) - np.min(measured_CH4))
            w_CO2 /= w_CH4
            w_CH4 = 1.
        loss = w_CO2*CO2_loss_value + w_CH4*CH4_loss_value

        if self.rate_penalty > 0:
            meas_CO2_rate = compute_rate(t, measured_CO2)
            pred_CO2_rate = compute_rate(t, predicted_CO2)
            meas_CH4_rate = compute_rate(t, measured_CH4)
            pred_CH4_rate = compute_rate(t, predicted_CH4)
            
            rate_loss_co2 = Loss(pred_CO2_rate, meas_CO2_rate)
            rate_loss_ch4 = Loss(pred_CH4_rate, meas_CH4_rate)
            
            rate_loss_value_co2 = rate_loss_co2.compute('mse')
            rate_loss_value_ch4 = rate_loss_ch4.compute('mse')
            
            pnlty = self.rate_penalty*(w_CO2*rate_loss_value_co2 + w_CH4*rate_loss_value_ch4)
            print('Loss', loss, 'penalty', pnlty)
            loss += pnlty
        return loss
    
    def plot(self):
        if not self.last_call:
            return
        
        predicted_CO2, predicted_CH4 = self.last_call
        days = self.replica.incubation['days']
        measured_CO2 = self.replica.incubation['CO2']
        measured_CH4 = self.replica.incubation['CH4']

        plt.figure()
        plt.plot(days, measured_CO2, 'rx')
        plt.plot(days, measured_CH4, 'bx')

        plt.plot(days, predicted_CO2, 'k-')
        plt.plot(days, predicted_CH4, 'k--')

        plt.title(str(self.replica))
        plt.show()
    
    def __str__(self):
        return f'fit to replica {self.replica}'
    
    

