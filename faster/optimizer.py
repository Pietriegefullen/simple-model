import gc
import os
import json
import traceback
import numpy as np
import scipy.optimize
import matplotlib.pyplot as plt
from datetime import datetime

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

def algo_kwargs(method):
    if method == 'PSO':
        return {'c1': .5,
                'c2': .3,
                'w': .9,
                'particles': 50,
                'iterations': 25}
    
    elif method == 'gradient':
        return {'method': 'L-BFGS-B',
                'iterations': 200}
    
    elif method == 'differential_evolution':
        return {'strategy': 'rand1bin',
                'updating': 'deferred',#'immediate',
                'popsize': 10,
                'workers': -1,
                'tol': 1e-4,
                'init': 'sobol',
                'polish': False,
                'recombination': .7, # CR
                'mutation': (.5,1.)  # F
                }
    
    elif method == 'direct' or method == 'dual_annealing':
        return {}
    elif method == 'COBYLA' or method == 'Powell':
        return {}
    
    else:
        raise NotImplementedError()

class Algorithm():
    def __init__(self, algorithm, **kwargs):
        self.algorithm = algorithm
        self.kwargs = kwargs
    
    def minimize(self, model, replicas, log_co2 = True, log_ch4 = True,
                 fit_from = 0, fit_to = None,
                 loss_weight_CO2 = 1, loss_weight_CH4 = 1, 
                 parameter_range = None,
                 weighted_measurements = False):
        variables = model.parameters().variables()
        print(len(variables), 'variables before setting bounds') 
        if not parameter_range is None:
            for p in parameter_range:
                for v in variables:
                    if v.name == p.name and not p.high == p.low:
                        v.high = p.high
                        v.low = p.low
        variables = model.parameters().variables()
        lower_bounds = np.reshape([v.transform(v.lower()) for v in variables], (-1,))
        upper_bounds = np.reshape([v.transform(v.upper()) for v in variables], (-1,))
            
        suffix = ''
        if not fit_from == 0 or not fit_to is None:
            suffix = '_' + str(fit_from) + '-'+ str(fit_to) 
        
        print()
        rep = '_'.join([str(r) for r in replicas])
        co2_log = 'lin' if not log_co2 else 'log'
        ch4_log = 'lin' if not log_ch4 else 'log'
        str_log = f'CO2 ({co2_log}), CH4 ({ch4_log}) '
        print(f'minimizing {str_log}with {self.algorithm} for {rep}')
        print('model:')
        print(str(model))
        print()
        title = 'Variable Parameters:'
        title += '\n' + '='*len(title) + '\n'
        sorted_params = sorted(variables, key = lambda x: x.name)
        print(title + '\n'.join([f'{i+1:3d}) ' + str(p) for i, p in enumerate(sorted_params)]))
    
        if not variables:
            raise Exception('Model has no variable parameters.')
              
        replica_obj = [ReplicaObjective(replica, model, log_co2 = log_co2, log_ch4 = log_ch4,
                                        fit_from = fit_from, fit_to = fit_to,
                                        loss_weight_CO2 = loss_weight_CO2,
                                        loss_weight_CH4 = loss_weight_CH4,
                                        weighted_measurements = weighted_measurements)
                       for replica in replicas]
        
        if self.algorithm == 'PSO':
            import pyswarms as ps
            if log:
                raise NotImplementedError()
            objective = ParticleObjective(replica_obj, variables, model)
            bounds = (np.array(lower_bounds), np.array(upper_bounds))
            
            particles = self.kwargs['particles']
            iterations = self.kwargs['iterations']
            options = {'c1': self.kwargs['c1'],
                       'c2': self.kwargs['c2'],
                       'w': self.kwargs['w']}
            
            optimizer = ps.single.GlobalBestPSO(n_particles=particles, 
                                                dimensions=len(variables),
                                                options=options,
                                                bounds = bounds)
            _ = optimizer.optimize(objective,
                                    iters = iterations,
                                    n_processes = None)
        
        elif self.algorithm == 'gradient':
            objective = Objective(replica_obj, variables, model, log = log)

            method = self.kwargs['method']
            iterations = self.kwargs['iterations']
            x0 = np.reshape([v.transform(v.value) for v in variables], (-1,))
            bounds = list(zip(lower_bounds, upper_bounds))
            _ = scipy.optimize.minimize(objective,
                                        x0,
                                        bounds = bounds,
                                        method = method,
                                        options = {'maxiter': iterations})
        
        elif self.algorithm == 'direct':
            objective = Objective(replica_obj, variables, model, log = log)
            bounds = list(zip(lower_bounds, upper_bounds))
            _ = scipy.optimize.direct(objective, bounds = bounds)
                
        elif self.algorithm == 'dual_annealing':
            objective = Objective(replica_obj, variables, model, log = log)
            bounds = list(zip(lower_bounds, upper_bounds))
            _ = scipy.optimize.dual_annealing(objective, bounds = bounds)
            
        elif self.algorithm == 'differential_evolution':
            strategy = self.kwargs['strategy']
            updating = self.kwargs['updating']
            
            objective = Objective(replica_obj, variables, model, suffix)
            objective.generation = 1
            bounds = list(zip(lower_bounds, upper_bounds))
            _ = scipy.optimize.differential_evolution(objective,
                                                      bounds = bounds,
                                                      strategy = strategy,
                                                      updating = updating, 
                                                      callback = objective.get_callback())

        elif self.algorithm == 'COBYLA' or self.algorithm == 'Powell':
            objective = Objective(replica_obj, variables, model, suffix)
            x0 = np.reshape([v.transform(v.value) for v in variables], (-1,))
            bounds = list(zip(lower_bounds, upper_bounds))
            _ = scipy.optimize.minimize(objective, x0, method = self.algorithm, bounds = bounds)
            
        else:
            raise NotImplementedError()
            
        _, optimal_parameters = objective.best_call()
        _ = [variable.set(value) for variable, value in zip(variables, optimal_parameters)]
        
        return objective.best_call()

def target_directory_path(replica_objectives, model_type, suffix = ''):
    str_model_type = '_' + model_type + '_'
    timestamp = datetime.now().strftime('%Y-%m-%d--%H-%M-%S')
    name = 'fit_' + '_'.join([str(r.replica)
                              for r in replica_objectives]) + str_model_type + timestamp
    cp_path = os.path.join(USER_VARIABLES.LOG_DIRECTORY + suffix, name) 
    return cp_path

class Objective():
    def __init__(self, replica_objectives, variables, model, suffix = ''):
        self.model = model
        self.replica_objectives = replica_objectives
        self.variables = variables
        self.generation = None
        self._call_count = 0
        self._best_call = None
        model_type = model.model_type()
        self.cp_path = target_directory_path(replica_objectives, model_type, suffix)
        if not os.path.isdir(self.cp_path):
            os.makedirs(self.cp_path)

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
            self._best_call = (total_loss, parameter_dict)
            
            for ro in self.replica_objectives:
                predicted_CO2_on_measured, predicted_CH4_on_measured = ro.last_call
                
                used_days = ro.days[ro.used_indices]
                replica_days, replica_CO2 = ro.replica.CO2()
                _, replica_CH4 = ro.replica.CH4()
                print(used_days)
                print(replica_days)
                input()
                for t in used_days:
                    print('used', t)
                    print(np.nonzero(replica_days == t)[0])
                used_indices = np.squeeze([np.nonzero(replica_days == t)[0] for t in used_days])
                used_CO2 = replica_CO2[used_indices]
                used_CH4 = replica_CH4[used_indices]
                
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

class ParticleObjective(Objective):
    def __call__(self, particle_model_parameter_values):
        particle_fitnesses = []
        for transformed_parameter_values in particle_model_parameter_values:
            parameter_values = self.inverse_transform(transformed_parameter_values)
            _ = [v.set(p) for v, p in zip(self.variables, np.squeeze(parameter_values))]
            losses = [obj() for obj in self.replica_objectives]
            total_loss = np.sum(losses)
            particle_fitnesses.append(total_loss)
        return np.array(particle_fitnesses)

class Loss():
    def __init__(self, predicted, measured):
        self.predicted = predicted
        self.measured = measured
    
    def MSE(self, weights = None):
        if not weights is None:
            return np.mean(weights*(self.predicted - self.measured)**2)
        return np.mean((self.predicted - self.measured)**2)
    
    def RMSE(self):
        return np.sqrt(self.MSE())
    
    def R2(self):
        raise NotImplementedError()

class ReplicaObjective():
    def __init__(self, replica, model, log_co2 = False, log_ch4 = False, 
                 fit_from = 0, fit_to = None,
                 loss_weight_CO2 = 1, loss_weight_CH4 = 1,
                 weighted_measurements = False):
        self.model = model
        self.replica = replica
        self.last_call = None
        self.log_co2 = log_co2
        self.log_ch4 = log_ch4
        self.w_CO2 = loss_weight_CO2
        self.w_CH4 = loss_weight_CH4
        self.weighted_measurements = weighted_measurements
        self._weights = None
        
        if not fit_to is None and fit_from >= fit_to:
            raise ValueError('"from" value >= "to" value')
        self._fit_from = fit_from
        
        if not fit_to is None and fit_to <= 0:
            raise Exception('"to" value must be > 0')
        self._fit_to = fit_to
        
        self.days = self.replica.incubation['days']
        
        self.used_indices = np.nonzero(self.days >= self._fit_from)[0]
        if not self._fit_to is None:
            self.used_indices = np.intersect1d(self.used_indices, np.nonzero(self.days <= self._fit_to)[0])
        
        
    def __call__(self):        
        try:
            results = self.model.predict(self.replica, 
                                         self.days[self.used_indices], 
                                         quiet = True)
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
        if 0 in t:
            measured_CO2 = measured_CO2[t != 0]
            predicted_CO2 = predicted_CO2[t != 0]
            w = w[t != 0]

        if self.log_co2:
            measured_CO2 = np.log(measured_CO2)
            predicted_CO2 = np.log(predicted_CO2)
            measured_CO2 = np.where(np.isfinite(measured_CO2), measured_CO2, -12)
            predicted_CO2 = np.where(np.isfinite(predicted_CO2), predicted_CO2, -12)
        CO2_loss = Loss(predicted_CO2, measured_CO2).MSE(w)
        
        measured_CH4 = self.replica.incubation['CH4'][self.used_indices]
        if 0 in t:
            measured_CH4 = measured_CH4[t != 0]
            predicted_CH4 = predicted_CH4[t != 0]
            
        if self.log_ch4:
            measured_CH4 = np.log(measured_CH4)
            predicted_CH4 = np.log(predicted_CH4)
            measured_CH4 = np.where(np.isfinite(measured_CH4), measured_CH4, -12)
            predicted_CH4 = np.where(np.isfinite(predicted_CH4), predicted_CH4, -12)
        CH4_loss = Loss(predicted_CH4, measured_CH4).MSE(w)


        loss = self.w_CO2*CO2_loss + self.w_CH4*CH4_loss

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
    
    
    
