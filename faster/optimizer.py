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
        return {'strategy': 'best1bin',
                'updating': 'immediate',
                'workers': -1,
                'recombination': .3, # CR
                'mutation': (.3,.8)  # F
                }
    
    elif method == 'direct' or method == 'dual_annealing':
        return {}

    else:
        raise NotImplementedError()

class Algorithm():
    def __init__(self, algorithm, **kwargs):
        self.algorithm = algorithm
        self.kwargs = kwargs
    
    def minimize(self, model, replicas, log = False):
        variables = model.parameters().variables()
        lower_bounds = np.reshape([v.transform(v.lower()) for v in variables], (-1,))
        upper_bounds = np.reshape([v.transform(v.upper()) for v in variables], (-1,))
        
        print()
        rep = '_'.join([str(r) for r in replicas])
        str_log = '' if not log else 'log transformed '
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
              
        replica_obj = [ReplicaObjective(replica, model, log = log)
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
            
            objective = Objective(replica_obj, variables, model, log = log)
            bounds = list(zip(lower_bounds, upper_bounds))
            _ = scipy.optimize.differential_evolution(objective,
                                                      bounds = bounds,
                                                      strategy = strategy,
                                                      updating = updating, 
                                                      callback = objective.get_callback())

        else:
            raise NotImplementedError()
            
        _, optimal_parameters = objective.best_call()
        _ = [variable.set(value) for variable, value in zip(variables, optimal_parameters)]
        
        return objective.best_call()

def target_directory_path(replica_objectives, log, model_type):
    str_log = '' if not log else '_log'
    str_model_type = '_' + model_type + '_'
    timestamp = datetime.now().strftime('%Y-%m-%d--%H-%M-%S')
    name = 'fit_' + '_'.join([str(r.replica)
                              for r in replica_objectives]) + str_log + str_model_type + timestamp
    cp_path = os.path.join(USER_VARIABLES.LOG_DIRECTORY, name) 
    return cp_path

class Objective():
    def __init__(self, replica_objectives, variables, model, log = False):
        self.model = model
        self.replica_objectives = replica_objectives
        self.variables = variables
        self._call_count = 0
        self._best_call = None
        self.log = log
        model_type = model.model_type()
        self.cp_path = target_directory_path(replica_objectives, log, model_type)
        if not os.path.isdir(self.cp_path):
            os.makedirs(self.cp_path)

    def transform(self, parameter_values):
        transformed_parameter_values = [v.transform(p)
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
        if self._call_count % 5 == 0:
            print('.', end = '', flush = True)
        parameter_values = self.inverse_transform(transformed_parameter_values)
        _ = [v.set(p) for v, p in zip(self.variables, np.squeeze(parameter_values))]
        total_loss = sum([obj() for obj in self.replica_objectives])   
        if np.isnan(total_loss):
            return 999.

        if not self._best_call or total_loss < self.best_call()[0]:
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

        gc.collect()
        return total_loss
  
    def file_name(self, total_loss):
        return f'call_{self._call_count:03d}_loss_{total_loss:.2f}'

    def get_callback(self):
        def callback(intermediate_result):
            print('callback')
            if not self._best_call is None:
                best_parameters = intermediate_result.x
                objective = intermediate_result.fun
                current_loss = objective(best_parameters)
                if current_loss < self._best_call[0]:
                    pop = intermediate_result.population
                    pop_file = os.path.join(self.cp_path, self.file_name(current_loss) + '.npy')
                    print('saving', pop_file)
                    with open(pop_file, 'w') as pf:
                        np.save(pop, pf)

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
    
    def RMSE(self):
        return np.sqrt(np.mean((self.predicted - self.measured)**2))


class ReplicaObjective():
    def __init__(self, replica, model, log = False):
        self.model = model
        self.replica = replica
        self.last_call = None
        self.log = log
        
    def __call__(self):        
        days = self.replica.incubation['days']
        try:
            results = self.model.predict(self.replica, days, quiet = True)
        except KeyboardInterrupt:
            while True:
                inp = input('continue ? [Y/n]> ')
                if inp == '' or inp == 'y':
                    return np.nan # continue
                elif inp == 'n':
                    raise Exception()
        except:
            # This also catches timeout Exception or KeyboardInterrupt, should LSODA get stuck
            return np.nan
        
        _, predicted_CO2 = results['CO2']
        measured_CO2 = self.replica.incubation['CO2']
       
        if self.log:
            measured_CO2 = np.log(measured_CO2)
            predicted_CO2 = np.log(predicted_CO2)
            
            measured_CO2 = np.where(np.isfinite(measured_CO2), measured_CO2, -12)
            predicted_CO2 = np.where(np.isfinite(predicted_CO2), predicted_CO2, -12)

        CO2_loss = Loss(predicted_CO2, measured_CO2).RMSE()
        
        _, predicted_CH4 = results['CH4']
        measured_CH4 = self.replica.incubation['CH4']
       
        if self.log:
            measured_CH4 = np.log(measured_CH4)
            predicted_CH4 = np.log(predicted_CH4)
            
            measured_CH4 = np.where(np.isfinite(measured_CH4), measured_CH4, -12)
            predicted_CH4 = np.where(np.isfinite(predicted_CH4), predicted_CH4, -12)

 
        CH4_loss = Loss(predicted_CH4, measured_CH4).RMSE()
        
        loss = CO2_loss + CH4_loss
        
        self.last_call = predicted_CO2, predicted_CH4

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
    
    
    
