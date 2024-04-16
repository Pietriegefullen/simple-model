import os
import json
import numpy as np
import scipy.optimize
import matplotlib.pyplot as plt
from datetime import datetime

import USER_VARIABLES

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
                'updating': 'immediate'}
    
    elif method == 'direct' or method == 'dual_annealing':
        return {}

    else:
        raise NotImplementedError()

class Algorithm():
    def __init__(self, algorithm, **kwargs):
        self.algorithm = algorithm
        self.kwargs = kwargs
    
    def minimize(self, model, replicas):
        variables = model.parameters().variables()
        lower_bounds = np.reshape([v.transform(v.lower()) for v in variables], (-1,))
        upper_bounds = np.reshape([v.transform(v.upper()) for v in variables], (-1,))
        
        print()
        rep = '_'.join([str(r) for r in replicas])
        print(f'minimizing with {self.algorithm} for {rep}')
        print('model:')
        print(str(model))
        print()
        title = 'Variable Parameters:'
        title += '\n' + '='*len(title) + '\n'
        sorted_params = sorted(variables, key = lambda x: x.name)
        print(title + '\n'.join([f'{i+1:3d}) ' + str(p) for i, p in enumerate(sorted_params)]))
    
        if not variables:
            raise Exception('Model has no variable parameters.')
              
        replica_obj = [ReplicaObjective(replica, model) for replica in replicas]
        
        if self.algorithm == 'PSO':
            import pyswarms as ps
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
            objective = Objective(replica_obj, variables, model)

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
            objective = Objective(replica_obj, variables, model)
            bounds = list(zip(lower_bounds, upper_bounds))
            _ = scipy.optimize.direct(objective, bounds = bounds)
                
        elif self.algorithm == 'dual_annealing':
            objective = Objective(replica_obj, variables, model)
            bounds = list(zip(lower_bounds, upper_bounds))
            _ = scipy.optimize.dual_annealing(objective, bounds = bounds)
            
        elif self.algorithm == 'differential_evolution':
            strategy = self.kwargs['strategy']
            updating = self.kwargs['updating']
            
            objective = Objective(replica_obj, variables, model)
            bounds = list(zip(lower_bounds, upper_bounds))
            _ = scipy.optimize.differential_evolution(objective,
                                                      bounds = bounds,
                                                      strategy = strategy,
                                                      updating = updating)
        else:
            raise NotImplementedError()
            
        _, optimal_parameters = objective.best_call()
        _ = [variable.set(value) for variable, value in zip(variables, optimal_parameters)]
        
        return objective.best_call()

class Objective():
    def __init__(self, replica_objectives, variables, model):
        self.model = model
        self.replica_objectives = replica_objectives
        self._calls = []
        self.variables = variables
        self._call_count = 0
        timestamp = datetime.now().strftime('%Y-%m-%d--%H-%M-%S')
        name = 'fit_' + '_'.join([str(r.replica)
                                  for r in replica_objectives]) + '_' + timestamp
        self.cp_path = os.path.join(USER_VARIABLES.LOG_DIRECTORY, name)
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
        parameter_values = self.inverse_transform(transformed_parameter_values)
        _ = [v.set(p) for v, p in zip(self.variables, np.squeeze(parameter_values))]
        total_loss = sum([obj() for obj in self.replica_objectives])   

        if not self._calls or total_loss < self.best_call()[0]:
            print('calls', f'{self._call_count:6d}', 'best total loss', total_loss)
            parameter_dict = self.model.parameters().get_config()
            #parameter_dict = {var.name:p
            #                for var, p in zip(self.variables, parameter_values)}
            replica_objectives = self.replica_objectives
            if not isinstance(replica_objectives, list):
                replica_objectives = [replica_objectives]
            file_name = f'call_{self._call_count:03d}_loss_{total_loss:.2f}'
            checkpoint_file = os.path.join(self.cp_path, file_name)
            with open(checkpoint_file, 'w') as cf:
                json.dump(parameter_dict, cf, indent = 4)
            self._calls.append((total_loss, parameter_dict))

        return total_loss
    
    def best_call(self):
        sorted_by_loss = sorted(self._calls, key = lambda x: x[0])
        return sorted_by_loss[0]
    
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
            self._calls.append((total_loss, parameter_values))
        return np.array(particle_fitnesses)

class Loss():
    def __init__(self, predicted, measured):
        self.predicted = predicted
        self.measured = measured
    
    def RMSE(self):
        return np.sqrt(np.mean((self.predicted - self.measured)**2))


class ReplicaObjective():
    def __init__(self, replica, model):
        self.model = model
        self.replica = replica
        self.last_call = None
        
    def __call__(self):        
        days = self.replica.incubation['days']
        results = self.model.predict(self.replica, days, quiet = True)
        
        _, predicted_CO2 = zip(*results['CO2'])
        measured_CO2 = self.replica.incubation['CO2']
        
        CO2_loss = Loss(predicted_CO2, measured_CO2).RMSE()
        
        _, predicted_CH4 = zip(*results['CH4'])
        measured_CH4 = self.replica.incubation['CH4']
        
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
    
    
    
