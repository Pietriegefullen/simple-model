import os

import numpy as np
import scipy.integrate
import matplotlib.pyplot as plt
import json

import system
import optimizer
import pathways


OPTIMIZATION_ALGORITHM = 'PSO' #'dual_annealing' #'differential_evolution' #'direct' # 'gradient' # 'PSO'


class Model():
    def __init__(self, pwys):
        self.system_state_log = ModelRun()
        self.model_parameters = ModelParameters()
        
        pathway_classes = [pathways.pathway_by_name(p) if isinstance(p, str) else p 
                           for p in pwys]
        self._unbuilt_contributing_pathways = pathway_classes
        self.contributing_pathways = None
        self.build(quiet = True)
        
    def build(self, quiet = False):
        self.contributing_pathways = [p(self.model_parameters) 
                                      for p in self._unbuilt_contributing_pathways]
        for p in self.contributing_pathways:
            p.inject_logger(self.system_state_log)
        
        # to initialize model parameters used in initial state
        _ = system.initial_state(None, self.model_parameters)
        
    def __call__(self, t, S):
        S = np.where(S < 1e-40, 0, S)
        
        dSj_dt = np.stack([pathway_j(t, S)
                            for pathway_j in self.contributing_pathways], axis = -1)
        
        dS_dt = np.sum(dSj_dt, axis = -1)
        dS_dt = np.clip(dS_dt, -S, np.inf) # don't let pools become negative
        return dS_dt
    
    def fit(self, replicas):
        if not isinstance(replicas, list):
            replicas = [replicas]
            
        algo = optimizer.Algorithm(OPTIMIZATION_ALGORITHM, 
                                   **optimizer.algo_kwargs(OPTIMIZATION_ALGORITHM))
        return algo.minimize(self, replicas)
        
        
    def predict(self, replica, t = None, quiet = False):
        if t is None:
            t = replica['days']
        self.build(quiet = quiet)
        S0 = system.initial_state(replica, self.parameters())
        self.parameters().check()
        self.system_state_log.reset()
        
        solver_result = scipy.integrate.solve_ivp(self, (0, max(t)),
                                                  S0, 
                                                  t_eval = t,
                                                  method = 'LSODA',
                                                  max_step = 10,
                                                  first_step = 1e-6)
        
        for t, S in zip(t, np.transpose(solver_result.y)):
            for Si, pool_name in zip(S, system.SYSTEM):
                self.system_state_log.log(pool_name, t, Si)
        
        return self.system_state_log
    
    def parameters(self):
        return self.model_parameters
    
    def __str__(self):
        model_string = f'Model with {len(self.contributing_pathways)} Pathways:\n'
        model_string += len(model_string)*'=' + '\n'
        model_string += '\n'.join([str(p) for p in self.contributing_pathways])
        model_string += '\n' + 'Parameters:\n' + '='*len('Parameters') + '\n'
        model_string += str(self.model_parameters)
        model_string += '\n'.join( [str(p) for p in self.contributing_pathways])
        return model_string

    def save(self, target_directory, file_name):
        cfg = {'pathways': [p.__class__.__name__ 
                            for p in self.contributing_pathways],
               'parameters': self.parameters().get_config()}
        if not os.path.isdir(target_directory):
            os.makedirs(target_directory)
        with open(os.path.join(target_directory, file_name + '.json'), 'w') as df:
            json.dump(cfg, df, indent = 4)
            
    def load(self, file):
        with open(file, 'r') as df:
            cfg = json.load(df)
        
        self.__init__(cfg['pathways'])
        self.model_parameters.set(cfg['parameters'])

class LogTransform():
    def transform(self, value): return np.log(value)
    def inverse(self, value): return np.exp(value)
    
class IdentityTransform():
    def transform(self, value): return value
    def inverse(self, value): return value


class Parameter():
    def __init__(self, name, value = np.nan, range = None, scale = 'log'):
        self.name = name
        self.value = value
        
        self.low = None
        self.high = None
        self.scale = scale
        self.transformer = None
        if not range is None:
            self.low = range[0]
            self.high = range[1]
    
    def lower(self):
        if not self.is_variable():
            raise Exception()
        return self.low
    
    def upper(self):
        if not self.is_variable():
            raise Exception()
        return self.high
    
    def constant(self, value):
        self.value = value
        self.low = None
        self.high = None
        return self
    
    def variable(self, value, range):
        self.value = value
        self.low = range[0]
        self.high = range[1]
        return self
    
    def set(self, p):
        if isinstance(p, Parameter):
            self.value = p.value
            self.low = p.low
            self.high = p.high
            self.scale = p.scale
        elif isinstance(p, (int, float)):
            self.value = float(p)
        else:
            raise NotImplementedError(str(p))
    
    def is_unset(self):
        return np.isnan(self.value)
                        
    def is_variable(self):
        return not self.is_unset() and not self.low is None and not self.high is None
    
    def get_transform(self):
        if self.transformer is None:
            if self.scale == 'linear':
                self.transformer = IdentityTransform()
            
            elif self.is_variable() and self.scale == 'log':
                self.transformer = LogTransform()
            
            else:
                raise NotImplementedError()
        return self.transformer
    
    def transform(self, value):
        tf = self.get_transform()
        return tf.transform(value)
    
    def inverse_transform(self, value):
        tf = self.get_transform()
        return tf.inverse(value)
        
    def __str__(self):
        var = ''
        if self.is_variable():
            var = f'  ({self.low:.3g}, {self.high:.3g})   {self.scale}'
        return f'{self.name} = {self.value:.3g}' + var
    
    def __float__(self):
        return float(self.value)

    def __add__(self, other):
        return self.value + float(other)
    
    def __sub__(self, other):
        return self.value - float(other)
    
    def __rsub__(self, other):
        return float(other) - self.value
    
    def __mul__(self, other):
        return self.value*float(other)
    
    def __rmul__(self, other):
        return float(other)*self.value
    
    def __div__(self, other):
        return self.value/float(other)
    
    def __truediv__(self, other):
        return self.__div__(other)
    
    def __rdiv__(self, other):
        return float(other)/self.value
    
    
def default_model_parameters(model_parameters = None):
    p = [
         Parameter('Hydrolysis_v_max', 1, [1e-8, 1]),
         Parameter('Hydrolysis_Kmb', 800, [0.0005, 2000]),
         
         Parameter('Ferm_v_max',       4, [0.001, 5]),
         #Variable('Ferm_Kmb',       890, [0.0005, 2000]),
         Parameter('Ferm_Km',        780, [0.0005, 1000]),
         Parameter('Ferm_inhibition',  7, [0.001, 20]),
         Parameter('Ferm_CUE',        .3, [0, 1], 'linear'),
         
         Parameter('death_rate',  8.3e-5),
         
         Parameter('Hydro_Km_CO2',    77, [.0005, 1000]),
         Parameter('Hydro_v_max',    .24, [0.003, 1.]),
         Parameter('Hydro_CUE',       .3, [0, 1], 'linear'),
         Parameter('Hydro_Km_H2',     77, [.0005, 1000]),
         
         Parameter('Homo_Km_H2',      77, [0.0005, 1000]),
         Parameter('Homo_Km_CO2',     77, [0.0005, 1000]),
         Parameter('Homo_v_max',      .5, [0.005, 1.]),
         Parameter('Homo_CUE',        .3, [0, 1], 'linear'),
         
         Parameter('Aceto_Km_Ac',     77, [0.0005, 1000]),
         Parameter('Ac_v_max',       .56, [0.005, 1.]),
         Parameter('Ac_CUE',          .3, [0, 1], 'linear'), 
         
         Parameter('Fe3_Km_Ac',       77, [0.0005, 1000]),
         Parameter('Fe3_Km_Fe3',     738, [0.0005, 1000]),
         Parameter('Fe3_v_max',    0.898, [0.002, 3.]), 
         Parameter('Fe3_CUE',        0.3, [0, 1], 'linear'),
         
         Parameter('Acetate',         50, [0, 100], 'linear'),
         Parameter('Fe3',             20, [0, 300], 'linear'),
         
         Parameter('M_Ferm',          .2, [1e-8, 0.5]),
         Parameter('M_Hydro',      .0025, [1e-8, 0.5]),
         Parameter('M_Homo',       .0001, [1e-8, 0.5]),
         Parameter('M_Ac',         .0101, [1e-8, 0.5]),
        ]
        
    if not model_parameters is None:
        return [par for par in p if par in model_parameters]
    return p
    
class ModelParameters():
    def __init__(self):
        self._parameters = {}
    
    def __getitem__(self, key):
        if not key in self._parameters:
            self._parameters[key] = Parameter(key)
        return self._parameters[key]
            
    def set(self, parameters):
        if isinstance(parameters, str) and parameters == 'default':
            self.set(default_model_parameters(self))
            
        elif isinstance(parameters, dict):
            for p, value in parameters.items():
                if not p in self._parameters:
                    raise Exception('Setting parameter that does not exist:', p)
                try:
                    self._parameters[p].set(value)
                except Exception as ex:

                    raise Exception(parameters)
                    
        elif isinstance(parameters, list):
            for p in parameters:
                if p.name in self._parameters:
                    self._parameters[p.name].set(p)
                else:
                    raise Exception('Setting parameter that does not exist: ' + str(p))
        elif isinstance(parameters, tuple):
            name, value = parameters
            self.set({name:value})
        else:
            raise NotImplementedError()
            
    def variables(self):
        return [p for p in self._parameters.values() if p.is_variable()]
    
    def unset(self):
        return [p for p in self._parameters.values() if p.is_unset()]
    
    def check(self):
        nan_pars = [v for v in self._parameters.values() if v.is_unset()]
        if len(nan_pars) > 0:
            raise Exception('Model parameters are NaN:\n' + '\n'.join([str(p) for p in nan_pars]) )
    
    def __contains__(self, other):
        compare = str(other)
        if hasattr(other, 'name'):
            compare = other.name
        return compare in self._parameters.keys()
    
    def get_config(self):
        return {p.name: p.value for p in self._parameters.values()}
    
    def __str__(self):
        title = 'Model Parameters:'
        title += '\n' + '='*len(title) + '\n'
        sorted_params = sorted(self._parameters.values(), key = lambda x: x.name)
        return title + '\n'.join([f'{i+1:3d}) ' + str(p) 
                                  for i,p in enumerate(sorted_params)])


class ModelRun():
    def __init__(self):
        self._log = {}
        
    def __eq__(self, other):
        return self._log == other._log
    
    def __getitem__(self, key):
        return self._log[key]
        
    def log(self, name, t, value):
        if not name in self._log:
            self._log[name] = []
        
        self._log[name].append((t,value))
        
    def reset(self):
        self._log = {}
        
    def plot(self, name, newfigure = True):
        if not isinstance(name, list):
            name = [name]
        
        for n in name:
            if not n in self._log:
                print(n + ' not logged')
            if newfigure:
                plt.figure()
            x, y = zip(*self._log[n])
            plt.plot(x, y, '-', label = n)
            plt.title(n)
        
    def __str__(self):
        run_string = 'Model run:'
        run_string += '\n' + '='*len(run_string) + '\n'
        run_string += '\n'.join([name + ' ' + str(self._log[name])
                          for name in sorted(self._log.keys())])
        return run_string
    
    