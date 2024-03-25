import numpy as np
import scipy.integrate
import matplotlib.pyplot as plt

import system
import chemistry
import CONSTANTS
import optimizer

HENRYS_LAW = system.henrys_law()

OPTIMIZATION_ALGORITHM = 'dual_annealing' #'differential_evolution' #'direct' # 'gradient' # 'PSO'

class Pathway():
    def __init__(self, microbe, educts, products):        
        self.educts = educts
        self.products = products
        self.microbe = microbe
        
        self.Km = np.sum(np.stack([system.vector(0, educt, educt['Km'])
                                   for educt in educts], axis = -1), axis = -1)
        self.Km += system.vector(0, microbe, microbe['Kmb'])
        
        stoich_vector = np.sum(np.stack([system.vector(0, subst, -subst['stoichiometry'])
                                         for subst in educts], axis = -1), axis = -1)
        stoich_vector += np.sum(np.stack([system.vector(0, subst, subst['stoichiometry'])
                                         for subst in products], axis = -1), axis = -1)
        self.stoichiometry = stoich_vector
        self.inhibition = system.vector(np.inf)
        for product in products:
            self.inhibition[system.index(product)] = product['inhibition']
            
        self.v_max = microbe['v_max']
        self.death_rate = system.vector(0, microbe, microbe['death_rate'])
        
        C_source = microbe['C_source']
        self.anabolism = system.vector(0)
        if not C_source is None:
            C_atoms = chemistry.C_atoms(C_source)
            C_source_stoich = [educt for educt in educts 
                               if str(educt) == C_source][0]['stoichiometry']
            CUE = microbe['CUE']
            anabolism_fraction = C_source_stoich*CUE/(1.-CUE + 1e-7)
            microbe_growth = system.vector(0, microbe, anabolism_fraction*C_atoms*CONSTANTS.MOLAR_MASS_C)
            C_source_reduction = system.vector(0, microbe['C_source'], anabolism_fraction)
            self.anabolism = microbe_growth - C_source_reduction
            
        self.pathway_vector = self.stoichiometry + self.anabolism
        
        self.microbe_index = system.index(microbe)
        self.use_thermodynamics = microbe['use_thermodynamics']
        
        if self.use_thermodynamics:
            self.deltaG_f = np.sum(np.stack([system.vector(0, str(subst), 
                                        chemistry.GIBBS_FORMATION[str(subst)])
                                        for subst in (educts + products)],
                                             axis = -1), axis = -1)
            self.deltaG_s = np.sum(self.stoichiometry*self.deltaG_f)

        
        self.state_logger = None
        
        
    def inject_logger(self, state_logger):
        self.state_logger = state_logger
    
    def log(self, name, t, value):
        if not self.state_logger is None:
            self.state_logger.log(self.__class__.__name__ + '_' + name, t, value)
        
    def thermodynamics(self, t, S):
        thermodynamic_factor = 1.
        if self.use_thermodynamics:
            R = CONSTANTS.GAS_CONSTANT
            T = 4. + CONSTANTS.KELVIN
            log_Q = system.vector(0)
            
            contributes = np.logical_and(self.stoichiometry != 0, S > 0)
            log_Q[contributes] = np.log(1e-6*S[contributes])
    
            deltaG_r = self.deltaG_s + R*T*np.sum(self.stoichiometry*log_Q)
            deltaG_rmin = chemistry.GIBBS_MINIMUM
            
            thermodynamic_factor = 1 - np.exp(np.minimum(0.,deltaG_r - deltaG_rmin)/(R*T))
            self.log('deltaG_r', t, deltaG_r)
            
        self.log('thermodynamic_factor', t, thermodynamic_factor)
        return thermodynamic_factor
    
    def __call__(self, t, S):
        biomass = S[self.microbe_index]
        biomass = np.clip(biomass, 1e-8,np.inf)
            
        dissolved_S = HENRYS_LAW*S

        eps = np.where(dissolved_S == 0, 1e-8, 0) # no effect, only to suppress warning of invalid value

        MM = np.where((self.Km + dissolved_S) == 0, 
                      1,
                      np.where(dissolved_S == 0,
                               0,
                               dissolved_S/(self.Km + dissolved_S + eps)))
        
        inhib = np.where(np.logical_or(dissolved_S == 0, (self.inhibition + dissolved_S) == 0),
                         1,
                         1 - dissolved_S/(self.inhibition + dissolved_S))
        inhib = np.where(self.inhibition == np.inf, 1, inhib)
        thermodynamic_factor = self.thermodynamics(t, S)
        
        MM_factor = np.prod(MM)
        inhib_factor = np.prod(inhib)
        v = self.v_max * MM_factor * inhib_factor * thermodynamic_factor

        dS_dt = biomass * v * self.pathway_vector - biomass * self.death_rate
        dS_dt = np.clip(dS_dt, -S, np.inf)        
        
        self.log('MM', t, MM_factor)
        self.log('inhib', t, inhib_factor)
        self.log('v', t, v)
        
        return np.reshape(dS_dt, (-1,))

    def __str__(self):
        educts = ' + '.join([str(s.stoichiometry) + ' ' + str(s) for s in self.educts])
        products = ' + '.join([str(s.stoichiometry) + ' ' + str(s) for s in self.products])
        pwy_string = f'{self.__class__.__name__: <12s}: {educts} -> {products}'    
        return pwy_string


class Model():
    def __init__(self, contributing_pathways):
        self.system_state_log = ModelRun()
        self.model_parameters = ModelParameters()
        self._unbuilt_contributing_pathways = contributing_pathways
        self.contributing_pathways = None
        self.build(quiet = True)
        
    def __call__(self, t, S):
        S = np.where(S < 1e-40, 0, S)
        
        dSj_dt = np.stack([pathway_j(t, S)
                            for pathway_j in self.contributing_pathways], axis = -1)
        
        dS_dt = np.sum(dSj_dt, axis = -1)
        dS_dt = np.clip(dS_dt, -S, np.inf) # don't let pools become negative
        return dS_dt
    
    def build(self, quiet = False):
        self.contributing_pathways = [p(self.model_parameters) 
                                      for p in self._unbuilt_contributing_pathways]
        for p in self.contributing_pathways:
            p.inject_logger(self.system_state_log)
        
        # to initialize model parameters used in initial state
        _ = system.initial_state(None, self.model_parameters)
        
    def fit(self, replicas):
        if not isinstance(replicas, list):
            replicas = [replicas]
            
        algo = optimizer.Algorithm(OPTIMIZATION_ALGORITHM, 
                                   **optimizer.algo_kwargs(OPTIMIZATION_ALGORITHM))
        algo.minimize(self, replicas)
        
        
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
    
    def __str__(self):
        title = 'Model Parameters:'
        title += '\n' + '='*len(title) + '\n'
        sorted_params = sorted(self._parameters.values(), key = lambda x: x.name)
        return title + '\n'.join([f'{i+1:3d}) ' + str(p) 
                                  for i,p in enumerate(sorted_params)])


class ModelRun():
    def __init__(self):
        self._log = {}
        
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
    
    