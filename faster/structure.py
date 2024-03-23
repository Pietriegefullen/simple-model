import numpy as np
import scipy.integrate
import matplotlib.pyplot as plt

import system
import chemistry
import CONSTANTS

HENRYS_LAW = system.henrys_law()

# TODO: logging values: log only a single final result per t
# TODO: build model after setting parameters/before each run (if parameters changed)

class Parameter(float):
    def __new__(cls, name, value):
        if value is None:
            value = np.nan
        return float.__new__(cls, float(value))
    
    def __init__(self, name, value):
        self._name = name
        self._value = value
        
    def __str__(self):
        return f'{self.__class__.__name__: <9s}({self._name[:15]: <15s}, {self._value:7.2g})'

    def value(self):
        return float(self)
        
class Variable(Parameter):
    def __new__(cls, name, bounds, value, scale):
        return super().__new__(cls, name, value)
    
    def __init__(self, name, bounds, default, scale = 'log'):
        super().__init__(name, default)
        self._low = bounds[0]
        self._high = bounds[1]
        self._scale = scale
    
    def __str__(self):
        return f'{self.__class__.__name__: <9s}({self._name[:15]: <15s}, {self._value:7.2g}, ({self._low: >7.2g}, {self._high: >7.2g}), {self._scale})'
        
class ModelParameters():
    def __init__(self):
        self._parameters = {}
                            
    def add(self, name, default_value = None, bounds = None, scale = 'linear'):
        if name in self._parameters:
            p = self._parameters[name]
            
            if not p.value is None:
                # parameter is already defined with a value (and bounds)
                if not default_value is None and not p.value == default_value:
                    raise Exception(f'Model parameter with name {name} is already defined with values.')
                else:
                    return p
            
            else:
                # parameter is defined without value and bounds
                if not default_value is None:
                    p.value = default_value
                    if not bounds is None:
                        p.low = bounds[0]
                        p.high = bounds[1]
                        p.scale = scale
                    return p
            
        if not bounds is None and default_value is None:
            raise Exception('Model parameter requires default value if bounds are specified.')
        
        if not scale == 'linear' and bounds is None:
            raise Exception('The "scale" argument only takes effect if bounds are set.s')
        
        if bounds is None:
            p = Parameter(name, default_value)
            
        else:
            p = Variable(name, bounds, default_value, scale)
            
        self._parameters[name] = p
        return p
    
    def check(self):
        nan_pars = [v for v in self._parameters.values() if np.isnan(v)]
        if len(nan_pars) > 0:
            raise Exception('Model parameters are NaN:\n' + '\n'.join([str(p) for p in nan_pars]) )
    
    def __getitem__(self, key):
        return self._parameters[key]
    
    def __contains__(self, other):
        compare = str(other)
        if hasattr(other, 'name'):
            compare = other.name
        return compare in self._parameters.keys()
    
    def __str__(self):
        return '\n'.join([str(p) for p in self._parameters.values()])


class Model():
    def __init__(self, contributing_pathways):
        self.system_state_log = ModelRun()
        self.model_parameters = ModelParameters()
        
        self.contributing_pathways = [p(self.model_parameters) 
                                      for p in contributing_pathways]
        for p in self.contributing_pathways:
            p.inject_logger(self.system_state_log)
        
    def __call__(self, t, S):
        print(t)
        S = np.where(S < 1e-40, 0, S)
        
        for Si, pool_name in zip(S, system.SYSTEM):
            self.system_state_log.log(pool_name, t, Si)
        
        dSj_dt = np.stack([pathway_j(t, S)
                            for pathway_j in self.contributing_pathways], axis = -1)
        
        dS_dt = np.sum(dSj_dt, axis = -1)
        dS_dt = np.clip(dS_dt, -S, np.inf) # don't let pools become negative
        return dS_dt
    
    def fit(self, replicas):
        if not isinstance(replicas, list):
            replicas = [replicas]
        
        objective = Objective([ReplicaObjective(replica, self) 
                               for replica in replicas])
        
        # get order of parameters (for each ReplicaObjective!)
        # get initial_values
        # get bounds
        # 
        
        raise NotImplementedError()
        # minimize the objective function
        # return optimal parameters
        
        
    def predict(self, t, replica):
        S0 = system.initial_state(replica, self.parameters())
        self.parameters().check()
        self.system_state_log.reset()
        solver_result = scipy.integrate.solve_ivp(self, (0, max(t)),
                                                  S0, t_eval = t,
                                                  method = 'LSODA',
                                                  first_step = 1e-6)
        return self.system_state_log
    
    def parameters(self):
        return self.model_parameters
    
    def __str__(self):
        model_string = f'Model with {len(self.contributing_pathways)} Pathways:\n'
        model_string += len(model_string)*'=' + '\n'
        model_string += '\n'.join([str(p) for p in self.contributing_pathways])
        model_string += '\n' + 'Parameters:\n' + '='*len('Parameters') + '\n'
        model_string += str(self.model_parameters)
        return model_string

class Objective():
    def __init__(self, replica_objectives):
        self.replica_objectives = replica_objectives
        
    def __call__(self, parameter_values):
        total_loss = sum([obj(parameter_values) for obj in self.replica_objectives])   
        return total_loss
    
class ReplicaObjective():
    def __init__(self, replica, model):
        self.model = model
        self.replica = replica
        
    def __call__(self, parameter_values):
        model_parameters = {} # TODO: transform the parameter values!
        raise NotImplementedError()
        self.model.set_parameters(model_parameters)
        
        days = self.replica.incubation['days']
        results = self.model.predict(days, self.replica)
        
        _, predicted_CO2 = zip(*results['CO2'])
        measured_CO2 = self.replica.incubation['CO2']
        
        _, predicted_CH4 = zip(*results['CH4'])
        measured_CH4 = self.replica.incubation['CH4']
        
        loss = np.mean((predicted_CO2 - measured_CO2)**2) + np.mean((predicted_CH4 - measured_CH4)**2)
        
        return loss
    
class ModelRun():
    def __init__(self):
        self._log = {}
        
    def log(self, name, t, value):
        if not name in self._log:
            self._log[name] = []
        
        self._log[name].append((t,value))
        
    def reset(self):
        self._log = {}
        
    def compute_R2(self, replica):
        raise NotImplemented()
        
    def plot(self, name, newfigure = True):
        if not isinstance(name, list):
            name = [name]
        
        for n in name:
            if not n in self._log:
                print(n + ' not logged')
            if newfigure:
                plt.figure()
            x, y = zip(*self._log[n])
            plt.plot(x, y, 'x-', label = n)
            plt.title(n)
        
    def __str__(self):
        run_string = 'Model run:'
        run_string += '\n' + '='*len(run_string) + '\n'
        run_string += '\n'.join([name + ' ' + str(self._log[name])
                          for name in sorted(self._log.keys())])
        return run_string
    
class Pathway():
    def __init__(self, microbe, educts, products):        
        self.educts = educts
        self.products = products
        self.microbe = microbe
        
        self.Km = np.sum(np.stack([system.vector(0, educt, educt['Km'])
                                   for educt in educts], axis = -1), axis = -1)
        
        stoich_vector = np.sum(np.stack([system.vector(0, subst, subst['stoichiometry'])
                                         for subst in (educts + products)], axis = -1), axis = -1)
        self.stoichiometry = stoich_vector
        self.inhibition = system.vector(np.inf)
        for product in products:
            self.inhibition[system.index(product)] = product['inhibition']
            
        self.v_max = microbe['v_max']
        self.death_rate = system.vector(0, microbe, microbe['death_rate'])
        
        CUE = microbe['CUE']
        C_source = microbe['C_source']
        self.anabolism = system.vector(0)
        if not C_source is None:
            C_atoms = chemistry.C_atoms(C_source)
            growth = system.vector(0, microbe, CUE/(1-CUE)*C_atoms*CONSTANTS.MOLAR_MASS_C)
            self.anabolism = growth - system.vector(0, microbe['C_source'],CUE/(1-CUE))
            
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
            self.log('deltaG_r', t, deltaG_r)
            deltaG_rmin = chemistry.GIBBS_MINIMUM
            
            thermodynamic_factor = 1 - np.exp(np.minimum(0.,deltaG_r - deltaG_rmin)/(R*T))
        
        self.log('thermodynamic_factor', t, thermodynamic_factor)
        return thermodynamic_factor
    
    def __call__(self, t, S):        
        biomass = S[self.microbe_index]
        biomass = np.clip(biomass, 1e-8,np.inf)
        
        pathway_vector = self.stoichiometry + self.anabolism
    
        dissolved_S = HENRYS_LAW*S
        
        MM = np.where((self.Km + dissolved_S) == 0, 
                      1,
                      np.where(dissolved_S == 0,
                               0,
                               dissolved_S/(self.Km + dissolved_S + 1e-7)))
        
        inhib = np.where(np.logical_or(dissolved_S == 0, (self.inhibition + dissolved_S) == 0),
                         1,
                         1 - dissolved_S/(self.inhibition + dissolved_S + 1e-7))
        inhib = np.where(self.inhibition == np.inf, 1, inhib)
        thermodynamic_factor = self.thermodynamics(t, S)
        
        MM_factor = np.prod(MM)
        inhib_factor = np.prod(inhib)
        v = self.v_max * MM_factor * inhib_factor * thermodynamic_factor
        
        self.log('MM', t, MM_factor)
        self.log('inhib', t, inhib_factor)
        self.log('v', t, v)
        
        dS_dt = biomass * v * pathway_vector - biomass * self.death_rate
        dS_dt = np.clip(dS_dt, -S, np.inf)        
        
        return np.reshape(dS_dt, (-1,))

    def __str__(self):
        educts = ' + '.join([str(s.stoichiometry) + ' ' + str(s) for s in self.educts])
        products = ' + '.join([str(s.stoichiometry) + ' ' + str(s) for s in self.products])
        return f'{self.__class__.__name__: <12s}: {educts} -> {products}'
    
class Microbe():
    def __init__(self, name, v_max, 
                 death_rate = 0, 
                 Kmb = 0, 
                 CUE = 0, 
                 C_source = None, 
                 use_thermodynamics = True):
        
        if not CUE == 0:
            assert C_source is not None
            
        self.name = name
        self.death_rate = death_rate
        self.v_max = v_max
        self.Kmb = Kmb
        self.CUE = CUE
        self.C_source = C_source
        self.use_thermodynamics = use_thermodynamics
        
    def __getitem__(self, key):
        return getattr(self,key)
        
    def __str__(self):
        return self.name

class Substance():
    def __init__(self, stoichiometry, name, 
                 Km = None, 
                 inhibition = None):
        assert not (Km is not None and inhibition is not None)
        
        self.name = name
        self.stoichiometry = stoichiometry
        self.Km = Km if not Km is None else 0
        self.inhibition = np.inf if inhibition is None else inhibition
    
    def __getitem__(self, key):
        return getattr(self,key)
    
    def __str__(self):
        return self.name

class Hydrolysis(Pathway):
    def __init__(self, model_parameters):
        educts = [Substance(1, 'C')]
        products = [Substance(1, 'DOC')]
        microbe = Microbe(name = 'M_Ferm',
                          v_max = model_parameters.add('Hydrolysis_v_max', 1, [1e-8, 1]),
                          Kmb = model_parameters.add('Hydrolysis_Kmb', 800, [0.0005, 2000]),
                          use_thermodynamics = False)
        super().__init__(microbe, educts, products)
    
class Fermentation(Pathway):
    def __init__(self, model_parameters):
        educts = [Substance(6, 'DOC', Km = model_parameters.add('Ferm_Km', 780, [0.0005, 1000]))]
        products = [Substance(3.5, 'Acetate', inhibition = model_parameters.add('Ferm_inhibition', 7, [0.001, 20])),
                    Substance(3, 'CO2'),
                    Substance(6, 'H2')]
        microbe = Microbe(name = 'M_Ferm',
                          v_max = model_parameters.add('Ferm_v_max', 4, [0.001, 5]),
                          Kmb = model_parameters.add('Ferm_Kmb', 890, [0.0005, 2000]),
                          CUE = model_parameters.add('Ferm_CUE', .3, [0, 1]),
                          death_rate = model_parameters.add('death_rate', 8.3e-5),
                          C_source = 'DOC',
                          use_thermodynamics = False)
        super().__init__(microbe, educts, products)
        

class Hydro(Pathway):
     def __init__(self, model_parameters):
        educts = [Substance(4, 'H2', Km = model_parameters.add('Hydro_Km_H2', 77, [.0005, 1000])),
                  Substance(1, 'CO2', Km = model_parameters.add('Hydro_Km_CO2', 77, [.0005, 1000]))]
        products = [Substance(1, 'CH4'),
                    Substance(2, 'H2O')]
        microbe = Microbe(name = 'M_Hydro',
                          v_max = model_parameters.add('Hydro_v_max', .24, [0.003, 1.]),
                          CUE = model_parameters.add('Hydro_CUE', .3, [0, 1]),
                          death_rate = model_parameters.add('death_rate'),
                          C_source = 'CO2')
        super().__init__(microbe, educts, products)

class Homo(Pathway):
    def __init__(self, model_parameters):
        educts = [Substance(4, 'H2', Km = model_parameters.add('Homo_Km_H2', 77, [0.0005, 1000])),
                  Substance(2, 'CO2', Km = model_parameters.add('Homo_Km_CO2', 77, [0.0005, 1000]))]
        products = [Substance(1, 'Acetate'),
                    Substance(2, 'H2O')]
        microbe = Microbe(name = 'M_Homo',
                          v_max = model_parameters.add('Homo_v_max', .5, [0.005, 1.]),
                          CUE = model_parameters.add('Homo_CUE', .3, [0, 1]),
                          death_rate = model_parameters.add('death_rate'),
                          C_source = 'CO2')
        super().__init__(microbe, educts, products)

class Aceto(Pathway):
    def __init__(self, model_parameters):
        educts = [Substance(1, 'Acetate', Km = model_parameters.add('Aceto_Km_Ac', 77, [0.0005, 1000]))]
        products = [Substance(1, 'CH4'),
                    Substance(1, 'CO2')]
        microbe = Microbe(name = 'M_Ac',
                          v_max = model_parameters.add('Ac_v_max', .56, [0.005, 1.]),
                          CUE = model_parameters.add('Ac_CUE', .3, [0, 1]),
                          death_rate = model_parameters.add('death_rate'),
                          C_source = 'Acetate')
        super().__init__(microbe, educts, products)

class Fe3(Pathway):
    def __init__(self, model_parameters):
        educts = [Substance(1, 'Acetate', Km = model_parameters.add('Fe3_Km_Ac', 77, [0.0005, 1000])),
                  Substance(1, 'H2O'),
                  Substance(8, 'Fe3', Km = model_parameters.add('Fe3_Km_Fe3', 738, [0.0005, 1000]))]
        products = [Substance(8, 'Fe2'),
                    Substance(2, 'CO2')]
        microbe = Microbe(name = 'M_Fe3',
                          v_max = model_parameters.add('Fe3_v_max', 0.898, [0.002, 3.]),
                          CUE = model_parameters.add('Fe3_CUE', 0.3, [0, 1]),
                          death_rate = model_parameters.add('death_rate'),
                          C_source = 'Acetate')
        super().__init__(microbe, educts, products)

def get_pathways(model_type):
    basic = [Hydrolysis,
             Fermentation,
             Hydro,
             Aceto]
    if model_type == 'complex':
        return basic + [Homo]
    elif model_type == 'simple':
        return basic
    else:
        raise NotImplementedError()

if __name__ == '__main__':
    # TODO: woher kommen die Werte für die initial pools?
    import data
    d = data.get_data_before_carex()
        
    model_type = 'simple'
    chosen_pathways = get_pathways(model_type)
    model = Model(chosen_pathways)
    
    mp = model.parameters()
    
    replica = d['13514']
    #model.fit(replica)
    print(model)


    t = list(range(1500))
    results = model.predict(t, replica)
    results.plot(['Hydrolysis_v', 'Hydrolysis_inhib', 'M_Hydro', 'M_Homo'])
    
    
    
    
    