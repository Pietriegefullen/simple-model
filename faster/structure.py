import numpy as np
import scipy.integrate
import matplotlib.pyplot as plt

import system
import chemistry
import CONSTANTS

HENRYS_LAW = system.henrys_law()

# TODO: logging values: log only a single final result per t, which?
# TODO: build model after setting parameters/before each run (if parameters changed)
# TODO: Parameter and Variable are not working! test and do it properly!
#       standard classes, only adding a __float__
#       make really clear which return value (log or linear) is used!!!
# TODO: parallel coordinates plot for optimization runs
# TODO: throw unused microbes out of the system!
# TODO: using debug mode, compare henry's vector and pathway vectors!
#       => compare model parameter values!!!


class Parameter():
    def __init__(self, name, value = np.nan):
        self._value = float(value)
        self._name = name

    def set(self, value):
        self._value = value
        
    def __str__(self):
        return f'{self.__class__.__name__: <9s}({self._name[:15]: <15s}, {self._value:7.2g})'

    def __float__(self):
        if self._value is None:
            raise Exception(f'Model parameter {self._name} not set.')
        return float(self._value)

    def value(self):
        return float(self)
    
    def variable(self):
        return False
    
    def __eq__(self, other):
        if not isinstance(other, Parameter):
            return False
        return self._name == other._name and self._value == other._value
    
    def __add__(self, other):
        return self._value + float(other)
    
    def __sub__(self, other):
        return self._value - float(other)
    
    def __rsub__(self, other):
        return float(other) - self._value
    
    def __mul__(self, other):
        return self._value*float(other)
    
    def __div__(self, other):
        return self._value/float(other)
    
    def __truediv__(self, other):
        return self.__div__(other)
    
    def __rdiv__(self, other):
        return float(other)/self._value
    
    
class Constant(Parameter):
    pass
    
class Variable(Parameter):
    def __init__(self, name, default, bounds, scale = 'log'):
        super().__init__(name, default)
        if scale == 'log' and bounds[0] <= 0:
            raise Exception(f'Lower bound of {name} cannot be <= 0 if log scale is used.')
        if scale == 'log' and bounds[1] <= 0:
            raise Exception(f'Upper bound of {name} cannot be <= 0 if log scale is used.')
        self._low = float(bounds[0])
        self._high = float(bounds[1])
        self._scale = scale
    
    def variable(self):
        return True
    
    def __float__(self):
        return float(self.value())
    
    def value(self):
        return self._value
        
    def upper(self):
        return self._high
    
    def lower(self):
        return self._low
    
    def transformed(self):
        if self._scale == 'linear':
            return self
        elif self._scale == 'log':
            p = Variable(self._name, 
                          np.log(self._value),
                         [np.log(self._low), np.log(self._high)],
                          scale = 'linear')
            return p
        else:
            raise NotImplementedError()
    
    def set(self, value):
        if not value <= self._high or not value >= self._low:
            raise Exception(f'Value for {self._name} not within bounds.')
        self._value = value
    
    def __eq__(self, other):
        if not isinstance(other, Variable):
            return False
        val_eq = self._name == other._name and self._value == other._value
        bds_eq = self._low == other._low and self._high == other._high
        return val_eq and bds_eq
    
    def __str__(self):
        return f'{self.__class__.__name__: <9s}({self._name[:15]: <15s}, {self._value:7.2g}, ({self._low: >7.2g}, {self._high: >7.2g}), {self._scale})'
        
def default_model_parameters():
    p = [
         Variable('Hydrolysis_v_max', 1, [1e-8, 1]),
         Variable('Hydrolysis_Kmb', 800, [0.0005, 2000]),
         
         Variable('Ferm_v_max',       4, [0.001, 5]),
         #Variable('Ferm_Kmb',       890, [0.0005, 2000]),
         Variable('Ferm_Km',        780, [0.0005, 1000]),
         Variable('Ferm_inhibition',  7, [0.001, 20]),
         Variable('Ferm_CUE',        .3, [0, 1], 'linear'),
         
         Constant('death_rate',  8.3e-5),
         
         Variable('Hydro_Km_CO2',    77, [.0005, 1000]),
         Variable('Hydro_v_max',    .24, [0.003, 1.]),
         Variable('Hydro_CUE',       .3, [0, 1], 'linear'),
         Variable('Hydro_Km_H2',     77, [.0005, 1000]),
         
         Variable('Homo_Km_H2',      77, [0.0005, 1000]),
         Variable('Homo_Km_CO2',     77, [0.0005, 1000]),
         Variable('Homo_v_max',      .5, [0.005, 1.]),
         Variable('Homo_CUE',        .3, [0, 1], 'linear'),
         
         Variable('Aceto_Km_Ac',     77, [0.0005, 1000]),
         Variable('Ac_v_max',       .56, [0.005, 1.]),
         Variable('Ac_CUE',          .3, [0, 1], 'linear'), 
         
         Variable('Fe3_Km_Ac',       77, [0.0005, 1000]),
         Variable('Fe3_Km_Fe3',     738, [0.0005, 1000]),
         Variable('Fe3_v_max',    0.898, [0.002, 3.]), 
         Variable('Fe3_CUE',        0.3, [0, 1], 'linear'),
         
         Variable('Acetate',         50, [0, 100], 'linear'),
         Variable('Fe3',             20, [0, 300], 'linear'),
         
         Variable('M_Ferm',          .2, [1e-8, 0.5]),
         Variable('M_Hydro',      .0025, [1e-8, 0.5]),
         Variable('M_Homo',       .0001, [1e-8, 0.5]),
         Variable('M_Ac',         .0101, [1e-8, 0.5]),
        ]
        
    return p
    
class ModelParameters():
    def __init__(self):
        self._parameters = {}
    
    def __getitem__(self, key):
        if not key in self._parameters:
            self._parameters[key] = Parameter(key)
        return float(self._parameters.__getitem__(key))
            
    def set(self, parameters):
        if isinstance(parameters, dict):
            parameters = [Constant(k, v) for k,v in parameters.items()]
        if not isinstance(parameters, list):
            parameters = [parameters]
        for p in parameters:
            self._parameters[p._name] = p
    
    def variables(self):
        return [p for p in self._parameters.values() if p.variable()]
    
    def check(self):
        nan_pars = [v for v in self._parameters.values() if np.isnan(float(v))]
        if len(nan_pars) > 0:
            raise Exception('Model parameters are NaN:\n' + '\n'.join([str(p) for p in nan_pars]) )
    
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
        self._unbuilt_contributing_pathways = contributing_pathways
        self.contributing_pathways = None
        self.build()
        
    def build(self):
        self.contributing_pathways = [p(self.model_parameters) 
                                      for p in self._unbuilt_contributing_pathways]
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
        
        variables = [v.transformed() for v in self.model_parameters.variables()]
        variable_names = [v._name for v in variables]
        x0 = np.reshape([v.value() for v in variables], (-1,))
        lower_bound = np.reshape([v.lower() for v in variables], (-1,))
        upper_bound = np.reshape([v.upper() for v in variables], (-1,))
        
        objective = Objective([ReplicaObjective(replica, self, variable_names) 
                               for replica in replicas])
        
        loss = objective(x0)
        print(loss)
        input()
        
        import pyswarms as ps
        import OPTIMIZATION_PARAMETERS
        options = OPTIMIZATION_PARAMETERS.PSO_PARAMETERS['options']
        particles = OPTIMIZATION_PARAMETERS.PSO_PARAMETERS['particles']
        iters = OPTIMIZATION_PARAMETERS.PSO_PARAMETERS['iterations']
        
        # Call instance of GlobalBestPSO
        optimizer = ps.single.GlobalBestPSO(n_particles=particles, 
                                            dimensions=len(variables),
                                            options=options,
                                            bounds = (np.array(lower_bound), 
                                                      np.array(upper_bound)))
        _, optimal_parameter_values = optimizer.optimize(objective,
                                                        iters = iters,
                                                        n_processes = 8)
        self.model_parameters.set(variable_names, optimal_parameter_values)
        return self.model_parameters
        
        
    def predict(self, t, replica):
        print('building for prediction')
        self.build()
        S0 = system.initial_state(replica, self.parameters())
        self.parameters().check()
        self.system_state_log.reset()
        print('predicting for')
        print(str(self.parameters()))
        solver_result = scipy.integrate.solve_ivp(self, (0, max(t)),
                                                  S0, 
                                                  t_eval = t,
                                                  method = 'LSODA',
                                                  first_step = 1e-6)
        
        # TODO: appending, not overwriting!
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
        return model_string

class Objective():
    def __init__(self, replica_objectives):
        self.replica_objectives = replica_objectives
        
    def __call__(self, parameter_values):
        total_loss = sum([obj(parameter_values) for obj in self.replica_objectives])   
        return total_loss
    
class ReplicaObjective():
    def __init__(self, replica, model, parameter_names):
        self.model = model
        self.replica = replica
        self.parameter_names = parameter_names
        
    def __call__(self, parameter_values):
        self.model.parameters().set([Parameter(n,v) 
                                     for n,v in zip(self.parameter_names, 
                                                    parameter_values)])
        self.model.build()
        
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
        
    def __getitem__(self, key):
        return self._log[key]
        
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
    
    
def print_array(arr, title = '', columns = None):
    np.set_printoptions(precision = 4,
                        suppress=True)
    if not len(arr.shape) == 2:
        arr = np.reshape(arr, (-1,1))
    arr_str = np.array2string(arr).replace('[','').replace(']','')
    spl = arr_str.split('\n')
    cells = [[c.strip() for c in row.split(' ') if len(c) > 0] for row in spl  ]

    label_width = 10
    labeled_str = '\n'.join(list([f'{p:<{label_width}}' + ''.join([f'{c: >{label_width}}' for c in row])
                                  for p, row in zip(system.SYSTEM, cells)]))

    print('')
    if not title == '':
        print(title)
        print('='*len(title))

    if columns is None:
        columns = list()
    print(' '*label_width + ' ' + '     '.join([f'{c:5}' for c in columns]))

    print(labeled_str)


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
        
        CUE = microbe['CUE']
        C_source = microbe['C_source']
        self.anabolism = system.vector(0)
        if not C_source is None:
            C_atoms = chemistry.C_atoms(C_source)
            growth = system.vector(0, microbe, CUE/(1.-CUE)*C_atoms*CONSTANTS.MOLAR_MASS_C)
            self.anabolism = growth - system.vector(0, microbe['C_source'],CUE/(1.-CUE))
            
        self.microbe_index = system.index(microbe)
        self.use_thermodynamics = microbe['use_thermodynamics']
        
        if self.use_thermodynamics:
            self.deltaG_f = np.sum(np.stack([system.vector(0, str(subst), 
                                        chemistry.GIBBS_FORMATION[str(subst)])
                                        for subst in (educts + products)],
                                             axis = -1), axis = -1)
            self.deltaG_s = np.sum(self.stoichiometry*self.deltaG_f)

        
        self.state_logger = None
        
        
        print_matrix = np.concatenate([np.reshape(self.stoichiometry + self.anabolism, (-1,1)),
                                       np.reshape(HENRYS_LAW, (-1, 1)),
                                       np.reshape(self.Km, (-1, 1)),
                                       np.reshape(self.inhibition, (-1, 1)),
                                       np.reshape(self.death_rate, (-1,1))], axis = 1)

        print('')
        print('building: ', microbe['name'], 'pathway')
        print('===========' + '='*len(microbe['name']) + '========')
        print_array(print_matrix, columns = ['pathway', 'henry', 'Km', 'inhib', 'grow'])
        input()
        

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
        
        pathway_vector = self.stoichiometry + self.anabolism
    
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

        dS_dt = biomass * v * pathway_vector - biomass * self.death_rate
        dS_dt = np.clip(dS_dt, -S, np.inf)        
        
        self.log('MM', t, MM_factor)
        self.log('inhib', t, inhib_factor)
        self.log('v', t, v)
        
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
                          v_max = model_parameters['Hydrolysis_v_max'],
                          Kmb = model_parameters['Hydrolysis_Kmb'],
                          use_thermodynamics = False)
        super().__init__(microbe, educts, products)
    
class Fermentation(Pathway):
    def __init__(self, model_parameters):
        educts = [Substance(6, 'DOC', 
                            Km = model_parameters['Ferm_Km'])]
        products = [Substance(3.5, 'Acetate', 
                              inhibition = model_parameters['Ferm_inhibition']),
                    Substance(3, 'CO2'),
                    Substance(6, 'H2')]
        microbe = Microbe(name = 'M_Ferm',
                          v_max = model_parameters['Ferm_v_max'],
                          CUE = model_parameters['Ferm_CUE'],
                          death_rate = model_parameters['death_rate'],
                          C_source = 'DOC',
                          use_thermodynamics = False)
        super().__init__(microbe, educts, products)
        

class Hydro(Pathway):
     def __init__(self, model_parameters):
        educts = [Substance(4, 'H2', 
                            Km = model_parameters['Hydro_Km_H2']),
                  Substance(1, 'CO2', 
                            Km = model_parameters['Hydro_Km_CO2'])]
        products = [Substance(1, 'CH4'),
                    Substance(2, 'H2O')]
        microbe = Microbe(name = 'M_Hydro',
                          v_max = model_parameters['Hydro_v_max'],
                          CUE = model_parameters['Hydro_CUE'],
                          death_rate = model_parameters['death_rate'],
                          C_source = 'CO2')
        super().__init__(microbe, educts, products)

class Homo(Pathway):
    def __init__(self, model_parameters):
        educts = [Substance(4, 'H2', 
                            Km = model_parameters['Homo_Km_H2']),
                  Substance(2, 'CO2', 
                            Km = model_parameters['Homo_Km_CO2'])]
        products = [Substance(1, 'Acetate'),
                    Substance(2, 'H2O')]
        microbe = Microbe(name = 'M_Homo',
                          v_max = model_parameters['Homo_v_max'],
                          CUE = model_parameters['Homo_CUE'],
                          death_rate = model_parameters['death_rate'],
                          C_source = 'CO2')
        super().__init__(microbe, educts, products)

class Aceto(Pathway):
    def __init__(self, model_parameters):
        educts = [Substance(1, 'Acetate', 
                            Km = model_parameters['Aceto_Km_Ac'])]
        products = [Substance(1, 'CH4'),
                    Substance(1, 'CO2')]
        microbe = Microbe(name = 'M_Ac',
                          v_max = model_parameters['Ac_v_max'],
                          CUE = model_parameters['Ac_CUE'],
                          death_rate = model_parameters['death_rate'],
                          C_source = 'Acetate')
        super().__init__(microbe, educts, products)

class Fe3(Pathway):
    def __init__(self, model_parameters):
        educts = [Substance(1, 'Acetate', 
                            Km = model_parameters['Fe3_Km_Ac']),
                  Substance(1, 'H2O'),
                  Substance(8, 'Fe3', 
                            Km = model_parameters['Fe3_Km_Fe3'])]
        products = [Substance(8, 'Fe2'),
                    Substance(2, 'CO2')]
        microbe = Microbe(name = 'M_Fe3',
                          v_max = model_parameters['Fe3_v_max'],
                          CUE = model_parameters['Fe3_CUE'],
                          death_rate = model_parameters['death_rate'],
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
    import data
    d = data.get_data_before_carex()
        
    model_type = 'simple'
    chosen_pathways = get_pathways(model_type)
    model = Model(chosen_pathways)
    
    
    par = {
    "death_rate": 8.33e-05,
    "Acetate": 1,
    "temperature": 4.0,
    "C": 2546.5533333333337,
    "DOC": 50.93106666666667,
    "pH": 3.95,
    "weight": 11.82,
    "water": 4.0,
    "H2O": 222033.74024716797,
    "M_Fe3": 0.15377556552732402,
    "M_Ferm": 0.29386173195040044,
    "M_Hydro": 0.43281383907529236,
    "M_Homo": 0.2574949776769526,
    "Hydrolysis_v_max": 0.6166340111649226,
    "Ferm_v_max": 1.070276371909682,
    "Fe3_v_max": 1.3804244562902253,
    "Homo_v_max": 0.9318093189492231,
    "Hydro_v_max": 0.7064582813317815,
    "Ac_v_max": 0.4599047701351146,
    "Hydrolysis_Kmb": 288.99678942466437,
    "Aceto_Km_Ac": 145.99705636830586,
    "Km_Homo_CO2": 376.8013896720814,
    "Km_Homo_H2": 688.3240608121672,
    "Hydro_Km_CO2": 661.7562751340953,
    "Hydro_Km_H2": 497.8934720994153,
    "Km_Fe3_Fe3": 173.34626557916957,
    "Km_Fe3_Acetate": 637.4030208609411,
    "Ferm_Km": 160.15461453008587,
    "Ferm_inhibition": 4.643075236732733,
    "Fe3": 81.99097055433658,
    "M_Ac": 0.014042559314258995,
    "Ferm_CUE": 0.30944032735284144,
    "CUE_Fe3": 0.012291327263939777,
    "Ac_CUE": 0.5724913805190271,
    "Homo_CUE": 0.4988932700684491,
    "Hydro_CUE": 0.5054549655151662
}
    
    model.parameters().set(par)
    
    replica = d['13514']
    #mp = model.fit(replica)
    #print(mp)


    results = model.predict(range(1500), replica)
    results.plot(['CH4', 'CO2'])
    plt.show()
    
    
    
    