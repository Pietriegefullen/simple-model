import numpy as np
import scipy.integrate
import matplotlib.pyplot as plt

import system
import chemistry
import CONSTANTS

HENRYS_LAW = system.henrys_law()

# TODO: parallel coordinates plot for optimization runs
# TODO: throw unused microbes out of the system!
# TODO: model parameters: setting during optimization (keep variable)
#       setting after build.
# TODO: structure run results larger than compare_models!!!
# TODO: shorter, cleaner output
# TODO: PSO might get stuck for some parameters b/c integration takes forever?
# TODO: encapsulate optimization.

class Parameter():
    def __init__(self, name, value = np.nan, range = None, scale = 'log'):
        self.name = name
        self.value = value
        
        self.low = None
        self.high = None
        self.scale = scale
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
    
    def transform(self):
        if self.scale == 'linear':
            return IdentityTransform()
        
        elif self.is_variable() and self.scale == 'log':
            return LogTransform()
        
        else:
            raise NotImplementedError()
        
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


class Model():
    def __init__(self, contributing_pathways):
        self.system_state_log = ModelRun()
        self.model_parameters = ModelParameters()
        self._unbuilt_contributing_pathways = contributing_pathways
        self.contributing_pathways = None
        self.build(quiet = True)
        
    def build(self, quiet = False):
        self.contributing_pathways = [p(self.model_parameters) 
                                      for p in self._unbuilt_contributing_pathways]
        for p in self.contributing_pathways:
            p.inject_logger(self.system_state_log)
        
        # to initialize model parameters used in initial state
        _ = system.initial_state(None, self.model_parameters)
        
        if not quiet:
            for p in self.contributing_pathways:
                #print('')
                print(p.__class__.__name__)
                #print('='*len(p.__class__.__name__))
                #print_matrix = np.concatenate([np.reshape(p.pathway_vector, (-1,1)),
                #                       np.reshape(HENRYS_LAW, (-1, 1)),
                #                       np.reshape(p.Km, (-1, 1)),
                #                       np.reshape(p.inhibition, (-1, 1)),
                #                       np.reshape(p.death_rate, (-1,1))], axis = 1)
#
#                print_array(print_matrix, columns = ['pathway', 'henry', 'Km', 'inhib', 'grow'])        

        
    def __call__(self, t, S):
        if False:
            print(f't = {t}')
            for v, n in zip(S, system.SYSTEM):
                print(f'{n[:15]:15}  {v}')
            input()
        S = np.where(S < 1e-40, 0, S)
        
        dSj_dt = np.stack([pathway_j(t, S)
                            for pathway_j in self.contributing_pathways], axis = -1)
        
        dS_dt = np.sum(dSj_dt, axis = -1)
        dS_dt = np.clip(dS_dt, -S, np.inf) # don't let pools become negative
        return dS_dt
    
    def fit(self, replicas):
        if not isinstance(replicas, list):
            replicas = [replicas]
        
        variables = self.model_parameters.variables()

        variable_names = [v.name for v in variables]
        x0 = np.reshape([v.transform().transform(v.value) 
                         for v in variables], (-1,))
        lower_bound = np.reshape([v.transform().transform(v.lower()) 
                                  for v in variables], (-1,))
        upper_bound = np.reshape([v.transform().transform(v.upper())
                                  for v in variables], (-1,))
        objective = Objective([ReplicaObjective(replica, self, variable_names) 
                               for replica in replicas],
                              [v.transform() for v in variables])
        
        import OPTIMIZATION_PARAMETERS
        method_name = OPTIMIZATION_PARAMETERS.GRADIENT_PARAMETERS['method']
        
        initial_guess_bounds = list(zip(lower_bound, upper_bound))
        _ = scipy.optimize.minimize(objective,
                                    x0,
                                    bounds = initial_guess_bounds,
                                    **OPTIMIZATION_PARAMETERS.GRADIENT_PARAMETERS)

        
        #import pyswarms as ps
        #import OPTIMIZATION_PARAMETERS
        #options = OPTIMIZATION_PARAMETERS.PSO_PARAMETERS['options']
        #particles = OPTIMIZATION_PARAMETERS.PSO_PARAMETERS['particles']
        #iters = OPTIMIZATION_PARAMETERS.PSO_PARAMETERS['iterations']
                
        ## Call instance of GlobalBestPSO
        #optimizer = ps.single.GlobalBestPSO(n_particles=particles, 
        #                                    dimensions=len(variables),
        #                                    options=options,
        #                                    bounds = (np.array(lower_bound), 
        #                                              np.array(upper_bound)))
        #optimizer.optimize(objective, iters = iters, n_processes = 8)

        optimal_parameter_values = objective.best_call()
        
        self.model_parameters.set(dict(zip(variable_names, optimal_parameter_values)))
        return self.model_parameters
        
        
    def predict(self, t, replica, quiet = False):
        self.build(quiet = quiet)
        S0 = system.initial_state(replica, self.parameters())
        self.parameters().check()
        self.system_state_log.reset()
        if not quiet:
            print('predicting for')
            print(str(self.parameters()))
            print(self)
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

class Objective():
    def __init__(self, replica_objectives, transforms):
        self.replica_objectives = replica_objectives
        self._calls = []
        self.transforms = transforms
    
    def transform(self, parameter_values):
        transformed_parameter_values = [t.transform(v) 
                                        for t, v in zip(self.transforms, parameter_values)]
        return np.squeeze(transformed_parameter_values)
    
    def inverse_transform(self, transformed_parameter_values):
        try:
            parameter_values = [t.inverse(v) 
                            for t, v in zip(self.transforms, transformed_parameter_values)]
        except Exception as ex:
            raise Exception(str(transformed_parameter_values))
        return np.squeeze(parameter_values)
    
    def __call__(self, transformed_parameter_values):
        parameter_values = self.inverse_transform(transformed_parameter_values)
        total_loss = sum([obj(parameter_values) 
                          for obj in self.replica_objectives])   
        self._calls.append((total_loss, parameter_values))
        print(total_loss)
        return total_loss
    
    def best_call(self):
        sorted_by_loss = sorted(self._calls)
        best = sorted_by_loss[0]
        best_loss, best_transformed_parameter_values = best
        best_parameter_values = self.inverse_transform(best_transformed_parameter_values)
        return best_parameter_values
    
    def __str__(self):
        return 'Objective function: sum of loss from\n' + '\n'.join([str(s) 
                                                         for s in self.replica_objectives])

class ParticleObjective(Objective):
        
    def __call__(self, particle_model_parameter_values):
            particle_fitnesses = []
            for transformed_parameter_values in particle_model_parameter_values:
                parameter_values = self.inverse_transform(transformed_parameter_values)
                losses = [obj(parameter_values) 
                          for obj in self.replica_objectives]
                total_loss = np.sum(losses)
                particle_fitnesses.append(total_loss)
            self._calls.append((total_loss, particle_fitnesses))
            return np.array(particle_fitnesses)

class Loss():
    def __init__(self, predicted, measured):
        self.predicted = predicted
        self.measured = measured
    
    def RMSE(self):
        return np.sqrt(np.mean((self.predicted - self.measured)**2))

class ReplicaObjective():
    def __init__(self, replica, model, parameter_names):
        self.model = model
        self.replica = replica
        self.parameter_names = parameter_names
        
    def __call__(self, parameter_values):
        print('called', parameter_values)
        
        self.model.parameters().set({n:v for n,v in zip(self.parameter_names, 
                                                        np.squeeze(parameter_values))})        
        days = self.replica.incubation['days']
        results = self.model.predict(days, self.replica, quiet = True)
        
        _, predicted_CO2 = zip(*results['CO2'])
        measured_CO2 = self.replica.incubation['CO2']
        
        CO2_loss = Loss(predicted_CO2, measured_CO2).RMSE()
        
        _, predicted_CH4 = zip(*results['CH4'])
        measured_CH4 = self.replica.incubation['CH4']
        
        CH4_loss = Loss(predicted_CH4, measured_CH4).RMSE()
        
        loss = CO2_loss + CH4_loss
        
        plt.figure()
        plt.plot(self.replica['days'], measured_CO2, 'rx')
        plt.plot(self.replica['days'], predicted_CO2, 'k-')
        plt.show()
        
        return loss
    
    def __str__(self):
        return f'fit to replica {self.replica}'
    
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
    
        #pwy_string += '\n' + pathway_formatter(self.microbe, self.educts, self.products)
    
        return pwy_string
    
def pathway_formatter(microbe, educts, products):
    pwy = ''
    pwy += 'parameters for' + microbe['name'] + 'pathway' + '\n'
    pwy += '===============' + '='*len(microbe['name']) + '=========' + '\n'
    for k,v in microbe.get_config().items():
        if k == 'name': continue
        pwy += f'   {k:10} {v}' + '\n'


    pwy += 'educts:' + '\n'
    for i,educt in enumerate(educts):
        cnt = f'{i+1:2d})'
        for k,v in educt.get_config().items():
            if not k == 'name':
                cnt = '   '
            pwy += f'{cnt}   {k:10} {v}' + '\n'


    pwy += 'products:' + '\n'
    for i,product in enumerate(products):
        cnt = f'{i+1:2d})'
        for k,v in product.get_config().items():
            if not k == 'name':
                cnt = '   '
            pwy += f'{cnt}   {k:10} {v}' + '\n'
            
    return pwy
    
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
    
    def get_config(self):
        return {'name': self.name,
                'death_rate': self.death_rate,
                'v_max':self.v_max,
                'Kmb': self.Kmb,
                'CUE': self.CUE,
                'C_source': self.C_source,
                'use_thermodynamics': self.use_thermodynamics}

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

    def get_config(self):
        return {'name': self.name,
                'stoichiometry': self.stoichiometry,
                'Km': self.Km,
                'inhibition': self.inhibition}


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
        return basic + [Homo,
                        Fe3]
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
    #"Acetate": 1,
    #"temperature": 4.0,
    #"C": 2546.5533333333337,
    #"DOC": 50.93106666666667,
    #"pH": 3.95,
    #"weight": 11.82,
    #"water": 4.0,
    #"H2O": 222033.74024716797,
    #"M_Fe3": 0.15377556552732402,
    #"M_Ferm": 0.29386173195040044,
    #"M_Hydro": 0.43281383907529236,
    #"M_Homo": 0.2574949776769526,
    "Hydrolysis_v_max": 0.6166340111649226,
    "Ferm_v_max": 1.070276371909682,
    #"Fe3_v_max": 1.3804244562902253,
    #"Homo_v_max": 0.9318093189492231,
    "Hydro_v_max": 0.7064582813317815,
    "Ac_v_max": 0.4599047701351146,
    "Hydrolysis_Kmb": 288.99678942466437,
    "Aceto_Km_Ac": 145.99705636830586,
    #"Km_Homo_CO2": 376.8013896720814,
    #"Km_Homo_H2": 688.3240608121672,
    "Hydro_Km_CO2": 661.7562751340953,
    "Hydro_Km_H2": 497.8934720994153,
    #"Km_Fe3_Fe3": 173.34626557916957,
    #"Km_Fe3_Acetate": 637.4030208609411,
    "Ferm_Km": 160.15461453008587,
    "Ferm_inhibition": 4.643075236732733,
    #"Fe3": 81.99097055433658,
    #"M_Ac": 0.014042559314258995,
    "Ferm_CUE": 0.30944032735284144,
    #"CUE_Fe3": 0.012291327263939777,
    "Ac_CUE": 0.5724913805190271,
    #"Homo_CUE": 0.4988932700684491,
    "Hydro_CUE": 0.5054549655151662
}
 
    model.parameters().set('default')
    model.parameters().set(par)

    replica = d['13514']
    #mp = model.fit(replica)
    
    results = model.predict(range(1500), replica)
    #print(results['CO2'][-1])
    results.plot(['CH4', 'CO2'], newfigure = False)
    plt.show()
    
    
    
    