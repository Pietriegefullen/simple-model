import os
import numpy as np
import matplotlib.pyplot as plt
import json
from abc import ABC, abstractmethod

def default_model_parameters(normalize_parameters = False):
    p = [
         Parameter('Hydrolysis_v_max', .0083, [1e-8, 0.1], normalize = normalize_parameters),
         Parameter('Hydrolysis_Kmb', 20, [1e-10, 200], normalize = normalize_parameters),
         
         Parameter('Ferm_v_max',     .525, [0.0001, 5], normalize = normalize_parameters),
         #Variable('Ferm_Kmb',       890, [0.0005, 2000]),
         Parameter('Ferm_Km',        83, [0.0001, 100], normalize = normalize_parameters), #833
         Parameter('Ferm_inhibition', 16, [0.001, 200], normalize = normalize_parameters),
         Parameter('Ferm_CUE',        .5, [0, 1], 'linear', normalize = normalize_parameters),
         
         Parameter('death_rate',  8.3e-5, normalize = normalize_parameters),
         
         Parameter('Hydro_Km_CO2',   500, [.0005, 1000], normalize = normalize_parameters),
         Parameter('Hydro_v_max',    .17, [0.003, 1.], normalize = normalize_parameters),
         Parameter('Hydro_CUE',       .5, [0, 1], 'linear', normalize = normalize_parameters),
         Parameter('Hydro_Km_H2',    500, [.0005, 1000], normalize = normalize_parameters),
         
         Parameter('Homo_Km_H2',     500, [0.0005, 1000], normalize = normalize_parameters),
         Parameter('Homo_Km_CO2',    500, [0.0005, 1000], normalize = normalize_parameters),
         Parameter('Homo_v_max',      .5, [0.005, 1.], normalize = normalize_parameters),
         Parameter('Homo_CUE',        .5, [0, 1], 'linear', normalize = normalize_parameters),
         
         Parameter('Aceto_Km_Ac',    25, [0.5, 500], normalize = normalize_parameters), #166
         Parameter('Ac_v_max',       .0083, [0.005, 1.], normalize = normalize_parameters),
         Parameter('Ac_CUE',          .5, [0, 1], 'linear', normalize = normalize_parameters), 
         
         Parameter('Fe3_Km_Ac',      50, [0.0005, 1000], normalize = normalize_parameters),#500
         Parameter('Fe3_Km_Fe3',     500, [0.0005, 1000], normalize = normalize_parameters),
         Parameter('Fe3_v_max',      1.5, [0.2, 3.], normalize = normalize_parameters), 
         Parameter('Fe3_CUE',        0.5, [0, 1], 'linear', normalize = normalize_parameters),
         
         Parameter('Acetate',         20, [0, 100], 'linear', normalize = normalize_parameters),
         Parameter('Fe3',            150, [0, 300], 'linear', normalize = normalize_parameters),
         
         Parameter('M_Ferm',         .42, [1e-8, 5], normalize = normalize_parameters),
         Parameter('M_Hydro',       .083, [1e-8, 5], normalize = normalize_parameters),
         Parameter('M_Fe3',          .25, [1e-8, 5], normalize = normalize_parameters),
         Parameter('M_Homo',         .25, [1e-8, 5], normalize = normalize_parameters),
         Parameter('M_Ac',         .0033, [1e-8, 5], normalize = normalize_parameters),
        ]
        
    return p

class Transform(ABC):
    def __init__(self, transform = None):
        self._transf = IdentityTransform() if transform is None else transform
        
    def transform(self, value):
        value = self._transf.transform(value)
        return self._transform(value)
        
    @abstractmethod
    def _transform(self, value):
        raise NotImplementedError()
    
    def inverse(self, value):
        inv = self._inverse(value)
        return self._transf.inverse(inv)
    
    @abstractmethod
    def _inverse(self, value):
        raise NotImplementedError()
        
class IdentityTransform():
    def transform(self, value): return value
    def inverse(self, value): return value

class LogTransform(Transform):
    def _transform(self, value): return np.log(value)
    def _inverse(self, value): return np.exp(value)
    
class Normalization(Transform):
    def __init__(self, transform, low, high, target = (0,1)):
        super().__init__(transform)
        self.low = low
        self.high = high
        self.target_low = target[0]
        self.target_high = target[1]

    def _transform(self, value):
        return self.target_low + (value - self.low)/(self.high - self.low)*self.target_high
    
    def _inverse(self, value):
        return (value - self.target_low)/self.target_high*(self.high - self.low) + self.low

class ModelParameters():
    def __init__(self, d = None):
        if isinstance(d, dict):
            d =  {k:(v if isinstance(v, Parameter) else Parameter(k,v))
                for k,v in d.items()}
        
        elif d is None:
            d = {}
        else:
            raise NotImplementedError()
        self._parameters = d
    
    def search_space(self):
        default_ranges = default_model_parameters()
        reductions = []
        for p in self.variables():
            ind = next((i for i, dp in enumerate(default_ranges) if dp.name == p.name), None)
            p_def = default_ranges[ind]
            p_range = p.transform(p.high) - p.transform(p.low)
            p_default_range = p.transform(p_def.high) - p.transform(p_def.low)
            w_i = p_range/p_default_range
            if w_i == 0:
                print(p.name, p.scale, p.transformer, p_range, p_default_range, w_i)
            reductions.append(w_i)

        volume = np.prod(reductions)
        dimension = len(reductions)
        average_reduction = np.power(volume, 1/dimension)            
        return dimension, volume, average_reduction
    
    def as_dict(self):
        return self._parameters
    
    def __iter__(self):
        return iter(self._parameters.values())
    
    def __getitem__(self, key):
        if not key in self._parameters:
            self._parameters[key] = Parameter(key)
        return self._parameters[key]
            
    def set(self, parameters, normalized = False):
        if isinstance(parameters, str) and parameters == 'default':
            self.set(default_model_parameters(normalize_parameters = normalized))
            
        elif isinstance(parameters, dict):
            for p, value in parameters.items():
                if not p in self._parameters:
                    print('Ignoring parameter that does not exist in model:', p)
                else:
                    try:
                        self._parameters[p].set(value)
                    except Exception as ex:
                        raise Exception(parameters)
                    
        elif isinstance(parameters, list):
            for p in parameters:
                if p.name in self._parameters:
                    self._parameters[p.name].set(p)
                else:
                    print('Ignoring parameter that does not exist in model:', p)

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


class Parameter():
    def __init__(self, name, value = np.nan, range = None, scale = 'log', normalize = False):
        self.name = name
        self.value = value
        
        self.low = None
        self.high = None
        self.scale = scale
        self.normalize = normalize
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
            self.normalize = p.normalize
        elif isinstance(p, (int, float)):
            if self.is_variable() and (not p <= self.high or not p >= self.low):
                print(f'WARNING: Setting {self.name} to {p} is out of bounds.')
            self.value = float(p)
        else:
            raise NotImplementedError(str(p))
    
    def is_unset(self):
        return np.isnan(self.value)
                        
    def is_variable(self):
        if self.low == self.high:
            return False
        return not self.is_unset() and not self.low is None and not self.high is None
    
    def get_transform(self):
        if self.transformer is None:
            if self.scale == 'linear':
                if self.is_variable() and self.normalize:
                    tf = Normalization(None, self.lower(), self.upper())
                else:
                    tf = IdentityTransform()
                self.transformer = tf
            
            elif self.is_variable() and self.scale == 'log':
                tf = LogTransform()
                if self.is_variable() and self.normalize:
                    tf = Normalization(tf, tf.transform(self.lower()), tf.transform(self.upper()))
                self.transformer = tf 
            
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
            var = f'  ({self.low:8.3g}, {self.high:8.3g})   {self.scale}   norm.: {self.normalize}'
        return f'{self.name[:15]:15} = {self.value:8.3g}' + var
    
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
    
def boxplots(loaded_parameters, save_target = None):
    parameter_names = list(loaded_parameters['simple'].keys()) + list(loaded_parameters['complex'].keys()) 
    parameter_names = list(set(parameter_names)) # unique names

    parameter_groups = {'pools': [],
                        'microbes': [],
                        'CUE': [],
                        'Kmb': [],
                        'Km': [],
                        'v_max': []}
    for p in parameter_names:
        if p == 'death_rate': 
            continue
        if not '_' in p:
            parameter_groups['pools'].append(p)
        elif p.startswith('M_'):
            parameter_groups['microbes'].append(p)
        elif 'CUE' in p:
            parameter_groups['CUE'].append(p)
        elif 'Kmb' in p:
            parameter_groups['Kmb'].append(p)
        elif 'Km' in p:
            parameter_groups['Km'].append(p)
        elif 'v_max' in p:
            parameter_groups['v_max'].append(p)
            
    d = .25
    for group_name, group in parameter_groups.items():
        plt.figure()
        i = 0
        plt.title(group_name)
        x_tick_labels = []
        x_tick_pos = []
        for par_name in group:
            pos_simple = 2*i + 1.5 - d
            pos_complex = 2*i + 1.5 + d
            tick_pos = 2*i + 1.5
            i += 1
            
            parameter_name = par_name.replace(group_name, '').replace('_', ' ')
            x_tick_labels.append(parameter_name)
            x_tick_pos.append(tick_pos)
            
            simple_data = [np.nan]
            if par_name in loaded_parameters['simple']:
                simple_data = loaded_parameters['simple'][par_name]
            complex_data = loaded_parameters['complex'][par_name]

            print(f'found {len(simple_data):2d} parameter values for {par_name} (simple)')
            print(f'found {len(complex_data):2d} parameter values for {par_name} (complex)')
            box_data = [simple_data, complex_data]
            plt.boxplot(box_data, 
                        positions = [pos_simple, pos_complex],
                        widths = 2*d)
            
        plt.xticks(rotation=90)
        plt.xticks(x_tick_pos, x_tick_labels)
        if not group_name == 'CUE':# and not group_name == 'pools':
            plt.yscale('log')

        if not save_target is None:
            figure_name = group_name
            plt.savefig(os.path.join(save_target, figure_name) + '.svg')

    if save_target is None:
        plt.show()

def sort_dict(d):
    s =  {}
    for k in sorted(d.keys()):
        ss = d[k]
        if isinstance(ss, dict):
            ss = sort_dict(ss)
        s[k] = ss
    return s


def R2plot():
    import model
    import data
    
    r2_values = {}

    kd = data.get_data_before_day()
    _, found = load_best()
    total = len([r for s in found.values() for r in s.values()])
    counter = 0
    for sample_name, repl_dict in found.items():
        for repl, (best_loss, best_parameters) in repl_dict.items():
            counter += 1
            print(f'{counter/total*100:.2f}%')
            model_type = 'simple' if 'simple' in repl else 'complex'
            loaded_model = model.Model(model.get_pathways(model_type))
            loaded_model.parameters().set(best_parameters)
  
            sample = kd[sample_name]
            val_replica = None
            for s in sample.leave_one_out_split():
                fit_repl = '/'.join(sorted([r.replica_number for r in s['fit']]))
                if fit_repl in repl:
                    val_replica = s['val']
                    break
            if val_replica is None:
                print('could not identify validation replica for', sample_name, repl)
                continue

            model_run = loaded_model.predict(val_replica)

            co2_r2 = model_run['R2']['CO2']
            ch4_r2 = model_run['R2']['CH4']
            r2_values[sample_name + '/' + val_replica.replica_number] = (co2_r2, ch4_r2)

            #plt.figure()
            #model_run.plot(['CO2', 'CH4'], newfigure = False)
            #val_replica.plot(log = False, newfigure = False)
            #ax = plt.gca()
            #t = ax.get_title()
            #plt.title(str(val_replica) + ' validation' )
            #plt.show()

    plt.figure()
    all_r2 = np.array([v for v in r2_values.values()])
    plt.plot(all_r2[:,0], all_r2[:,1], 'k.')
    plt.xlabel('R2 CO2')
    plt.ylabel('R2 CH4')
    plt.ylim([0,1])
    plt.xlim([0,1])
    plt.show()
    
def load_best():
    import os
    import model
    import USER_VARIABLES
    all_parameters = {  'simple': {},
                        'complex': {}}
    found = {}
    result_source = USER_VARIABLES.LOG_DIRECTORY
    for f in os.listdir(result_source):
        parameter_source = os.path.join(result_source, f)
        if not os.path.isdir(parameter_source) or not f.startswith('fit'):
            continue

        if not 'log' in f:
            continue

        try:
            best_loss, best_parameters = model.get_best_loss_parameters(parameter_source)
            model_type = 'complex' if 'complex' in f else 'simple'
        
            replicas = [s for s in f.split('log')[0].replace('fit_','').split('_') if not s == '']
            sample = replicas[0][:4]
            repl = '/'.join([r[-1] for r in sorted(replicas)]) + ' ' + model_type
            if not sample in found:
                found[sample] = {}
            if not repl in found[sample]:
                found[sample][repl] = None

    
            if found[sample][repl] is None or best_loss < found[sample][repl][0]:
                found[sample][repl] = (best_loss, best_parameters)
            
            else:
                continue
                
        except:
            print('no parameters found in', parameter_source)
            continue

        for k, p in best_parameters.items():
            if not k in all_parameters[model_type]:
                all_parameters[model_type][k] = []
            all_parameters[model_type][k].append(p)

    found = sort_dict(found)

    return all_parameters, found

if __name__ == '__main__':

    R2plot()
    1/0
    all_parameters, found = load_best()

    boxplots(all_parameters)

