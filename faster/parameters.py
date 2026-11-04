import os
import numpy as np
import matplotlib.pyplot as plt
import json
from abc import ABC, abstractmethod
import hashing
import checkpoint
import system
import traceback
import USER_VARIABLES    


BOOLEAN_TYPES = (bool, np.bool_)

def load_parameters(init_config, return_run_config = False, source_directory = None):
    print('determine initial parameters')
    init_config = {k: v for k,v in init_config.items() if not v is None}

    if 'parameters' in init_config:
        initial_parameters = init_config['parameters']
        
        if 'range' in init_config:
            parameter_range = init_config['range']
            _ = [parameter_range[p.name].set(p) for p in initial_parameters]
            initial_parameters = parameter_range
        
        return initial_parameters
    
    if 'default' in init_config and init_config['default'] or len(init_config) == 0 or (len(init_config) == 1 and 'normalized' in init_config):
        normalized = init_config['normalized'] if 'normalized' in init_config else False
        default_range = ModelParameters()
        return default_range

    if 'estimate' in init_config and init_config['estimate']:
        normalized = init_config['normalized'] if 'normalized' in init_config else False
        estimated_range = ModelParameters({
            p.name: p for p in estimated_model_parameters(normalized)})
        return estimated_range

    if 'file' in init_config and not init_config['file'] is None:
        loaded_parameters, loss, run_config = load_parameter_file(init_config['file'])
        if 'range' in init_config and init_config['range'] == 'default':
            default_range = ModelParameters()
            _ = [default_range[p.name].set(p.value) for p in loaded_parameters
                 if p in default_range]
            for p in loaded_parameters:
                if not p.name in default_range:
                    default_range[p.name].set(p)
            loaded_parameters = default_range

        elif len(init_config) == 1:
            pass
        else:
            raise NotImplementedError()
        if return_run_config:
            return loaded_parameters, run_config
        return loaded_parameters

    print('loading existing checkpoints')
    candidate_files = get_candidates(init_config, source_directory = source_directory)
    selected_candidates = select_candidates(candidate_files, init_config)
    pl = '' if len(selected_candidates) == 1 else 's'
    print(f'selected {len(selected_candidates)} candidate{pl}')
    if len(selected_candidates) == 0:
        raise Exception()
    
    loaded_range, _, _, file_path = selected_candidates[0]
    if len(selected_candidates) > 1:
        selected_parameters, _, _, file_path = zip(*selected_candidates[1:])
        for par, pth in zip(selected_parameters, file_path):
            #print(os.path.split(pth)[-1])
            loaded_range.extend_range_by_value(par, ignore_constants = True)
    return loaded_range

def load_parameter_file(file_path):
    exs = []
    for attempt in [_load_annotated_parameter_file]:
        try:
            parameters = attempt(file_path)
            return parameters
        except Exception as ex:
            #exs.append((ex, file_path))
            exs.append((ex, traceback.format_exc()))
            continue

    for x in exs:
        print(x[1])
        print('====')
    raise Exception('Could not load parameters from checkpoint file. ' + os.path.split(file_path)[-1])

def load_file(file_path):
    with open(file_path, 'r') as pf:
        file_content = json.load(pf)
    return file_content


def set_initial_pool(parameters):
    pool_names = system.SYSTEM
    for pool in system.SYSTEM:
        if not pool in parameters:
            parameters[pool].constant(0)
            

def _load_annotated_parameter_file(file_path):
    file_content = load_file(file_path)
    assert 'parameters' in file_content
    assert 'total_loss' in file_content
    assert 'run_config' in file_content
    parameters = ModelParameters(file_content['parameters'])
    loss = float(file_content['total_loss'])
    return parameters, loss, file_content['run_config']

def get_candidates(init_config, source_directory = None):
    if source_directory is None:
        source_directory = USER_VARIABLES.RESULTS_DIRECTORY
    checks = []
    def check(criterion, value, mode = 'eq'):
        def check_function(run_config):
            run_config_value = checkpoint.get_value(run_config, criterion)
            if mode == 'eq':
                return str(run_config_value) == str(value)
            elif mode == 'in':
                return str(value) in str(run_config_value)
        check_function.__name__ = 'check_' + str(criterion) + '_' + str(value)
        return check_function
    
    if 'replica' in init_config:
        replica = str(init_config['replica'])
        replica = [c for c in replica if c in '0123456789' ]
        assert len(replica) == 5
        init_config['sample'] = int(replica[:4])
        init_config['validation_replica'] = int(replica[-1:])
        del init_config['replica']

    if 'run_ID' in init_config:
        checks.append(check(('run', hashing.build_run_id), {'run', init_config}))

    if 'model' in init_config:
        checks.append(check(('model',hashing.build_model_id), init_config['model']))

    if 'sample' in init_config:
        checks.append(check(('chosen', 'sample'), init_config['sample']))
        
    if 'validation_replica' in init_config:
        checks.append(check(('chosen', 'validation_replica'), init_config['validation_replica']))
        
    if 'source' in init_config:
        checks.append(check('source', init_config['source'], 'in'))

    candidates = []
    for root, dirs, files in os.walk(source_directory):
        for file in files:
            file_path = os.path.join(root, file)
            if any([s.startswith('.') for s in file_path.split(os.sep)]):
                continue
            try:
                checkpoint_data = load_file(file_path)
            except (json.JSONDecodeError, UnicodeDecodeError):
                continue
            check_data = checkpoint_data.copy()
            if not 'run_config' in check_data:
                continue
            check_data['run_config']['source'] = file_path
            if not all([ c(check_data['run_config']) for c in checks]): 
                continue
            
            parameters = ModelParameters(checkpoint_data['parameters'])
            
            loss = float(checkpoint_data['total_loss'])
            candidates.append((parameters, loss, checkpoint_data['run_config'], file_path))
    return candidates

def select_candidates(candidates, init_config):
    if 'best_N' in init_config and not init_config['best_N'] is None:
        candidates = sorted(candidates, key = lambda c: c[1])
        return candidates[:int(init_config['best_N'])]
    
    if 'worst' in init_config:
        candidates = sorted(candidates, key = lambda c: -c[1])
        return candidates[:int(init_config['worst'])]
    
    return candidates

def provided_parameters(replica):
    import system
    p = ModelParameters()
    _ = system.initial_state(replica, p)
    return p

class Transform(ABC):
    def __init__(self, transform = None):
        self._transf = IdentityTransform() if transform in ('linear', None) else transform
    
    def __call__(self, value):
        return self.transform(value)
        
    def transform(self, value):
        value = self._transf.transform(value)
        value = self._transform(value)
        return value
        
    @abstractmethod
    def _transform(self, value):
        raise NotImplementedError()
    
    def inverse(self, value):
        inv = self._inverse(value)
        return self._transf.inverse(inv)
    
    @abstractmethod
    def _inverse(self, value):
        raise NotImplementedError()
    
    def __str__(self):
        return self.__class__.__name__
    
    def short(self):
        return self.__class__.__name__
            
    def operator(self):
        l, r = self._transf.operator()
        return self.short() + '(' + l, r + ')'
    
class IdentityTransform():
    def __call__(self, value):
        return self.transform(value)
    def transform(self, value): return value
    def inverse(self, value): return value
    def operator(self): return '', ''
    def __str__(self): return ''

class LogTransform(Transform):
    def __init__(self, transform = None):
        super().__init__(transform)
        
    def _transform(self, value):
        with np.errstate(invalid='ignore', divide = 'ignore'):
            result = np.log(value)
        return result
    
    def _inverse(self, value): 
        return np.exp(value)
    
    def short(self):
        return 'log'
    
class Normalization(Transform):
    def __init__(self, transform, low, high, target = (0,1)):
        super().__init__(transform)
        self.low = low
        self.high = high
        self.target_low = target[0]
        self.target_high = target[1]

    def _transform(self, value):
        target_range = self.target_high - self.target_low
        result = self.target_low + (value - self.low)/(self.high - self.low)*target_range
        return result
    
    def _inverse(self, value):
        target_range = self.target_high - self.target_low
        return (value - self.target_low)/target_range*(self.high - self.low) + self.low
    
    def short(self):
        return 'norm'
    
    def operator(self):
        l, r, = super().operator()
        r = f'),{self.low:.2g}, {self.high:.2g}' + r[1:]
        return l, r

class ModelParameters():
    def __init__(self, d = None):
        if isinstance(d, ModelParameters):
            d = d.as_dict()

        if isinstance(d, list):
            d = {p.name: p for p in d}

        if isinstance(d, dict):
            parameter_values = d
            d = {}
            for k, v in parameter_values.items():
                if isinstance(v, Parameter):
                    p = v
                elif isinstance(v, dict):
                    p = Parameter.from_config(v)
                elif isinstance(v, (int, float) + BOOLEAN_TYPES):
                    p = Parameter(k, v)
                else:
                    raise TypeError(f'Unsupported value for parameter {k}: {type(v).__name__}')
                d[k] = p
        
        elif d is None:
            d = {}

        else:
            print(d)
            raise NotImplementedError()
        self._parameters = d
    
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
            raise NotImplementedError()
            
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
                    print('Ignoring parameter that does not exist in model:', p.name)

        elif isinstance(parameters, tuple):
            name, value = parameters
            self.set({name:value})
            
        elif isinstance(parameters, ModelParameters):
            self.set([p for p in parameters])
            
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
    
    def get_values(self):
        return {p.name: p.value for p in self._parameters.values()}
    
    def get_config(self, only_range = False):
        return {p.name: p.get_config(only_range = only_range) 
                for p in self._parameters.values()}
    
    def __str__(self):
        title = 'Model Parameters:'
        title += '\n' + '='*len(title) + '\n'
        sorted_params = sorted(self._parameters.values(), key = lambda x: x.name)
        s = title + '\n'.join([f'{i+1:3d}) ' + str(p) 
                                  for i,p in enumerate(sorted_params)])
        return s

    def extend_range_by_value(self, parameters, ignore_constants = True):
        if isinstance(parameters, Parameter):
            parameters = [parameters]
        for p in parameters:
            self[p.name].extend_range_by_value(p, ignore_constants)

class Parameter():
    def __init__(self, name, value = np.nan, range = None, scale = None, normalize = False):
        self.name = name
        self._value = value
        
        self.options = None
        self.low = None
        self.high = None
        
        if isinstance(range, set):
            self.options = range
        
        elif isinstance(range, (tuple, list)):
            self.low = range[0]
            self.high = range[1]

        # Standalone parameters (for example in tests and ad-hoc tools) do
        # not belong to a versioned model catalogue.  Preserve the historic
        # logarithmic default for ranged numeric parameters; versioned model
        # definitions always pass their scale explicitly.
        if not range is None and scale is None:
            scale = 'log'

        self.scale = scale
        self.normalize = normalize
        self.transformer = None
        
    @property
    def value(self):
        if not isinstance(self._value, (type(None), int, float) + BOOLEAN_TYPES):
            raise ValueError()
        return self._value
        
    @value.setter
    def value(self, value):
        if not isinstance(value, (int, float) + BOOLEAN_TYPES):
            raise ValueError()
        self._value = value

    @classmethod
    def from_config(cls, cfg):
        if not 'value' in cfg:
            cfg['value'] = np.nan
        p = Parameter(cfg['name'], cfg['value'], cfg['range'], cfg['scale'], cfg['normalize'])
        return p
    
    def get_config(self, only_range = False):
        rng = None
        if self.is_variable():
            if not self.options is None:
                rng = frozenset(self.options)
            else:
                rng = (self.low, self.high)
        cfg = {
                'name': self.name,
                'range': rng,
                'scale': self.scale,
                'normalize': self.normalize,
                    }
        if not only_range:
            cfg['value'] = self.value
        return cfg
    
    def lower(self):
        if not self.is_variable():
            raise Exception()
        if self.low is None:
            return np.min(self.options)
        return self.low
    
    def upper(self):
        if not self.is_variable():
            raise Exception()
        if self.high is None:
            return np.max(self.options)
        return self.high
    
    def constant(self, value):
        self.value = value
        self.low = None
        self.high = None
        self.options = None
        return self
    
    def variable(self, value, range):
        self.value = value
        if isinstance(range, set):
            self.low = None
            self.high = None
            self.options = range
        elif isinstance(range, tuple):
            self.low = range[0]
            self.high = range[1]
            self.options = None
        return self
    
    def make_unset(self):
        self.value = np.nan

    def set(self, p):
        if isinstance(p, Parameter):
            self.value = p.value
            self.scale = p.scale
            self.normalize = p.normalize
            self.low = p.low
            self.high = p.high
            self.options = p.options

        elif isinstance(p, (int, float) + BOOLEAN_TYPES):
            if self.is_variable() and np.isfinite(p):
                out_of_range = False
                if self.options is None:
                    # An exact endpoint in transformed coordinates can map
                    # back through exp/log normalisation a few ULPs outside
                    # its stored physical bound.  Canonicalise only that
                    # round-off spillover; do not hide a real violation.
                    tolerance = (64 * np.finfo(float).eps
                                 * max(abs(self.low), abs(self.high), np.finfo(float).tiny))
                    if p < self.low and self.low - p <= tolerance:
                        p = self.low
                    elif p > self.high and p - self.high <= tolerance:
                        p = self.high
                    elif p < self.low or p > self.high:
                        out_of_range = True
                elif isinstance(self.options, set) and not p in self.options:
                    out_of_range = True

                if out_of_range:
                    print(f'WARNING: Setting {self.name} to {p} is out of bounds.')
            self.value = float(p)

        else:
            raise NotImplementedError(str(p))

    def extend_range_by_value(self, p, ignore_constants = False):
        if self.is_variable() or ignore_constants:
            value = p.value
            
            if not self.options is None:
                self.options.add(value)
            
            else:
                if self.high is None:
                    self.high = self.value
                    
                if self.low is None:
                    self.low = self.value
                
                
                if value > self.high:
                    self.high = value
        
                if value < self.low:
                    self.low = value

    def is_unset(self):
        return self.value == float('nan') or self.value is None or np.isnan(self.value)

    def is_variable(self):
        if not self.options is None and len(self.options) > 1:
            return True

        elif isinstance(self.options, set) and len(self.options) <= 1:
            return False

        if (not self.low is None) and np.isclose(self.low,self.high):
            return False

        return not self.is_unset() and not self.low is None and not self.high is None
    
    def is_constant(self):
        return not self.is_variable() and not self.is_unset()
    
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
            var = f'  ({self.low:8.3g}, {self.high:8.3g})   {self.scale:6s}   norm.: {self.normalize}'
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

