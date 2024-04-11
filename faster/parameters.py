import numpy as np

def default_model_parameters(model_parameters = None):
    p = [
         Parameter('DOC_per_TOC',    6.),

         Parameter('Hydrolysis_v_max', .83, [1e-8, 1]),
         Parameter('Hydrolysis_Kmb', 333, [0.0005, 2000]),
         
         Parameter('Ferm_v_max',     2.5, [0.001, 5]),
         #Variable('Ferm_Kmb',       890, [0.0005, 2000]),
         Parameter('Ferm_Km',        833, [0.0005, 1000]),
         Parameter('Ferm_inhibition', 16, [0.001, 20]),
         Parameter('Ferm_CUE',        .16, [0, 1], 'linear'),
         
         Parameter('death_rate',  8.3e-5),
         
         Parameter('Hydro_Km_CO2',   833, [.0005, 1000]),
         Parameter('Hydro_v_max',    .17, [0.003, 1.]),
         Parameter('Hydro_CUE',       .16, [0, 1], 'linear'),
         Parameter('Hydro_Km_H2',     833, [.0005, 1000]),
         
         Parameter('Homo_Km_H2',      500, [0.0005, 1000]),
         Parameter('Homo_Km_CO2',     500, [0.0005, 1000]),
         Parameter('Homo_v_max',      .5, [0.005, 1.]),
         Parameter('Homo_CUE',        .5, [0, 1], 'linear'),
         
         Parameter('Aceto_Km_Ac',     166, [0.0005, 1000]),
         Parameter('Ac_v_max',       .83, [0.005, 1.]),
         Parameter('Ac_CUE',          .5, [0, 1], 'linear'), 
         
         Parameter('Fe3_Km_Ac',      500, [0.0005, 1000]),
         Parameter('Fe3_Km_Fe3',     500, [0.0005, 1000]),
         Parameter('Fe3_v_max',      1.5, [0.002, 3.]), 
         Parameter('Fe3_CUE',        0.5, [0, 1], 'linear'),
         
         Parameter('Acetate',         1, [0, 100], 'linear'),
         Parameter('Fe3',            150, [0, 300], 'linear'),
         
         Parameter('M_Ferm',         .42, [1e-8, 0.5]),
         Parameter('M_Hydro',       .083, [1e-8, 0.5]),
         Parameter('M_Fe3',          .25, [1e-8, 0.5]),
         Parameter('M_Homo',         .25, [1e-8, 0.5]),
         Parameter('M_Ac',         .0033, [1e-8, 0.5]),
        ]
        
    if not model_parameters is None:
        return [par for par in p if par in model_parameters]
    return p

class LogTransform():
    def transform(self, value): return np.log(value)
    def inverse(self, value): return np.exp(value)
    
class IdentityTransform():
    def transform(self, value): return value
    def inverse(self, value): return value


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
            if self.is_variable() and (not p <= self.high or not p >= self.low):
                raise Exception(f'Setting {self.name} to {p} is out of bounds.')
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
    
