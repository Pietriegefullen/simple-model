import numpy as np

import chemistry

SYSTEM = [
            'C',
            'DOC',
            'CH4',
            'CO2',
            'Fe3',
            'M_Ac',
            'M_Ferm',
            'M_Fe3',
            'M_Hydro',
            'M_Homo',
            'Acetate',
            'H2',
            'Fe2',
            'H2O',
            ]



def index(name):
    return SYSTEM.index(str(name))

def stoichiometry(educts, products):
    
    return stoichiometry_vector

def henrys_law():
    henrys_values = []
    for pool in SYSTEM:
        value = chemistry.henrys_law(pool)
        henrys_values.append(value)
    henrys_law_vector = np.reshape(henrys_values, (-1,))
    return henrys_law_vector

def vector(full, pool = None, value = None):
    _vector = np.full((len(SYSTEM), ), float(full))
    if (pool is None) != (value is None):
        raise Exception('Require pool AND value or both None')
    if not pool is None:
        _vector[index(pool)] = value
    return _vector

def initial_state(replica, model_parameters):
    S0 = np.zeros((len(SYSTEM),))
    S0 += vector(0, 'C', replica.initial_C())
    S0 += vector(0, 'DOC', replica.initial_DOC())
    S0 += vector(0, 'H2O', replica.water_content)
    
    model_parameters.add('Fe3', 20, [0, 300])
    model_parameters.add('M_Ferm', .2, [1e-8, 0.5])
    model_parameters.add('M_Hydro', .0025, [1e-8, 0.5])
    model_parameters.add('M_Homo', .0001, [1e-8, 0.5])
    model_parameters.add('M_Ac', .0101, [1e-8, 0.5])
    model_parameters.add('Acetate', 50, [0, 100])
    
    for pool in SYSTEM:
        if pool in model_parameters:
            init = model_parameters[pool].value()
            v = vector(0, pool, init)
            S0 += v
    
    return S0