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
    
    if not replica is None:
        model_parameters['C'].constant(replica.initial_C())
        DOC = replica.initial_C()*.02*model_parameters['DOC_per_TOC']
        model_parameters['DOC'].constant(DOC)
        model_parameters['H2O'].constant(replica.water_content)
    
    else:
        model_parameters['C']
        model_parameters['DOC']
        model_parameters['H2O']
    
    model_parameters['Fe3']
    model_parameters['M_Ferm']
    model_parameters['M_Hydro']
    model_parameters['M_Homo']
    model_parameters['M_Ac']
    #model_parameters['M_Fe3']
    model_parameters['Acetate']
    
    for pool in SYSTEM:
        if pool in model_parameters:
            init = model_parameters[pool]
            v = vector(0, pool, init)
            S0 += v
    
    return S0
