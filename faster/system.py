import numpy as np

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
            'pH',
            'H2O',
            'weight',
            'water'
            ]

def index(name):
    return SYSTEM.index(name)

def stoichiometry(educts, products):
    
    return stoichiometry_vector

def henrys_law():
    henrys_values = []
    for pool in SYSTEM:
        value = 
        henrys_values.append(value)
    henrys_law_vector = np.reshape(henrys_values, (-1,))
    return henrys_law_vector

def vector(pool, value):
    _vector = np.zeros((len(SYSTEM), ))
    _vector[index(pool)] = value
    return _vector
