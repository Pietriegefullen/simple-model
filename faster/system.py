import numpy as np

SYSTEM = [
            'TOC',
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

MICROBES = np.array([pool_name.startswith('M_') for pool_name in SYSTEM])

def index(name):
    return SYSTEM.index(str(name))

def vector(full, pool = None, value = None):
    _vector = np.full((len(SYSTEM), ), float(full))
    if (pool is None) != (value is None):
        raise Exception('Require pool AND value or both None')
    if not pool is None:
        _vector[index(pool)] = value
    return _vector

def initial_state(replica, model_parameters, initial_mean_days = 0, reference = 'gdw'):
    S0 = np.zeros((len(SYSTEM),))

    if not replica is None:
        conversion = 1.
        if reference == 'L': 
            L_per_gdw = (replica.water_content*1e-3)/replica.dry_weight # L/gdw
            conversion = L_per_gdw

        model_parameters['TOC'].constant(replica.initial_TOC()/conversion)
        model_parameters['DOC'].constant(replica.initial_DOC()/conversion)
        model_parameters['H2O'].constant(replica.initial_H2O()/conversion)

    else:
        model_parameters['TOC']
        model_parameters['DOC']
        model_parameters['H2O']

    for pool in SYSTEM:
        S0 += vector(0, pool, model_parameters[pool])
    
    if not replica is None and str(replica.sample.sample_name).startswith('2'):
        # inocculated sample
        S0[index('Fe3')] = 0
        S0[index('Acetate')] = 0
        S0[index('DOC')] = 0
    
    elif not replica is None:
        d = int(np.nonzero(replica.CO2()[0] >= initial_mean_days)[0][0])
        initial_mean = replica.CO2()[1][0]
        if d > 0:
            initial_mean = np.mean(replica.CO2()[1][:d])
        S0[index('CO2')] = initial_mean/conversion
        
        d = int(np.nonzero(replica.CH4()[0] >= initial_mean_days)[0][0])
        initial_mean = replica.CH4()[1][0]
        if d > 0:
            initial_mean = np.mean(replica.CH4()[1][:d])
        S0[index('CH4')] = initial_mean/conversion
        
    return S0
