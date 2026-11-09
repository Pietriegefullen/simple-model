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

def initial_state(replica, model_parameters, initial_mean_days = 0):
    S0 = np.zeros((len(SYSTEM),))
    
    if not replica is None:
        gdw_per_L = replica.dry_weight/replica.water_content*1e-3 # gdw/L
        H2O_per_L = replica.initial_H2O()*gdw_per_L # micromol L-1
        TOC_per_L = replica.initial_TOC()*gdw_per_L
        DOC_per_L = replica.initial_DOC()*gdw_per_L

        model_parameters['TOC'].constant(TOC_per_L)
        model_parameters['DOC'].constant(DOC_per_L)
        model_parameters['H2O'].constant(H2O_per_L)
    
    else:
        model_parameters['TOC']
        model_parameters['DOC']
        model_parameters['H2O']

    # in parameters, microbes are in micromol L-1.
    # here, convert to milligram per gram dry weight
    for pool in SYSTEM:
        if pool in model_parameters:
            init = model_parameters[pool] # micromol L-1
            v = vector(0, pool, init)
            S0 += v
        else:
            model_parameters[pool].constant(0)
    
    if not replica is None and str(replica.sample.sample_name).startswith('2'):
        # inocculated sample
        S0[index('Fe3')] = model_parameters['Fe3'].constant(0)
        S0[index('Acetate')] = model_parameters['Acetate'].constant(0)
        S0[index('DOC')] = model_parameters['DOC'].constant(0)
    
    elif not replica is None:
        d = int(np.nonzero(replica.CO2()[0] >= initial_mean_days)[0][0])
        initial_mean = replica.CO2()[1][0]
        if d > 0:
            initial_mean = np.mean(replica.CO2()[1][:d])
        S0[index('CO2')] = initial_mean*gdw_per_L
        
        d = int(np.nonzero(replica.CH4()[0] >= initial_mean_days)[0][0])
        initial_mean = replica.CH4()[1][0]
        if d > 0:
            initial_mean = np.mean(replica.CH4()[1][:d])
        S0[index('CH4')] = initial_mean*gdw_per_L
        
    return S0

