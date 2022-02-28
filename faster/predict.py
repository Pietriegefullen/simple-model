import scipy.integrate
import traceback
import numpy as np

# run the model (IVP, initial value problem) for a set of initial values.

# allow two (several?) run modes: one basic, to be used for integration,
# and another that returns additional information, not just the pool values.

import pathways
import model
from ORDER import POOL_ORDER, pool_index

METHOD = 'Radau' #'LSODA' #'BDF'

def predictor(t_eval,
              model_parameters,
              chosen_pathways,
              verbose = False,
              mark = None,
              ask = False,
              extended_output = None):
    
    environment = {} 
    if 'temperature' in model_parameters:
        environment.update({'temperature': model_parameters['temperature']})
                            
    if 'pH' in model_parameters:
        environment.update({'pH': model_parameters['pH']})
    
    defined_pathways = [pathway_definer(model_parameters) for pathway_definer in chosen_pathways]
    right_hand_side = model.builder(defined_pathways, environment)

    initial_system_state = np.array([model_parameters[name]
                                     if name in model_parameters else 0.
                                     for name in POOL_ORDER])

    if verbose:
        print('')
        print('running model:')
        print('=============')
        print('pathways:')
        for p in chosen_pathways:
            print('   ' + p.__name__)
           #print(f'   {p.__name__:s}')
        print('')

        print('initial system:')
        mark = list() if mark is None else mark
        for i, pool in enumerate(POOL_ORDER):
            star = '*' if pool in mark else ''
            print(f'   {pool:10s}{star:2s} {initial_system_state[i]:g}')
        print('')

        print('model parameters:')
        for k, v in model_parameters.items():
            star = '*' if k in mark else ''
            print(f'   {k[:20]:20}{star:2s} {v}')
        print('')

    if ask:
        input('Start?')

    extended_results = {}

    try:
        solver_result = scipy.integrate.solve_ivp(right_hand_side,
                                                  (0, max(t_eval)),
                                                  initial_system_state,
                                                  method = METHOD,
                                                  t_eval = t_eval)#,
                                                  # atol = 1e-100,
                                                  # rtol = 1e-1,
                                                  # max_step = 10)
        pool_results = solver_result.y
  
        # extended output
        if not extended_output is None:
            extended_model = model.builder(defined_pathways, environment, extended_output)
            extended_results = {}
            for t, system_state in zip(t_eval, np.transpose(pool_results)):
                extended_values = extended_model(t,system_state)
                for k,v in extended_values.items():
                    if not k in extended_results:
                        extended_results[k] = list()
                    extended_results[k].append(v)
                                           
    except Exception:
        print(traceback.format_exc())
        print('EXCEPTION IN SOLVER')
        n_pools = len(POOL_ORDER)
        pool_results = np.empty((n_pools, t_eval.size))
        pool_results[:] = np.nan

    n_pools, n_days = pool_results.shape
    nan_array = np.empty((n_pools, t_eval.size - n_days))
    nan_array[:] = np.nan
    padded_pools = np.concatenate([pool_results,
                                   nan_array], axis = -1)

    pool_dict = dict(zip(POOL_ORDER, padded_pools))
    
    pool_dict.update(extended_results)
    
    return pool_dict
