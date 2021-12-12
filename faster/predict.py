import scipy.integrate
import traceback
import numpy as np

# run the model (IVP, initial value problem) for a set of initial values.

# allow two (several?) run modes: one basic, to be used for integration,
# and another that returns additional information, not just the pool values.

import pathways
import model
import INTEGRATION_PARAMETERS
from ORDER import POOL_ORDER, pool_index

counter = 0

def predictor(t_eval, model_parameters, all_pathways, initial_system_state, verbose = False, mark = None):
    # global counter
    # counter += 1
    # print(counter)

    built_pathways = [pathway_builder(model_parameters) for pathway_builder in all_pathways]
    right_hand_side = model.builder(model_parameters, built_pathways)

    if verbose:
        print('running model:')
        print('=============')
        print('pathways:')
        for p in all_pathways:
            print(f'   {p.__name__}')
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

    try:
        solver_result = scipy.integrate.solve_ivp(right_hand_side,
                                                  (0, max(t_eval)),
                                                  initial_system_state,
                                                  method = INTEGRATION_PARAMETERS.METHOD,
                                                  t_eval = t_eval)
                                                  # atol = 1e-100,
                                                  # rtol = 1e-1,
                                                  # max_step = 10)

    except Exception:
        print(traceback.format_exc())

    n_pools, n_days = solver_result.y.shape
    nan_array = np.empty((n_pools, t_eval.size - n_days))
    nan_array[:] = np.nan
    padded_pools = np.concatenate([solver_result.y,
                                   nan_array], axis = -1)

    pool_dict = dict(zip(POOL_ORDER, padded_pools))

    return pool_dict
