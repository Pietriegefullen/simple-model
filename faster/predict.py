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

def predictor(t_eval, model_parameters, all_pathways, initial_system_state):

    built_pathways = [pathway_builder(model_parameters) for pathway_builder in all_pathways]
    right_hand_side = model.builder(model_parameters, built_pathways)

    print('running model:')
    print('=============')
    print('pathways:')
    for p in all_pathways:
        print(f'   {p.__name__}')
    print('')

    print('initial system:')
    for i, pool in enumerate(POOL_ORDER):
        print(f'   {pool:10s} {initial_system_state[i]:g}')
    print('')

    print('model parameters:')
    for k, v in model_parameters.items():
        print(f'   {k[:15]:15} {v}')
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


def run(model_parameters):
    chosen_pathways = [pathways.Ferm_help,
                       pathways.Ferm,
                       pathways.Fe3,
                       pathways.Hydro,
                       pathways.Homo,
                       pathways.Ac
                       ]

    initial_system_state = np.array([model_parameters[name + '_init']
                                     if name + '_init' in model_parameters else 0.
                                     for name in POOL_ORDER])

    all_days = np.arange(4500)                                                 # Days to make predictions for
    pool_value_dict = predictor(all_days,                                      # vorhersagen erstellen
                                model_parameters,
                                chosen_pathways,
                                initial_system_state)

    return pool_value_dict
