import numpy as np
import matplotlib.pyplot as plt

import predict
import pathways
import plot
import data
import optimizer
from ORDER import POOL_ORDER

# what are the required functionalities?
# run for a single sample? run for arbitrary parameters?
# run for several samples?
def run(all_days, model_parameters):

    chosen_pathways = [pathways.Ferm_help,
                       pathways.Ferm,
                       pathways.Fe3,
                       pathways.Hydro,
                       pathways.Homo,
                       pathways.Ac
                       ]

    initial_system_state = np.array([model_parameters[name]
                                     if name in model_parameters else 0.
                                     for name in POOL_ORDER])

    pool_value_dict = predict.predictor(all_days,                                      # vorhersagen erstellen
                                        model_parameters,
                                        chosen_pathways,
                                        initial_system_state)

    return pool_value_dict


def run_specimen(specimen_index):

    model_parameters = data.model_parameters_from_data(specimen_index, site = 'all')

    model_parameters.update({'M_Fe3':           0.2,
                            'M_Ferm':           0.7,
                            'M_Hydro':          0.0025,
                            'M_Homo':           0.0001,
                            'M_Ac':             0.001,

                            'Sensenmann':       8.33e-5,

                            'Vmax_Fe3':         1.709,
                            'Vmax_help_Ferm':   0.177,
                            'Vmax_Ferm':        0.606,
                            'Vmax_Homo':        0.85,
                            'Vmax_Hydro':       0.423,
                            'Vmax_Ac':          0.159,

                            'Kmb_help_Ferm':    1.0,
                            'Inhibition_Ferm':  1.0,

                            'Fe3':              3.645})

    all_days = np.arange(4500)                                                 # Days to make predictions for
    pool_value_dict = run(all_days, model_parameters)

    plot.all_pools(pool_value_dict, all_days)


def fit_specimen(specimen_index, site):

    chosen_pathways = [pathways.Ferm_help,
                       pathways.Ferm,
                       pathways.Fe3,
                       pathways.Hydro,
                       pathways.Homo,
                       pathways.Ac
                       ]

    fixed_parameters = {'Sensenmann':           8.33e-5}

    optimal_parameters = optimizer.fit_specimen(specimen_index,
                                                site,
                                                chosen_pathways,
                                                fixed_parameters,
                                                algo = 'gradient')

    model_parameters = data.model_parameters_from_data(specimen_index, site = 'all')
    model_parameters.update({'M_Fe3':           0.2,
                            'M_Ferm':           0.7,
                            'M_Hydro':          0.0025,
                            'M_Homo':           0.0001,
                            'M_Ac':             0.001,

                            'Sensenmann':       8.33e-5,

                            'Vmax_Fe3':         1.709,
                            'Vmax_help_Ferm':   0.177,
                            'Vmax_Ferm':        0.606,
                            'Vmax_Homo':        0.85,
                            'Vmax_Hydro':       0.423,
                            'Vmax_Ac':          0.159,

                            'Kmb_help_Ferm':    1.0,
                            'Inhibition_Ferm':  1.0,

                            'Fe3':              3.645})
    model_parameters.update(optimal_parameters)

    all_days = np.arange(4500)
    pool_value_dict = run(all_days, model_parameters)

    # plot.all_pools(pool_value_dict, all_days)

    measured_data_dict = data.specimen_data(specimen_index, site)

    days = measured_data_dict['measured_time']
    CO2 =  measured_data_dict['CO2']
    CH4 =  measured_data_dict['CH4']

    plot.fit(days, CO2, pool_value_dict['CO2'], all_days)
    plot.fit(days, CH4, pool_value_dict['CH4'], all_days)

    plt.show()

if __name__ == '__main__':
    fit_specimen(9, 'all')


    # TODO:
    # extract extra dicts after run
    # do all the plotting
    # optimization tolerances
    # explicit separation of model parameters, changeables, ...
    # callback print during optimization?
    # save optimization checkpoints
    # fix model errors
