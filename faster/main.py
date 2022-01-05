import numpy as np
import matplotlib.pyplot as plt

import predict
import pathways
import plot
import data
import optimizer
from ORDER import POOL_ORDER
import OPTIMIZATION_PARAMETERS

def run_model(specimen_index, site):

    model_parameters = pathways.default_model_parameters(specimen_index, site)

    all_days = np.arange(4500)     # Days to make predictions for

    chosen_pathways = [
                       pathways.Ferm_help,
                       pathways.Ferm,
                       pathways.Fe3,
                       pathways.Hydro,
                       pathways.Homo,
                       pathways.Ac
                       ]

    pool_value_dict = predict.predictor(all_days,
                                        model_parameters,
                                        chosen_pathways,
                                        verbose = True)


    measured_data = data.specimen_data(specimen_index, site)
    plot.all_pools(pool_value_dict, all_days, measured_data)


def fit_model(specimen_index, site):

    chosen_pathways = [pathways.Ferm_help,
                       pathways.Ferm,
                       pathways.Fe3,
                       pathways.Hydro,
                       pathways.Homo,
                       pathways.Ac
                       ]

    fixed_parameters = pathways.default_model_parameters(specimen_index)

    optimal_parameters = optimizer.fit_specimen(specimen_index,
                                                site,
                                                chosen_pathways,
                                                fixed_parameters,
                                                algo = OPTIMIZATION_PARAMETERS.ALGORITHM)

    model_parameters = fixed_parameters
    model_parameters.update(data.model_parameters_from_data(specimen_index, site = site))

    model_parameters.update(optimal_parameters)

    all_days = np.arange(4500)

    pool_value_dict = predict.predictor(all_days,
                                        model_parameters,
                                        chosen_pathways,
                                        verbose = True)

    measured_data_dict = data.specimen_data(specimen_index, site)

    days = measured_data_dict['measured_time']
    CO2 =  measured_data_dict['CO2']
    CH4 =  measured_data_dict['CH4']

    plot.fit(days, CO2, pool_value_dict['CO2'], all_days)
    plot.fit(days, CH4, pool_value_dict['CH4'], all_days)

    plt.show()

if __name__ == '__main__':
    fit_model(9, 'all')
    #run_model(9, site = 'all')

    # TODO: print setup, then ask for confirmation

    # TODO:
    # extract extra dicts after run
    # do all the plotting
    # optimization tolerances
    # callback print during optimization?
    # save optimization checkpoints
    # fix model errors
    #
    # scaling for the optimization algorithm!
    #
