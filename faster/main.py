import os
import numpy as np
import matplotlib.pyplot as plt

import argparse
from datetime import datetime
import json

import predict
import pathways
import plot
import data
import optimizer
from ORDER import POOL_ORDER
import OPTIMIZATION_PARAMETERS
import USER_VARIABLES

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

    start_time = datetime.now()    

    chosen_pathways = [pathways.Ferm_help,
                       pathways.Ferm,
                       pathways.Fe3,
                       #pathways.Hydro,
                       #pathways.Homo,
                       pathways.Ac
                       ]

    fixed_parameters = pathways.default_model_parameters(specimen_index)

    optimal_parameters = optimizer.fit_specimen(specimen_index,
                                                site,
                                                chosen_pathways,
                                                fixed_parameters,
                                                algo = OPTIMIZATION_PARAMETERS.ALGORITHM)

    model_parameters = dict(fixed_parameters)
    model_parameters.update(optimal_parameters)
    
    file_name = '_'.join([start_time.strftime('%Y-%m-%d_%H-%M-%S'),
                          'specimen',
                          str(specimen_index),
                          'site',
                          site]) + '.json'
    parameter_file = os.path.join(USER_VARIABLES.LOG_DIRECTORY, file_name)
    with open(parameter_file, 'w') as pf:
        json.dump(model_parameters, pf,indent = 4)

    all_days = np.arange(4500)

    pool_value_dict = predict.predictor(all_days,
                                        model_parameters,
                                        chosen_pathways,
                                        verbose = True)

    measured_data_dict = data.specimen_data(specimen_index, site)

    days = measured_data_dict['measured_time']
    CO2 =  measured_data_dict['CO2']
    CH4 =  measured_data_dict['CH4']

    plot.all_pools(pool_value_dict, all_days)
    plot.fit(days, CO2, pool_value_dict['CO2'], all_days)
    plot.fit(days, CH4, pool_value_dict['CH4'], all_days)

    plt.show()

my_parser = argparse.ArgumentParser(description='simple model main')

my_parser.add_argument('-w', '-workers',
                       metavar='workers',
                       type=int,
                       help='number of workers used by optimization',
                       default = 8)

args = my_parser.parse_args()
OPTIMIZATION_PARAMETERS.WORKERS = args.w

if __name__ == '__main__':
    # 9 ist die probe die ich normalerweise hab
    fit_model(17, 'all')
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
