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

def load_model_parameters(file):
    
    file_path = os.path.join(USER_VARIABLES.LOG_DIRECTORY, file)
    if not file_path.endswith('.json'):
        file_path += '.json'
        
    with open(file_path, 'r') as pf:
        model_parameters = json.load(pf)

    return model_parameters

def load_and_plot(file):
    
    splitparts = file.split('_')
    specimen_index = int(splitparts[splitparts.index('specimen')+1])
    site = splitparts[splitparts.index('site')+1]
    
    all_days = np.arange(4500)     # Days to make predictions for

    model_parameters = load_model_parameters(file)
    pool_value_dict = run_model(model_parameters, all_days)
    
    measured_data = data.specimen_data(specimen_index, site)
    plot.all_pools(pool_value_dict, all_days, measured_data)

def run_and_plot(specimen_index, site):
    
    all_days = np.arange(4500)     # Days to make predictions for

    model_parameters = pathways.default_model_parameters(specimen_index, site)
    
    pool_value_dict = run_model(model_parameters, all_days)
    
    measured_data = data.specimen_data(specimen_index, site)
    plot.all_pools(pool_value_dict, all_days, measured_data)


    
def run_model(model_parameters, all_days):

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


    return pool_value_dict

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
    
    timestamp = start_time.strftime('%Y-%m-%d_%H-%M-%S')
    file_name = '_'.join([timestamp,
                          'specimen',
                          str(specimen_index),
                          'site',
                          site]) + '.json'
    parameter_file = os.path.join(USER_VARIABLES.LOG_DIRECTORY, file_name)
    with open(parameter_file, 'w') as pf:
        json.dump(model_parameters, pf,indent = 4)
        
    pathway_file_name = '_'.join([timestamp,
                                  'pathways']) + '.txt'
    pathway_file = os.path.join(USER_VARIABLES.LOG_DIRECTORY, pathway_file_name)
    with open(pathway_file, 'w') as pf:
        for p in chosen_pathways:
            pf.write(p.__name__ + '\n')

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
    #load_and_plot('2022-01-21_10-32-02_specimen_17_site_all')
    #fit_model(17, 'all')
    
    run_and_plot("13782", site = 'all')

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
