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
    specimen_index = str(splitparts[splitparts.index('specimen')+1])
    site = splitparts[splitparts.index('site')+1]

    pathway_file = os.path.join(USER_VARIABLES.LOG_DIRECTORY,
                                file + '_pathways.txt')
    with open(pathway_file, 'r') as pwf:
        used_pathways = [line.replace('\n', '') for line in pwf.readlines()]

    all_days = np.arange(4500)     # Days to make predictions for

    default_model_parameters = pathways.default_model_parameters(specimen_index, site)
    model_parameters = load_model_parameters(file)
    
    for default_key, default_value in default_model_parameters.items():
        if not default_key in model_parameters.keys():
            print(f'Could not find {default_key} in model parameters. Using default value.')
            model_parameters[default_key] = default_value

    pool_value_dict = run_model(model_parameters, all_days, pathway_names = used_pathways)

    measured_data = data.specimen_data(str(specimen_index), site)
    plot.all_pools(pool_value_dict, all_days, measured_data)

def run_and_plot(specimen_index, site, extended_output = None, pathway_names = None):

    all_days = np.arange(4500)     # Days to make predictions for

    model_parameters = pathways.default_model_parameters(specimen_index, site)

    pool_value_dict = run_model(model_parameters, all_days, extended_output, pathway_names)

    measured_data = data.specimen_data(specimen_index, site)
    plot.all_pools(pool_value_dict, all_days, measured_data)

def run_model(model_parameters, all_days, extended_output = None, pathway_names = None):
    all_pathways = [
                       pathways.Ferm_help,
                       pathways.Ferm,
                       pathways.Fe3,
                       pathways.Hydro,
                       pathways.Homo,
                       pathways.Ac
                       ]
    
    if pathway_names is None:
        chosen_pathways = all_pathways
        
    else:    
        chosen_pathways = [pw for pw in all_pathways if pw.__name__ in pathway_names]    


    pool_value_dict = predict.predictor(all_days,
                                        model_parameters,
                                        chosen_pathways,
                                        verbose = True,
                                        extended_output = extended_output)


    return pool_value_dict

def save_model(specimen_index, site, model_parameters, prefix = ''):
    save_time = datetime.now()
    timestamp = save_time.strftime('%Y-%m-%d_%H-%M-%S')
    file_name = '_'.join([prefix,
                          timestamp,
                          'specimen',
                          str(specimen_index),
                          'site',
                          site]) + '.json'
    parameter_file = os.path.join(USER_VARIABLES.LOG_DIRECTORY, file_name)

    if not os.path.isdir(USER_VARIABLES.LOG_DIRECTORY):
        os.makedirs(USER_VARIABLES.LOG_DIRECTORY)

    with open(parameter_file, 'w') as pf:
        json.dump(model_parameters, pf,indent = 4)
    
    return parameter_file

def fit_model(specimen_index, site, pathway_names = None):
    all_pathways = [
                       pathways.Ferm_help,
                       pathways.Ferm,
                       pathways.Fe3,
                       pathways.Hydro,
                       pathways.Homo,
                       pathways.Ac
                       ]
    
    if pathway_names is None:
        chosen_pathways = all_pathways
        
    else:    
        chosen_pathways = [pw for pw in all_pathways if pw.__name__ in pathway_names]   
    fixed_parameters = pathways.default_model_parameters(specimen_index)

    optimal_parameters = optimizer.fit_specimen(specimen_index,
                                                site,
                                                chosen_pathways,
                                                fixed_parameters,
                                                algo = OPTIMIZATION_PARAMETERS.ALGORITHM)

    model_parameters = dict(fixed_parameters)
    model_parameters.update(optimal_parameters)

    parameter_file = save_model(specimen_index, site, model_parameters)
    pathway_file = parameter_file.replace('.json', '_pathways.txt')
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
    load_and_plot('_2022-02-25_19-09-04_specimen_13560_site_all')

    speciemen_identifier = "13560"
#=============================================================================
    # run_and_plot(speciemen_identifier, 
    #              site = 'all', 
    #              extended_output = ['deltaGr',
    #                                 'deltaCO2',
    #                                 'deltaCH4',
    #                                 'thermo',
    #                                 'MM',
    #                                 'v',
    #                                 'deltaGs',
    #                                 'inhibition'],
    #               pathway_names = [
    #                                 'Ferm_help',
    #                                 'Ferm',
    #                                 'Fe3',
    #                                 # 'Ac',
    #                                 # 'Hydro',
    #                                 # 'Homo'
    #                                 ])
#=============================================================================
    fit_model(speciemen_identifier, 
             site = 'all', 
             pathway_names = [
                                'Ferm_help',
                                'Ferm',
                                'Fe3',
                                # 'Ac',
                                # 'Hydro',
                                # 'Homo'
                                ])

"""
'13510', '13511', '13512', '13520', '13521', '13530', '13531', '13670', '13671', '13672',
 '13690', '13691', '13692', '13700', '13701', '13702', '13720', '13721', '13722', '13730',
 '13731', '13732', '13740', '13741', '13742', '13750', '13751', '13752', '13770', '13771',
 '13772', '13780', '13781', '13782', '13540', '13542', '13550', '13551', '13560', '13562',
 '13571', '13572', '13580', '13581', '13590', '13591', '13600', '13602', '13610', '13612',
 '13620', '13622', '13630', '13631', '13640', '13641', '13650', '13651', '13652', '13661',
 '13662', '13680', '13681', '13682', '13710', '13711', '13712', '13760', '13761', '13762',
 '13790', '13791', '13792', '13800', '13801', '13802'

    
    Rep_No_CH4 = [13560, 13562, 13580, 13581, 13590, 13591, 13600, 13602, 13622, 13641]

"""



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
