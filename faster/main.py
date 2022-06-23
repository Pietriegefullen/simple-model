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
from data import load_matlab
import math
import random 
import seaborn as sns



# print('sample', all_the_samples[0])

# for sample_numper in all_the_samples:
#     print('sample numer', sample_numper)

def load_model_parameters(file):

    file_path = os.path.join(USER_VARIABLES.LOG_DIRECTORY, file)
    if not file_path.endswith('.json'):
        file_path += '.json'

    with open(file_path, 'r') as pf:
        model_parameters = json.load(pf)

    return model_parameters

def load_and_plot(file, extended_output = None, which_plot = 'measured_and_modelled', show = True):

    figures = []    

    splitparts = file.split('_')
    specimen_index = str(splitparts[splitparts.index('specimen')+1])
    site = splitparts[splitparts.index('site')+1]

    pathway_file = os.path.join(USER_VARIABLES.LOG_DIRECTORY,
                                file + '_pathways.txt')
    with open(pathway_file, 'r') as pwf:
        used_pathways = [line.replace('\n', '') for line in pwf.readlines()]

    all_days = np.arange(4500)     # Days to make predictions for

    default_model_parameters = {}
    if not site == 'None': # um das durchschnittsmodell zu plotten
        default_model_parameters = pathways.default_model_parameters(specimen_index, site)
    model_parameters = load_model_parameters(file)
    
    for default_key, default_value in default_model_parameters.items():
        if not default_key in model_parameters.keys():
            print(f'Could not find {default_key} in model parameters. Using default value.')
            model_parameters[default_key] = default_value
    
    pool_value_dict = run_model(model_parameters, 
                                all_days, 
                                extended_output = extended_output,
                                pathway_names = used_pathways)

    measured_data = None# um das durchschnittsmodell zu plotten
    if not site == 'None':# um das durchschnittsmodell zu plotten
        measured_data = data.specimen_data(str(specimen_index), site)
    
    if which_plot == 'measured_and_modelled':
        # plot both model and data:
        figures = plot.all_pools(pool_value_dict, all_days, specimen_index, 
                                 measured_data, show = show)
        
    elif which_plot == 'only_data':
        days = measured_data['measured_time']
        for name, measurements in measured_data.items():
            if 'time' in name: 
                continue
            if not np.array(measurements).size == len(days):
                continue
            plot.plot_pool(name, measurements, days, 'xr')
        
        figures = plot.plot_pool(specimen_index, measured_data[1], all_days)

    return figures


def run_and_plot(specimen_index, site, extended_output = None, pathway_names = None):

    all_days = np.arange(4500)     # Days to make predictions for

    model_parameters = pathways.default_model_parameters(specimen_index, site)

    pool_value_dict = run_model(model_parameters, all_days, extended_output, pathway_names)

    measured_data = data.specimen_data(specimen_index, site)
#    plot.all_pools({'CO2':pool_value_dict['CO2']}, all_days, specimen_index, measured_data)
    plot.all_pools(pool_value_dict, all_days, specimen_index, measured_data)

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
    notplot = True
    
    all_pathways = [
                       pathways.Ferm_help,
                       pathways.Ferm,
                       #pathways.Fe3,
                       #pathways.Homo,
                       ]
    
    if site == 'all':
        all_pathways += [pathways.Ac,
                        # pathways.Hydro
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

    if notplot:
        return

    all_days = np.arange(4500)

    pool_value_dict = predict.predictor(all_days,
                                        model_parameters,
                                        chosen_pathways,
                                        verbose = True)

    measured_data_dict = data.specimen_data(specimen_index, site)

    days = measured_data_dict['measured_time']
    CO2 =  measured_data_dict['CO2']
    CH4 =  measured_data_dict['CH4']

    plot.all_pools(pool_value_dict, all_days, specimen_index,measured_data_dict)
    plot.fit(days, CO2, pool_value_dict['CO2'], all_days)
    plot.fit(days, CH4, pool_value_dict['CH4'], all_days)

    plt.show()





def evaluate_loss(file):
    splitparts = file.split('_')
    specimen_index = str(splitparts[splitparts.index('specimen')+1])
    site = splitparts[splitparts.index('site')+1]

    pathway_file = os.path.join(USER_VARIABLES.LOG_DIRECTORY,
                                file + '_pathways.txt')
    with open(pathway_file, 'r') as pwf:
        used_pathways = [line.replace('\n', '') for line in pwf.readlines()]

    model_parameters = load_model_parameters(file)
    
    all_pathways = [
                       pathways.Ferm_help,
                       pathways.Ferm,
                       #pathways.Fe3,
                       #pathways.Hydro,
                       #pathways.Homo,
                       pathways.Ac
                       ]
    
    if used_pathways is None:
        chosen_pathways = all_pathways
        
    else:    
        chosen_pathways = [pw for pw in all_pathways if pw.__name__ in used_pathways]  
    
    loss = optimizer.evaluate_loss(specimen_index, site, chosen_pathways, model_parameters)
    #anzeigen vom loss einer bestimmten Probe
    print('the loss of',file,'is', loss)
    return (specimen_index, loss)


my_parser = argparse.ArgumentParser(description='simple model main')

my_parser.add_argument('-w', '-workers',
                       metavar='workers',
                       type=int,
                       help='number of workers used by optimization',
                       default = 8)

args = my_parser.parse_args()
OPTIMIZATION_PARAMETERS.WORKERS = args.w




def mean_model_and_datapoints():
    
    replica_list_No_CH4, superdata_No_CH4_vor_Impfung,superdata_No_CH4_nach_Impfung,superdata_after_No_CH4,superdata_bevor_No_CH4, superdata_No_CH4,  superdata, replica_list, superdata_carex, superdata_Kuru, superdata_Sam, replica_list_Kuru, replica_list_Sam,superdata_2021_all, replica_list_superdata_2021_all, superdata_ohne_Fe3, Rep_ohne_Fe3,superdata_mit_Fe3, Rep_mit_Fe3 = load_matlab()

    #site = superdata_2021_all
    #site = superdata_Kuru
    site = superdata_Sam
    
    plotted_figures = load_and_plot('Mean_Parameters_specimen_000000_site_None', 
                                    which_plot = 'measured_and_modelled',
                                    show = False)
    
    CO2_plot = plotted_figures['CO2']
    CH4_plot = plotted_figures['CH4']
    

    Random_day_values = dict()
    
    all_measured_days = [int(day) for specimen in site.values() for day in specimen['measured_time']]

    # extract for all specimens the days for which most specimens were measured
    consider_more_than = 10
    unique, counts = np.unique(all_measured_days, return_counts = True)    
    sorted_by_freq = reversed(sorted(zip(unique, counts), key = lambda x: x[1]))
    only_most_freq = sorted(list(filter(lambda x: x[1] >= consider_more_than, sorted_by_freq)))
    chosen_days, _ = zip(*only_most_freq)

    # from all specimens that have been measured on the chosen day, 
    # extract the CO2 and CH4 measurements
    boxplot_dict = {'CO2':{},
                    'CH4':{}}
    CO2_boxplot_data = {}
    CH4_boxplot_data = {}
    for day in chosen_days:
        CO2_boxplot_data[day] = []
        CH4_boxplot_data[day] = []
        for specimen, specimen_data in site.items():
            int_measured_time = [int(d) for d in specimen_data['measured_time']]
            if not day in int_measured_time:
                continue
            day_index = int_measured_time.index(day)
            CO2_boxplot_data[day].append(specimen_data['CO2'][day_index])
            CH4_boxplot_data[day].append(specimen_data['CH4'][day_index])
       
    
    plt.figure('CO2')    
    plt.boxplot(list(CO2_boxplot_data.values()), 
                showmeans=True,
                positions = list(CO2_boxplot_data.keys()))
    
    plt.figure('CH4')    
    plt.boxplot(list(CH4_boxplot_data.values()), 
                showmeans=True,
                positions = list(CH4_boxplot_data.keys()))


    
   
      
    
    




#'13510', '13511','13512','13520',  '13521', '13530', '13531', '13670', '13671',
#                  '13672', '13690', '13691', '13692', '13700',
#
                  
all_the_samples =[ '13701', '13702', '13720', '13721',
                  '13722', '13730', '13731', '13732', '13740', '13741', '13742', '13750', '13751', 
                  '13752', '13770', '13771','13772', '13780', '13781', '13782', '13540', '13542', 
                  '13550', '13551', '13560', '13562','13571', '13572', '13580', '13581', '13590', 
                  '13591', '13600', '13602', '13610', '13612','13620', '13622', '13630', '13631', 
                  '13640', '13641', '13650', '13651', '13652', '13661','13662', '13680', '13681',
                  '13682', '13710', '13711', '13712', '13760', '13761', '13762','13790', '13791', 
                  '13792', '13800', '13801', '13802']

#all_the_samples =[ '13692']

# automatically assign a site to each specimen.
all_samples_and_sites = []
for sample, sites in data.specimen_sites(all_the_samples).items():
    chosen_site = 'all'
    if 'No-CH4' in sites:
        chosen_site = 'No-CH4'
    all_samples_and_sites.append((sample, chosen_site))
    
    
# laden aller log files um den loss zu bekommen    
all_files =[]
for filename in os.listdir(USER_VARIABLES.LOG_DIRECTORY):
     
     if filename.endswith(".json") : 
         all_files.append(filename.replace('.json', ''))
             
         
#print(all_files)         
         

if __name__ == '__main__':
    
    # Mean_Parameters_specimen_000000_site_None # durchschnittliche Parameterwerte
    mean_model_and_datapoints() 
    
    file = ('Mean_Parameters_specimen_000000_site_None')
    
    #anzeigen vom Loss einer bestimmten Probe
    #evaluate_loss(file)

    
    # erstellen der Loss liste aller fits:
    
    # all_losses = []
    # for file_0 in all_files:
    #     loss_tuple = evaluate_loss(file_0)
    #     all_losses.append(loss_tuple)
    # print(all_losses)
    
    # if not os.path.isdir(USER_VARIABLES.LOG_DIRECTORY):
    #     os.makedirs(USER_VARIABLES.LOG_DIRECTORY)
    
    # file_name= ''.join('all_losses'+'.txt' ) 
        
    # parameter_file = os.path.join(USER_VARIABLES.LOG_DIRECTORY, file_name)

    # with open(parameter_file, 'w') as pf:
    #     json.dump(all_losses, pf,indent = 4)
        
    # plot the mean model with sample datapoints   
     
        
        
        
# # 1369 ist die Probe auf der meine Annahmen basieren
#     load_and_plot(file, 
#                   extended_output = ['deltaGr',
#                                     'deltaCO2',
#                                     'deltaCH4',
#                                     'thermo',
#                                     'MM',
#                                     'v',
#                                     'deltaGs',
#                                     'inhibition',
#                                     'logQ',
#                                     'deltaH2',
#                                     'logFe3',
#                                     'logQH2O',
#                                     'dissH2O'],
#                                      which_plot = 'measured_and_modelled') # 'only_data' oder 'measured_and_modelled




    for sample_number, site_name in all_samples_and_sites: 
        #print('sample_number', sample_number)
        #print('site_name', site_name)
            
        speciemen_identifier = sample_number
      #  print(sample_number)
    #=============================================================================
        # run_and_plot(speciemen_identifier, 
        #               site = 'all', 
        #               extended_output = ['deltaGr',
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
        #                                 'Ac',
        #                                 'Hydro',
        #                                 'Homo'
        #                                 ])
    #=============================================================================
# # #=============================================================================
    
#         fit_model(speciemen_identifier, 
#                        site = site_name, 
#                       pathway_names = [
#                                            'Ferm_help',
#                                            'Ferm',
#                                            #'Fe3',
#                                            'Ac',
#                                            #'Hydro',
#                                            #'Homo'
#                                            ])
# # #=============================================================================
    
"""
  '13510',
  '13511', 
  '13512',
  '13520',
  '13521',
  '13530',
  '13531',
  '13670',
  '13671',
  '13672',
  '13690',
  '13691',
  '13692',
  '13700',
  '13701',
  '13702',
  '13720',
  '13721',
  '13722',
  '13730',
  '13731',
  '13732', 
  '13740',
  '13741',
  '13742',
  '13750',
  '13751',
  '13752',
  '13770',
  '13771',
  '13772',
  '13780',
  '13781',
  '13782',
  '13540',
  '13542',
  '13550',
  '13551',
  '13560',
  '13562',
  '13571', 
  '13572',
  '13580',
  '13581',
  '13590',
  '13591',
  '13600',
  '13602',
  '13610',
  '13612',
  '13620',
  '13622',
  '13630',
  '13631',
  '13640',
  '13641',
  '13650',
  '13651',
  '13652',
  '13661',
  '13662',
  '13680',
  '13681',
  '13682',
  '13710',
  '13711',
  '13712',
  '13760',
  '13761',
  '13762',
  '13790',
  '13791',
  '13792',
  '13800',
  '13801',
  '13802'
 
 
 Rep_No_CH4 = [13560, 13562, 13580, 13581, 13590, 13591, 13600, 13602, 13622, 13641]
 
"""
    #
  #  
    
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
