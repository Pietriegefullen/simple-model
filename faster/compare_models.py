# -*- coding: utf-8 -*-
"""
Created on Sat Feb 24 13:20:56 2024

@author: Lara
"""
import os
import matplotlib.pyplot as plt
import numpy as np

import data
import optimizer
import OPTIMIZATION_PARAMETERS
import pathways
import predict
import main
import USER_VARIABLES


COLORS = {'CO2': 'r',
          'CH4': 'b'}

before = 1500 # TODO: test!
save_dir = 'before_' + str(before)

# TODO: unclear which parameters are fixed and which optimized.#
#       change s.t. either are given explicitly, terminology "default parameters", not "fixed"
# TODO: cleanup, collect settings
# TODO: daten beschneiden vor dem fit? -> delete after day ????
# TODO: same axis formats for scatter (only major tick labels)
# TODO: plot all you find in folder.

# TODO: compute R2 (Or other measures of goondess of fit)
# TODO: remove unused parameters (faster optimization?)
# TODO: log final loss value 
# TODO: use standard error of regression (=RMSE?) => value independent of number of samples

def build_replica_groups(sample_numbers):
    replica_groups = {}
    
    for sample in sample_numbers:
        assert len(sample) == 5
        sample_id, replica = sample[:-1], sample[-1]
        if not sample_id in replica_groups:
            replica_groups[sample_id] = []
        replica_groups[sample_id].append(sample)
    
    return replica_groups

model_pathways = {'complex': [pathways.Ferm,
                                  pathways.Ferm_help,
                                  pathways.Fe3,
                                  pathways.Hydro,
                                  pathways.Homo,
                                  pathways.Ac],
                      'simple': [pathways.Ferm,
                                 pathways.Ferm_help,
                                  pathways.Hydro,
                                  pathways.Ac],}

def file_exists(replica_name, model_type, save_dir = None):
    file = None
    
    file_ending = '_'.join(['specimen',
                            replica_name,
                            'site',
                            'all']) + '.json'
    if save_dir is None:
        save_dir = USER_VARIABLES.LOG_DIRECTORY
    else:
        parent, _ = os.path.split(USER_VARIABLES.LOG_DIRECTORY)
        save_dir = os.path.join(parent, save_dir)
    
   
    if not os.path.isdir(save_dir):
        return None
    
    found = False
    for f in os.listdir(save_dir):
        if f.endswith(file_ending) and f.startswith(model_type) and 'leave-one-out' in f:
            file = os.path.join(save_dir, f)
            found = True
            break
    if not found:
        return None
    return file

def fit_specimens():
    
    superdata = data.load_matlab(['superdata_2021_all'])
    sample_numbers = superdata.keys()
    all_specimen_groups = list(build_replica_groups(sample_numbers).values())
    all_specimen_groups = [g for g in all_specimen_groups if len(g) == 3]
    
    goodness = {}
    for specimen_replicas in all_specimen_groups:
        
        for replica in range(3):
            validation_replica = specimen_replicas[replica]
            training_replicas = [specimen_replicas[(replica+1)%3],
                                 specimen_replicas[(replica+2)%3]]
                        
            replica_default_parameters = {}
            for sn in training_replicas:
                replica_default_parameters[sn] = pathways.default_model_parameters(sn, 'all')
            
            for k,v in replica_default_parameters[training_replicas[0]].items():
                v2 = replica_default_parameters[training_replicas[1]][k]
                assert v == v2
            
            for model_type in ['simple', 'complex']:
                
                if file_exists('_'.join(training_replicas), model_type, save_dir):
                    print('skipping', model_type, 'for', validation_replica)
                    continue
                                
                print('fitting ' + model_type + 'to' + '_'.join(training_replicas))
                optimal_parameters = optimizer.fit_specimen(training_replicas, 
                                               'all',
                                               model_pathways[model_type], 
                                               replica_default_parameters[training_replicas[0]], 
                                               OPTIMIZATION_PARAMETERS.ALGORITHM,
                                               before = before)
                
                # save model parameters and pathways
                model_parameters = replica_default_parameters[training_replicas[0]]
                model_parameters.update(optimal_parameters)
                parameter_file =  main.save_model('_'.join(training_replicas), 
                                        'all',
                                        model_parameters, 
                                        prefix = model_type + '_leave-one-out',
                                        save_dir = save_dir)
                pathway_file = parameter_file.replace('.json', '_pathways.txt')
                with open(pathway_file, 'w') as pf:
                    for p in model_pathways[model_type]:
                        pf.write(p.__name__ + '\n')
                
def load_fitted():
    plot_specimens = ['13691']
    
    goodness = {}
    
    superdata = data.load_matlab(['superdata_2021_all'])
    sample_numbers = superdata.keys()
    all_specimen_groups = list(build_replica_groups(sample_numbers).values())
    all_specimen_groups = [g for g in all_specimen_groups if len(g) == 3]
    
    save_name = {}
    for specimen_replicas in all_specimen_groups:
        for replica in range(3):
            validation_replica = specimen_replicas[replica]
            training_replicas = [specimen_replicas[(replica+1)%3],
                                 specimen_replicas[(replica+2)%3]]
            save_name[validation_replica] = '_'.join(training_replicas)

    for specimen_number in sample_numbers:
        if not specimen_number in plot_specimens and not specimen_number[:4] in plot_specimens:
            continue
        
        for model_type in ['simple', 'complex']:
            
            if not specimen_number in save_name:
                print('specimen not saved:', specimen_number)
                continue
            specimen_name = save_name[specimen_number]
            file = file_exists(specimen_name, model_type, save_dir)
            if file is None:
                print('Could not find saved parameters for specimen ', specimen_number, ' (' + model_type + ')')
                continue

            default_model_parameters = pathways.default_model_parameters(specimen_number, 'all')
            model_parameters = main.load_model_parameters(file)
            print(specimen_number, model_type, 'loaded model parameters')
            
            for default_key, default_value in default_model_parameters.items():
                if not default_key in model_parameters.keys():
                    model_parameters[default_key] = default_value            
                    
            # run model with optimal parameters on validation replica
            validation_replica_data = data.specimen_data(specimen_number, 'all')
            
            if not before is None:
                days_before = []
                for d in validation_replica_data['measured_time']:
                    if d < before:
                        days_before.append(d)
                    else:
                        break
                validation_replica_data['measured_time'] = np.array(days_before)
                validation_replica_data['CO2'] = validation_replica_data['CO2'][:len(days_before)]
                validation_replica_data['CH4'] = validation_replica_data['CH4'][:len(days_before)]
            
            
            measurement_days = validation_replica_data['measured_time']
            print(specimen_number, model_type, 'predicting')
            pools = predict.predictor(measurement_days,
                                      model_parameters,
                                      model_pathways[model_type])
            
            model_results = {'specimen_number': specimen_number,
                             'model_type': model_type,
                             'optimal_parameters': model_parameters,
                             'days': measurement_days,
                             'measured_CO2': validation_replica_data['CO2'],
                             'measured_CH4': validation_replica_data['CH4'],
                             'CO2': pools['CO2'],
                             'CH4': pools['CH4'],
                             'CO2_R2': None,
                             'CH4_R2': None}
            
            if not specimen_number in goodness:
                goodness[specimen_number] = {}
            goodness[specimen_number][model_type] = model_results
            
        if not specimen_number in goodness:
            print(specimen_number, ' results missing')
            continue
        
        if not 'complex' in goodness[specimen_number] or not 'simple' in goodness[specimen_number]:
            print(specimen_number, 'missing model')
            continue
        
        print(specimen_number, ' plotting pools')
        plot_time_series(goodness[specimen_number])
                
        print(specimen_number, ' plotting correlation')
        for model_type, model_results in goodness[specimen_number].items():
            plot_scatter(model_results, specimen_number, model_type)
        
    return goodness

def plot_time_series(specimen_results):
    for pool in ['CO2', 'CH4']:
        fig = plot_pool(specimen_results['complex'], pool, 
                        model_line = 'k-')
        _   = plot_pool(specimen_results['simple'], pool, 
                        model_line = 'k--',
                        fig = fig,
                        plot_measured = False)
        
def plot_pool(model_results, pool ,model_line = 'k-', fig = None, plot_measured = True):
    if fig is None:
        fig = plt.figure()
    else:    
        fig = plt.figure(fig.number)
        
    if plot_measured:
        plt.plot(model_results['days'],
                 model_results['measured_' + pool], COLORS[pool]+'x',
                 label = 'measured '+ pool)
    plt.plot(model_results['days'],
             model_results[pool], model_line,
             label = 'modelled ' + pool + ' ('+model_results['model_type']+')')
    plt.title(' '.join([pool, 
                        model_results['specimen_number']]))
    ax = plt.gca()
    ax.legend()
    ax.set_ylim([0,40])
    plt.xlabel('day')
    plt.ylabel(pool)
    return fig

def plot_goodness():
    plt.close('all')
    goodness = load_fitted()

    for sample_number, models in goodness.items():
        for model_type, model_results in models.items():
            plot_scatter(model_results, sample_number, model_type)
        
        for pool in ['CO2', 'CH4']:
            if not 'complex' in goodness[sample_number] or not 'simple' in goodness[sample_number]:
                continue
            fig = plot_pool(goodness[sample_number]['complex'], pool, 
                            model_line = 'k-')
            _   = plot_pool(goodness[sample_number]['simple'], pool, 
                            model_line = 'k--',
                            fig = fig,
                            plot_measured = False)
                    
   
def plot_scatter(model_results, specimen_number, model_type):

    plt.figure()
    
    max_value = 0
    for gas in ['CO2', 'CH4']:    
        measured = model_results['measured_' + gas]
        modelled = model_results[gas]
        max_value = max(max_value, max(max(measured), max(modelled)))
        
        plt.plot(measured, modelled, COLORS[gas] + '.', label = gas)
        
    plt.plot([0, max_value], [0, max_value], 'k-', linewidth=1) # plot diagonal
    plt.xlabel('measured')
    plt.xlabel('modelled')
    
    ax = plt.gca()
    ax.legend()
    ax.set_xlim([0, max_value])
    ax.set_ylim([0, max_value])
    ax.set_aspect('equal', 'box')
    plt.title(' '.join([specimen_number, model_type]))




if __name__ == '__main__':
    #fit_specimens()
    load_fitted()



