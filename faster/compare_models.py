import sys
import os
import traceback
import matplotlib.pyplot as plt
from datetime import datetime
import json
import numpy as np

import model
import data
import USER_VARIABLES

def fit_sample(sample_name, val_replica_number, model_type, 
               log_co2 = False, log_ch4 = False, confirm = False,
               fit_from = 0, fit_to = None, normalized_parameters = False,
               loss_weight_CO2 = 1, loss_weight_CH4 = 1, 
               loss_function_co2 = 'mse', 
               loss_function_ch4 = 'mse',
               parameter_range = None,
               local_search = False, 
               initial_parameters = None,
               weighted_measurements = False,
               rate_penalty = 0,
               normalized = False,
               initial_mean_days = 0,
               fit_mode = 'split'):
    
    d = data.get_data_before_day()
    sample = d[sample_name]
    
    print('fitting')
    print(str(sample))
    print(val_replica_number)
    print(model_type)
    print('from', fit_from, 'to', fit_to)
    print('local', local_search)
    
    fit_replicas = None
    try:
        if fit_mode == 'split':
            splits = sample.leave_one_out_split()
        elif fit_mode == 'single':
            splits = [{'fit':replica, 'val': replica}
                      for replica in sample.replicas]
        else:
            raise NotImplementedError()
        
        for s in splits:
            if int(s['val'].replica_number) == int(val_replica_number):
                fit_replicas = s['fit']
                break
        if fit_replicas is None:
            raise Exception('Split could not be determined')
            
    except Exception as ex:
        print('skipping', str(sample), str(ex))
        return
    
    chosen_pathways = model.get_pathways(model_type)
    pathway_model = model.Model(chosen_pathways)
    pathway_model.parameters().set('default', normalized = normalized_parameters)

    if not initial_parameters is None:
        for k,v in initial_parameters.items():
            pathway_model.parameters()[k].set(float(v))
    
    if not 'Fe3' in chosen_pathways:
        pathway_model.parameters()['Fe3'].constant(0)
    if not 'Homo' in chosen_pathways:
        pathway_model.parameters()['M_Homo'].constant(0)

    best_loss, _ = pathway_model.fit(fit_replicas, log_co2 = log_co2, log_ch4 = log_ch4, 
                                     fit_from = fit_from, fit_to = fit_to,
                                     loss_weight_CO2 = loss_weight_CO2, 
                                     loss_weight_CH4 = loss_weight_CH4,
                                     loss_function_co2 = loss_function_co2,
                                     loss_function_ch4 = loss_function_ch4,
                                     parameter_range = parameter_range,
                                     algorithm = None if not local_search else 'Powell',
                                     weighted_measurements = weighted_measurements,
                                     rate_penalty = rate_penalty,
                                     normalized = normalized,
                                     suffix = fit_mode,
                                     initial_mean_days = initial_mean_days)
    
if __name__ == '__main__':
    from main_file import load_parameter_range, load_fitted_parameters
    default_sample = 1351
    default_val_replica_number = 4
    default_model_type = 'complex'
    
    default_fit_mode = 'split' #'split'
    
    default_initial_mean_days = 0
            
    fit_from = 0
    fit_to = None
    
    normalized_parameters = True
    loss_weight_CO2 = 1.0
    loss_weight_CH4 = 1.0
    loss_function_co2 = 'mse'
    loss_function_ch4 = 'mse'
    normalized = True
    narrower_range = True
    best_N = 8
    local_search = False
    weighted_measurements = False
    model_type = default_model_type
    rate_penalty = 0.#1000
    
    log_co2 = False
    log_ch4 = True
   
    hasargs = len(sys.argv) > 1
    
    sample = sys.argv[1] if hasargs else default_sample
    
    if not hasargs:
        val_replica_number = default_val_replica_number

    else:        
        val_replica_number = int(sys.argv[2])
        
        #if val_replica_number < 4 or val_replica_number > 6:
        #    raise Exception()
            
        if 'narrow' in sys.argv:
            narrower_range = True
            best_N = int(sys.argv[sys.argv.index('narrow')+1])
            
        else:
            narrower_range = False
    
        if  'from' in sys.argv:
            idx = sys.argv.index('from')
            fit_from = int(sys.argv[idx+1])
    
        if  'to' in sys.argv:
            idx = sys.argv.index('to')
            fit_to = int(sys.argv[idx+1])
    
    suffix = ''
    if not fit_from == 0 or not fit_from is None:
        s_fit_to = str(None) if fit_to is None else str(int(fit_to))
        suffix = '_' + str(int(fit_from)) + '-' + s_fit_to
    if suffix == '_0-None':
        suffix = ''
        
        
    fit_mode = default_fit_mode
    if (hasargs and 'split' in sys.argv) or default_fit_mode == 'split':
        fit_mode = 'split'
    elif (hasargs and 'single' in sys.argv) or default_fit_mode == 'single':
        fit_mode = 'single'
        suffix = '_'.join([suffix, 'single'])
        
    initial_mean_days = default_initial_mean_days
    if hasargs and 'init' in sys.argv:
        idx = sys.argv.index('init')
        initial_mean_days = int(sys.argv[idx+1])
        
    loaded_range = None
    replica_name = str(val_replica_number)
    if narrower_range:
        try:
            loaded_range = load_parameter_range(sample, 
                                                replica_name, 
                                                model_type,
                                                best_N = best_N,
                                                suffix = suffix)

        except Exception as ex:
            if 'single result file' in str(ex):
                pass
            else:
                raise ex

    if hasargs and 'simple' in sys.argv:
        raise Exception('Using simple model? Why?')
    
    if not hasargs:
        pass
        
    elif 'local' in sys.argv:
        local_search = True
    
    else:
        local_search = False
        
   
    best_parameters = None
    if local_search:
        best_parameters = load_fitted_parameters(sample, 
                                                 replica_name, 
                                                 model_type,
                                                 best = None)
    if best_parameters is None and not loaded_range is None:
        best_parameters = loaded_range.as_dict()
        
    
    fit_sample(sample, val_replica_number, model_type, log_co2 = log_co2, log_ch4 = log_ch4, confirm = hasargs,
               fit_from = fit_from,
               fit_to = fit_to,
               normalized_parameters = normalized_parameters,
               loss_weight_CO2 = loss_weight_CO2,
               loss_weight_CH4 = loss_weight_CH4,
               loss_function_co2 = loss_function_co2,
               loss_function_ch4 = loss_function_ch4,
               parameter_range = loaded_range,
               local_search = local_search,
               initial_parameters = best_parameters,
               weighted_measurements = weighted_measurements,
               rate_penalty = rate_penalty,
               normalized = normalized,
               initial_mean_days = initial_mean_days,
               fit_mode = fit_mode)

