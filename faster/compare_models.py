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

# TODO: should initial acetate and initial Fe3 be linear or log?
# TODO: check matlab data and compare. why different?
# TODO: compute measures of fit (Whose responsibility?)
# TODO: load model, then continue optimizing with gradient method?
# TODO: remove unused pools/variables, e.g. initial Fe3 in simple model.
# TODO: 1354-6 has no carex but usable only up to day 1309?
# TODO: timeout for solve_ivp in case parameter combination causes too small time steps?
def fit(include_samples = None, exclude_samples = None):
    dataset = data.get_data_before_carex()
    target_directory = USER_VARIABLES.LOG_DIRECTORY
        
    for sample in dataset.samples:
        skip_sample = False # used for KeyboardInterrupt
        if include_samples and not str(sample.sample_name) in include_samples:
            print('skipping', str(sample))
            continue
        
        if exclude_samples and str(sample.sample_name) in exclude_samples:
            print('skipping', str(sample))
            continue
        
        try:
            splits = sample.leave_one_out_split()
            
        except Exception as ex:
            print('skipping', str(sample), str(ex))
            continue
        
        for split in splits:
            if skip_sample: break
                
            fit_replicas = split['fit']
            validation_replica = split['val']
            
            for model_type in ['simple', 'complex']:
                if skip_sample: break
                print()
                print(f'fitting {model_type} to {str(sample)}')
                
                chosen_pathways = model.get_pathways(model_type)
                pathway_model = model.Model(chosen_pathways)
                pathway_model.parameters().set('default')
                
                if not 'Fe3' in chosen_pathways:
                    pathway_model.parameters()['Fe3'].constant(0)
                if not 'Homo' in chosen_pathways:
                    pathway_model.parameters()['M_Homo'].constant(0)

                try:
                    best_loss, _ = pathway_model.fit(fit_replicas)
                    
                    fit_repl = '_'.join([str(r) for r in fit_replicas])
                    now = datetime.now()
                    timestamp = '_' + now.strftime('%Y-%m-%d_%H-%M-%S')
                    file_name = f'fit_{fit_repl}_{model_type}_{best_loss:.3g}_from_{fit_from:d}_to_{str(fit_to)}' + timestamp
                    pathway_model.save(target_directory, file_name)
                    
                    # plot fit over each fit replica
                    for replica in fit_replicas:
                        fit_results = pathway_model.predict(replica)
                        plt.close('all')
                        fit_results.plot(['CH4', 'CO2'], newfigure = False)
                        replica.plot()
                        figure_name = file_name.replace(fit_repl, str(replica)) + '.svg'
                        plt.savefig(os.path.join(target_directory, figure_name))

                    # plot fit over validation replica
                    val_results = pathway_model.predict(validation_replica)
                    plt.close('all')
                    val_results.plot(['CH4', 'CO2'], newfigure = False)
                    validation_replica.plot()
                    figure_name = file_name + '.svg'
                    val_repl = str(validation_replica)
                    plt.savefig(os.path.join(target_directory, 
                                             figure_name.replace('fit_' + fit_repl, 
                                                                 'val_' + val_repl)))
                
                except KeyboardInterrupt as ex:
                    while True:
                        i = input('> ')
                        if i == 'q':
                            return
                        break
                    print('skipping', str(sample))
                    skip_sample = True
                
                except Exception as ex:
                    print()
                    print('while trying to fit ', str(sample))
                    print(traceback.format_exc())
                    print()
   
def fit_sample(sample_name, split_number, model_type, log_co2 = False, log_ch4 = False, confirm = False,
               fit_from = 0, fit_to = None, normalized_parameters = False,
               loss_weight_CO2 = 1, loss_weight_CH4 = 1):
    target_directory = USER_VARIABLES.LOG_DIRECTORY
    d = data.get_data_before_day()
    sample = d[sample_name]
    
    try:
        splits = sample.leave_one_out_split()
        fit_replicas = splits[split_number]['fit']
        
    except Exception as ex:
        print('skipping', str(sample), str(ex))
        return
    
    chosen_pathways = model.get_pathways(model_type)
    pathway_model = model.Model(chosen_pathways)
    pathway_model.parameters().set('default', normalized = normalized_parameters)
    
    if not 'Fe3' in chosen_pathways:
        pathway_model.parameters()['Fe3'].constant(0)
    if not 'Homo' in chosen_pathways:
        pathway_model.parameters()['M_Homo'].constant(0)

    best_loss, _ = pathway_model.fit(fit_replicas, log_co2 = log_co2, log_ch4 = log_ch4, 
                                     fit_from = fit_from, fit_to = fit_to,
                                     loss_weight_CO2 = loss_weight_CO2, 
                                     loss_weight_CH4 = loss_weight_CH4)
    
if __name__ == '__main__':
    default_sample = 1367
    default_model_type = 'complex'
    fit_from = 0
    fit_to = None
    normalized_parameters = True
    loss_weight_CO2 = 0.1
    loss_weight_CH4 = 1.0
    
    log_co2 = False
    log_ch4 = True
    # None means no fit at all!
   
    hasargs = False
    if len(sys.argv) == 1:
        sample = default_sample
    else:
        hasargs = True
        sample = sys.argv[1]
        
    split = 0
    if '0' in sys.argv:
        split = 0
    elif '1' in sys.argv:
        split = 1
    elif '2' in sys.argv:
        split = 2
        
    #log = True
    #if 'log' in sys.argv:
    #    log = True
    #elif 'lin' in sys.argv:
    #    log = False
    #    raise Exception('Non-log fitting. Check what you are doing.')
        
    model_type = 'simple'
    if 'complex' in sys.argv:
        model_type = 'complex'
    elif 'simple' in sys.argv:
        model_type = 'simple'
    else:
        model_type = default_model_type
    
    if 'from' in sys.argv:
        idx = sys.argv.index('from')
        fit_from = int(sys.argv[idx+1])

    if 'to' in sys.argv:
        idx = sys.argv.index('to')
        fit_to = int(sys.argv[idx+1])
    
    fit_sample(sample, split, model_type, log_co2 = log_co2, log_ch4 = log_ch4, confirm = hasargs,
               fit_from = fit_from,
               fit_to = fit_to,
               normalized_parameters = normalized_parameters,
               loss_weight_CO2 = loss_weight_CO2,
               loss_weight_CH4 = loss_weight_CH4)
