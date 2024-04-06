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
                    file_name = f'fit_{fit_repl}_{model_type}_{best_loss:.3g}' + timestamp
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
    
    
def boxplots():
    target_directory = USER_VARIABLES.LOG_DIRECTORY
    
    loaded_parameters = {'simple': {}, 
                         'complex': {}}
    for f in os.listdir(target_directory):
        file = os.path.join(target_directory, f)
        if not os.path.isfile(file) or not f.endswith('.json'):
            continue
        if not f.startswith('fit_'):
            continue
        
        model_type = 'simple' if 'simple' in f else 'complex'
        
        with open(file, 'r') as pf:
            parameters = json.load(pf)['parameters']
            
        for p, v in parameters.items():
            if not p in loaded_parameters[model_type]:
                loaded_parameters[model_type][p] = []
            loaded_parameters[model_type][p].append(v)
        
    
    parameter_names = list(loaded_parameters['simple'].keys()) + list(loaded_parameters['complex'].keys()) 
    parameter_names = list(set(parameter_names)) # unique names

    parameter_groups = {'pools': [],
                        'microbes': [],
                        'CUE': [],
                        'Kmb': [],
                        'Km': [],
                        'v_max': []}
    for p in parameter_names:
        if p == 'death_rate': 
            continue
        if not '_' in p:
            parameter_groups['pools'].append(p)
        elif p.startswith('M_'):
            parameter_groups['microbes'].append(p)
        elif 'CUE' in p:
            parameter_groups['CUE'].append(p)
        elif 'Kmb' in p:
            parameter_groups['Kmb'].append(p)
        elif 'Km' in p:
            parameter_groups['Km'].append(p)
        elif 'v_max' in p:
            parameter_groups['v_max'].append(p)
            
    d = .25
    for group_name, group in parameter_groups.items():
        plt.figure()
        i = 0
        plt.title(group_name)
        x_tick_labels = []
        x_tick_pos = []
        for par_name in group:
            pos_simple = 2*i + 1.5 - d
            pos_complex = 2*i + 1.5 + d
            tick_pos = 2*i + 1.5
            i += 1
            
            parameter_name = par_name.replace(group_name, '').replace('_', ' ')
            x_tick_labels.append(parameter_name)
            x_tick_pos.append(tick_pos)
            
            simple_data = np.nan
            if par_name in loaded_parameters['simple']:
                simple_data = [v for v in loaded_parameters['simple'][par_name]]
            complex_data = [v for v in loaded_parameters['complex'][par_name]]

            box_data = [simple_data, 
                        complex_data]
            plt.boxplot(box_data, 
                        positions = [pos_simple, pos_complex],
                        widths = 2*d)
            
        plt.xticks(rotation=90)
        plt.xticks(x_tick_pos, x_tick_labels)
        if not group_name == 'CUE':# and not group_name == 'pools':
            plt.yscale('log')

    plt.show()    
    
if __name__ == '__main__':
    #boxplots()
    fit()