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
    
if __name__ == '__main__':
    fit(['1358'])
