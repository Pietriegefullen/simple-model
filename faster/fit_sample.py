import os

import model
import data
import optimizer
import parameters
import checkpoint
import hashing
import numpy as np

# TODO: handle few usable sample points!!!
# TODO: in Objective, determine t_start, t_end for all loss contributions
#       predict only as necessary.
# TODO: complete init_config specification for fits using new checkpoints
# TODO: loss normalization is not equal for loaded checkpoints?

# TODO: model id is GJT, not AJ???
# TODO: checkpoint stores only if better in target directory!
#       => allows -1 workers.

# TODO: plot all replicas or single according to fit_mode.

def run(run_config, initial_parameters, **kwargs):
    chosen = run_config['chosen']
    objective_config = run_config['objective']
    algo_config = {chosen['algorithm']: run_config['algo']}
    init_config = {'range': parameters.ModelParameters(run_config['range']),
                   'parameters': initial_parameters}
    pathway_model, objective, run_log = fit(chosen, objective_config, algo_config, init_config, 
                                          store_checkpoints = False,
                                          minimize = False,
                                          **kwargs)
    return pathway_model, objective, run_log

def fit(chosen, objective_config, algo_config, init_config, 
        store_checkpoints = True,
        verbose_callback = False,
        minimize = True, 
        cp_target = None):
    # get sample from dataset
    dataset = data.get_data_before_day()
    sample = dataset[chosen['sample']]
    split = sample.get_split(chosen['validation_replica'], 
                             chosen['fit_mode'])
    
    legacy_path = None
    if 'parameters' in init_config:
        initial_parameters = init_config['parameters']
        
        if 'range' in init_config:
            parameter_range = init_config['range']
            _ = [parameter_range[p.name].set(p) for p in initial_parameters]
            initial_parameters = parameter_range
            
    else:
        if not 'file' in init_config and (not 'model' in init_config or init_config['model'] is None):
            init_config['model'] = model_id
        
        legacy_file = init_config['file']
        legacy_path = None if legacy_file is None else os.path.split(legacy_file)[0]
        initial_parameters = parameters.load_parameters(init_config)
    
    # override model parameters
    for p_name, p_value in chosen['parameter_override'].items():
        initial_parameters[p_name].constant(p_value)
    
    # build model and configure parameters
    pathway_model = model.configure_model(chosen['pathways'],
                                          chosen['normalized_parameters'], 
                                          initial_parameters)

    model_id = hashing.build_model_id(pathway_model.get_config(only_structure = True))

    # select optimiser
    algo = optimizer.get(chosen['algorithm'])
    algo.configure(algo_config[chosen['algorithm']])
    
    total_objective = optimizer.build_objective_function(pathway_model,
                                                         split['fit'],
                                                         objective_config)

    # make replica-provided parameters nan
    for name in ['H2O', 'CH4', 'CO2', 'TOC', 'DOC']:
        if name in initial_parameters:
            del initial_parameters._parameters[name] 

    run_config = {'model': pathway_model.get_config(only_structure = True),
                 'chosen': chosen,
                 'objective': objective_config,
                 'algo': algo.get_config(),
                 'range': initial_parameters.get_config(only_range = True)
                 }
    if not legacy_path is None:
        run_config.update({'legacy': legacy_path, 'legacy_file': legacy_file})
    
    if store_checkpoints:
        total_objective.add_callback(checkpoint.CheckpointCallback(run_config, 
                                                                   keep_only_n = 10,
                                                                   verbose = verbose_callback,
                                                                   target = cp_target))
    total_objective.add_callback(checkpoint.PrintCallback(run_config))

    if not minimize:
        print(total_objective)
        total_objective(initial_parameters, transformed = False)
        print()

    else:    
        algo.minimize(total_objective, initial_parameters)
    
    run_log = total_objective.model().system_state_log

    return pathway_model, total_objective, run_log

if __name__ == '__main__':
    import user_input
    
    args = user_input.parse_args_fit()

    chosen = {
            'sample':                   args.sample,
            'validation_replica':       args.validation_replica, 
            
            't_start':                  args.t_start,
            't_end':                    args.t_end,
            
            'fit_mode':                 'single' if args.single else 'split',
            'pathways':                 ['Hydrolysis',
                                         'Fermentation',
                                         'Hydro',
                                         'Aceto',
                                         'Homo',
                                         'Fe3'],
            'parameter_override':       args.override,  # e.g. 'Homo_thermodynamics': False
            'normalized_parameters':    True,
            'algorithm':                'powell' if args.local else 'differential_evolution',
            }

    objective_config = {
            'loss_weight':              {'CO2': 1.,
                                         'CH4': 1.},
            'reduction':                {'CO2': 'mse',
                                         'CH4': 'mse'},
            'transform':                {'CO2': ['normalize'],
                                         'CH4': ['log', 'normalize']},
            }
    
    algo_config = {
            'differential_evolution':   {'workers' : 1}, # empty dict uses default
            'powell':                   {}
                }
    
    init_config = {
            'default':                  args.default,
            'best_N':                   None if args.default else args.best,
            'sample':                   None if args.default else chosen['sample'],
            'validation_replica':       None if args.default else chosen['validation_replica'],
            'model':                    None,
            'run_ID':                   None,
            'file':                     None,
    }
    
    for omitted_pathway in args.omit:
        chosen['pathways'].remove(omitted_pathway)
    
    fit(chosen, objective_config, algo_config, init_config)
    
    
