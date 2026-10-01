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

def run(run_config, initial_parameters, **kwargs):
    chosen = run_config['chosen']
    objective_config = run_config['objective']
    algo_config = {chosen['algorithm']: run_config['algo']}
    init_config = {'range': parameters.ModelParameters(run_config['range']),
                   'parameters': initial_parameters}
    run_log = fit(chosen, objective_config, algo_config, init_config, 
                  store_checkpoints = False,
                  minimize = False,
                  **kwargs)
    return run_log

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
    # build model and configure parameters
    pathway_model = model.Model(chosen['pathways'])
    pathway_model.parameters().set('default', normalized = chosen['normalized_parameters'])
    model_id = hashing.build_model_id(pathway_model.get_config(only_structure = True))
    
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
        pathway_model.parameters()[p_name].constant(p_value)
    
    # select optimiser
    algo = optimizer.get(chosen['algorithm'])
    algo.configure(algo_config[chosen['algorithm']])
    
    # build objective function
    replica_objectives = []
    for replica in split['fit']:
        replica_objective = optimizer.Objective(pathway_model, replica)
        for pool in ['CO2', 'CH4']:
            
            # make replica-provided parameters nan
            for name in ['H2O', 'CH4', 'CO2', 'TOC', 'DOC']:
                if name in initial_parameters:
                    del initial_parameters._parameters[name] 
                 
            
            tf = parameters.IdentityTransform()
            for t in objective_config['transform'][pool]:
                if t == 'normalize':
                    if pool == 'CO2':
                        _,pool_values = replica.CO2()
                    elif pool == 'CH4':
                        _,pool_values = replica.CH4()
                    else:
                        raise NotImplementedError()
                    values = tf.transform(pool_values)
                    finite_values = values[np.isfinite(values)]
                    replica_low = np.min(finite_values)
                    replica_high = np.max(finite_values)
                    tf = parameters.Normalization(tf, replica_low, replica_high)
                elif t == 'log':
                    tf = parameters.LogTransform(tf)
            
            pool_loss = optimizer.get_loss_function(pool, 
                                                    objective_config['reduction'][pool], 
                                                    tf,
                                                    t_start = chosen['t_start'],
                                                    t_end = chosen['t_end'])
            replica_objective.add_loss(pool_loss, objective_config['loss_weight'][pool])
        replica_objectives.append(replica_objective)

    run_config= {'model': pathway_model.get_config(only_structure = True),
                 'chosen': chosen,
                 'objective': objective_config,
                 'algo': algo.get_config(),
                 'range': initial_parameters.get_config(only_range = True)
                 }
    if not legacy_path is None:
        run_config['legacy'] = legacy_path
        run_config['legacy_file'] = legacy_file
    total_objective = sum(replica_objectives)
    
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
    return run_log

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
    
    
