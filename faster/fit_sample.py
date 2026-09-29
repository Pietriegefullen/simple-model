import os

import model
import data
import optimizer
import parameters
import checkpoint
import hashing

# TODO: for converting, derive configuration from checkpoint file and folder
# TODO: is a timestamp meaningful? only for CP, not folder
# TODO: parse initial parameter range (narrow <best_N>, or other criteria!)
# TODO: parse loss configuration (weights, other reduction functions)
# TODO: handle few usable sample points!!!
# TODO: warn if loaded parameters have incompatible origin -> input()
# TODO: in Objective, determine t_start, t_end for all loss contributions
#       predict only as necessary.
# TODO: save hyperparameters with every plot (how?) -> maintain origin: model version, ...
# TODO: if sample has only two replicas, fit_mode split is equivalent to single.
#       => only the meaning of validation changes.

def fit(chosen, objective_config, algo_config, init_config, 
        verbose_callback = False,
        convert = False):

    # get sample from dataset
    dataset = data.get_data_before_day()
    sample = dataset[chosen['sample']]
    split = sample.get_split(chosen['validation_replica'], 
                             chosen['fit_mode'])
    
    # build model and configure parameters
    pathway_model = model.Model(chosen['pathways'])
    pathway_model.parameters().set('default', normalized = chosen['normalized_parameters'])
    
    model_id = hashing.build_model_id(pathway_model.get_config(only_structure = True))
    
    if not 'file' in init_config and (not 'model' in init_config or init_config['model'] is None):
        init_config['model'] = model_id
    
    legacy_path = None if not 'file' in init_config else os.path.split(init_config['file'])[0]
    initial_parameters = parameters.load_parameters(init_config)
    pathway_model.parameters().set(initial_parameters)

    
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
            transform = parameters.get_transform(objective_config['transform'][pool])
            pool_loss = optimizer.get_loss_function(pool, 
                                                    objective_config['reduction'][pool], 
                                                    transform,
                                                    t_start = chosen['t_start'],
                                                    t_end = chosen['t_end'])
            replica_objective.add_loss(pool_loss, objective_config['loss_weight'][pool])
        replica_objectives.append(replica_objective)

    run_config= {'model': pathway_model.get_config(only_structure = True),
                 'chosen': chosen,
                 'objective': objective_config,
                 'algo': algo.get_config(),
                 'initial': pathway_model.parameters().get_config()
                 }
    if not legacy_path is None:
        run_config['legacy'] = legacy_path
    total_objective = sum(replica_objectives)
    total_objective.add_callback(checkpoint.CheckpointCallback(run_config, 
                                                               keep_only_n = 10,
                                                               verbose = verbose_callback))
    total_objective.add_callback(checkpoint.PrintCallback(run_config))
    
    if convert:
        total_objective(initial_parameters, transformed = False)
        return
    
    algo.minimize(total_objective, initial_parameters)

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
            'differential_evolution':   {}, # empty dict uses default
            'powell':                   {}
                }
    
    init_config = {
            'best_N':                   3,
            'sample':                   chosen['sample'],
            'validation_replica':       chosen['validation_replica'],
            'model':                    None,
            'run_ID':                   None,
            'file':                     None,
    }
    
    for omitted_pathway in args.omit:
        chosen['pathways'].remove(omitted_pathway)
    
    fit(chosen, objective_config, algo_config, init_config)
    
    