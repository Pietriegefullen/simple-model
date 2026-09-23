import argparse

import model
import data
import optimizer
import parameters
import checkpoint
import hashing


parser = argparse.ArgumentParser(
                    prog='fit_sample',
                    description='Fits a model to replica data.',
                    epilog='')
parser.add_argument('sample', type = int)
parser.add_argument('validation_replica', type = int)

# omit pathways
parser.add_argument('--omit', nargs = '+', default = [])

# override parameters/switches

args = parser.parse_args()

# TODO: read args into config

# TODO: __str__ for Objective/Loss

# TODO: finish init config!
# TODO: warn if loaded parameters have incompatible origin -> input()
# TODO: loaded parameter range sets all to variable?

# TODO: adding two objectives returns Objective, not Addable? => Objective IS Addable.
#TODO: set model variable/constant parameters , initial values LATER
# TODO: save hyperparameters with every plot (how?) -> maintain origin: model version, ...

chosen = {
            'sample':                   1351,
            'validation_replica':       4, 
            
            't_start':                  None,
            't_end':                    None,
            
            'fit_mode':                 'split', # 'single' or 'split'
            'pathways':                 ['Hydrolysis',
                                         'Fermentation',
                                         'Hydro',
                                         'Aceto',
                                         'Homo',
                                         'Fe3'],
            
            # set specified parameters to constant value
            'parameter_override':       {},  # e.g. 'Homo_thermodynamics': False
            
            'normalized_parameters':    True,
            'algorithm':                'differential_evolution',
            
            'loss_weight':              {'CO2': 1.,
                                         'CH4': 1.},
            'reduction':                {'CO2': 'mse',
                                         'CH4': 'mse'},
            'transform':                {'CO2': ['normalize'],
                                         'CH4': ['log', 'normalize']},
            }

algo_config = {'differential_evolution': {}, # empty dict uses default
               'powell':                 {}
               }

chosen['sample'] = args.sample
chosen['validation_replica'] = args.validation_replica
for omitted_pathway in args.omit:
    chosen['pathways'].remove(omitted_pathway)
    


# get sample from dataset
dataset = data.get_data_before_day()
sample = dataset[chosen['sample']]
split = sample.get_split(chosen['validation_replica'], 
                         chosen['fit_mode'])

# build model and configure parameters
pathway_model = model.Model(chosen['pathways'])
pathway_model.parameters().set('default', normalized = chosen['normalized_parameters'])

model_id = hashing.build_model_id(pathway_model.get_config(only_structure = True))

init_config = {
                'best_N':               3,
                'sample':               chosen['sample'],
                'validation_replica':   chosen['validation_replica'],
                'model':                model_id,
                'run_ID':               None,
                
}
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
        transform = parameters.get_transform(chosen['transform'][pool])
        pool_loss = optimizer.get_loss_function(pool, 
                                                chosen['reduction'][pool], 
                                                transform,
                                                t_start = chosen['t_start'],
                                                t_end = chosen['t_end'])
        replica_objective.add_loss(pool_loss, chosen['loss_weight'][pool])
    replica_objectives.append(replica_objective)

run_config= {'model': pathway_model.get_config(only_structure = True),
             'chosen': chosen,
             'algo': algo.get_config(),
             'initial': pathway_model.parameters().get_config()}
total_objective = sum(replica_objectives)
total_objective.add_callback(checkpoint.CheckpointCallback(run_config))
total_objective.add_callback(checkpoint.PrintCallback(run_config))

algo.minimize(total_objective, initial_parameters)
