import os

import argparse
import fit_sample
import user_input
from USER_VARIABLES import PROJECT_DIRECTORY as ROOT

# TODO: for converted checkpoints, the run ID is identical
#       regardless of parameter ranges.
#       but Parameter ranges should be identical for the same run?
#       for legacy checkpoints, use default ranges?

# TODO: make sure run ID is identical for EQUIVALENT definitions.
#       algo_config is not used for any ID.
#       should t_start, t_end be used for objective ID?

parser = argparse.ArgumentParser(
                    prog='convert checkpoint',
                    description='',
                    epilog='')
parser.add_argument('path', nargs = '+') # override
parser.add_argument('--exclude', nargs = '+') # override
parser.add_argument('--p', nargs = '+') # override

args = parser.parse_args()
p = args.path
x = args.exclude
override = args.p

def get_checkpoint_files(path, criteria = None, exclude = None):
    file_list = []
    for d in os.listdir(path):
        if d == 'results': continue
        if d.startswith('.') or d.startswith('__'): continue
    
        folder = os.path.join(path, d)
        if not os.path.isdir(folder): 
            continue

        if d.startswith('fit_'):
            for f in os.listdir(folder):
                file_path = os.path.join(path, d, f)
                      
                pos = criteria is None or all([c in file_path for c in criteria])
                neg = exclude is None or not any([c in file_path for c in exclude])
                
                if pos and neg:
                    file_list.append(file_path)

        else:
            file_list += get_checkpoint_files(folder, criteria, exclude)
    return file_list
    
candidates = get_checkpoint_files(ROOT, p, x)
pl = "" if len(candidates) == 1 else "s"
print(f'Found {len(candidates):d} candidate{pl}.')

input()

for candidate in candidates:
    setattr(args, 't_start', None)
    setattr(args, 't_end', None)
    setattr(args, 'override', {})
    setattr(args, 'local', False)
    setattr(args, 'omit', [])

    series, folder_name = os.path.split(os.path.split(candidate)[0])
    replicas = []
    for token in folder_name.replace('fit_', '').split('_'):
        if not user_input.numeric(token):
            break
        replicas.append(str(token))
    
    if len(replicas) == 1:
        sample = replicas[0][:-1]
        validation_replica = replicas[0][-1]
        single = True
        
    else:
        samples = {r[:-1] for r in replicas}
        assert len(samples) == 1
        sample = samples.pop()
        fit_replicas = [r[-1] for r in replicas]
        validation_replica = next(r for r in '456' if r not in fit_replicas)
        single = False
        
    setattr(args, 'sample', replicas[0][:-1])
    setattr(args, 'validation_replica', replicas[0][-1])
    setattr(args, 'single', single)

    if not 'complex' in folder_name:
        raise NotImplementedError()
    
    if '_0-' in series or 'penalty' in series:
        raise NotImplementedError()
    
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
            'file':                     candidate,
            'range':                    'default',
    }
    
    for omitted_pathway in args.omit:
        chosen['pathways'].remove(omitted_pathway)
    
    fit_sample.fit(chosen, objective_config, algo_config, init_config, 
                   verbose_callback = True,
                   convert = True)
    