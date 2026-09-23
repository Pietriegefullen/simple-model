import os
import json

import hashing

from USER_VARIABLES import RESULTS_DIRECTORY as CP_ROOT

checkpoint_plain = [('model', hashing.build_model_id), 
                    ('chosen', 'sample'), 
                    ('chosen', 'validation_replica'), 
                    ('chosen', 'fit_mode')
                    # run ID is always appended
                    ]
# configure callback:
#   keep only N best, if N is None, keep all
#   store location

def get_value(d, ks, hash_function = None):
    if not isinstance(ks, (list, tuple)):
        ks = [ks]
    if len(ks) == 2 and callable(ks[-1]):
        ks, hash_function = ks
        return hash_function(d[ks])
    
    if len(ks) > 1:
        hash_function = None
        return get_value(d[ks[0]], ks[1:], hash_function)
    
    else:
        return d[ks[0]]
    
    raise Exception('Should not happen.')

class Callback():
    def __init__(self, run_config):
        self.objective = None
        self.target_directory = None
        self.run_config = run_config
        self.run_id = hashing.build_run_id(run_config)

    def set_objective(self, objective):
        self.objective = objective
        model_id = hashing.build_model_id(objective.model().get_config())
        self.target_directory = os.path.join(CP_ROOT, self.run_dir())

    def run_dir(self):
        plain = []
        for keys in checkpoint_plain:
            value = get_value(self.run_config, keys)
            #value = str(self.run_config['chosen'][key])
            plain.append(value)
        plain.append(self.run_id)
        return '_'.join([str(p) for p in plain])

class PrintCallback(Callback):
    def __init__(self, run_config):
        super().__init__(run_config)
        self.print_r2 = False
        
    def __call__(self):
        args, kwargs, loss_value = self.objective.last_call()
        _, _, best_loss = self.objective.best_call()
        cnt = self.objective.call_count()
        if self.print_r2:
            print(self.objective)
            print(f'call {cnt:6d}: loss value {loss_value:8.3g}, best loss: {best_loss:8.3g}')
        else:
            print(f'\rcall {cnt:6d}: loss value {loss_value:8.3g}, best loss: {best_loss:8.3g}', end = '')

class CheckpointCallback(Callback):
    def __init__(self, run_config, keep_only_n = None):
        super().__init__(run_config)
        self.keep_only_n = keep_only_n

    def cleanup(self, save_dir):
        if not self.keep_only_n is None:
            raise NotImplementedError()
            # sort by loss value
            # keep only 
            # skip config file
        
    def __call__(self):
        cp_transformed_parameters, _, last_loss = self.objective.last_call()
        _,_, best_loss = self.objective.best_call()
        
        if not last_loss == best_loss:
            return

        if not os.path.isdir(self.target_directory):
            os.makedirs(self.target_directory)

        checkpoint_parameters = self.objective.model().parameters().get_config()
        
        checkpoint_data = {
                            'parameters': checkpoint_parameters,
                            'total_loss': last_loss,
                            'run_config': self.run_config,
                           }
       
        f_loss = f'{last_loss:.3f}'.replace('.', '')
        if len(f_loss) > 8:
            # overflow
            f_loss = '9'*8
        else:
            while len(f_loss) < 8:
                f_loss = '0' + f_loss
        
        file_name = self.run_dir() + f'_loss-{f_loss}'
        checkpoint_file = os.path.join(self.target_directory, file_name)
        with open(checkpoint_file, 'w') as cf:
            json.dump(checkpoint_data, cf, indent = 4)

        self.cleanup(self.target_directory)

