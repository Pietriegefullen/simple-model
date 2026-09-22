import os
import hashlib
import json

import wn

from USER_VARIABLES import RESULTS_DIRECTORY as CP_ROOT

ALPHABET = 'abcdefghijklmnopqrstuvwxyz'.upper()

checkpoint_plain = ['sample', 'validation_replica', 'fit_mode']

def freeze(obj):
    if isinstance(obj, dict):
        return tuple(sorted((k, freeze(v)) for k, v in obj.items()))
    if isinstance(obj, (list, tuple)):
        return tuple(freeze(v) for v in obj)
    if isinstance(obj, set):
        return tuple(freeze(v) for v in obj)
    return obj

def compute_hash(config):
    assert isinstance(config, dict)
    digest = hash(tuple(sorted(freeze(config))))
    return digest

def build_id(config, group_length = 3, length = 2):
    digest = compute_hash(config)
    
    number = abs(digest)

    result = []
    base = len(ALPHABET)

    while number:
        number, remainder = divmod(number, base)
        result.append(ALPHABET[remainder])

    s_hash = "".join(reversed(result)).zfill(length*group_length)[-length*group_length:]
    
    identifier = "-".join(
                        s_hash[i:i+group_length] 
                        for i in range(0, len(s_hash), group_length)
                    )
    return identifier

# configure callback:
#   keep only N best, if N is None, keep all
#   store location

class Callback():
    def __init__(self, run_config):
        self.objective = None
        self.target_directory = None
        self.run_config = run_config
        self.run_id = build_id(run_config, 3, 4)

    def set_objective(self, objective):
        self.objective = objective
        model_id = build_id(objective.model().get_config(), 3, 2)
        self.target_directory = os.path.join(CP_ROOT, self.run_dir())

    def run_dir(self):
        plain = [self.run_id]
        for key in checkpoint_plain:
            value = str(self.run_config['chosen'][key])
            plain.append(value)
        return '_'.join(plain)

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
            config_file = os.path.join(self.target_directory,
                                        '_'.join([self.run_id, 'config']))
            with open(config_file, 'w') as cf:
                json.dump(self.run_config, cf, indent = 4)

        variables = self.objective.model().parameters().variables()
        checkpoint_parameters = [var.inverse_transform(p) 
                            for var, p in zip(variables, cp_transformed_parameters)]
        
        checkpoint_data = {
                            'parameters': checkpoint_parameters,
                            'total_loss': last_loss,
                           }
       
        file_name = self.run_dir() + f'_loss-{round(last_loss*1000):06d}'
        checkpoint_file = os.path.join(self.target_directory, file_name)
        with open(checkpoint_file, 'w') as cf:
            #json.dump(checkpoint_data, cf, indent = 4)
            print('dumping:', checkpoint_file)
            input()

        self.cleanup(self.target_directory)

if __name__ == '__main__':
    d = {'b': 456, 'a': 123, 'c': set([1,2,4])}
    print(build_id(d))
