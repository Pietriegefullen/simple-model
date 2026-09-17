import hashlib
import json

import wn

from USER_VARIABLES import RESULTS_DIRECTORY as CP_ROOT

ALPHABET = 'abcdefghijklmnopqrstuvwxyz'.upper()

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

class PrintCallback():
    def __init__(self):
        self.print_r2 = False
        
    def __call__(self, objective):
        args, kwargs, loss_value = objective.last_call()
        _, _, best_loss = objective.best_call()
        cnt = objective.call_count()
        if self.print_r2:
            print(objective)
            print(f'call {cnt:6d}: loss value {loss_value:8.3g}, best loss: {best_loss:8.3g}')
        else:
            print(f'\rcall {cnt:6d}: loss value {loss_value:8.3g}, best loss: {best_loss:8.3g}', end = '')

class CheckpointCallback():
    def __init__(self, keep_only_n = None):
        self.keep_only_n = keep_only_n
        
        
    def checkpoint_id(self, objective):
        model_id = build_id(objective.model().get_config(), 3, 2)
        
        # plain:
        # sample_number-validation_replica
        # fit_mode
        
        # loss_value w/o decimal point -> careful formatting!
    
    def cleanup(self, save_dir):
        if not self.keep_only_n is None:
            raise NotImplementedError()
            # sort by loss value
            # keep only 
        
        
    def __call__(self, objective):
        cp_transformed_parameters, _, last_loss = objective.last_call()
        _,_, best_loss = objective.best_call()
        
        if not last_loss == best_loss:
            return
        
        model_id = 'model-' + build_id(objective.model().get_config(), 3,2)
                    
        # optimisation ID: also including parameters (which are const, which variable!, bounds)
        
        # run ID: initial parameters!
        
        # run ID is unique for
        # sample/replica
        # fit_mode
        
        validation_replica = #TODO: get from where?
        
        # hashed: fit_period (t_start, t_end)
                
        
        save_dir = os.path.join(CP_PATH, run_id)
        if not os.path.isdir(save_dir):
            os.makedirs(save_dir)
        
        file_name = '_'.join([sample_name, ])
        file_path = os.path.join(save_dir, file_name)
        
        # get configuration
        # store parameters
        #   with which metadata?
        
        variables = objective.model().parameters().variables()
        checkpoint_parameters = [var.inverse_transform(p) 
                            for var, p in zip(variables, cp_transformed_parameters)]
        
        checkpoint_data = {
                            'parameters': checkpoint_parameters,
                            'total_loss': last_loss,
                            'model': model_config,
                           }
        
        #with open(checkpoint_file, 'w') as cf:
        #    json.dump(checkpoint_data, cf, indent = 4)

        #self.cleanup(save_dir)

if __name__ == '__main__':
    d = {'b': 456, 'a': 123, 'c': set([1,2,4])}
    print(build_id(d))