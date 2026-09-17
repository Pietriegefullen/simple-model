import hashlib
import json

from USER_VARIABLES import RESULTS_DIRECTORY as CP_ROOT

# RUN hash
# model hash
# 

ALPHABET = 'abcdefghijklmnopqrstuvwxyz'.upper()

def freeze(obj):
    if isinstance(obj, dict):
        return tuple(sorted((k, freeze(v)) for k, v in obj.items()))
    if isinstance(obj, (list, tuple)):
        return tuple(freeze(v) for v in obj)
    if isinstance(obj, set):
        return tuple(freeze(v) for v in obj)
    return obj

def build_id(config, group_length = 4, length = 12):
    assert isinstance(config, dict)
    print(freeze(config))
    data = json.dumps(freeze(config), sort_keys=True, separators=(",", ":"))
    #s_hash = hashlib.sha256(data.encode()).hexdigest()
   
    digest = hash(tuple(sorted(freeze(config))))
    
    number = abs(digest)

    result = []
    base = len(ALPHABET)

    while number:
        number, remainder = divmod(number, base)
        result.append(ALPHABET[remainder])

    s_hash = "".join(reversed(result)).zfill(length)[-length:]
    
    identifier = "-".join(
                        s_hash[i:i+3] for i in range(0, len(s_hash), group_length)
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
    def __call__(self, objective):
        
        # model ID: model version -> manually?
        #   low priority, expecting little change.
            
        # pathway ID: from pathway configurations
        
        # optimisation ID: also including parameters (which are const, which variable!, bounds)
        
        # run ID: initial parameters!
        
        # run ID is unique for
        # sample/replica
        # fit_mode
        
        # fit_period (t_start, t_end)
        
        # don't use hash for simple (bool) configuration such as fit_mode
        
        run_id = build_id()
        
        save_dir = os.path.join(CP_PATH, run_id)
        
        if not os.path.isdir(save_dir):
            os.makedirs(save_dir)
        
        # determine save path
        # get configuration
        # store parameters
        
        pass

if __name__ == '__main__':
    d = {'b': 456, 'a': 123, 'c': set([1,2,4])}
    print(build_id(d))