import os
import json

import hashing
import parameters

from USER_VARIABLES import RESULTS_DIRECTORY as CP_ROOT

checkpoint_plain = [
                    ('objective', hashing.build_loss_id),
                    ('model', hashing.build_model_id), 
                    ('chosen', 'sample'), 
                    ('chosen', 'validation_replica'), 
                    ('chosen', 'fit_mode')
                    # run ID is always appended
                    # checkpoint ID
                    ]

fill_length = {'fit_mode': 6}

def fill_value(value, key):
    if not key in fill_length:
        return value
    length = fill_length[key]
    if len(value) > length:
        return value[:length]
    while len(value) < length:
        value += '_'
    return value
        
    
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
        self._target_directory = None
        self.run_config = run_config
        if not run_config is None and 'legacy' in run_config:
            self.run_id = hashing.build_run_id({'legacy': run_config['legacy']})
        else:
            self.run_id = None if run_config is None else hashing.build_run_id(run_config)

    def set_objective(self, objective):
        self.objective = objective
        
    def target_directory(self):
        if self._target_directory is None and not self.run_id is None:
            self._target_directory = os.path.join(CP_ROOT, self.run_dir_name())
        return self._target_directory

    def run_dir_name(self):
        plain = []
        
        for keys in checkpoint_plain:
            value = get_value(self.run_config, keys)
            value = fill_value(str(value), keys[-1])
            plain.append(value)
        plain.append(self.run_id)
        return '_'.join([p for p in plain])

class SetAllConstant(Callback):
    def __init__(self):
        super().__init__(None)
        
    def __call__(self):
        for var in self.objective.model().parameters().variables():
            var.constant(var.value)

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
    def __init__(self, run_config, keep_only_n = None, verbose = False, legacy = False):
        super().__init__(run_config)
        assert keep_only_n > 0
        self.keep_only_n = keep_only_n
        self.verbose = False
        self.legacy = legacy
        
        # set range of legacy parameters to default range
        if self.legacy:
            default_parameters = parameters.default_model_parameters()
            parameter_names = self.run_config()['initial'].keys()
            for p_name in parameter_names:
                p = parameters.Parameter.from_config(self.run_config['initial'][p_name])
                if not p.is_variable(): continue
                self.run_config['initial'][p_name] = default_parameters[p_name].get_config()
        

    def cleanup(self, save_dir):
        if self.keep_only_n is None:
            return
        files = [os.path.join(save_dir, f) for f in os.listdir(save_dir)]
        if len(files) > self.keep_only_n:            
            all_files = sorted([(parameters.load_parameter_file(f)[1], f)
                                for f in files])

            for _, f in all_files[self.keep_only_n:]:
                os.remove(f)
        
    def __call__(self):
        cp_transformed_parameters, _, last_loss = self.objective.last_call()
        _,_, best_loss = self.objective.best_call()
        
        if not last_loss == best_loss:
            return

        if not os.path.isdir(self.target_directory()):
            os.makedirs(self.target_directory())

        checkpoint_parameters = self.objective.model().parameters().get_values()
        
        checkpoint_data = {
                            'parameters': checkpoint_parameters,
                            'total_loss': last_loss,
                            'run_config': self.run_config,
                           }
       
        cp_id = hashing.build_checkpoint_id(checkpoint_data)
        
        f_loss = f'{last_loss:.3f}'.replace('.', '')
        if len(f_loss) > 8:
            # overflow
            f_loss = '9'*8
        else:
            while len(f_loss) < 8:
                f_loss = '0' + f_loss
        
        file_name = self.run_dir_name() + f'_loss-{f_loss}_' + cp_id
        checkpoint_file = os.path.join(self.target_directory(), file_name)
        with open(checkpoint_file, 'w') as cf:
            json.dump(checkpoint_data, cf, indent = 4)
            
        #if self.verbose:
        print()
        print('Saved checkpoint ', checkpoint_file)
        self.cleanup(self.target_directory())

