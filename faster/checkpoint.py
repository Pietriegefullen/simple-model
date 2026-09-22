import os
import json

import hashing

from USER_VARIABLES import RESULTS_DIRECTORY as CP_ROOT

checkpoint_plain = ['sample', 'validation_replica', 'fit_mode']
# configure callback:
#   keep only N best, if N is None, keep all
#   store location

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

        checkpoint_parameters = self.objective.model().parameters().get_config()
        
        checkpoint_data = {
                            'parameters': checkpoint_parameters,
                            'total_loss': last_loss,
                           }
       
        file_name = self.run_dir() + f'_loss-{round(last_loss*1000):06d}'
        checkpoint_file = os.path.join(self.target_directory, file_name)
        with open(checkpoint_file, 'w') as cf:
            json.dump(checkpoint_data, cf, indent = 4)

        self.cleanup(self.target_directory)

