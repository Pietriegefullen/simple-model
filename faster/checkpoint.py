import os
import json
import fcntl
import tempfile
from contextlib import contextmanager

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
            plain.append(value)
        return '_'.join([str(p) for p in plain])

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
    def __init__(self, run_config, keep_only_n = None, verbose = False, legacy = False,
                 target = None):
        super().__init__(run_config)
        assert keep_only_n is None or keep_only_n > 0
        self.keep_only_n = keep_only_n
        self.verbose = verbose
        self.legacy = legacy
        if not target is None:
            self._target_directory = target
        
        # set range of legacy parameters to default range
        if self.legacy:
            default_parameters = parameters.default_model_parameters()
            parameter_names = self.run_config()['initial'].keys()
            for p_name in parameter_names:
                p = parameters.Parameter.from_config(self.run_config['initial'][p_name])
                if not p.is_variable(): continue
                self.run_config['initial'][p_name] = default_parameters[p_name].get_config()
        

    def checkpoints(self, save_dir):
        checkpoints = []
        for entry in os.scandir(save_dir):
            if not entry.is_file() or entry.name.startswith('.'):
                continue
            try:
                with open(entry.path) as checkpoint_file:
                    loss = json.load(checkpoint_file)['total_loss']
            except (json.JSONDecodeError, KeyError, OSError):
                continue
            checkpoints.append((loss, entry.path))
        return sorted(checkpoints)

    @contextmanager
    def locked_directory(self):
        save_dir = self.target_directory()
        os.makedirs(save_dir, exist_ok = True)
        lock_file = os.path.join(save_dir, '.checkpoint.lock')
        with open(lock_file, 'a') as lock:
            fcntl.flock(lock, fcntl.LOCK_EX)
            try:
                yield save_dir
            finally:
                fcntl.flock(lock, fcntl.LOCK_UN)

    def cleanup(self, save_dir):
        if self.keep_only_n is None:
            return
        for _, checkpoint_file in self.checkpoints(save_dir)[self.keep_only_n:]:
            os.remove(checkpoint_file)
        
    def __call__(self):
        cp_transformed_parameters, _, last_loss = self.objective.last_call()
        with self.locked_directory() as save_dir:
            checkpoints = self.checkpoints(save_dir)
            if (self.keep_only_n is not None
                    and len(checkpoints) >= self.keep_only_n
                    and last_loss >= checkpoints[-1][0]):
                return

            checkpoint_data = {
                                'parameters': self.objective.model().parameters().get_values(),
                                'total_loss': last_loss,
                                'run_config': self.run_config,
                               }
            cp_id = hashing.build_checkpoint_id(checkpoint_data)

            f_loss = f'{last_loss:.6f}'.replace('.', '')
            if len(f_loss) > 8:
                f_loss = '9'*8
            else:
                f_loss = f_loss.zfill(8)

            file_name = '_'.join([self.run_dir_name(), f'loss-{f_loss}',
                                  self.run_id, cp_id])
            checkpoint_file = os.path.join(save_dir, file_name)
            file_descriptor, temporary_file = tempfile.mkstemp(
                dir = save_dir, prefix = '.checkpoint-', text = True)
            try:
                with os.fdopen(file_descriptor, 'w') as checkpoint:
                    json.dump(checkpoint_data, checkpoint, indent = 4)
                os.replace(temporary_file, checkpoint_file)
            finally:
                if os.path.exists(temporary_file):
                    os.remove(temporary_file)

            self.cleanup(save_dir)

        if self.verbose:
            print(f'\nSaved checkpoint {checkpoint_file}')
