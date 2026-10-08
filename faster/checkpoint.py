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
    def __init__(self, run_config, initial_best_loss = None,
                 initial_best_r2 = None):
        super().__init__(run_config)
        # A fit can resume in a directory that already contains results.  The
        # stored loss is an incumbent too, even though it is not part of this
        # Objective's call log.
        self.initial_best_loss = initial_best_loss
        self.initial_best_r2 = initial_best_r2
        
    def __call__(self):
        args, kwargs, loss_value = self.objective.last_call()
        cnt = self.objective.call_count()
        generation = self.objective.generation()
        generation_text = 'n/a' if generation is None else str(generation)
        best_loss = self.objective.best_loss()
        best_r2 = self.objective.best_R2()
        if (self.initial_best_loss is not None
                and self.initial_best_loss <= best_loss):
            best_loss = self.initial_best_loss
            # Checkpoints do not retain R², so do not pair a historical loss
            # with R² from a different evaluation.  ``fit`` recreates the
            # checkpoint's run log before supplying this value.
            best_r2 = self.initial_best_r2
        r2_text = 'n/a' if best_r2 is None else f'{best_r2:.2f}'
        print(f'\rgeneration {generation_text}: call {cnt:6d}: '
              f'loss value {loss_value:8.3g}, '
              f'best loss: {best_loss:8.3g}, best R2 = {r2_text}', end = '')

class CheckpointCallback(Callback):
    def __init__(self, run_config, keep_only_n = None, verbose = False, legacy = False,
                 target = None, overwrite_existing = False):
        super().__init__(run_config)
        assert keep_only_n is None or keep_only_n > 0
        self.keep_only_n = keep_only_n
        self.verbose = verbose
        self.legacy = legacy
        self.overwrite_existing = overwrite_existing
        self._existing_checkpoint_paths = None
        self._existing_best_loss = None
        self._current_run_best_loss = None
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

    def prepare_existing_checkpoints(self, checkpoints):
        """Snapshot or explicitly discard checkpoints that predate this run."""
        if self._existing_checkpoint_paths is not None:
            return checkpoints

        self._existing_checkpoint_paths = {path for _, path in checkpoints}
        if checkpoints:
            self._existing_best_loss = checkpoints[0][0]

        if not self.overwrite_existing:
            return checkpoints

        for _, checkpoint_file in checkpoints:
            os.remove(checkpoint_file)
        return []

    def load_best_existing_checkpoint(self):
        """Return the best valid checkpoint already in this run directory.

        Initialising this before minimisation makes the existing result the
        incumbent for both checkpoint retention and live progress output.
        ``--overwrite-checkpoints`` deliberately keeps its documented
        behaviour of starting with an empty checkpoint directory.
        """
        with self.locked_directory() as save_dir:
            checkpoints = self.prepare_existing_checkpoints(
                self.checkpoints(save_dir))
            if not checkpoints:
                return None

            _, checkpoint_file = checkpoints[0]
            try:
                with open(checkpoint_file) as handle:
                    return checkpoint_file, json.load(handle)
            except (json.JSONDecodeError, OSError):
                # A file can become unreadable after the initial scan only if
                # another process has modified it outside this lock.  Treat it
                # as unavailable rather than preventing the fit from running.
                return None

    def existing_checkpoint_is_better(self):
        return (not self.overwrite_existing
                and self._existing_best_loss is not None
                and self._current_run_best_loss is not None
                and self._existing_best_loss < self._current_run_best_loss)
        
    def __call__(self):
        cp_transformed_parameters, _, last_loss = self.objective.last_call()
        if (self._current_run_best_loss is None
                or last_loss < self._current_run_best_loss):
            self._current_run_best_loss = last_loss

        with self.locked_directory() as save_dir:
            checkpoints = self.prepare_existing_checkpoints(self.checkpoints(save_dir))
            if self.existing_checkpoint_is_better():
                return

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
