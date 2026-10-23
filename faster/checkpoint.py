import os
import json
import fcntl
import tempfile
import base64
import pickle
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
        self.run_config = (None if run_config is None
                           else hashing.normalize_run_config(run_config))

        self.run_id = None if self.run_config is None else hashing.build_run_id(self.run_config)

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
                 initial_best_r2 = None, header = None, footer = None, carriage_return = True,
                 show_penalty_count = False):
        super().__init__(run_config)
        # A fit can resume in a directory that already contains results.  The
        # stored loss is an incumbent too, even though it is not part of this
        # Objective's call log.
        self.initial_best_loss = initial_best_loss
        self.initial_best_r2 = initial_best_r2
        self.header = header
        self.footer = footer
        self.carriage_return = carriage_return
        self.show_penalty_count = show_penalty_count
        
    def __call__(self):
        args, kwargs, loss_value = self.objective.last_call()
        cnt = getattr(self.objective, 'evaluation_count',
                      self.objective.call_count)()
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
        header = '' if self.header is None else self.header + ' '
        cret = '\r' if self.carriage_return else ''
        footer = '' if self.footer is None else (' ' + self.footer)
        penalty_text = ''
        if self.show_penalty_count:
            penalty_count = getattr(self.objective, 'penalty_count', lambda: 0)()
            penalty_text = f', penalties {penalty_count}/{cnt}'
        print(f'{cret}{header}generation {generation_text}: call {cnt:6d}: '
              f'loss value {loss_value:8.3g}, '
              f'best loss: {best_loss:8.3g}, best R2 = {r2_text}'
              + penalty_text + footer, end = '')


class CallLimitInterrupt(Exception):
    pass


class StalledEvaluationInterrupt(CallLimitInterrupt):
    """No initial-sampling worker completed an evaluation in time."""

    def __init__(self, timeout):
        self.timeout = timeout
        super().__init__(f'No objective evaluation completed within {timeout:g} seconds')


class CMAStateCheckpoint:
    """Persist the latest CMA-ES progress state for later resumption.

    Ordinary checkpoints are deliberately ranked by objective value.  That is
    useful for selecting a model, but is not an optimiser checkpoint: CMA-ES
    also needs its covariance matrix, evolution paths and random-number state.
    This companion file is instead replaced as CMA-ES advances.  A state that
    was interrupted within a population includes its pending samples and
    completed losses, so it remains internally coherent.
    It lives beside the ordinary checkpoints, starts with a dot so older
    checkpoint scanners ignore it, and is atomically replaced under the same
    directory lock.

    The opaque payload is pickle data produced by the maintained ``cma``
    package.  It is local run state, not an interchange format; the run ID and
    a format version prevent accidentally loading it for a different fit.
    """

    FILE_NAME = '.cma-es-state.json'
    VERSION = 1

    def __init__(self, run_config, target=None):
        self.run_config = hashing.normalize_run_config(run_config)
        self.run_id = hashing.build_run_id(self.run_config)
        self._target_directory = target

    def target_directory(self):
        if self._target_directory is None:
            # Match Callback's default result location without requiring a
            # Callback instance (which also owns an Objective).
            plain = [get_value(self.run_config, keys) for keys in checkpoint_plain]
            self._target_directory = os.path.join(CP_ROOT, '_'.join(map(str, plain)))
        return self._target_directory

    @contextmanager
    def _locked_directory(self):
        save_dir = self.target_directory()
        os.makedirs(save_dir, exist_ok=True)
        lock_file = os.path.join(save_dir, '.checkpoint.lock')
        with open(lock_file, 'a') as lock:
            fcntl.flock(lock, fcntl.LOCK_EX)
            try:
                yield save_dir
            finally:
                fcntl.flock(lock, fcntl.LOCK_UN)

    def save(self, state):
        """Atomically store an opaque, JSON-compatible CMA state payload."""
        encoded_state = base64.b64encode(
            pickle.dumps(state, protocol=pickle.HIGHEST_PROTOCOL)).decode('ascii')
        document = {
            'version': self.VERSION,
            'run_id': self.run_id,
            'state': encoded_state,
        }
        with self._locked_directory() as save_dir:
            destination = os.path.join(save_dir, self.FILE_NAME)
            descriptor, temporary_file = tempfile.mkstemp(
                dir=save_dir, prefix='.cma-es-state-', text=True)
            try:
                with os.fdopen(descriptor, 'w') as handle:
                    json.dump(document, handle, sort_keys=True)
                os.replace(temporary_file, destination)
            finally:
                if os.path.exists(temporary_file):
                    os.remove(temporary_file)

    def load(self):
        """Return a compatible saved state, or ``None`` when none is usable."""
        with self._locked_directory() as save_dir:
            source = os.path.join(save_dir, self.FILE_NAME)
            try:
                with open(source) as handle:
                    document = json.load(handle)
                if (document.get('version') != self.VERSION
                        or document.get('run_id') != self.run_id):
                    return None
                return pickle.loads(base64.b64decode(document['state']))
            except (OSError, ValueError, KeyError, TypeError,
                    pickle.PickleError, EOFError):
                # A damaged state file must never make a valuable ordinary
                # result checkpoint unusable.  The caller can safely restart
                # from its normal x0 seed instead.
                return None

    def clear(self):
        """Discard this run's saved CMA state before a deliberate fresh run."""
        with self._locked_directory() as save_dir:
            try:
                os.remove(os.path.join(save_dir, self.FILE_NAME))
            except FileNotFoundError:
                pass


class PenaltyCounter:
    """A small process-safe counter for invalid objective evaluations.

    Unlike manager proxies, this remains reliable when an inner optimiser
    repeatedly sends independently pickled objectives through one reusable
    process pool.  It intentionally shares the checkpoint directory lock so
    that increments cannot race with checkpoint writes or cleanup.
    """

    FILE_NAME = '.penalty-count.json'

    def __init__(self, target):
        self.target = target

    @contextmanager
    def _locked_directory(self):
        os.makedirs(self.target, exist_ok=True)
        lock_file = os.path.join(self.target, '.checkpoint.lock')
        with open(lock_file, 'a') as lock:
            fcntl.flock(lock, fcntl.LOCK_EX)
            try:
                yield
            finally:
                fcntl.flock(lock, fcntl.LOCK_UN)

    def _read_unlocked(self):
        try:
            with open(os.path.join(self.target, self.FILE_NAME)) as handle:
                document = json.load(handle)
                return {
                    'count': int(document['count']),
                    'calls': int(document.get('calls', 0)),
                }
        except (OSError, ValueError, KeyError, TypeError):
            return {'count': 0, 'calls': 0}

    def reset(self):
        with self._locked_directory():
            destination = os.path.join(self.target, self.FILE_NAME)
            descriptor, temporary_file = tempfile.mkstemp(
                dir=self.target, prefix='.penalty-count-', text=True)
            try:
                with os.fdopen(descriptor, 'w') as handle:
                    json.dump({'count': 0, 'calls': 0}, handle)
                os.replace(temporary_file, destination)
            finally:
                if os.path.exists(temporary_file):
                    os.remove(temporary_file)

    def record(self, penalized):
        with self._locked_directory():
            current = self._read_unlocked()
            destination = os.path.join(self.target, self.FILE_NAME)
            descriptor, temporary_file = tempfile.mkstemp(
                dir=self.target, prefix='.penalty-count-', text=True)
            try:
                with os.fdopen(descriptor, 'w') as handle:
                    json.dump({
                        'count': current['count'] + int(penalized),
                        'calls': current['calls'] + 1,
                    }, handle)
                os.replace(temporary_file, destination)
            finally:
                if os.path.exists(temporary_file):
                    os.remove(temporary_file)

    def value(self):
        with self._locked_directory():
            return self._read_unlocked()['count']

    def call_count(self):
        with self._locked_directory():
            return self._read_unlocked()['calls']

    def increment(self):
        """Compatibility shorthand for recording one penalised call."""
        self.record(True)


class Interrupt(Callback):
    def __init__(self, call_limit):
        super().__init__(None)
        self.call_limit = call_limit

    def __call__(self):
        if self.objective.call_count() >= self.call_limit:
            raise CallLimitInterrupt('Call limit')


class CheckpointCallback(Callback):
    def __init__(self, run_config, keep_only_n = None, verbose = False, legacy = False,
                 target = None, overwrite_existing = False,
                 ignore_existing = False):
        super().__init__(run_config)
        assert keep_only_n is None or keep_only_n > 0
        self.keep_only_n = keep_only_n
        self.verbose = verbose
        self.legacy = legacy
        self.overwrite_existing = overwrite_existing
        self.ignore_existing = ignore_existing
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

    def cleanup(self, save_dir, checkpoints=None):
        if self.keep_only_n is None:
            return
        if checkpoints is None:
            checkpoints = self.checkpoints(save_dir)
        for _, checkpoint_file in checkpoints[self.keep_only_n:]:
            os.remove(checkpoint_file)

    def prepare_existing_checkpoints(self, checkpoints):
        """Snapshot or explicitly discard checkpoints that predate this run."""
        if self._existing_checkpoint_paths is not None:
            return [checkpoint for checkpoint in checkpoints
                    if checkpoint[1] not in self._existing_checkpoint_paths]

        self._existing_checkpoint_paths = {path for _, path in checkpoints}
        if self.ignore_existing:
            return []

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

        Once a normal resumable run has started, a retained pre-existing
        checkpoint is still a valid incumbent.  In particular, a new failed
        or worse evaluation must not make that incumbent disappear from a
        post-optimisation lookup.
        """
        with self.locked_directory() as save_dir:
            all_checkpoints = self.checkpoints(save_dir)
            if self._existing_checkpoint_paths is None:
                checkpoints = self.prepare_existing_checkpoints(all_checkpoints)
            elif self.ignore_existing or self.overwrite_existing:
                checkpoints = [item for item in all_checkpoints
                               if item[1] not in self._existing_checkpoint_paths]
            else:
                checkpoints = all_checkpoints
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

    def checkpoint_availability(self):
        """Summarise why a checkpoint lookup may not have found a result.

        This is diagnostic-only: it neither repairs nor removes files.  It is
        used when a nested fit must keep going without a durable inner result,
        so the warning distinguishes "no evaluation" from filesystem and
        format problems.
        """
        with self.locked_directory() as save_dir:
            counts = {
                'readable': 0,
                'hidden': 0,
                'non_file': 0,
                'unreadable': 0,
                'invalid_json': 0,
                'missing_total_loss': 0,
            }
            for entry in os.scandir(save_dir):
                if entry.name.startswith('.'):
                    counts['hidden'] += 1
                    continue
                if not entry.is_file():
                    counts['non_file'] += 1
                    continue
                try:
                    with open(entry.path) as checkpoint_file:
                        document = json.load(checkpoint_file)
                except OSError:
                    counts['unreadable'] += 1
                    continue
                except json.JSONDecodeError:
                    counts['invalid_json'] += 1
                    continue
                if not isinstance(document, dict) or 'total_loss' not in document:
                    counts['missing_total_loss'] += 1
                    continue
                counts['readable'] += 1
            return ', '.join(f'{name}={count}' for name, count in counts.items()
                             if count) or 'directory is empty'

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
            # Some compound objectives have useful per-evaluation state that
            # is not represented by a single model's parameter collection.
            # Let those objectives persist that state without making ordinary
            # single-model checkpoints any larger or changing their format.
            metadata = getattr(self.objective, 'checkpoint_metadata', None)
            if callable(metadata):
                checkpoint_data.update(metadata())
            checkpoint_data = hashing.normalize_checkpoint_config(checkpoint_data)
            cp_id = hashing.build_checkpoint_id(checkpoint_data)
            last_loss = checkpoint_data['total_loss']

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
                    json.dump(checkpoint_data, checkpoint, indent = 4,
                              sort_keys = True)
                os.replace(temporary_file, checkpoint_file)
            finally:
                if os.path.exists(temporary_file):
                    os.remove(temporary_file)

            self.cleanup(save_dir, self.prepare_existing_checkpoints(
                self.checkpoints(save_dir)))

        if self.verbose:
            print(f'\nSaved checkpoint {checkpoint_file}')
