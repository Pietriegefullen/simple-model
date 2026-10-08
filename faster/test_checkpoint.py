from contextlib import redirect_stdout
import io
import json
import os
import tempfile
import unittest

import checkpoint
from test_hashing import checkpoint_with


class _Parameters:
    def get_values(self):
        return {'Hydro_thermodynamics': True, 'rate': 1.0}


class _Model:
    def parameters(self):
        return _Parameters()


class _Objective:
    def __init__(self, loss):
        self.loss = loss

    def last_call(self):
        return None, None, self.loss

    def model(self):
        return _Model()


class _PrintObjective:
    def __init__(self, calls):
        self.calls = calls
        self.index = 0
        self.generation_number = 1

    def last_call(self):
        return None, None, self.calls[self.index][0]

    def last_R2(self):
        return self.calls[self.index][1]

    def best_loss(self):
        return min(loss for loss, _ in self.calls[:self.index + 1])

    def best_R2(self):
        best_index = min(range(self.index + 1), key=lambda i: self.calls[i][0])
        return self.calls[best_index][1]

    def call_count(self):
        return self.index + 1

    def generation(self):
        return self.generation_number


def _write_checkpoint(path, loss):
    with open(path, 'w') as handle:
        json.dump({'total_loss': loss}, handle)


class ExistingCheckpointTest(unittest.TestCase):
    def _callback(self, directory, overwrite_existing=False):
        run_config = checkpoint_with()['run_config']
        callback = checkpoint.CheckpointCallback(
            run_config, target=directory, overwrite_existing=overwrite_existing)
        callback.set_objective(_Objective(2.0))
        return callback

    def test_better_existing_checkpoint_is_kept_by_default(self):
        with tempfile.TemporaryDirectory() as directory:
            old_checkpoint = os.path.join(directory, 'old-checkpoint')
            _write_checkpoint(old_checkpoint, 1.0)

            self._callback(directory)()

            self.assertEqual(set(os.listdir(directory)),
                             {'old-checkpoint', '.checkpoint.lock'})

    def test_checkpointing_resumes_once_the_current_run_beats_existing(self):
        with tempfile.TemporaryDirectory() as directory:
            _write_checkpoint(os.path.join(directory, 'old-checkpoint'), 1.0)
            callback = self._callback(directory)

            callback()
            callback.objective.loss = 0.5
            callback()

            checkpoints = [name for name in os.listdir(directory)
                           if not name.startswith('.')]
            self.assertEqual(len(checkpoints), 2)

    def test_loads_lowest_loss_existing_checkpoint(self):
        with tempfile.TemporaryDirectory() as directory:
            _write_checkpoint(os.path.join(directory, 'worse-checkpoint'), 2.0)
            best_checkpoint = os.path.join(directory, 'best-checkpoint')
            _write_checkpoint(best_checkpoint, 1.0)

            loaded = self._callback(directory).load_best_existing_checkpoint()

            self.assertEqual(loaded[0], best_checkpoint)
            self.assertEqual(loaded[1]['total_loss'], 1.0)

    def test_overwrite_existing_checkpoints_replaces_them(self):
        with tempfile.TemporaryDirectory() as directory:
            old_checkpoint = os.path.join(directory, 'old-checkpoint')
            _write_checkpoint(old_checkpoint, 1.0)

            self._callback(directory, overwrite_existing=True)()

            checkpoints = [name for name in os.listdir(directory)
                           if not name.startswith('.')]
            self.assertEqual(len(checkpoints), 1)
            self.assertNotEqual(checkpoints[0], 'old-checkpoint')

    def test_print_callback_keeps_r2_from_the_best_loss_evaluation(self):
        objective = _PrintObjective(((3.0, 0.1), (4.0, 0.9), (2.0, 0.2)))
        callback = checkpoint.PrintCallback(None)
        callback.set_objective(objective)

        output = io.StringIO()
        with redirect_stdout(output):
            callback()
            objective.index = 1
            callback()
            objective.index = 2
            callback()

        self.assertIn('best loss:        2, best R2 = 0.20', output.getvalue())

    def test_print_callback_includes_generation_count(self):
        objective = _PrintObjective(((3.0, 0.1),))
        objective.generation_number = 4
        callback = checkpoint.PrintCallback(None)
        callback.set_objective(objective)

        output = io.StringIO()
        with redirect_stdout(output):
            callback()

        self.assertIn('generation 4: call', output.getvalue())

    def test_print_callback_reports_existing_best_r2_until_it_is_beaten(self):
        objective = _PrintObjective(((2.0, 0.9), (0.5, 0.2)))
        callback = checkpoint.PrintCallback(
            None, initial_best_loss=1.0, initial_best_r2=0.7)
        callback.set_objective(objective)

        output = io.StringIO()
        with redirect_stdout(output):
            callback()
            objective.index = 1
            callback()

        self.assertIn('best loss:        1, best R2 = 0.70', output.getvalue())
        self.assertIn('best loss:      0.5, best R2 = 0.20', output.getvalue())


if __name__ == '__main__':
    unittest.main()
