import threading
import tempfile
from contextlib import redirect_stdout
import io
import unittest
from unittest import mock

import numpy as np

import optimizer
import checkpoint


class ValueObjective(optimizer.Addable):
    def _call(self, value, **kwargs):
        return value


class MetricValueObjective(optimizer.Addable):
    def _call(self, value, **kwargs):
        loss, residual, total = value
        self._last_mse = residual, total
        return loss

    def last_weighted_mse(self):
        return self._last_mse


class CaptureBestLoss:
    def set_objective(self, objective):
        self.objective = objective
        self.values = []

    def __call__(self):
        self.values.append(self.objective.best_loss())


class SharedValue:
    def __init__(self, value):
        self.value = value


class SamplingVariable:
    def __init__(self, lower, upper):
        self._lower = lower
        self._upper = upper
        self.value = (lower + upper) / 2

    def lower(self):
        return self._lower

    def upper(self):
        return self._upper

    def transform(self, value):
        return value

    def inverse_transform(self, value):
        return value


class SamplingParameters:
    def __init__(self):
        self._variables = [SamplingVariable(-1, 1), SamplingVariable(0, 2)]

    def variables(self):
        return self._variables


class SamplingModel:
    def __init__(self):
        self._parameters = SamplingParameters()

    def parameters(self):
        return self._parameters


class SamplingObjective:
    def __init__(self):
        self._model = SamplingModel()
        self.evaluated = []
        self.selected = None

    def model(self):
        return self._model

    def __call__(self, candidate):
        candidate = np.asarray(candidate)
        self.evaluated.append(candidate.copy())
        return np.sum((candidate - np.array([0.2, 1.3])) ** 2)

    def set_parameters(self, variables, values, transformed=False):
        self.selected = np.asarray(values).copy()
        self.selected_was_transformed = transformed


class InterruptingSamplingObjective(SamplingObjective):
    def __init__(self, interrupt_after):
        super().__init__()
        self.interrupt_after = interrupt_after

    def __call__(self, candidate):
        if len(self.evaluated) >= self.interrupt_after:
            raise KeyboardInterrupt()
        return super().__call__(candidate)


class _NoVariables:
    def variables(self):
        return []


class _NoParameterModel:
    def parameters(self):
        return _NoVariables()

    def predict(self, replica, t_eval):
        return {}


class _InvalidLoss:
    def t_eval(self, replica):
        return np.array([1.0])

    def __call__(self, replica, run_log):
        raise optimizer.NonFiniteLoss('mse(CH4) returned inf')


class _FailedIntegrationModel(_NoParameterModel):
    def predict(self, replica, t_eval):
        raise optimizer.IntegrationFailure('Integration exceeded 5 seconds')


class GlobalBestLossTest(unittest.TestCase):
    def test_callbacks_share_the_best_loss_across_workers(self):
        counter = SharedValue(0)
        lock = threading.Lock()
        best_loss = SharedValue(float('inf'))

        first_worker = ValueObjective()
        second_worker = ValueObjective()
        for worker in (first_worker, second_worker):
            worker.set_global_call_counter(counter, lock, best_loss)

        first_callback = CaptureBestLoss()
        second_callback = CaptureBestLoss()
        first_worker.add_callback(first_callback)
        second_worker.add_callback(second_callback)

        first_worker(5.0)
        second_worker(3.0)
        first_worker(4.0)

        self.assertEqual(first_callback.values, [5.0, 3.0])
        self.assertEqual(second_callback.values, [3.0])
        self.assertEqual(first_worker.best_loss(), 3.0)
        self.assertEqual(second_worker.best_loss(), 3.0)

    def test_global_best_r2_stays_paired_with_global_best_loss(self):
        counter = SharedValue(0)
        lock = threading.Lock()
        best_loss = SharedValue(float('inf'))
        best_r2 = SharedValue(float('nan'))

        first_worker = MetricValueObjective()
        second_worker = MetricValueObjective()
        for worker in (first_worker, second_worker):
            worker.set_global_call_counter(counter, lock, best_loss, best_r2)

        first_worker((5.0, 5.0, 10.0))
        second_worker((3.0, 1.0, 10.0))
        first_worker((4.0, 0.0, 10.0))

        self.assertEqual(first_worker.best_loss(), 3.0)
        self.assertAlmostEqual(first_worker.best_R2(), 0.9)
        self.assertAlmostEqual(second_worker.best_R2(), 0.9)


class InitialSamplingTest(unittest.TestCase):
    def test_evaluates_exact_population_without_differential_evolution(self):
        objective = SamplingObjective()
        algorithm = optimizer.InitialSampling()
        algorithm.configure({'init': 'sobol', 'samples': 7, 'seed': 11})

        with mock.patch.object(
                optimizer.scipy.optimize, 'differential_evolution',
                side_effect=AssertionError('differential evolution must not run')):
            algorithm._minimize(objective, None)

        self.assertEqual(len(objective.evaluated), 7)
        losses = [np.sum((candidate - np.array([0.2, 1.3])) ** 2)
                  for candidate in objective.evaluated]
        self.assertTrue(np.array_equal(
            objective.selected, objective.evaluated[np.argmin(losses)]))
        self.assertTrue(objective.selected_was_transformed)
        self.assertEqual(algorithm.best_fun, min(losses))


class DifferentialEvolutionBudgetTest(unittest.TestCase):
    def test_ignores_and_warns_for_unsupported_configuration_keys(self):
        algorithm = optimizer.DifferentialEvolution()

        with self.assertWarnsRegex(RuntimeWarning, 'ignored unsupported.*name'):
            algorithm.configure({'name': 'explore', 'tol': 0.01})

        self.assertNotIn('name', algorithm.get_config())
        self.assertEqual(algorithm.get_config()['tol'], 0.01)

    def test_generation_callback_stops_at_the_requested_evaluation_budget(self):
        algorithm = optimizer.DifferentialEvolution()
        algorithm._max_evaluations = 5
        algorithm._objective = type(
            'Objective', (), {'call_count': lambda self: 5,
                              'set_generation': lambda self, value: None})()
        algorithm._shared_generation = None
        algorithm.generation = 1

        self.assertTrue(algorithm.generation_counter())


class CMAEvolutionStrategyTest(unittest.TestCase):
    def test_respects_budget_and_leaves_objective_at_best_candidate(self):
        objective = SamplingObjective()
        algorithm = optimizer.CMAEvolutionStrategy()
        algorithm.configure({'seed': 5, 'popsize': 8, 'sigma0': 0.35,
                             'max_evaluations': 40})

        algorithm._minimize(objective, None)

        self.assertEqual(len(objective.evaluated), 40)
        losses = [np.sum((candidate - np.array([0.2, 1.3])) ** 2)
                  for candidate in objective.evaluated]
        self.assertTrue(np.array_equal(
            objective.selected, objective.evaluated[np.argmin(losses)]))
        self.assertTrue(objective.selected_was_transformed)

    def test_allows_a_budget_smaller_than_one_population(self):
        objective = SamplingObjective()
        algorithm = optimizer.CMAEvolutionStrategy()
        algorithm.configure({'seed': 2, 'popsize': 8, 'max_evaluations': 3})

        algorithm._minimize(objective, None)

        self.assertEqual(len(objective.evaluated), 3)

    def test_continues_from_a_generation_state_checkpoint(self):
        """A resumed CMA run has the same next populations as an intact run."""
        with tempfile.TemporaryDirectory() as directory:
            state = checkpoint.CMAStateCheckpoint(
                {'model': {}, 'chosen': {}, 'objective': {}, 'algo': {}, 'range': {}},
                target=directory)
            first = SamplingObjective()
            interrupted = optimizer.CMAEvolutionStrategy()
            interrupted.configure({'seed': 11, 'popsize': 8, 'sigma0': 0.3,
                                   'max_evaluations': 16})
            interrupted.set_state_checkpoint(state)
            interrupted._minimize(first, None)

            continued = SamplingObjective()
            resumed = optimizer.CMAEvolutionStrategy()
            resumed.configure({'seed': 11, 'popsize': 8, 'sigma0': 0.3,
                               'max_evaluations': 32})
            resumed.set_state_checkpoint(state)
            resumed._minimize(continued, None)

            whole = SamplingObjective()
            uninterrupted = optimizer.CMAEvolutionStrategy()
            uninterrupted.configure({'seed': 11, 'popsize': 8, 'sigma0': 0.3,
                                     'max_evaluations': 32})
            uninterrupted._minimize(whole, None)

        combined = np.vstack(first.evaluated + continued.evaluated)
        self.assertEqual(combined.shape, np.vstack(whole.evaluated).shape)
        self.assertTrue(np.array_equal(combined, np.vstack(whole.evaluated)))

    def test_interrupt_within_a_population_preserves_completed_candidates(self):
        with tempfile.TemporaryDirectory() as directory:
            state = checkpoint.CMAStateCheckpoint(
                {'model': {}, 'chosen': {}, 'objective': {}, 'algo': {}, 'range': {}},
                target=directory)
            interrupted_objective = InterruptingSamplingObjective(2)
            interrupted = optimizer.CMAEvolutionStrategy()
            interrupted.configure({'seed': 17, 'popsize': 8, 'max_evaluations': 32})
            interrupted.set_state_checkpoint(state)
            with self.assertRaises(KeyboardInterrupt):
                interrupted._minimize(interrupted_objective, None)

            resumed_objective = SamplingObjective()
            resumed = optimizer.CMAEvolutionStrategy()
            resumed.configure({'seed': 17, 'popsize': 8, 'max_evaluations': 32})
            resumed.set_state_checkpoint(state)
            resumed._minimize(resumed_objective, None)

            whole = SamplingObjective()
            uninterrupted = optimizer.CMAEvolutionStrategy()
            uninterrupted.configure({'seed': 17, 'popsize': 8, 'max_evaluations': 32})
            uninterrupted._minimize(whole, None)

        combined = np.vstack(interrupted_objective.evaluated + resumed_objective.evaluated)
        self.assertTrue(np.array_equal(combined, np.vstack(whole.evaluated)))


class CombinedObjectiveTest(unittest.TestCase):
    def test_combined_objective_delegates_parameter_updates_to_a_leaf(self):
        receiver = mock.Mock()
        combined = optimizer.Addable(receiver, None)

        combined.set_parameters(['parameter'], [0.5], transformed=True)

        receiver.set_parameters.assert_called_once_with(
            ['parameter'], [0.5], True)


class NonFiniteLossTest(unittest.TestCase):
    def test_empty_or_infinite_mse_is_classified_as_non_finite_loss(self):
        loss = optimizer.Loss('CH4', optimizer.mse)
        loss.get_values = lambda replica, run_log: (np.array([]), np.array([]))
        with self.assertRaisesRegex(optimizer.NonFiniteLoss, 'no usable'):
            loss(None, None)

        loss.get_values = lambda replica, run_log: (
            np.array([1e308]), np.array([-1e308]))
        with self.assertRaisesRegex(optimizer.NonFiniteLoss, 'returned inf'):
            loss(None, None)

    def test_invalid_single_replica_receives_a_finite_penalty(self):
        objective = optimizer.Objective(_NoParameterModel(), object())
        objective.add_loss(_InvalidLoss(), 1.0)

        output = io.StringIO()
        with redirect_stdout(output):
            value = objective._call([], transformed=False)

        self.assertEqual(value, optimizer.INVALID_OBJECTIVE_PENALTY)
        self.assertIsNone(objective.last_R2())
        self.assertIn('WARNING: Replica', output.getvalue())

    def test_failed_integration_receives_a_finite_penalty(self):
        objective = optimizer.Objective(_FailedIntegrationModel(), object())
        objective.add_loss(_InvalidLoss(), 1.0)

        value = objective._call([], transformed=False)

        self.assertEqual(value, optimizer.INVALID_OBJECTIVE_PENALTY)
        self.assertIsNone(objective.last_R2())


if __name__ == '__main__':
    unittest.main()
