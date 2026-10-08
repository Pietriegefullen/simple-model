import threading
from contextlib import redirect_stdout
import io
import unittest
from unittest import mock

import numpy as np

import optimizer


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


if __name__ == '__main__':
    unittest.main()
