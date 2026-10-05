import threading
import unittest

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


if __name__ == '__main__':
    unittest.main()
