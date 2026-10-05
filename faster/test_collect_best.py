import unittest

import collect_best


class SampleYLimitsTest(unittest.TestCase):
    def test_limits_include_every_replica_and_log_ch4_excludes_zero(self):
        class Replica:
            def __init__(self, co2, ch4):
                self._co2 = co2
                self._ch4 = ch4

            def CO2(self):
                return [0, 1], self._co2

            def CH4(self):
                return [0, 1], self._ch4

        class Sample:
            replicas = [
                Replica([1.0, 4.0], [0.0, 1e-3]),
                Replica([2.0, 10.0], [1e-4, 1e-2]),
            ]

        limits = collect_best.sample_y_limits(Sample())
        self.assertLess(limits['CO2'][0], 1.0)
        self.assertGreater(limits['CO2'][1], 10.0)
        self.assertGreater(limits['CH4'][0], 0.0)
        self.assertLess(limits['CH4'][0], 1e-4)
        self.assertGreater(limits['CH4'][1], 1e-2)


class CheckpointFilterTest(unittest.TestCase):
    def test_filters_are_combined(self):
        class Checkpoint:
            def __init__(self, sample, replica, model_id, fit_mode, loss_id):
                self.replica = (sample, replica)
                self.model_id = model_id
                self.fit_mode = fit_mode
                self.loss_id = loss_id

        checkpoints = [
            Checkpoint('1351', '4', 'model-A', 'single', 'loss-1'),
            Checkpoint('1351', '5', 'model-A', 'split', 'loss-1'),
            Checkpoint('1352', '4', 'model-B', 'split', 'loss-2'),
        ]
        selected = collect_best.filter_checkpoints(
            checkpoints, samples=['1351'], model_ids=['model-A'],
            fit_modes=['split'], replicas=['5'], loss_ids=['loss-1'])
        self.assertEqual(selected, [checkpoints[1]])


if __name__ == '__main__':
    unittest.main()
