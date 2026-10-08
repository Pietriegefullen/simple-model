import unittest

import fit_status


class _Replica:
    def __init__(self, number, co2, ch4):
        self.replica_number = str(number)
        self._co2 = co2
        self._ch4 = ch4

    def CO2(self):
        return self._co2

    def CH4(self):
        return self._ch4


class _Sample:
    sample_name = '1000'

    def __init__(self):
        self.replicas = [
            _Replica('4', ([0., 1., 2.], [2., 4., 6.]),
                     ([0., 1., 2.], [0., 1., 4.])),
            _Replica('5', ([0., 1., 2.], [3., 6., 9.]),
                     ([0., 1., 2.], [0., 2., 8.])),
        ]
        for replica in self.replicas:
            replica.sample = self

    def get_split(self, validation_replica, fit_mode):
        validation = next(replica for replica in self.replicas
                          if replica.replica_number == str(validation_replica))
        return {'fit': [validation] if fit_mode == 'single' else
                [replica for replica in self.replicas if replica is not validation]}

    def has_replicas(self):
        return len(self.replicas)


class _Dataset:
    def __init__(self):
        self.samples = [_Sample()]

    def __getitem__(self, sample):
        return next(item for item in self.samples if item.sample_name == str(sample))


class FitStatusTest(unittest.TestCase):
    def setUp(self):
        self.dataset = _Dataset()
        self.objective = {
            'loss_weight': {'CO2': 1., 'CH4': 1.},
            'transform': {'CO2': ['normalize'], 'CH4': ['log', 'normalize']},
        }

    def test_observed_variance_uses_fitted_replicas_and_objective_transforms(self):
        variance = fit_status.observed_variance(
            self.dataset,
            {'sample': '1000', 'validation_replica': '4', 'fit_mode': 'single'},
            self.objective)
        self.assertAlmostEqual(variance, 0.3125)

    def test_report_marks_missing_replicas_and_counts_coverage(self):
        result = fit_status.FitResult('1000', '4', 'single', 'A', .875, None)
        report = fit_status.format_report([result], self.dataset)
        self.assertIn('1000        ', report)
        self.assertIn('0.875', report)
        self.assertIn('coverage: 1/2 fitted', report)


if __name__ == '__main__':
    unittest.main()
