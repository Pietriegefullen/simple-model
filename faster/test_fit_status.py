import sys
import unittest
from pathlib import Path
from unittest import mock

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

    def test_report_highlights_values_below_the_default_threshold(self):
        results = [
            fit_status.FitResult('1000', '4', 'single', 'A', .875, None),
            fit_status.FitResult('1000', '5', 'single', 'A', .799, None),
        ]
        report = fit_status.format_report(results, self.dataset, color=True)
        self.assertIn('  1000', report)
        self.assertIn(' 0.875', report)
        self.assertIn(f'{fit_status.HIGHLIGHT_STYLE} 0.799', report)
        self.assertNotIn('0.799!', report)
        self.assertIn('R² < 0.800', report)
        self.assertIn('coverage: 2/2 fitted', report)
        self.assertIn('sample      4      5', report)
        self.assertNotIn('sample      1', report)

    def test_report_accepts_a_different_highlight_threshold(self):
        result = fit_status.FitResult('1000', '4', 'single', 'A', .875, None)
        report = fit_status.format_report(
            [result], self.dataset, highlight_threshold=.9, color=True)
        self.assertIn('R² < 0.900', report)
        self.assertIn(f'{fit_status.HIGHLIGHT_STYLE} 0.875', report)

    def test_report_omits_samples_with_requested_prefixes(self):
        omitted_sample = _Sample()
        omitted_sample.sample_name = '2000'
        self.dataset.samples.append(omitted_sample)
        result = fit_status.FitResult('1000', '4', 'single', 'A', .875, None)

        report = fit_status.format_report(
            [result], self.dataset, omit_sample_prefixes=('2', '3'))

        self.assertIn('1000', report)
        self.assertNotIn('2000', report)
        self.assertIn('coverage: 1/2 fitted', report)

    def test_best_results_filters_on_a_partial_model_id(self):
        checkpoint = {
            'parameters': {},
            'total_loss': 0.1,
            'run_config': {
                'model': {'name': 'matching'},
                'chosen': {
                    'sample': '1000', 'validation_replica': '4',
                    'fit_mode': 'single'},
                'objective': self.objective,
            },
        }
        excluded = {
            **checkpoint,
            'run_config': {**checkpoint['run_config'], 'model': {'name': 'other'}},
        }
        with mock.patch.object(
                fit_status, 'checkpoint_files',
                return_value=[(Path('matching.json'), checkpoint),
                              (Path('other.json'), excluded)]), \
                mock.patch.object(fit_status.hashing, 'build_model_id', side_effect=[
                    'model-3NWE-JDSY', 'model-FA75-GSEM']), \
                mock.patch.object(fit_status, 'observed_variance', return_value=1.0):
            results, errors = fit_status.best_results(
                'unused', self.dataset, model_ids=['3NWE'])

        self.assertEqual(errors, [])
        self.assertEqual(len(results), 1)
        self.assertEqual(results[0].path, Path('matching.json'))

    def test_model_id_option_replaces_variant_option(self):
        with mock.patch.object(sys, 'argv', [
                'fit_status.py', '--model-id', '3NWE', '--model-id', 'JDSY']):
            args = fit_status.parse_args()

        self.assertEqual(args.model_ids, ['3NWE', 'JDSY'])
        self.assertFalse(hasattr(args, 'variants'))


if __name__ == '__main__':
    unittest.main()
