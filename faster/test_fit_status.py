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
        self.assertIn(f'{fit_status.HIGHLIGHT_STYLE}  0.799', report)
        self.assertNotIn('0.799!', report)
        self.assertIn('R² < 0.800', report)
        self.assertIn('coverage: 2/2 fitted', report)
        self.assertIn('sample      4      5', report)
        self.assertNotIn('sample      1', report)

    def test_r2_cells_reserve_space_for_a_negative_sign(self):
        positive = fit_status.format_r2(.1234, .8)
        negative = fit_status.format_r2(-.1234, .8)

        self.assertEqual(positive, '  0.123')
        self.assertEqual(negative, ' -0.123')
        self.assertEqual(len(positive), len(negative))
        self.assertEqual(positive.index('.'), negative.index('.'))

    def test_report_accepts_a_different_highlight_threshold(self):
        result = fit_status.FitResult('1000', '4', 'single', 'A', .875, None)
        report = fit_status.format_report(
            [result], self.dataset, highlight_threshold=.9, color=True)
        self.assertIn('R² < 0.900', report)
        self.assertIn(f'{fit_status.HIGHLIGHT_STYLE}  0.875', report)

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

    def test_multi_results_use_every_replica_score_from_best_global_checkpoint(self):
        def multi_checkpoint(loss, r2_values):
            return {
                'parameters': {},
                'total_loss': loss,
                'run_config': {
                    'model': {'name': 'joint'},
                    'chosen': {'sample': 0, 'validation_replica': 0,
                               'fit_mode': 'multi'},
                    'objective': self.objective,
                    'multi_replica': {'replicas': ['10004', '10005']},
                },
                'replicas': [
                    {'sample': '1000', 'replica_number': replica,
                     'replica': f'1000{replica}', 'r2': r2}
                    for replica, r2 in r2_values.items()
                ],
            }

        best_global = multi_checkpoint(.1, {'4': .71, '5': .72})
        better_local_scores = multi_checkpoint(.2, {'4': .99, '5': .98})
        with mock.patch.object(
                fit_status, 'checkpoint_files',
                return_value=[(Path('best.json'), best_global),
                              (Path('other.json'), better_local_scores)]), \
                mock.patch.object(fit_status.hashing, 'build_model_id',
                                  return_value='model-JOINT'), \
                mock.patch.object(fit_status.hashing, 'build_run_id',
                                  return_value='run-JOINT'):
            results, errors = fit_status.best_results('unused', self.dataset)

        self.assertEqual(errors, [])
        self.assertEqual(
            {(result.replica, result.total_r2, result.path) for result in results},
            {('4', .71, Path('best.json')), ('5', .72, Path('best.json'))})

    def test_polished_results_follow_only_the_best_global_checkpoint(self):
        def outer(loss, identifier):
            return {
                'parameters': {}, 'total_loss': loss,
                'run_config': {
                    'model': {'name': 'joint'},
                    'chosen': {'sample': 0, 'validation_replica': 0,
                               'fit_mode': 'multi'},
                    'objective': self.objective,
                    'multi_replica': {'replicas': ['10004']},
                },
                'replicas': [{'sample': '1000', 'replica_number': '4'}],
            }

        def polished(loss):
            return {
                'parameters': {}, 'total_loss': loss,
                'run_config': {
                    'model': {'name': 'joint'},
                    'chosen': {'sample': '1000', 'validation_replica': '4',
                               'fit_mode': 'single'},
                    'objective': self.objective,
                },
            }

        root = Path('unused')
        best_root = root / 'joint_polish' / 'cp-BEST'
        other_root = root / 'joint_polish' / 'cp-OTHER'

        def files(directory):
            directory = Path(directory)
            if directory == root:
                return [(Path('outer_cp-BEST'), outer(.1, 'BEST')),
                        (Path('outer_cp-OTHER'), outer(.2, 'OTHER')),
                        (Path('single-fit'), polished(.4))]
            if directory == best_root:
                return [(Path('best-polish-1'), polished(.3)),
                        (Path('best-polish-2'), polished(.2))]
            if directory == other_root:
                return [(Path('other-polish'), polished(.01))]
            return []

        with mock.patch.object(fit_status, 'checkpoint_files', side_effect=files), \
                mock.patch.object(Path, 'is_dir', return_value=True), \
                mock.patch.object(fit_status.hashing, 'build_model_id',
                                  side_effect=['model-POLISHED', 'model-POLISHED',
                                               'model-SINGLE']), \
                mock.patch.object(fit_status.hashing, 'build_run_id',
                                  return_value='run-JOINT'), \
                mock.patch.object(fit_status, 'observed_variance', return_value=1.0):
            results, errors = fit_status.polished_results(
                root, self.dataset, single_model_ids=['SINGLE'])

        self.assertEqual(errors, [])
        self.assertEqual(len(results), 2)
        by_mode = {result.fit_mode: result for result in results}
        self.assertAlmostEqual(by_mode['single'].total_r2, .6)
        self.assertEqual(by_mode['single'].path, Path('single-fit'))
        self.assertAlmostEqual(by_mode['polished'].total_r2, .8)
        self.assertEqual(by_mode['polished'].path, Path('best-polish-2'))

    def test_polished_comparison_places_single_and_polished_values_side_by_side(self):
        report = fit_status.format_polished_comparison([
            fit_status.FitResult('1000', '4', 'single', 'single', .61, None),
            fit_status.FitResult('1000', '4', 'polished', 'polished', .82, None),
        ], self.dataset, highlight_threshold=.8)

        self.assertIn('single polished', report)
        self.assertIn('0.610   0.820', report)

    def test_polished_comparison_blanks_replicas_outside_the_multi_fit(self):
        report = fit_status.format_polished_comparison([
            fit_status.FitResult('1000', '4', 'single', 'single', .61, None),
            fit_status.FitResult('1000', '4', 'polished', 'polished', .82, None),
        ], self.dataset, highlight_threshold=.8, polished_scope={('1000', '4')})

        row = next(line for line in report.splitlines() if line.startswith('  1000'))
        self.assertIn('0.610', row)
        self.assertNotIn('--', row)

    def test_model_id_option_replaces_variant_option(self):
        with mock.patch.object(sys, 'argv', [
                'fit_status.py', '--model-id', '3NWE', '--model-id', 'JDSY']):
            args = fit_status.parse_args()

        self.assertEqual(args.model_ids, ['3NWE', 'JDSY'])
        self.assertFalse(hasattr(args, 'variants'))


if __name__ == '__main__':
    unittest.main()
