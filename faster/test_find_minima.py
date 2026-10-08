import tempfile
import unittest
from pathlib import Path
from unittest import mock

import find_minima


class ClusteringTest(unittest.TestCase):
    def test_rms_distance_is_dimension_normalized(self):
        self.assertAlmostEqual(find_minima.rms_distance([0, 0], [0.1, 0.1]), .1)

    def test_connected_components_uses_transitive_single_linkage(self):
        labels = find_minima.connected_components([[0.0], [.04], [.08], [.30]], .05)
        self.assertEqual(labels, [1, 1, 1, 2])

    def test_loss_cutoff_never_discards_saved_candidates(self):
        candidates = [
            {'status': 'complete', 'total_loss': 1.0},
            {'status': 'complete', 'total_loss': 1.05},
            {'status': 'complete', 'total_loss': 1.2},
            {'status': 'failed'},
        ]
        accepted, threshold = find_minima.accepted_candidates(candidates, .05, 0.)
        self.assertAlmostEqual(threshold, 1.05)
        self.assertEqual(len(accepted), 2)
        self.assertEqual(len(candidates), 4)

    def test_clustering_labels_are_written_without_comparing_array_dicts(self):
        candidates = [
            {'status': 'complete', 'total_loss': 1., 'parameter_coordinates': [0.],
             'prediction_signature': [1.], 'search': {'local_success': True}},
            {'status': 'complete', 'total_loss': 2., 'parameter_coordinates': [.8],
             'prediction_signature': [2.], 'search': {'local_success': True}},
        ]
        accepted, _ = find_minima.cluster_candidates(candidates, .1, 0., .1, .1)
        self.assertEqual(len(accepted), 1)
        self.assertTrue(candidates[0]['accepted'])
        self.assertFalse(candidates[1]['accepted'])

    def test_main_preserves_a_checkpoint_for_every_restart(self):
        def candidate(args, restart, seed):
            return {
                'schema': 'test', 'status': 'complete', 'restart': restart, 'seed': seed,
                'parameters': {'rate': .1 + restart}, 'total_loss': 1. + .01 * restart,
                'run_config': {},
                'search': {'local_success': True, 'global_loss': 1.,
                           'global_evaluations': 2, 'local_evaluations': 3},
                'fit_r2': .8, 'validation_loss': 1., 'validation_r2': .7,
                'parameter_coordinates': [.1 + .01 * restart],
                'parameter_names': ['rate'], 'prediction_signature': [1.],
                'trajectories': {
                    'CO2': {'time': [1], 'measured': [1], 'predicted': [1]},
                    'CH4': {'time': [1], 'measured': [1], 'predicted': [1]},
                },
            }

        with tempfile.TemporaryDirectory() as directory:
            with mock.patch.object(find_minima, 'fitted_candidate', side_effect=candidate):
                result = find_minima.main([
                    '1366', '6', '--restarts', '2', '--output-root', directory,
                    '--run-name', 'test-run', '--loss-relative', '.1',
                ])
            self.assertEqual(len(list((Path(result) / 'candidates').glob('restart-*.json'))), 2)
            self.assertTrue((Path(result) / 'basin_summary.csv').is_file())


if __name__ == '__main__':
    unittest.main()
