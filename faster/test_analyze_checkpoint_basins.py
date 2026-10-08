import unittest

import analyze_checkpoint_basins as basins


class CheckpointBasinTest(unittest.TestCase):
    def test_scaled_coordinates_use_log_parameter_range(self):
        candidates = [
            {'parameters': {'Hydrolysis_v_max': 1e-8}},
            {'parameters': {'Hydrolysis_v_max': 1e-1}},
        ]
        names, coordinates = basins.scaled_coordinates(candidates)
        index = names.index('Hydrolysis_v_max')
        self.assertAlmostEqual(coordinates[0, index], 0.)
        self.assertAlmostEqual(coordinates[1, index], 1.)

    def test_loss_envelope_retains_all_good_candidates(self):
        candidates = [{'loss': 1.}, {'loss': 1.05}, {'loss': 1.2}]
        accepted, threshold = basins.loss_envelope(candidates, .05, 0.)
        self.assertAlmostEqual(threshold, 1.05)
        self.assertEqual(accepted, candidates[:2])

    def test_connected_components_identifies_separate_parameter_basin(self):
        labels = basins.connected_components([[0.], [.04], [.3]], .05)
        self.assertEqual(labels, [1, 1, 2])


if __name__ == '__main__':
    unittest.main()
