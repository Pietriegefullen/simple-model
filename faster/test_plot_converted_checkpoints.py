import unittest

import plot_converted_checkpoints as converted_plots


class ConvertedCheckpointPlotTest(unittest.TestCase):
    def test_safe_name_preserves_readable_checkpoint_identity(self):
        self.assertEqual(converted_plots.safe_name('cp A:test.json'), 'cp-A-test.json')

    def test_replica_metadata_is_required(self):
        document = {'parameters': {}, 'unit_conversion': {
            'target_units': {'M_*': 'micromol/L'}}}
        with self.assertRaises(ValueError):
            converted_plots.replica_for(document, {})


if __name__ == '__main__':
    unittest.main()
