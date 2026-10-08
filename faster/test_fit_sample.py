from contextlib import redirect_stdout
import io
import unittest

import fit_sample
import hashing
import parameters
from test_hashing import checkpoint_with


class RequestedCheckpointTest(unittest.TestCase):
    def setUp(self):
        self.checkpoint = checkpoint_with()
        self.chosen = self.checkpoint['run_config']['chosen'].copy()
        self.model = self.checkpoint['run_config']['model']

    def test_compatible_checkpoint_needs_no_confirmation(self):
        warnings = fit_sample.checkpoint_compatibility_warnings(
            self.checkpoint, self.chosen, self.model)

        self.assertEqual(warnings, [])

    def test_incompatible_metadata_is_reported_for_confirmation(self):
        chosen = self.chosen.copy()
        chosen['sample'] = 9999

        warnings = fit_sample.checkpoint_compatibility_warnings(
            self.checkpoint, chosen, self.model)

        self.assertTrue(any(warning.startswith('sample:') for warning in warnings))
        with redirect_stdout(io.StringIO()):
            with self.assertRaisesRegex(RuntimeError, 'not confirmed'):
                fit_sample.confirm_checkpoint_compatibility(
                    'checkpoint.json', warnings, input_function=lambda _: 'no')
            fit_sample.confirm_checkpoint_compatibility(
                'checkpoint.json', warnings, input_function=lambda _: 'yes')

    def test_checkpoint_values_must_be_inside_current_range(self):
        fit_range = parameters.ModelParameters({
            'rate': parameters.Parameter('rate', 0.5, [0.0, 1.0], 'linear'),
        })
        checkpoint = {'parameters': {'rate': 1.5}}

        errors = fit_sample.checkpoint_range_errors(checkpoint, fit_range)

        self.assertEqual(errors, ['rate=1.5 is outside [0, 1]'])


if __name__ == '__main__':
    unittest.main()
