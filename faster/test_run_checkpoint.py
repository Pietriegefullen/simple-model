import unittest
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest import mock

import run_checkpoint
from test_hashing import checkpoint_with


class ReplayConfigurationTest(unittest.TestCase):
    def test_saved_configuration_is_copied_before_overrides(self):
        checkpoint = checkpoint_with()
        original_loader = run_checkpoint._requested_checkpoint
        original_resolver = run_checkpoint.resolve_checkpoint
        run_checkpoint._requested_checkpoint = lambda _: checkpoint
        run_checkpoint.resolve_checkpoint = lambda _: Path('checkpoint.json')
        try:
            config, loaded_parameters = run_checkpoint.replay_configuration(
                'checkpoint.json', {'rate': 0.75}, t_start=5, t_end=40)
        finally:
            run_checkpoint._requested_checkpoint = original_loader
            run_checkpoint.resolve_checkpoint = original_resolver

        self.assertEqual(config['chosen']['parameter_override'],
                         {'Hydro_thermodynamics': True, 'rate': 0.75})
        self.assertEqual(config['chosen']['t_start'], 5)
        self.assertEqual(config['chosen']['t_end'], 40)
        self.assertEqual(loaded_parameters['rate'].value, 1.0)
        self.assertNotIn('rate', checkpoint['run_config']['chosen']['parameter_override'])
        self.assertIsNone(checkpoint['run_config']['chosen']['t_start'])

    def test_cli_parses_parameter_and_time_overrides(self):
        args = run_checkpoint.parse_args([
            'checkpoint.json', '--output-dir', 'report',
            '--p', 'rate', '0.75', '--t', '0', '40'])

        self.assertEqual(args.parameter_overrides, {'rate': 0.75})
        self.assertIsNone(args.t_start)
        self.assertEqual(args.t_end, 40)

    def test_resolves_unique_partial_checkpoint_id(self):
        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            checkpoint = root / 'fit_loss-1_cp-ABCD-EFGH'
            checkpoint.touch()

            resolved = run_checkpoint.resolve_checkpoint('abcd', root)

        self.assertEqual(resolved, checkpoint)

    def test_rejects_an_ambiguous_partial_checkpoint_id(self):
        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / 'fit_loss-1_cp-ABCD-EFGH').touch()
            (root / 'fit_loss-2_cp-ABCD-IJKL').touch()

            with self.assertRaisesRegex(ValueError, 'ambiguous'):
                run_checkpoint.resolve_checkpoint('abcd', root)

    def test_output_directory_must_be_new_or_empty(self):
        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            output = root / 'report'

            self.assertEqual(run_checkpoint.prepare_output_directory(output), output)
            (output / 'existing.svg').touch()
            with self.assertRaisesRegex(ValueError, 'not empty'):
                run_checkpoint.prepare_output_directory(output)

    def test_relative_output_directory_is_rooted_at_the_project(self):
        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            with mock.patch.object(run_checkpoint, 'PROJECT_DIRECTORY', root):
                output = run_checkpoint.prepare_output_directory('reports/run-1')

            self.assertEqual(output, root / 'reports/run-1')

    def test_metadata_values_are_valid_json_and_keep_infinities_replayable(self):
        serialized = json.dumps(run_checkpoint._json_safe({
            'positive': float('inf'),
            'negative': float('-inf'),
            'undefined': float('nan'),
        }), allow_nan=False)

        self.assertEqual(json.loads(serialized), {
            'positive': 'Infinity', 'negative': '-Infinity', 'undefined': 'NaN'})


if __name__ == '__main__':
    unittest.main()
