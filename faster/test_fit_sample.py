from contextlib import redirect_stderr, redirect_stdout
import io
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest import mock

import fit_sample
import hashing
import user_input
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


class OptimizerArgumentTest(unittest.TestCase):

    def test_fit_mode_must_be_explicit(self):
        stderr = io.StringIO()
        with redirect_stderr(stderr), self.assertRaises(SystemExit) as error:
            user_input.parse_args_fit(['1351', '4'])

        self.assertEqual(error.exception.code, 2)
        self.assertIn('Specify --single', stderr.getvalue())

    def test_fit_mode_flags_are_mutually_exclusive(self):
        with redirect_stderr(io.StringIO()), self.assertRaises(SystemExit) as error:
            user_input.parse_args_fit(['1351', '4', '--single', '--split'])

        self.assertEqual(error.exception.code, 2)

    def test_each_fit_mode_can_be_selected(self):
        single = user_input.parse_args_fit(['1351', '4', '--single'])
        split = user_input.parse_args_fit(['1351', '4', '--split'])

        self.assertEqual((single.single, single.split), (True, False))
        self.assertEqual((split.single, split.split), (False, True))

    def test_checkpoint_override_resolves_partial_id_and_distinct_values(self):
        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            converted = root / 'v0.4.3'
            converted.mkdir()
            source = converted / 'single_run-ABC-DEF_from-AAAA-BBBB_cp-ABCD-EFGH.json'
            document = {
                'run_config': {'model': {'version': '0.3'}},
                'parameters': {'Ac_CUE': .1, 'Ferm_CUE': .2, 'Ferm_Km': 123.},
                'total_loss': .01,
                'unit_conversion': {'target_units': {'M_*': 'micromol/L'}},
                'converted_checkpoint': {'format': 'litre-normalized-v1',
                                         'replay_model_version': '0.4.3'},
            }
            source.write_text(json.dumps(document))
            with mock.patch.object(user_input, 'CONVERTED_CHECKPOINTS_DIRECTORY', root):
                names = ['Ac_CUE', 'Ferm_CUE', 'Ferm_Km', 'Hydro_v_max']
                args = user_input.parse_args_fit([
                    '1351', '4', '--single', '--p', '*_CUE', '*_Km*', 'efgh',
                    '--p', 'Hydro_v_max', '1e-2'])
                fixed = {}
                for group in args.p:
                    fixed.update(user_input.parse_override(group, names, model_version='0.4.3'))
                self.assertEqual(fixed, {**document['parameters'], 'Hydro_v_max': .01})
                with self.assertRaisesRegex(ValueError, 'no saved value for: Hydro_v_max'):
                    user_input.parse_override(['Hydro_v_max', 'efgh'], names, model_version='0.4.3')
            del document['converted_checkpoint'], document['unit_conversion']
            source.write_text(json.dumps(document))
            with self.assertRaisesRegex(ValueError, 'supplies v0.3 values'):
                user_input.parse_override(['*_CUE', str(source)], names, model_version='0.4.3')

    def test_fix_v_max_except_argument_is_parsed(self):
        args = user_input.parse_args_fit([
            '1351', '4', '--single', '--fix-v-max-except', 'Hydrolysis_v_max',
            '--fix-v-max-value', '0.25'])

        self.assertEqual(args.fix_v_max_except, ['Hydrolysis_v_max'])
        self.assertEqual(args.fix_v_max_value, 0.25)

    def test_fix_v_max_value_requires_v_max_selection(self):
        with (redirect_stdout(io.StringIO()), redirect_stderr(io.StringIO()),
              self.assertRaises(SystemExit)):
            user_input.parse_args_fit(['1351', '4', '--single', '--fix-v-max-value', '0.25'])

    def test_default_and_estimate_are_mutually_exclusive(self):
        with (redirect_stdout(io.StringIO()), redirect_stderr(io.StringIO()),
              self.assertRaises(SystemExit)):
            user_input.parse_args_fit(['1351', '4', '--single', '--default', '--estimate'])

    def test_converted_envelope_is_a_fresh_parameter_source(self):
        args = user_input.parse_args_fit(['1351', '4', '--single', '--converted-envelope'])

        self.assertTrue(args.converted_envelope)
        with (redirect_stdout(io.StringIO()), redirect_stderr(io.StringIO()),
              self.assertRaises(SystemExit)):
            user_input.parse_args_fit(
                ['1351', '4', '--single', '--default', '--converted-envelope'])

    def test_hypothesis_is_a_fresh_parameter_source(self):
        args = user_input.parse_args_fit(['1351', '4', '--single', '--hypothesis'])

        self.assertTrue(args.hypothesis)
        with (redirect_stdout(io.StringIO()), redirect_stderr(io.StringIO()),
              self.assertRaises(SystemExit)):
            user_input.parse_args_fit(['1351', '4', '--single', '--default', '--hypothesis'])

    def test_bounds_is_a_parameter_source(self):
        args = user_input.parse_args_fit(['1351', '4', '--single', '--bounds', 'ABCD'])

        self.assertEqual(args.bounds, 'ABCD')
        with (redirect_stdout(io.StringIO()), redirect_stderr(io.StringIO()),
              self.assertRaises(SystemExit)):
            user_input.parse_args_fit(['1351', '4', '--single', '--default', '--bounds', 'ABCD'])

    def test_cma_es_options_are_available_for_single_sample_fits(self):
        args = user_input.parse_args_fit([
            '1351', '4', '--single', '--algorithm', 'cma-es', '--cma-sigma', '0.1',
            '--cma-restarts', '1', '--cma-restart-strategy', 'bipop', '--funtol', '0.002', '--popsize', '12',
            '--max-evaluations', '300'])

        self.assertEqual(args.algorithm, 'cma-es')
        self.assertEqual(args.cma_sigma, 0.1)
        self.assertEqual(args.cma_restarts, 1)
        self.assertEqual(args.cma_restart_strategy, 'bipop')
        self.assertEqual(args.funtol, 0.002)
        self.assertEqual(args.popsize, 12)
        self.assertEqual(args.max_evaluations, 300)

    def test_gas_limits_are_opt_in(self):
        default = user_input.parse_args_fit(['1351', '4', '--single'])
        enabled = user_input.parse_args_fit(['1351', '4', '--single', '--gas-limits'])

        self.assertFalse(default.gas_limits)
        self.assertTrue(enabled.gas_limits)


if __name__ == '__main__':
    unittest.main()
