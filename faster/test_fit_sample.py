from contextlib import redirect_stderr, redirect_stdout
import io
import unittest

import fit_sample
import hashing
import parameters
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

    def test_checkpoint_values_must_be_inside_current_range(self):
        fit_range = parameters.ModelParameters({
            'rate': parameters.Parameter('rate', 0.5, [0.0, 1.0], 'linear'),
        })
        checkpoint = {'parameters': {'rate': 1.5}}

        errors = fit_sample.checkpoint_range_errors(checkpoint, fit_range)

        self.assertEqual(errors, ['rate=1.5 is outside [0, 1]'])


class OptimizerArgumentTest(unittest.TestCase):
    def test_fix_v_max_except_keeps_requested_parameter_variable(self):
        model_parameters = parameters.ModelParameters({
            parameter.name: parameter for parameter in parameters.default_model_parameters()})

        fit_sample.fix_v_max_except(model_parameters, ['Hydrolysis_v_max'])

        self.assertTrue(model_parameters['Hydrolysis_v_max'].is_variable())
        self.assertFalse(model_parameters['Ferm_v_max'].is_variable())
        self.assertTrue(model_parameters['Ferm_Km'].is_variable())

    def test_fix_v_max_except_sets_the_fixed_parameter_values(self):
        model_parameters = parameters.ModelParameters({
            parameter.name: parameter for parameter in parameters.default_model_parameters()})

        fit_sample.fix_v_max_except(
            model_parameters, ['Hydrolysis_v_max'], value=0.25)

        self.assertEqual(model_parameters['Ferm_v_max'].value, 0.25)
        self.assertFalse(model_parameters['Ferm_v_max'].is_variable())
        self.assertNotEqual(model_parameters['Hydrolysis_v_max'].value, 0.25)
        self.assertTrue(model_parameters['Hydrolysis_v_max'].is_variable())

    def test_fix_v_max_except_rejects_non_v_max_parameter(self):
        model_parameters = parameters.ModelParameters({
            parameter.name: parameter for parameter in parameters.default_model_parameters()})

        with self.assertRaisesRegex(ValueError, 'only accepts'):
            fit_sample.fix_v_max_except(model_parameters, ['Ferm_Km'])

    def test_fix_v_max_except_argument_is_parsed(self):
        args = user_input.parse_args_fit([
            '1351', '4', '--fix-v-max-except', 'Hydrolysis_v_max',
            '--fix-v-max-value', '0.25'])

        self.assertEqual(args.fix_v_max_except, ['Hydrolysis_v_max'])
        self.assertEqual(args.fix_v_max_value, 0.25)

    def test_fix_v_max_value_requires_v_max_selection(self):
        with (redirect_stdout(io.StringIO()), redirect_stderr(io.StringIO()),
              self.assertRaises(SystemExit)):
            user_input.parse_args_fit(['1351', '4', '--fix-v-max-value', '0.25'])

    def test_estimate_uses_estimated_parameter_catalogue(self):
        args = user_input.parse_args_fit(['1351', '4', '--estimate'])

        self.assertTrue(args.estimate)
        initial = parameters.load_parameters({'estimate': args.estimate})
        estimated = {
            parameter.name: parameter for parameter in parameters.estimated_model_parameters()}

        self.assertEqual(initial['Ferm_Km'].value, estimated['Ferm_Km'].value)
        self.assertEqual(initial['Ferm_Km'].high, estimated['Ferm_Km'].high)

    def test_default_and_estimate_are_mutually_exclusive(self):
        with (redirect_stdout(io.StringIO()), redirect_stderr(io.StringIO()),
              self.assertRaises(SystemExit)):
            user_input.parse_args_fit(['1351', '4', '--default', '--estimate'])

    def test_cma_es_options_are_available_for_single_sample_fits(self):
        args = user_input.parse_args_fit([
            '1351', '4', '--algorithm', 'cma-es', '--cma-sigma', '0.1',
            '--cma-restarts', '1', '--popsize', '12',
            '--max-evaluations', '300'])

        self.assertEqual(args.algorithm, 'cma-es')
        self.assertEqual(args.cma_sigma, 0.1)
        self.assertEqual(args.cma_restarts, 1)
        self.assertEqual(args.popsize, 12)
        self.assertEqual(args.max_evaluations, 300)


if __name__ == '__main__':
    unittest.main()
