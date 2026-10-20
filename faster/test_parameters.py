import unittest
import io
from contextlib import redirect_stdout

import numpy as np

import parameters


class NumPyBooleanParameterTest(unittest.TestCase):
    def test_numpy_boolean_scalars_are_accepted(self):
        value = np.bool_(False)
        model_parameters = parameters.ModelParameters({'flag': value})

        self.assertIs(model_parameters['flag'].value, value)

        parameter = parameters.Parameter('flag')
        parameter.value = value
        self.assertIs(parameter.value, value)

        model_parameters['flag'].set(value)
        self.assertEqual(model_parameters['flag'].value, 0.0)


class TransformedBoundaryTest(unittest.TestCase):
    def test_log_transform_endpoint_is_clamped_without_out_of_range_warning(self):
        parameter = parameters.Parameter('rate', 1.0, (0.2, 3.0), 'log', normalize=True)
        round_trip_upper = parameter.inverse_transform(parameter.transform(parameter.upper()))

        output = io.StringIO()
        with redirect_stdout(output):
            parameter.set(round_trip_upper)

        self.assertEqual(parameter.value, parameter.upper())
        self.assertEqual(output.getvalue(), '')

    def test_meaningfully_out_of_range_value_still_warns(self):
        parameter = parameters.Parameter('rate', 1.0, (0.2, 3.0), 'log', normalize=True)
        output = io.StringIO()
        with redirect_stdout(output):
            parameter.set(3.01)

        self.assertIn('out of bounds', output.getvalue())


if __name__ == '__main__':
    unittest.main()
