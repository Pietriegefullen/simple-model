import unittest

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


if __name__ == '__main__':
    unittest.main()
