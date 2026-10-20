import unittest

import numpy as np

import model
import optimizer
import system


class IntegrationGuardTest(unittest.TestCase):
    def setUp(self):
        self.initial_state = np.zeros(len(system.SYSTEM))
        self.times = np.array([0., 1.])
        self.limits = {'CO2': 1., 'CH4': 1e6}

    def test_co2_runaway_event_aborts_integration(self):
        def runaway_rhs(_time, _state):
            derivative = np.zeros(len(system.SYSTEM))
            derivative[system.index('CO2')] = 10.
            return derivative

        with self.assertRaisesRegex(optimizer.IntegrationFailure, 'CO2 exceeded'):
            model.integrate(runaway_rhs, self.times, self.initial_state, [], None,
                            self.limits)

if __name__ == '__main__':
    unittest.main()
