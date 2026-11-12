import unittest

import numpy as np

import parameters
import system
from model_v0_3 import model


class Replica:
    dry_weight = 2.0
    water_content = 100.0

    class Sample:
        sample_name = '1351'

    sample = Sample()

    def initial_TOC(self):
        return 10.0

    def initial_DOC(self):
        return 20.0

    def initial_H2O(self):
        return 30.0

    def CO2(self):
        return np.array([0.0, 2.0]), np.array([40.0, 60.0])

    def CH4(self):
        return np.array([0.0, 2.0]), np.array([50.0, 70.0])


class LegacyInitialStateTest(unittest.TestCase):
    def test_v03_retains_its_pre_litre_initialisation(self):
        values = parameters.ModelParameters({pool: 1.0 for pool in system.SYSTEM})

        state = model.initial_state(Replica(), values)

        self.assertEqual(state[system.index('TOC')], 10.0)
        self.assertEqual(state[system.index('DOC')], 20.0)
        self.assertEqual(state[system.index('H2O')], 30.0)
        self.assertEqual(state[system.index('CO2')], 40.0)
        self.assertEqual(state[system.index('CH4')], 50.0)

    def test_rebuild_preserves_loaded_checkpoint_values(self):
        pathway_model = model.build_model(['Hydrolysis', 'Hydro'], True)
        pathway_model.parameters()['Hydro_v_max'].constant(0.314159)

        pathway_model.build()

        self.assertEqual(pathway_model.parameters()['Hydro_v_max'].value,
                         0.314159)


if __name__ == '__main__':
    unittest.main()
