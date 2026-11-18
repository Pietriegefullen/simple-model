import unittest

import numpy as np

import plot_m0_vmax_surface as surface


class TestM0VmaxSurface(unittest.TestCase):
    def test_pathway_mapping(self):
        self.assertEqual(surface.pathway_parameters('hydro'),
                         ('Hydro', 'M_Hydro', 'Hydro_v_max'))

    def test_shared_biomass_pathway_is_explained(self):
        with self.assertRaisesRegex(ValueError, 'shares M_Ferm'):
            surface.pathway_parameters('Hydrolysis')

    def test_logarithmic_grid_contains_checkpoint(self):
        grid = surface.logarithmic_grid(1., (.01, 100.), 1., 5)
        np.testing.assert_allclose(grid, [.1, 10 ** -.5, 1., 10 ** .5, 10.])

    def test_logarithmic_grid_respects_bounds(self):
        grid = surface.logarithmic_grid(1., (.5, 2.), 2., 5)
        self.assertEqual(grid[0], .5)
        self.assertEqual(grid[-1], 2.)
