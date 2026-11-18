import unittest

import numpy as np

import fpca


class FunctionalPCATest(unittest.TestCase):
    def test_log_activity_clips_zero_rates_at_the_configured_floor(self):
        logged = fpca.log_activity([[0., 1e-9, 1e-3]], floor=1e-6)

        np.testing.assert_allclose(logged, [[np.log(1e-6), np.log(1e-6), np.log(1e-3)]])

    def test_log_activity_rejects_negative_rates(self):
        with self.assertRaisesRegex(ValueError, 'negative'):
            fpca.log_activity([[-1e-8]], floor=1e-6)

    def test_shared_component_ylim_contains_every_component_with_margin(self):
        limits = fpca.shared_component_ylim(np.array([
            [[-3., 2.]],
            [[-1., 5.]],
        ]))

        self.assertLess(limits[0], -3.)
        self.assertGreater(limits[1], 5.)

    def test_component_factor_argument_defaults_to_ten(self):
        self.assertEqual(fpca.parse_args([]).component_factor, 10.)

    def test_components_are_orthonormal_under_time_integral(self):
        time = np.linspace(0., 1., 11)
        signatures = np.zeros((3, len(fpca.PATHWAYS), len(time)))
        signatures[1, 0] = time
        signatures[2, 1] = time ** 2

        mean, basis, scores, explained = fpca.functional_pca(signatures, time, 2)

        self.assertEqual(mean.shape, (len(fpca.PATHWAYS), len(time)))
        self.assertEqual(scores.shape, (3, 2))
        self.assertAlmostEqual(explained.sum(), 1.)
        self.assertEqual(basis.shape, (2, len(fpca.PATHWAYS), len(time)))
        weights = fpca.trapezoid_weights(time)
        gram = np.einsum('ipt,jpt,t->ij', basis, basis, weights)
        np.testing.assert_allclose(gram, np.eye(2), atol=1e-10)

    def test_rejects_more_components_than_checkpoints(self):
        signatures = np.zeros((2, len(fpca.PATHWAYS), 3))
        with self.assertRaisesRegex(ValueError, 'components'):
            fpca.functional_pca(signatures, [0., 1., 2.], 3)

    def test_confirmation_accepts_yes(self):
        fpca.confirm(2, input_function=lambda _: 'yes')


if __name__ == '__main__':
    unittest.main()
