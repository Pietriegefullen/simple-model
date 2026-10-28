"""Parity tests for the opt-in Cython ODE right-hand-side backend."""

import unittest

import numpy as np

import cython_model
import data
import fit_replicas
import model
import optimizer


PATHWAYS = ['Hydrolysis', 'Fermentation', 'Hydro', 'Aceto', 'Homo', 'Fe3']


class CythonModelParityTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # Build the extension before starting integration workers: it makes a
        # missing compiler/dependency fail at an obvious point in this test.
        cython_model._rhs_module()
        cls.replica = data.get_data_before_day()[1367].replicas[0]

    def make_models(self):
        python_model = model.build_model(PATHWAYS, normalized=True)
        compiled_model = cython_model.build_model(PATHWAYS, normalized=True)
        return python_model, compiled_model

    def test_rhs_matches_python_for_representative_states(self):
        python_model, compiled_model = self.make_models()
        python_model.replica = self.replica
        compiled_model.replica = self.replica
        python_model.build()
        compiled_model.build()

        states = [
            np.full(14, 1.0),
            np.array([2.1, 0.4, 0.2, 0.3, 1.0, 0.1, 0.2, 0.3,
                      0.4, 0.5, 0.6, 0.7, 0.8, 0.9]),
            np.array([-1.0, 0.0, 0.1, 0.2, 2.0, 0.0, 0.3, 0.4,
                      0.5, 0.6, 0.7, 0.8, 0.9, 1.0]),
        ]
        for state in states:
            np.testing.assert_allclose(
                compiled_model(0.0, state), python_model(0.0, state),
                rtol=1e-12, atol=1e-12)

    def test_trajectory_matches_python(self):
        python_model, compiled_model = self.make_models()
        times = np.array([0.0, 1.0, 3.0, 7.0])
        python_run = python_model.predict(self.replica, t_eval=times)
        compiled_run = compiled_model.predict(self.replica, t_eval=times)

        for pool in ('CO2', 'CH4'):
            np.testing.assert_allclose(compiled_run[pool][0], python_run[pool][0])
            np.testing.assert_allclose(compiled_run[pool][1], python_run[pool][1],
                                       rtol=1e-10, atol=1e-11)

    def test_objective_loss_matches_python(self):
        python_model, compiled_model = self.make_models()
        config = fit_replicas.objective_config()
        python_objective = optimizer.build_objective_function(
            python_model, [self.replica], config)
        compiled_objective = optimizer.build_objective_function(
            compiled_model, [self.replica], config)

        python_values = np.asarray([
            variable.transform(variable.value)
            for variable in python_objective.variables()
        ])
        compiled_values = np.asarray([
            variable.transform(variable.value)
            for variable in compiled_objective.variables()
        ])
        self.assertAlmostEqual(
            compiled_objective(compiled_values), python_objective(python_values),
            places=10)

    def test_multi_fit_parser_selects_cython_backend(self):
        args = fit_replicas.parse_args(['--shared', 'Hydro_v_max', '--c'])
        self.assertTrue(args.use_cython)
        self.assertIsInstance(
            fit_replicas.build_model(PATHWAYS, normalized_parameters=True,
                                     use_cython=args.use_cython),
            cython_model.CythonModel)


if __name__ == '__main__':
    unittest.main()
