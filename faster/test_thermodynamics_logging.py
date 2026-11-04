import unittest

import numpy as np

import system
from parameters import ModelParameters
from model_run import ModelRun

from model_factory import build_model, pathway_by_name


class Replica:
    sample = type('Sample', (), {'sample_name': '1351'})()
    dry_weight = 1.
    total_gas_volume = 100.

    def __getitem__(self, key):
        if key == 'days':
            return np.array([0., 1., 2.])
        raise KeyError(key)

    def initial_TOC(self):
        return 1000.

    def initial_DOC(self):
        return 20.

    def initial_H2O(self):
        return 20000.

    def CO2(self):
        return np.array([0.]), np.array([0.])

    def CH4(self):
        return np.array([0.]), np.array([0.])


class ThermodynamicsLoggingTest(unittest.TestCase):
    def setUp(self):
        version = "0.4"
        parameters = ModelParameters(build_model(version, ['Fe3'], True).parameters())
        self.pathway = pathway_by_name(version, 'Fe3')(parameters)

    def state(self, *, acetate=20., fe3=150., fe2=1., co2=1., water=20000.):
        state = np.ones(len(system.SYSTEM))
        state[system.index('Acetate')] = acetate
        state[system.index('Fe3')] = fe3
        state[system.index('Fe2')] = fe2
        state[system.index('CO2')] = co2
        state[system.index('H2O')] = water
        return state

    def fe3_model(self):
        pathway_model = build_model(['Fe3'], normalized = True)
        pathway_model.model_parameters.set({
            'Fe3_Km_Ac': 50., 'Fe3_Km_Fe3': 500., 'Fe3_v_max': 1.5,
            'Fe3_CUE': .5, 'death_rate': 8.3e-5, 'Fe3_thermodynamics': True,
            'M_Fe3': .25, 'Fe3': 150., 'Acetate': 20.,
        })
        return pathway_model

    def test_delta_g_accepts_a_time_series(self):
        states = np.column_stack((
            self.state(fe2=1.),
            self.state(fe2=10.),
        ))

        expected = np.array([
            self.pathway.deltaG_r(states[:, 0]),
            self.pathway.deltaG_r(states[:, 1]),
        ])

        np.testing.assert_allclose(self.pathway.deltaG_r(states), expected)

    def test_delta_g_is_not_logged_during_pathway_evaluation(self):
        run_log = ModelRun()
        self.pathway.inject_logger(run_log)

        self.pathway.thermodynamics(0., self.state())

        self.assertNotIn('Fe3_deltaG_r', run_log.keys())
        self.assertNotIn('Fe3_thermodynamic_factor', run_log.keys())

    def test_model_logs_delta_g_at_accepted_times(self):
        pathway_model = build_model(['Fe3'])
        pathway_model.model_parameters.set({
            self.pathway.microbe.name: .25,
            'Fe3_Km_Ac': 50.,
            'Fe3_Km_Fe3': 500.,
            'Fe3_v_max': 1.5,
            'Fe3_CUE': .5,
            'death_rate': 8.3e-5,
            'Fe3_thermodynamics': True,
        })
        pathway_model.build()
        pathway_model.replica = Replica()
        pathway_model.system_state_log.reset(None)

        states = np.column_stack((self.state(fe2=1.), self.state(fe2=10.)))
        time = np.array([2., 7.])
        pathway_model.log_accepted_thermodynamics(time, states)

        dissolved_states = np.column_stack([
            pathway_model.compute_dissolved(state) for state in states.T
        ])

        logged_time, logged_delta_g = pathway_model.system_state_log['Fe3_deltaG_r']
        self.assertIs(logged_time, time)
        np.testing.assert_allclose(
            logged_delta_g,
            pathway_model.contributing_pathways[0].deltaG_r(dissolved_states))

    def test_prediction_logs_delta_g_only_at_t_eval(self):
        pathway_model = self.fe3_model()

        run_log = pathway_model.predict_with_logging(Replica())
        time, delta_g = run_log['Fe3_deltaG_r']
        pool_time, fe3 = run_log['Fe3']

        np.testing.assert_array_equal(time, np.array([0., 1., 2.]))
        np.testing.assert_array_equal(time, pool_time)
        states = np.vstack([run_log[pool][1] for pool in system.SYSTEM])
        dissolved_states = np.column_stack([
            pathway_model.compute_dissolved(state) for state in states.T
        ])
        np.testing.assert_allclose(
            delta_g,
            pathway_model.contributing_pathways[0].deltaG_r(dissolved_states))
        self.assertEqual(len(delta_g), len(fe3))

    def test_fast_prediction_omits_diagnostic_series(self):
        pathway_model = self.fe3_model()

        run_log = pathway_model.predict(Replica())

        self.assertNotIn('Fe3_deltaG_r', run_log.keys())
        self.assertNotIn('Fe3_MM', run_log.keys())
        self.assertNotIn('Fe3_v', run_log.keys())

    def test_prediction_keeps_kinetic_diagnostics_from_the_solver_process(self):
        pathway_model = self.fe3_model()

        run_log = pathway_model.predict_with_logging(Replica())
        pool_time, _ = run_log['Fe3']
        for name in ('Fe3_MM', 'Fe3_v'):
            time, values = run_log[name]
            np.testing.assert_array_equal(time, pool_time)
            self.assertEqual(len(values), len(pool_time))

    def test_prediction_logs_each_active_mm_substance_factor(self):
        pathway_model = self.fe3_model()

        run_log = pathway_model.predict_with_logging(Replica())
        time, combined = run_log['Fe3_MM']
        _, acetate = run_log['Fe3_MM_Acetate']
        component_time, fe3 = run_log['Fe3_MM_Fe3']

        np.testing.assert_array_equal(component_time, time)
        np.testing.assert_allclose(combined, acetate * fe3)

    def test_diagnostic_prediction_logs_all_pools_and_effective_km_values(self):
        pathway_model = self.fe3_model()

        run_log = pathway_model.predict_with_logging(Replica())
        time, _ = run_log['TOC']
        for pool in system.SYSTEM:
            dissolved_time, dissolved = run_log[f'{pool}_dissolved']
            np.testing.assert_array_equal(dissolved_time, time)
            self.assertEqual(len(dissolved), len(time))

        _, acetate_km = run_log['Fe3_Km_Acetate']
        _, fe3_km = run_log['Fe3_Km_Fe3']
        total_time, total_km = run_log['Fe3_Km_total']
        np.testing.assert_array_equal(total_time, time)
        np.testing.assert_allclose(total_km, acetate_km + fe3_km)
        self.assertIn('Fe3_reaction_rate', run_log.keys())
        self.assertIn('Fe3_rate_Acetate', run_log.keys())

    def test_diagnostic_replay_preserves_every_lean_pool_trajectory(self):
        """Logging is derived after integration and cannot alter prediction."""
        lean_model = self.fe3_model()
        logged_model = self.fe3_model()

        lean_log = lean_model.predict(Replica())
        diagnostic_log = logged_model.predict_with_logging(Replica())

        for pool in system.SYSTEM:
            lean_time, lean_values = lean_log[pool]
            diagnostic_time, diagnostic_values = diagnostic_log[pool]
            np.testing.assert_array_equal(diagnostic_time, lean_time)
            np.testing.assert_allclose(diagnostic_values, lean_values,
                                       rtol=0, atol=0)


if __name__ == '__main__':
    unittest.main()
