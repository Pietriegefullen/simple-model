import unittest
import json

import numpy as np

import CONSTANTS
import optimizer
import system
from parameters import ModelParameters
from model_run import ModelRun

from model_factory import build_model, pathway_by_name
from model_v0_4 import model as model_v0_4


class Replica:
    sample = type('Sample', (), {'sample_name': '1351'})()
    dry_weight = 1.
    water_content = 1.
    total_gas_volume = 100.
    incubation = {'days': np.array([0., 1., 2.])}

    def __getitem__(self, key):
        if key == 'days':
            return np.array([0., 1., 2.])
        if key == 'CO2':
            return self.CO2()[1]
        if key == 'CH4':
            return self.CH4()[1]
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


class ObjectiveReplica(Replica):
    def CO2(self):
        return np.array([0., 1., 2.]), np.array([1., 2., 3.])

    def CH4(self):
        return np.array([0., 1., 2.]), np.array([1., 2., 3.])


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

    def test_gas_runaway_guard_is_disabled_when_no_limits_are_requested(self):
        """The optimiser must see the real trajectory unless it opts in."""
        initial_state = np.zeros(len(system.SYSTEM))
        co2_index = system.index('CO2')

        def produces_co2(_time, _state):
            derivative = np.zeros(len(system.SYSTEM))
            derivative[co2_index] = 1.
            return derivative

        without_limits = []
        model_v0_4.integrate(
            produces_co2, np.array([0., 1.]), initial_state,
            without_limits, None, None)
        self.assertAlmostEqual(without_limits[0].y[co2_index, -1], 1.)

        with self.assertRaises(optimizer.IntegrationFailure):
            model_v0_4.integrate(
                produces_co2, np.array([0., 1.]), initial_state,
                [], None, {'CO2': .5, 'CH4': np.inf})

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

        run_log = pathway_model.predict_with_logging(
            Replica(), reference='L', microbe_unit='micromol')
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

    def test_prediction_keeps_the_candidate_parameters_for_checkpointing(self):
        pathway_model = self.fe3_model()
        pathway_model.parameters()['Fe3_Km_Ac'].set(75.)
        pathway_model.parameters()['Fe3_v_max'].set(2.5)

        pathway_model.predict(Replica())

        self.assertEqual(pathway_model.parameters()['Fe3_Km_Ac'].value, 75.)
        self.assertEqual(pathway_model.parameters()['Fe3_v_max'].value, 2.5)
        fe3 = pathway_model.contributing_pathways[0]
        self.assertEqual(fe3.Km[system.index('Acetate')], 75.)
        self.assertEqual(float(fe3.v_max), 2.5)

    def test_prediction_does_not_change_the_optimisation_variable_list(self):
        pathway_model = self.fe3_model()
        pathway_model.parameters()['TOC'].variable(100., (1., 1_000.))
        names_before = [parameter.name for parameter
                        in pathway_model.parameters().variables()]

        pathway_model.predict(Replica())

        names_after = [parameter.name for parameter
                       in pathway_model.parameters().variables()]
        self.assertEqual(names_after, names_before)

    def test_saved_parameter_values_replay_with_identical_loss_and_r2(self):
        model = self.fe3_model()
        model.parameters()['Fe3_Km_Ac'].set(75.)
        model.parameters()['Fe3_v_max'].set(2.5)
        objective_config = {
            'loss_weight': {'CO2': 1., 'CH4': 1.},
            'reduction': {'CO2': 'mse', 'CH4': 'mse'},
            'transform': {'CO2': ['normalize'], 'CH4': ['log', 'normalize']},
        }
        objective = optimizer.build_objective_function(
            model, [ObjectiveReplica()], objective_config)
        loss = objective(model.parameters(), transformed=False)
        r2 = objective.last_R2()

        # Match the JSON representation stored in an ordinary checkpoint.
        saved_values = json.loads(json.dumps(model.parameters().get_values()))
        replay_model = self.fe3_model()
        replay_objective = optimizer.build_objective_function(
            replay_model, [ObjectiveReplica()], objective_config)
        replay_loss = replay_objective(
            ModelParameters(saved_values), transformed=False)

        self.assertEqual(replay_loss, loss)
        self.assertEqual(replay_objective.last_R2(), r2)

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

    def test_predict_with_logging_matches_predict_in_every_requested_unit(self):
        """Diagnostics must not alter either trajectories or output units."""
        replica = Replica()
        for reference in ('L', 'gdw'):
            for microbe_unit in ('micromol', 'mg'):
                lean_log = self.fe3_model().predict(
                    replica, reference=reference, microbe_unit=microbe_unit)
                diagnostic_log = self.fe3_model().predict_with_logging(
                    replica, reference=reference, microbe_unit=microbe_unit)
                for pool in system.SYSTEM:
                    lean_time, lean_values = lean_log[pool]
                    diagnostic_time, diagnostic_values = diagnostic_log[pool]
                    np.testing.assert_array_equal(diagnostic_time, lean_time)
                    np.testing.assert_allclose(
                        diagnostic_values, lean_values, rtol=0, atol=0)

    def test_predict_output_conversions_are_consistent(self):
        """The public L/micromol output converts exactly to gdw/mg output."""
        replica = Replica()
        litre_log = self.fe3_model().predict(
            replica, reference='L', microbe_unit='micromol')
        gdw_log = self.fe3_model().predict(
            replica, reference='gdw', microbe_unit='mg')
        L_per_gdw = (replica.water_content * 1e-3) / replica.dry_weight

        for pool in system.SYSTEM:
            time, litre_values = litre_log[pool]
            gdw_time, gdw_values = gdw_log[pool]
            expected = litre_values * L_per_gdw
            if pool.startswith('M_'):
                expected *= CONSTANTS.MOLAR_MASS_BIOMASS
            np.testing.assert_array_equal(gdw_time, time)
            np.testing.assert_allclose(gdw_values, expected, rtol=1e-12, atol=0)

    def test_internal_initial_state_is_litre_based(self):
        """The v0.4 RHS must receive experimental pools in micromol/L."""
        replica = Replica()
        run_log = self.fe3_model().predict(
            replica, reference='L', microbe_unit='micromol')
        L_per_gdw = (replica.water_content * 1e-3) / replica.dry_weight

        self.assertEqual(run_log['TOC'][1][0],
                         replica.initial_TOC() / L_per_gdw)
        self.assertEqual(run_log['DOC'][1][0],
                         replica.initial_DOC() / L_per_gdw)

    def test_accepted_delta_g_is_independent_of_display_reference(self):
        """Thermodynamics must use reconstructed solver (µmol/L) states."""
        replica = Replica()
        litre_model = self.fe3_model()
        litre_log = litre_model.predict_with_logging(
            replica, reference='L', microbe_unit='micromol')
        gdw_model = self.fe3_model()
        gdw_log = gdw_model.predict_with_logging(
            replica, reference='gdw', microbe_unit='mg')

        time, delta_g_litre = litre_log['Fe3_deltaG_r']
        gdw_time, delta_g_gdw = gdw_log['Fe3_deltaG_r']
        np.testing.assert_array_equal(gdw_time, time)
        np.testing.assert_allclose(delta_g_gdw, delta_g_litre,
                                   rtol=1e-12, atol=0)

        # Independently reconstruct the internal state from the displayed
        # gdw/mg pools and verify the logged reaction energies.
        L_per_gdw = (replica.water_content * 1e-3) / replica.dry_weight
        states = np.vstack([gdw_log[pool][1] for pool in system.SYSTEM])
        states /= L_per_gdw
        for pool_index, pool in enumerate(system.SYSTEM):
            if pool.startswith('M_'):
                states[pool_index] /= CONSTANTS.MOLAR_MASS_BIOMASS
        dissolved_states = np.column_stack([
            gdw_model.compute_dissolved(state) for state in states.T
        ])
        expected = gdw_model.contributing_pathways[0].deltaG_r(dissolved_states)
        np.testing.assert_allclose(delta_g_gdw, expected, rtol=1e-12, atol=0)


if __name__ == '__main__':
    unittest.main()
