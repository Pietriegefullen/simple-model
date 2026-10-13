import unittest
from unittest import mock

import fit_all
import optimizer


class StagedFitHelpersTest(unittest.TestCase):
    def test_model_identifier_uses_the_current_model_builder(self):
        chosen = {
            'pathways': fit_all.PATHWAYS,
            'normalized_parameters': True,
            'parameter_override': {},
        }

        self.assertTrue(fit_all.model_identifier(chosen).startswith('model-'))

    def test_restart_config_selects_requested_survivor_counts(self):
        chosen = {
            'sample': '1351',
            'validation_replica': '4',
        }
        self.assertEqual(fit_all.stage_init_config(chosen, 'model-id', 0),
                         {'default': True})
        self.assertEqual(
            fit_all.stage_init_config(chosen, 'model-id', 1)['best_N'], 10)
        self.assertEqual(
            fit_all.stage_init_config(chosen, 'model-id', 2)['best_N'], 3)

    def test_split_budget_counts_all_fitted_replicas(self):
        class Sample:
            replicas = [object(), object(), object(), object()]

        objective_calls, replica_count = fit_all._objective_calls_for_model_budget(
            Sample(), 'split', 100)
        self.assertEqual((objective_calls, replica_count), (34, 3))

    def test_initialization_only_stage_uses_exact_objective_budget(self):
        class Sample:
            replicas = [object(), object(), object()]

        config, objective_calls, replica_count = fit_all.stage_algorithm_config(
            0, 101, Sample(), 'split', workers=1, seed=7,
            scipy_polish=False, initialization_only=True,
            initialization='halton')

        self.assertEqual((objective_calls, replica_count), (51, 2))
        self.assertEqual(config, {
            'initial_sampling': {
                'init': 'halton',
                'samples': 51,
                'seed': 7,
                'workers': 1,
            },
        })

    def test_stage_label_is_not_forwarded_to_differential_evolution(self):
        class Sample:
            replicas = [object(), object(), object()]

        config, _, _ = fit_all.stage_algorithm_config(
            0, 100, Sample(), 'split', workers=1, seed=7)
        algorithm = optimizer.DifferentialEvolution()
        with self.assertWarnsRegex(RuntimeWarning, 'ignored unsupported.*name'):
            algorithm.configure(config['differential_evolution'])

        self.assertNotIn('name', algorithm.get_config())

    def test_exploration_retries_until_the_r2_target_is_reached(self):
        class Sample:
            sample_name = '1351'
            replicas = [object(), object()]

        class Replica:
            replica_number = '4'

        class Objective:
            def __init__(self, r2):
                self.r2 = r2

            def best_R2(self):
                return self.r2

            def call_count(self):
                return 17

            def best_loss(self):
                return 0.1

        objectives = [Objective(0.2), Objective(0.6), Objective(0.7),
                      Objective(0.8)]
        with mock.patch.object(fit_all, 'model_identifier', return_value='model-id'), \
                mock.patch.object(fit_all, 'fit',
                                  side_effect=[(None, objective, None)
                                               for objective in objectives]):
            result = fit_all.run_one_fit(
                None, 'default', fit_all.PATHWAYS, {}, Sample(), Replica(),
                'single', (1, 1, 1), 1, 7, False, mock.Mock(),
                stage_repeats=(3, 1, 1), r2_target=0.5)

        self.assertEqual(
            [(entry['name'], entry['attempt']) for entry in result],
            [('explore', 1), ('explore', 2), ('refine_top_10', 1),
             ('polish_top_3', 1)])


if __name__ == '__main__':
    unittest.main()
