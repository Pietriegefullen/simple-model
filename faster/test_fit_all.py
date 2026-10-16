import unittest
import tempfile
from unittest import mock

import checkpoint
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

    def test_selected_samples_excludes_requested_samples(self):
        class Sample:
            def __init__(self, sample_name, has_replicas=True):
                self.sample_name = sample_name
                self._has_replicas = has_replicas

            def has_replicas(self):
                return self._has_replicas

        class Dataset:
            samples = [Sample('1351'), Sample('1352'), Sample('1353'),
                       Sample('1354', has_replicas=False)]

        selected = fit_all.selected_samples(
            Dataset(), requested_samples=['1351', '1352', '1353'],
            excluded_samples=['1352'])

        self.assertEqual([sample.sample_name for sample in selected],
                         ['1351', '1353'])

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

    def test_exploration_uses_exact_size_sobol_initialization_by_default(self):
        class Sample:
            replicas = [object(), object(), object()]

        config, objective_calls, replica_count = fit_all.stage_algorithm_config(
            0, 101, Sample(), 'split', workers=1, seed=7)

        self.assertEqual((objective_calls, replica_count), (51, 2))
        self.assertEqual(config, {
            'initial_sampling': {
                'init': 'sobol',
                'samples': 51,
                'seed': 7,
                'workers': 1,
            },
        })

    def test_exploration_can_opt_back_into_differential_evolution(self):
        class Sample:
            replicas = [object(), object(), object()]

        config, _, _ = fit_all.stage_algorithm_config(
            0, 101, Sample(), 'split', workers=1, seed=7,
            exploration_initializer=False)

        self.assertIn('differential_evolution', config)

    def test_exploration_configures_a_worker_stall_timeout(self):
        class Sample:
            replicas = [object(), object()]

        config, objective_calls, _ = fit_all.stage_algorithm_config(
            0, 100, Sample(), 'single', workers=1, seed=7,
            exploration_call_limit=23, exploration_stall_timeout=4)

        self.assertEqual(objective_calls, 23)
        self.assertEqual(config['initial_sampling']['stall_timeout'], 4)

    def test_stage_label_is_not_forwarded_to_differential_evolution(self):
        class Sample:
            replicas = [object(), object(), object()]

        config, _, _ = fit_all.stage_algorithm_config(
            1, 100, Sample(), 'split', workers=1, seed=7)
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

    def test_exploration_always_uses_the_call_limit(self):
        class Sample:
            sample_name = '1351'
            replicas = [object(), object()]

        class Replica:
            replica_number = '4'

        class Objective:
            def best_R2(self):
                return 0.6

            def call_count(self):
                return 17

            def best_loss(self):
                return 0.1

        with mock.patch.object(fit_all, 'model_identifier', return_value='model-id'), \
                mock.patch.object(fit_all, 'fit', return_value=(None, Objective(), None)) as fit:
            fit_all.run_one_fit(
                None, 'default', fit_all.PATHWAYS, {}, Sample(), Replica(),
                'single', (100, 1, 1), 1, 7, False, mock.Mock(),
                stage_repeats=(1, 1, 1), r2_target=0.5,
                explore_call_limit=23)

        self.assertEqual(fit.call_args_list[0].kwargs['call_limit'], 23)

    def test_confirmed_interrupt_skips_only_the_current_stage(self):
        class Sample:
            sample_name = '1351'
            replicas = [object(), object()]

        class Replica:
            replica_number = '4'

        class Objective:
            def best_R2(self):
                return 0.6

            def call_count(self):
                return 17

            def best_loss(self):
                return 0.1

        with mock.patch.object(fit_all, 'model_identifier', return_value='model-id'), \
                mock.patch.object(fit_all, 'fit', side_effect=[
                    KeyboardInterrupt(), (None, Objective(), None),
                    (None, Objective(), None)]) as fit, \
                mock.patch('builtins.input', return_value='yes'):
            result = fit_all.run_one_fit(
                None, 'default', fit_all.PATHWAYS, {}, Sample(), Replica(),
                'single', (1, 1, 1), 1, 7, False, mock.Mock(),
                stage_repeats=(1, 1, 1), r2_target=0.5)

        self.assertEqual(fit.call_count, 3)
        self.assertEqual(result[0]['status'], 'skipped_by_user')
        self.assertEqual([entry['name'] for entry in result[1:]],
                         ['refine_top_10', 'polish_top_3'])

    def test_call_limit_continues_to_the_following_stage(self):
        class Sample:
            sample_name = '1351'
            replicas = [object(), object()]

        class Replica:
            replica_number = '4'

        class Objective:
            def best_R2(self):
                return 0.6

            def call_count(self):
                return 17

            def best_loss(self):
                return 0.1

        with mock.patch.object(fit_all, 'model_identifier', return_value='model-id'), \
                mock.patch.object(fit_all, 'fit', side_effect=[
                    checkpoint.CallLimitInterrupt(), (None, Objective(), None),
                    (None, Objective(), None)]) as fit:
            result = fit_all.run_one_fit(
                None, 'default', fit_all.PATHWAYS, {}, Sample(), Replica(),
                'single', (100, 1, 1), 1, 7, False, mock.Mock(),
                stage_repeats=(1, 1, 1), r2_target=0.5,
                explore_call_limit=23)

        self.assertEqual(fit.call_count, 3)
        self.assertEqual(result[0]['status'], 'call_limit')
        self.assertEqual([entry['name'] for entry in result[1:]],
                         ['refine_top_10', 'polish_top_3'])

    def test_stall_timeout_continues_to_the_following_stage(self):
        class Sample:
            sample_name = '1351'
            replicas = [object(), object()]

        class Replica:
            replica_number = '4'

        class Objective:
            def best_R2(self):
                return 0.6

            def call_count(self):
                return 17

            def best_loss(self):
                return 0.1

        with mock.patch.object(fit_all, 'model_identifier', return_value='model-id'), \
                mock.patch.object(fit_all, 'fit', side_effect=[
                    checkpoint.StalledEvaluationInterrupt(4),
                    (None, Objective(), None), (None, Objective(), None)]) as fit:
            result = fit_all.run_one_fit(
                None, 'default', fit_all.PATHWAYS, {}, Sample(), Replica(),
                'single', (100, 1, 1), 1, 7, False, mock.Mock(),
                stage_repeats=(1, 1, 1), r2_target=0.5,
                explore_call_limit=23, explore_stall_timeout=4)

        self.assertEqual(fit.call_count, 3)
        self.assertEqual(result[0]['status'], 'stall_timeout')
        self.assertEqual(result[0]['timeout_seconds'], 4)

    def test_second_interrupt_during_confirmation_aborts_the_fit(self):
        class Sample:
            sample_name = '1351'
            replicas = [object(), object()]

        class Replica:
            replica_number = '4'

        with mock.patch.object(fit_all, 'model_identifier', return_value='model-id'), \
                mock.patch.object(fit_all, 'fit', side_effect=KeyboardInterrupt()), \
                mock.patch('builtins.input', side_effect=KeyboardInterrupt()):
            with self.assertRaises(KeyboardInterrupt):
                fit_all.run_one_fit(
                    None, 'default', fit_all.PATHWAYS, {}, Sample(), Replica(),
                    'single', (1, 1, 1), 1, 7, False, mock.Mock(),
                    stage_repeats=(1, 1, 1), r2_target=0.5)

    def test_good_existing_fit_is_reported_as_skipped(self):
        class Replica:
            replica_number = '4'

        class Sample:
            sample_name = '1351'
            replicas = [Replica()]

            def has_replicas(self):
                return True

        class Dataset:
            samples = [Sample()]

        with tempfile.TemporaryDirectory() as output, \
                mock.patch.object(fit_all.data, 'get_data_before_day',
                                  return_value=Dataset()), \
                mock.patch.object(
                    fit_all, 'run_one_fit',
                    side_effect=fit_all.ExistingFitGoodEnough(0.91)):
            args = fit_all.parse_args([
                '--output', output, '--variants', 'default',
                '--fit-modes', 'single'])
            summary = fit_all.run_batch(args)

        self.assertEqual(len(summary['skipped']), 1)
        self.assertEqual(summary['skipped'][0]['reason'],
                         'existing_fit_good_enough')
        self.assertEqual(summary['skipped'][0]['best_r2'], 0.91)


if __name__ == '__main__':
    unittest.main()
