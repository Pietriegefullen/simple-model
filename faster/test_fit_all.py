import unittest

import fit_all


class StagedFitHelpersTest(unittest.TestCase):
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


if __name__ == '__main__':
    unittest.main()
