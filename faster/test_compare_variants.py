import unittest
from unittest import mock

import compare_variants
import parameters


class Checkpoint:
    def __init__(self, model_id, loss, parameter_value=1.0):
        self.model_id = model_id
        self.loss = loss
        self.loss_id = 'loss-test'
        self.fit_mode = 'split'
        self.replica = ('1351', '4')
        self.parameters = parameters.ModelParameters(
            {'Homo_v_max': parameter_value})
        self.run_config = {
            'model': {'Hydro': {}, 'Homo': {}},
            'chosen': {'pathways': ['Hydro', 'Homo']},
            'range': {'Hydro_v_max': {}, 'Homo_v_max': {}},
            'objective': {},
        }


class ComparisonGroupsTest(unittest.TestCase):
    def test_groups_keep_the_lowest_loss_checkpoint_per_variant(self):
        first_a = Checkpoint('model-FA75-first', 2.0)
        best_a = Checkpoint('model-PIRL-best', 1.0)
        b = Checkpoint('model-EM33-b', 3.0)

        groups = compare_variants.comparison_groups((first_a, best_a, b))

        self.assertEqual(len(groups), 1)
        group = ('loss-test', 'split', ('1351', '4'))
        self.assertIs(groups[group]['A'], best_a)
        self.assertIs(groups[group]['B'], b)

    def test_target_is_separate_by_loss_mode_and_replica(self):
        group = ('loss-test', 'split', ('1351', '4'))
        target = compare_variants.comparison_target('/tmp/comparisons', group)
        self.assertEqual(
            str(target), '/tmp/comparisons/loss-test_split/1351-4')

    def test_groups_ignore_models_outside_the_requested_variants(self):
        a = Checkpoint('model-FA75-a', 1.0)
        b = Checkpoint('model-EM33-b', 1.0)
        unknown = Checkpoint('model-UNKNOWN-x', 1.0)

        groups = compare_variants.comparison_groups((a, b, unknown))

        group = ('loss-test', 'split', ('1351', '4'))
        self.assertEqual(set(groups[group]), {'A', 'B'})

    def test_c_is_derived_from_a_parameters_and_omits_hydro(self):
        a = Checkpoint('model-FA75-a', 1.0, parameter_value=2.0)
        b = Checkpoint('model-EM33-b', 1.0, parameter_value=3.0)
        a.parameters['Hydro_v_max'].set(7.0)
        specs = compare_variants.variant_run_specs({'A': a, 'B': b}, ('A', 'B'))
        self.assertEqual([spec.variant for spec in specs], ['A', 'B', 'C'])

        with mock.patch.object(compare_variants, 'run',
                               side_effect=lambda config, values: (config, values)):
            runs = compare_variants.run_variants(specs)

        c_config, c_values = runs['C']
        self.assertNotIn('Hydro', c_config['chosen']['pathways'])
        self.assertNotIn('Hydro', c_config['model'])
        self.assertNotIn('Hydro_v_max', c_config['range'])
        self.assertEqual(c_values.get_values(), {'Homo_v_max': 2.0})
        self.assertNotIn('Hydro_v_max', c_values.get_values())
        self.assertIsNot(c_values, a.parameters)
        self.assertEqual(runs['A'][1].get_values(), a.parameters.get_values())
        self.assertEqual(runs['B'][1].get_values(), b.parameters.get_values())


if __name__ == '__main__':
    unittest.main()
