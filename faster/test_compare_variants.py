import unittest
from unittest import mock

import compare_checkpoints
import compare_variants
import matplotlib.pyplot as plt
import parameters
import system


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
        specs = compare_variants.variant_run_specs(
            {'A': a, 'B': b}, ('A', 'B', 'C'))
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

    def test_d_is_derived_from_a_parameters_with_thermodynamics_disabled(self):
        a = Checkpoint('model-FA75-a', 1.0, parameter_value=2.0)
        a.parameters['Hydro_v_max'].set(7.0)
        for pathway in ('Hydro', 'Homo', 'Aceto', 'Fe3'):
            a.parameters[f'{pathway}_thermodynamics'].constant(True)

        specifications = compare_variants.variant_run_specs({'A': a}, ('A', 'D'))
        self.assertEqual([specification.variant for specification in specifications],
                         ['A', 'D'])
        d = specifications[-1]

        self.assertEqual(d.variant, 'D')
        self.assertEqual(d.run_config['chosen']['pathways'],
                         a.run_config['chosen']['pathways'])
        self.assertEqual(set(d.run_config['model']), set(a.run_config['model']))
        self.assertEqual(d.run_config['range'], a.run_config['range'])
        self.assertEqual(d.parameter_values.get_values()['Hydro_v_max'], 7.0)
        self.assertEqual(d.parameter_values.get_values()['Homo_v_max'], 2.0)
        disabled = {
            name: value for name, value in d.parameter_values.get_values().items()
            if name.endswith('_thermodynamics')
        }
        self.assertEqual(disabled, {
            'Hydro_thermodynamics': False,
            'Homo_thermodynamics': False,
            'Aceto_thermodynamics': False,
            'Fe3_thermodynamics': False,
        })
        self.assertEqual(d.run_config['chosen']['parameter_override'], disabled)
        self.assertEqual(
            {name: config['use_thermodynamics']
             for name, config in d.run_config['model'].items()},
            {'Hydro': False, 'Homo': False})

    def test_variant_selection_derives_only_the_requested_counterfactual(self):
        a = Checkpoint('model-FA75-a', 1.0)
        b = Checkpoint('model-EM33-b', 1.0)
        for pathway in ('Hydro', 'Homo', 'Aceto', 'Fe3'):
            a.parameters[f'{pathway}_thermodynamics'].constant(True)

        specifications = compare_variants.variant_run_specs(
            {'A': a, 'B': b}, ('A', 'D'))

        self.assertEqual([specification.variant for specification in specifications],
                         ['A', 'D'])


class H2ConsumptionTest(unittest.TestCase):
    def test_pathway_rate_includes_biomass_and_h2_stoichiometry(self):
        pathway_type = type('Homo', (), {})
        pathway = pathway_type()
        pathway.stoichiometry = [0.0] * len(system.SYSTEM)
        pathway.stoichiometry[system.index('H2')] = -4.0
        pathway_model = type('Model', (), {
            'contributing_pathways': [pathway],
        })()
        run_log = {
            'Homo_v': ([0, 1, 2], [.5, .5, .5]),
            'M_Homo': ([0, 1, 2], [1, 2, 3]),
        }

        time, consumption = compare_variants.pathway_h2_consumption_rate(
            pathway_model, run_log, 'Homo')

        self.assertEqual(time.tolist(), [0.0, 1.0, 2.0])
        self.assertEqual(consumption.tolist(), [2.0, 4.0, 6.0])


class ParameterComparisonTest(unittest.TestCase):
    def test_rate_product_rows_combine_initial_biomass_and_v_max(self):
        rows = compare_checkpoints.rate_product_rows([
            {'M_Hydro': 2., 'Hydro_v_max': .5},
            {'M_Hydro': 4., 'Hydro_v_max': .25},
        ], show_constant=True)

        row = next(row for row in rows if row['name'] == 'Hydro_M0_v_max')
        self.assertEqual(row['group'], 'Hydro')
        self.assertEqual(row['values'].tolist(), [1., 1.])
        self.assertIn('M(0)', row['label'])

    def test_rate_product_rows_respect_constant_row_filtering(self):
        rows = compare_checkpoints.rate_product_rows([
            {'M_Hydro': 2., 'Hydro_v_max': .5},
            {'M_Hydro': 4., 'Hydro_v_max': .25},
        ], show_constant=False)

        self.assertFalse(any(row['name'] == 'Hydro_M0_v_max' for row in rows))


class ABDifferenceTest(unittest.TestCase):
    def test_co2_difference_is_percent_relative_to_model_a(self):
        difference = compare_variants.relative_percent_difference(
            [100, 50], [125, 25])
        self.assertEqual(difference.tolist(), [25.0, -50.0])

    def test_methane_log_change_is_symmetric_for_equal_fold_changes(self):
        difference = compare_variants.methane_log_relative_change(
            [1, 1], [2, .5])
        self.assertAlmostEqual(difference[0], -difference[1])
        self.assertAlmostEqual(difference[0], 100 * __import__('math').log(2))

    def test_shared_limits_expand_only_the_outlying_replica(self):
        def detail(co2, ch4):
            return {'differences': {
                'CO2': ([0, 1], co2),
                'CH4': ([0, 1], ch4),
            }}

        typical = detail([1, -2], [3, -4])
        outlier = detail([2, -3], [5, -100])
        shared = compare_variants.shared_difference_limits([typical, outlier])

        self.assertEqual(
            compare_variants.replica_difference_limits(typical, shared)['CO2'],
            (-shared['CO2'], shared['CO2']))
        self.assertGreater(
            compare_variants.replica_difference_limits(outlier, shared)['CH4'][1],
            shared['CH4'])


class ABScoreScatterTest(unittest.TestCase):
    def tearDown(self):
        plt.close('all')

    def test_scatter_has_equal_axes_and_an_equality_diagonal(self):
        figure = compare_variants.plot_ab_r2_scatter([
            {'fit_mode': 'single', 'r2_a': .7, 'r2_b': .8},
            {'fit_mode': 'split', 'r2_a': .9, 'r2_b': .85},
        ])
        axis = figure.axes[0]

        self.assertEqual(axis.get_xlim(), axis.get_ylim())
        diagonal = axis.lines[0]
        self.assertEqual(diagonal.get_xdata().tolist(), diagonal.get_ydata().tolist())


if __name__ == '__main__':
    unittest.main()
