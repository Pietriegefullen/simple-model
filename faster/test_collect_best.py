import unittest

import collect_best
import matplotlib.pyplot as plt
import system
from chemistry import GIBBS_MINIMUM


class SampleYLimitsTest(unittest.TestCase):
    def test_limits_include_every_replica_and_log_ch4_excludes_zero(self):
        class Replica:
            def __init__(self, co2, ch4):
                self._co2 = co2
                self._ch4 = ch4

            def CO2(self):
                return [0, 1], self._co2

            def CH4(self):
                return [0, 1], self._ch4

        class Sample:
            replicas = [
                Replica([1.0, 4.0], [0.0, 1e-3]),
                Replica([2.0, 10.0], [1e-4, 1e-2]),
            ]

        limits = collect_best.sample_y_limits(Sample())
        self.assertLess(limits['CO2'][0], 1.0)
        self.assertGreater(limits['CO2'][1], 10.0)
        self.assertGreater(limits['CH4'][0], 0.0)
        self.assertLess(limits['CH4'][0], 1e-4)
        self.assertGreater(limits['CH4'][1], 1e-2)


class CheckpointFilterTest(unittest.TestCase):
    def test_filters_are_combined(self):
        class Checkpoint:
            def __init__(self, sample, replica, model_id, fit_mode, loss_id):
                self.replica = (sample, replica)
                self.model_id = model_id
                self.fit_mode = fit_mode
                self.loss_id = loss_id

        checkpoints = [
            Checkpoint('1351', '4', 'model-A', 'single', 'loss-1'),
            Checkpoint('1351', '5', 'model-A', 'split', 'loss-1'),
            Checkpoint('1352', '4', 'model-B', 'split', 'loss-2'),
        ]
        selected = collect_best.filter_checkpoints(
            checkpoints, samples=['1351'], model_ids=['model-A'],
            fit_modes=['split'], replicas=['5'], loss_ids=['loss-1'])
        self.assertEqual(selected, [checkpoints[1]])

    def test_partial_model_id_matches_the_generated_identifier(self):
        class Checkpoint:
            def __init__(self, model_id):
                self.replica = ('1351', '4')
                self.model_id = model_id
                self.fit_mode = 'split'
                self.loss_id = 'loss-1'

        matching = Checkpoint('model-3NWE-JDSY')
        checkpoints = [matching, Checkpoint('model-FA75-GSEM')]

        selected = collect_best.filter_checkpoints(
            checkpoints, model_ids=['3NWE'])
        self.assertEqual(selected, [matching])


class Fe3AcetoDiagnosticsTest(unittest.TestCase):
    def tearDown(self):
        plt.close('all')

    def test_diagnostics_combine_the_two_pathways_and_acetate_uptake(self):
        def pathway(name, acetate_coefficient):
            instance = type(name, (), {})()
            instance.pathway_vector = [0.0] * len(system.SYSTEM)
            instance.pathway_vector[system.index('Acetate')] = acetate_coefficient
            return instance

        pathway_model = type('Model', (), {
            'contributing_pathways': [
                pathway('Fe3', -1.5),
                pathway('Aceto', -1.25),
            ],
        })()
        run_log = {
            'Acetate': ([0, 1], [20, 16]),
            'M_Fe3': ([0, 1], [2, 3]),
            'M_Ac': ([0, 1], [4, 5]),
            'Fe3_deltaG_r': ([0, 1], [-20, -21]),
            'Aceto_deltaG_r': ([0, 1], [-22, -23]),
            'Fe3_v': ([0, 1], [.5, .5]),
            'Aceto_v': ([0, 1], [.2, .2]),
            'Fe3_MM': ([0, 1], [.3, .4]),
            'Aceto_MM': ([0, 1], [.5, .6]),
        }

        figures = dict(collect_best.fe3_aceto_diagnostic_figures(
            None, (pathway_model, None, run_log), suptitle='test'))

        self.assertEqual(set(figures), {
            'Fe3_Aceto_microbes',
            'Fe3_Aceto_deltaG_r',
            'Fe3_Aceto_Acetate_consumption',
            'Fe3_Aceto_MM',
        })
        acetate_axis = figures['Fe3_Aceto_Acetate_consumption'][1]
        self.assertEqual(acetate_axis.lines[0].get_ydata().tolist(), [1.5, 2.25])
        self.assertEqual(acetate_axis.lines[1].get_ydata().tolist(), [1.0, 1.25])
        acetate_pool_axis = acetate_axis.figure.axes[1]
        self.assertEqual(acetate_pool_axis.lines[0].get_ydata().tolist(), [20, 16])
        self.assertIn('Acetate', acetate_pool_axis.get_ylabel())
        self.assertLessEqual(acetate_axis.get_position().x1, .76)
        self.assertEqual(figures['Fe3_Aceto_deltaG_r'][1].lines[-1].get_linestyle(), '--')
        self.assertEqual(figures['Fe3_Aceto_MM'][1].get_ylim(), (0.0, 1.0))


class PathwayQuantityPlotTest(unittest.TestCase):
    def tearDown(self):
        plt.close('all')

    def test_delta_g_plot_marks_the_gibbs_minimum(self):
        _, axis = collect_best.plot_quantity(
            'Hydro_deltaG_r', [0, 1], [-20, -30], None, suptitle='test')

        self.assertEqual(len(axis.lines), 2)
        self.assertEqual(axis.lines[-1].get_linestyle(), '--')
        self.assertEqual(list(axis.lines[-1].get_ydata()),
                         [GIBBS_MINIMUM, GIBBS_MINIMUM])

    def test_inhibition_factor_uses_a_zero_to_one_axis(self):
        _, axis = collect_best.plot_quantity(
            'Fermentation_inhib', [0, 1], [.2, .8], None, suptitle='test')

        self.assertEqual(axis.get_ylim(), (0.0, 1.0))

    def test_rate_plot_includes_v_max_with_headroom(self):
        _, axis = collect_best.plot_quantity(
            'Hydro_v', [0, 1], [.2, .8], None, v_max=1., suptitle='test')

        self.assertEqual(axis.get_ylim(), (0.0, 1.05))

    def test_ch4_pathway_production_is_stacked_with_final_amounts(self):
        run_log = {
            'CH4': ([0, 1], [3, 6]),
            'Hydro_CH4 from M_Hydro': ([0, 1], [2, 2]),
            'Aceto_CH4 from M_Ac': ([0, 1], [1, 1]),
        }

        figure, axis = collect_best.ch4_pathway_production_figure(
            None, run_log, suptitle='test')

        self.assertEqual(axis.get_yscale(), 'linear')
        self.assertEqual(len(axis.collections), 2)
        self.assertEqual(axis.lines[0].get_ydata().tolist(), [0, 3])
        labels = [text.get_text() for text in axis.get_legend().get_texts()]
        self.assertTrue(any('2' in label for label in labels))
        self.assertTrue(any('1' in label for label in labels))
        self.assertIn('3', labels[-1])
        self.assertIsNotNone(figure)


if __name__ == '__main__':
    unittest.main()
