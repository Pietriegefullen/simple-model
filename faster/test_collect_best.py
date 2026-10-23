import json
import os
import tempfile
import unittest
from unittest import mock

import checkpoint_report as plot_best
import matplotlib.pyplot as plt
import numpy as np
import parameters
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

        limits = plot_best.sample_y_limits(Sample())
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
        selected = plot_best.filter_checkpoints(
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

        selected = plot_best.filter_checkpoints(
            checkpoints, model_ids=['3NWE'])
        self.assertEqual(selected, [matching])


class CheckpointCollectionTest(unittest.TestCase):
    def test_ignores_generated_report_artifacts(self):
        """Only generated checkpoint filenames are candidates for loading."""
        with tempfile.TemporaryDirectory() as directory:
            for filename in ('parameter_influence.png', 'sensitivity.csv'):
                with open(os.path.join(directory, filename), 'w') as handle:
                    handle.write('not a checkpoint')
            with open(os.path.join(directory, 'metadata.json'), 'w') as handle:
                json.dump({}, handle)

            with mock.patch.object(parameters, 'load_parameter_file') as load:
                self.assertEqual(plot_best.collect_all(directory), [])

            load.assert_not_called()


class BestCheckpointTest(unittest.TestCase):
    def test_loss_id_does_not_create_a_second_best_fit(self):
        class Checkpoint:
            is_multi_replica = False

            def __init__(self, identifier, loss, loss_id):
                self.id = identifier
                self.loss = loss
                self.loss_id = loss_id
                self.model_id = 'model-A'
                self.fit_mode = 'split'
                self.replica = ('1351', '4')

        retained = plot_best.choose_best([
            Checkpoint('older', 4., 'loss-A'),
            Checkpoint('better', 2., 'loss-B'),
        ])

        self.assertEqual([checkpoint.id for checkpoint in retained], ['better'])


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

        figures = dict(plot_best.fe3_aceto_diagnostic_figures(
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
        self.assertIn('mg', figures['Fe3_Aceto_microbes'][1].get_ylabel())

    def test_combined_delta_g_plots_include_all_and_hydro_homo_pathways(self):
        run_log = {
            f'{pathway}_deltaG_r': ([0, 1], [-20, -21])
            for pathway in plot_best.PATHWAY_PLOT_ORDER
        }

        figures = dict(plot_best.pathway_delta_g_figures(
            None, run_log, suptitle='test'))

        self.assertEqual(set(figures), {
            'All_pathways_deltaG_r', 'Hydro_Homo_deltaG_r',
        })
        all_axis = figures['All_pathways_deltaG_r'][1]
        self.assertEqual(len(all_axis.lines), len(plot_best.PATHWAY_PLOT_ORDER) + 1)
        self.assertEqual(
            [line.get_label() for line in all_axis.lines[:-1]],
            [rf'$\mathrm{{{pathway}\ pathway}}$'
             for pathway in plot_best.PATHWAY_PLOT_ORDER])
        hydro_homo_axis = figures['Hydro_Homo_deltaG_r'][1]
        self.assertEqual(len(hydro_homo_axis.lines), 3)


class PathwayQuantityPlotTest(unittest.TestCase):
    def tearDown(self):
        plt.close('all')

    def test_delta_g_plot_marks_the_gibbs_minimum(self):
        _, axis = plot_best.plot_quantity(
            'Hydro_deltaG_r', [0, 1], [-20, -30], None, suptitle='test')

        self.assertEqual(len(axis.lines), 2)
        self.assertEqual(axis.lines[-1].get_linestyle(), '--')
        self.assertEqual(list(axis.lines[-1].get_ydata()),
                         [1e-3*GIBBS_MINIMUM, 1e-3*GIBBS_MINIMUM])

    def test_inhibition_factor_uses_a_zero_to_one_axis(self):
        _, axis = plot_best.plot_quantity(
            'Fermentation_inhib', [0, 1], [.2, .8], None, suptitle='test')

        self.assertEqual(axis.get_ylim(), (0.0, 1.0))

    def test_rate_plot_includes_v_max_with_headroom(self):
        _, axis = plot_best.plot_quantity(
            'Hydro_v', [0, 1], [.2, .8], None, v_max=1., suptitle='test')

        self.assertEqual(axis.get_ylim(), (0.0, 1.05))

    def test_unit_constant_v_max_labels_m_pool_as_pathway_capacity(self):
        _, axis = plot_best.plot_quantity(
            'M_Hydro', [0, 1], [.2, .8], None, suptitle='test',
            capacity_pathway='Hydro')

        axis.figure.canvas.draw()
        self.assertEqual(axis.get_title(), 'Hydro pathway capacity (M_Hydro)')
        ylabel = axis.get_ylabel().replace('\\', '')
        self.assertIn('pathway capacity', ylabel)
        self.assertIn('mg/g dry weight', ylabel)
        self.assertTrue(any('$v_{\\max} = 1$ (constant)' in text.get_text()
                            for text in axis.texts))

    def test_only_a_fixed_unit_v_max_marks_m_as_capacity(self):
        class RunLog(dict):
            pass

        run_log = RunLog({
            'Hydro_v': ([0, 1], [.2, .8]),
            'M_Hydro': ([0, 1], [.2, .8]),
        })
        run_log._parameters = parameters.ModelParameters({
            'Hydro_v_max': parameters.Parameter('Hydro_v_max', 1.).constant(1.),
        })

        quantities = list(
            plot_best.pathway_plot_quantities(run_log, include_pathway=True))
        microbe_quantity = next(item for item in quantities if item[0] == 'M_Hydro')
        self.assertEqual(microbe_quantity[-1], 'Hydro')
        self.assertTrue(
            plot_best.pathway_uses_capacity_parameter(run_log, 'Hydro'))

        run_log._parameters['Hydro_v_max'].variable(1., (.5, 1.5))
        self.assertFalse(
            plot_best.pathway_uses_capacity_parameter(run_log, 'Hydro'))

    def test_microbe_pool_uses_mass_biomass_units(self):
        _, axis = plot_best.plot_quantity(
            'M_Hydro', [0, 1], [.2, .8], None, suptitle='test')

        ylabel = axis.get_ylabel().replace('\\', '')
        self.assertIn('microbial biomass', ylabel)
        self.assertIn('mg/g dry weight', ylabel)
        self.assertNotIn('μmol', ylabel)

    def test_dissolved_h2_follows_total_h2_in_pool_plot_order(self):
        class PathwayModel:
            def compute_dissolved(self, states):
                assert states.ndim == 1
                dissolved = states.copy()
                dissolved[system.index('H2')] *= .25
                return dissolved

        run_log = {
            'H2': ([0, 1], [4., 8.]),
            'H2O': ([0, 1], [10., 10.]),
        }

        quantities = list(plot_best.pathway_plot_quantities(
            run_log, pathway_model=PathwayModel()))
        self.assertEqual([quantity[0] for quantity in quantities],
                         ['H2', 'H2_dissolved', 'H2O'])
        dissolved_h2 = quantities[1]
        np.testing.assert_allclose(dissolved_h2[1], [0, 1])
        self.assertEqual(dissolved_h2[2].tolist(), [1., 2.])

    def test_dissolved_h2_has_a_clear_plot_title(self):
        _, axis = plot_best.plot_quantity(
            'H2_dissolved', [0, 1], [.1, .2], None, suptitle='test')

        self.assertEqual(axis.get_title(), 'Dissolved H2')

    def test_reaction_rate_uses_biomass_on_the_specific_rate_grid(self):
        run_log = {
            'Hydro_v': ([0., .5, 1.5, 2.], [.1, .2, .4, .5]),
            'M_Hydro': ([0., 1., 2.], [2., 4., 8.]),
        }

        pathway, time, actual_rate, maximum_rate = next(
            plot_best.pathway_reaction_rate_curves(
                run_log, {'Hydro': .5}))

        self.assertEqual(pathway, 'Hydro')
        np.testing.assert_allclose(time, [0., .5, 1.5, 2.])
        np.testing.assert_allclose(actual_rate, [.2, .6, 2.4, 4.])
        np.testing.assert_allclose(maximum_rate, [1., 1.5, 3., 4.])

    def test_reaction_rate_plot_marks_biomass_scaled_v_max_as_dashed(self):
        _, axis = plot_best.plot_reaction_rate(
            'Hydro', [0, 1], [.2, .8], [1., 2.], None, suptitle='test')

        self.assertEqual(len(axis.lines), 2)
        self.assertEqual(axis.lines[0].get_ydata().tolist(), [.2, .8])
        self.assertEqual(axis.lines[1].get_ydata().tolist(), [1., 2.])
        self.assertEqual(axis.lines[1].get_linestyle(), '--')

    def test_limitation_factor_plot_combines_kinetics_and_thermodynamics(self):
        run_log = {
            'Homo_v': ([0, 1], [.2, .3]),
            'Homo_MM': ([0, 1], [.4, .8]),
            'Homo_thermodynamic_factor': ([0, 1], [.7, .9]),
            'Homo_inhib': ([0, 1], [.5, .25]),
            'Hydro_v': ([0, 1], [.2, .3]),
            'Hydro_MM': ([0, 1], [.5, .6]),
        }

        figures = dict(plot_best.pathway_limitation_factor_figures(
            None, run_log, suptitle='test'))
        homo_axis = figures['Homo_limitation_factors'][1]
        hydro_axis = figures['Hydro_limitation_factors'][1]

        self.assertEqual(homo_axis.get_ylim(), (0.0, 1.05))
        self.assertEqual(homo_axis.lines[0].get_color(), 'tab:blue')
        self.assertEqual(homo_axis.lines[1].get_linestyle(), '--')
        self.assertEqual(homo_axis.lines[2].get_linestyle(), '-.')
        np.testing.assert_allclose(homo_axis.lines[3].get_ydata(), [.14, .18])
        self.assertTrue(any('not logged' in text.get_text()
                            for text in hydro_axis.texts))

    def test_ch4_pathway_production_is_stacked_with_final_amounts(self):
        run_log = {
            'CH4': ([0, 1], [3, 6]),
            'Hydro_CH4 from M_Hydro': ([0, 1], [2, 2]),
            'Aceto_CH4 from M_Ac': ([0, 1], [1, 1]),
        }

        figure, axis = plot_best.ch4_pathway_production_figure(
            None, run_log, suptitle='test')
        _, log_axis = plot_best.ch4_pathway_production_figure(
            None, run_log, suptitle='test', log_scale=True)

        self.assertEqual(axis.get_yscale(), 'linear')
        self.assertEqual(log_axis.get_yscale(), 'log')
        self.assertEqual(len(axis.collections), 2)
        self.assertEqual(axis.lines[0].get_ydata().tolist(), [0, 3])
        labels = [text.get_text() for text in axis.get_legend().get_texts()]
        self.assertTrue(any('2' in label for label in labels))
        self.assertTrue(any('1' in label for label in labels))
        self.assertIn('3', labels[-1])
        self.assertIsNotNone(figure)

    def test_co2_pathway_production_includes_all_producers_on_both_scales(self):
        def pathway(name, co2_coefficient):
            instance = type(name, (), {})()
            instance.pathway_vector = [0.] * len(system.SYSTEM)
            instance.pathway_vector[system.index('CO2')] = co2_coefficient
            return instance

        pathway_model = type('Model', (), {
            'contributing_pathways': [
                pathway('Fermentation', 3.),
                pathway('Aceto', 1.),
                pathway('Hydro', -1.),
            ],
        })()
        run_log = {
            'CO2': ([0, 1], [3., 6.]),
            'Fermentation_v': ([0, 1], [.5, .5]),
            'M_Ferm': ([0, 1], [2., 2.]),
            'Aceto_v': ([0, 1], [1., 1.]),
            'M_Ac': ([0, 1], [3., 3.]),
            'Hydro_v': ([0, 1], [.5, .5]),
            'M_Hydro': ([0, 1], [4., 4.]),
        }

        figure, axis = plot_best.co2_pathway_production_figure(
            None, pathway_model, run_log, suptitle='test')
        _, log_axis = plot_best.co2_pathway_production_figure(
            None, pathway_model, run_log, suptitle='test', log_scale=True)

        self.assertEqual(axis.get_yscale(), 'linear')
        self.assertEqual(log_axis.get_yscale(), 'log')
        self.assertEqual(len(axis.collections), 2)
        labels = [text.get_text() for text in axis.get_legend().get_texts()]
        self.assertTrue(any('Fermentation' in label for label in labels))
        self.assertTrue(any('Aceto' in label for label in labels))
        self.assertFalse(any('Hydro' in label for label in labels))
        self.assertIn('3', labels[-1])
        self.assertIsNotNone(figure)

    def test_co2_production_uses_accepted_states_not_trial_rate_logs(self):
        class Fermentation:
            pathway_vector = np.zeros(len(system.SYSTEM))

            def __call__(self, time, state):
                return system.vector(0., 'CO2', state[system.index('DOC')])

        Fermentation.pathway_vector[system.index('CO2')] = 1.
        pathway_model = type('Model', (), {
            'contributing_pathways': [Fermentation()],
            'compute_dissolved': staticmethod(lambda state: state),
        })()
        run_log = {
            pool: ([0., 1.], [1., 2.]) for pool in system.SYSTEM
        }
        run_log['CO2'] = ([0., 1.], [3., 4.5])
        # Deliberately incompatible trial-state diagnostics: these must not
        # influence the cumulative curve when full accepted states are logged.
        run_log['Fermentation_v'] = ([0., 1.], [100., 100.])
        run_log['M_Ferm'] = ([0., 1.], [100., 100.])

        _, axis = plot_best.co2_pathway_production_figure(
            None, pathway_model, run_log, suptitle='test')

        np.testing.assert_allclose(axis.collections[0].get_paths()[0].vertices[:, 1].max(),
                                   1.5)
        np.testing.assert_allclose(axis.lines[0].get_ydata(), [0., 1.5])


if __name__ == '__main__':
    unittest.main()
