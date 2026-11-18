import unittest
from tempfile import TemporaryDirectory
from unittest import mock

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt

import plot_checkpoint
from USER_VARIABLES import PLOT_DIRECTORY


matplotlib.rcParams['text.usetex'] = False


class Checkpoint:
    def __init__(self, identifier, *, model='model-X', replica='1351-4',
                 fit='split', loss='loss-X', total_loss=.123456):
        self.id = identifier
        self.model_id = model
        self.replica_name = replica
        self.fit_mode = fit
        self.loss_id = loss
        self.loss = total_loss


class ComparisonLabelTest(unittest.TestCase):
    def test_common_attributes_move_to_title_and_ids_stay_in_legend(self):
        title, labels = plot_checkpoint.comparison_labels([
            Checkpoint('ONE'), Checkpoint('TWO'),
        ])

        self.assertIn('model model-X', title)
        self.assertIn('sample 1351', title)
        self.assertIn('replica 4', title)
        self.assertEqual(labels, ['cp-ONE', 'cp-TWO'])

    def test_only_differing_attributes_are_added_to_each_legend_label(self):
        title, labels = plot_checkpoint.comparison_labels([
            Checkpoint('ONE', model='model-A', replica='1351-4'),
            Checkpoint('TWO', model='model-B', replica='1351-5'),
        ])

        self.assertIn('fit split', title)
        self.assertNotIn('model ', title)
        self.assertIn('sample 1351', title)
        self.assertIn('model model-A', labels[0])
        self.assertIn('replica 5', labels[1])
        self.assertNotIn('fit split', labels[0])


class CommandLineTest(unittest.TestCase):
    def test_single_plot_filename_places_replica_after_figure_number(self):
        checkpoint = Checkpoint('ABCD', replica='1351-4')

        self.assertEqual(
            plot_checkpoint.single_plot_filename(0, checkpoint),
            '00_1351-4_fit.png')

    def test_default_output_directory_is_grouped_below_plot_directory(self):
        checkpoint = Checkpoint('ABCD', replica='1351-4', fit='split')
        result = plot_checkpoint.CheckpointPlot(checkpoint, (None, None, {}))

        self.assertEqual(plot_checkpoint.default_output_directory([result]),
                         plot_checkpoint.Path(PLOT_DIRECTORY) / 'checkpoint' /
                         'split' / '1351-4' / 'cp-ABCD_loss-0.123456')

    def test_comparison_default_keeps_the_requested_checkpoint_order(self):
        first = plot_checkpoint.CheckpointPlot(
            Checkpoint('ABCD', replica='1351-4', fit='split'), (None, None, {}))
        second = plot_checkpoint.CheckpointPlot(
            Checkpoint('EFGH', replica='1351-4', fit='split'), (None, None, {}))

        self.assertEqual(plot_checkpoint.default_output_directory([first, second]),
                         plot_checkpoint.Path(PLOT_DIRECTORY) /
                         'checkpoint-comparison' / 'split' / '1351-4' /
                         'cp-ABCD_loss-0.123456__cp-EFGH_loss-0.123456')

    def test_explicit_output_directory_overrides_the_default(self):
        args = plot_checkpoint.parse_args(['ABCD', '--output-dir', 'chosen'])

        self.assertEqual(args.output_dir, 'chosen')

    def test_best_restart_is_an_alternative_to_best(self):
        args = plot_checkpoint.parse_args(['--best-restart'])

        self.assertTrue(args.best_restart)
        with self.assertRaises(SystemExit):
            plot_checkpoint.parse_args(['--best', '--best-restart'])

    def test_named_comparison_uses_a_comparison_directory(self):
        args = plot_checkpoint.parse_args([
            'ABCD', '--compare', 'D', '*_thermodynamics', 'False'])

        self.assertEqual(
            args.comparisons,
            [plot_checkpoint.ComparisonSpec(
                'D', (('*_thermodynamics', False),))])
        result = plot_checkpoint.CheckpointPlot(
            Checkpoint('ABCD', replica='1351-4', fit='split'), (None, None, {}))
        replay = plot_checkpoint.CheckpointPlot(
            Checkpoint('ABCD', replica='1351-4', fit='split'), (None, None, {}), 'D')
        self.assertEqual(plot_checkpoint.default_output_directory(
                             [result, replay], args.comparisons),
                         plot_checkpoint.Path(PLOT_DIRECTORY) /
                         'checkpoint-comparison' / 'split' / '1351-4' /
                         'cp-ABCD_loss-0.123456__D')

    def test_compare_requires_complete_parameter_value_pairs(self):
        with self.assertRaises(ValueError):
            plot_checkpoint.parse_args(['ABCD', '--compare', 'D', 'Hydro_v_max'])

    def test_existing_comparison_figure_is_skipped_unless_overwritten(self):
        with TemporaryDirectory() as temporary:
            path = plot_checkpoint.Path(temporary) / '00_fit.png'
            path.touch()
            first, _ = plt.subplots()
            with mock.patch.object(
                    plot_checkpoint, 'comparison_figures',
                    return_value={'fit': (first, None)}):
                self.assertEqual(plot_checkpoint.save_comparison(
                    [object()], temporary, dataset={}, overwrite=False), [])

            second, _ = plt.subplots()
            with mock.patch.object(
                    plot_checkpoint, 'comparison_figures',
                    return_value={'fit': (second, None)}):
                written = plot_checkpoint.save_comparison(
                    [object()], temporary, dataset={}, overwrite=True)

        self.assertEqual(written, [path])


class BestRestartCheckpointTest(unittest.TestCase):
    class Checkpoint:
        is_multi_replica = False

        def __init__(self, identifier, loss, *, run_id='run-A', restart=None,
                     replica=('1351', '4')):
            self.id = identifier
            self.loss = loss
            self.run_id = run_id
            self.model_id = 'model-A'
            self.fit_mode = 'single'
            self.replica = replica
            self.metadata = ({} if restart is None
                             else {'cma_restart': {'index': restart}})

    def test_uses_ordinary_best_when_no_restarts_were_saved(self):
        selected = plot_checkpoint.best_restart_checkpoints([
            self.Checkpoint('older', 3.), self.Checkpoint('best', 1.),
        ])

        self.assertEqual([checkpoint.id for checkpoint in selected], ['best'])

    def test_uses_best_usable_checkpoint_from_each_restart(self):
        selected = plot_checkpoint.best_restart_checkpoints([
            self.Checkpoint('ordinary-best', .1),
            self.Checkpoint('restart-0-old', 4., restart=0),
            self.Checkpoint('restart-0-best', 2., restart=0),
            self.Checkpoint('restart-1-best', 3., restart=1),
            self.Checkpoint('restart-2-penalty', 1e100, restart=2),
        ])

        self.assertEqual(
            {checkpoint.id for checkpoint in selected},
            {'restart-0-best', 'restart-1-best'})

    def test_keeps_restart_numbers_from_independent_runs_separate(self):
        selected = plot_checkpoint.best_restart_checkpoints([
            self.Checkpoint('first-run', 1., run_id='run-A', restart=0),
            self.Checkpoint('second-run', 2., run_id='run-B', restart=0),
        ])

        self.assertEqual(
            {checkpoint.id for checkpoint in selected}, {'first-run', 'second-run'})


class ComparisonReplayTest(unittest.TestCase):
    def test_wildcards_expand_against_checkpoint_parameters(self):
        class Parameters:
            def get_values(self):
                return {
                    'Hydro_thermodynamics': True,
                    'Ferm_thermodynamics': True,
                    'Hydro_v_max': 1,
                }

        checkpoint = type('Checkpoint', (), {
            'parameters': Parameters(), 'id': 'ONE'})()
        specification = plot_checkpoint.ComparisonSpec(
            'D', (('*_thermodynamics', False),))

        self.assertEqual(plot_checkpoint._expanded_overrides(checkpoint, specification), {
            'Ferm_thermodynamics': False,
            'Hydro_thermodynamics': False,
        })

    def test_custom_variant_name_distinguishes_replayed_curve(self):
        original = plot_checkpoint.CheckpointPlot(
            Checkpoint('ONE', model='model-A'), (None, None, {}))
        replay = plot_checkpoint.CheckpointPlot(
            Checkpoint('ONE', model='model-A'), (None, None, {}), 'D')

        title, labels = plot_checkpoint.comparison_labels([original, replay])

        self.assertNotIn('model ', title)
        self.assertEqual(labels, [
            'cp-ONE | model model-A',
            'cp-ONE | model D',
        ])


class ComparisonQuantityTest(unittest.TestCase):
    def tearDown(self):
        plt.close('all')

    def test_matching_logged_quantities_share_one_figure_and_legend(self):
        first = plot_checkpoint.CheckpointPlot(
            Checkpoint('ONE'), (None, None, {
                'TOC': ([0, 1], [3., 2.]),
                'Hydro_MM': ([0, 1], [.4, .8]),
            }))
        second = plot_checkpoint.CheckpointPlot(
            Checkpoint('TWO'), (None, None, {
                'TOC': ([0, 1], [4., 1.]),
                'Hydro_MM': ([0, 1], [.5, .6]),
            }))

        figures = plot_checkpoint._comparison_quantities([first, second])
        axis = figures['TOC'].axes[0]

        self.assertEqual(len(axis.lines), 2)
        self.assertEqual([text.get_text() for text in axis.get_legend().texts],
                         ['cp-ONE', 'cp-TWO'])
        self.assertEqual(figures['Hydro_MM'].axes[0].get_ylim(), (0.0, 1.0))

    def test_single_result_diagnostic_keeps_its_legend(self):
        result = plot_checkpoint.CheckpointPlot(
            Checkpoint('ONE'), (None, None, {'TOC': ([0, 1], [3., 2.])}))

        figure = plot_checkpoint._comparison_quantities([result])['TOC']

        self.assertEqual([text.get_text() for text in figure.axes[0].get_legend().texts],
                         ['cp-ONE'])

    def test_split_fit_names_are_available_for_interactive_viewers(self):
        result = plot_checkpoint.CheckpointPlot(
            Checkpoint('ONE'), (None, None, {}))

        names = plot_checkpoint.comparison_plot_names([result], split_fit=True)

        self.assertEqual(names[:2], ['fit_CO2', 'fit_CH4'])

    def test_limitation_factors_share_one_pathway_figure(self):
        first = plot_checkpoint.CheckpointPlot(
            Checkpoint('ONE'), (None, None, {
                'Hydro_v': ([0, 1], [.2, .3]),
                'Hydro_MM': ([0, 1], [.4, .8]),
            }))
        second = plot_checkpoint.CheckpointPlot(
            Checkpoint('TWO'), (None, None, {
                'Hydro_v': ([0, 1], [.2, .3]),
                'Hydro_MM': ([0, 1], [.5, .6]),
            }))

        figures = plot_checkpoint._comparison_limitation_factors([first, second])
        axis = figures['Hydro_limitation_factors'].axes[0]

        # One MM line and one all-factor line per checkpoint, plus reference.
        self.assertEqual(len(axis.lines), 5)
        self.assertIn('cp-ONE', [text.get_text() for text in axis.get_legend().texts])

    def test_reaction_rate_signature_replaces_per_pathway_comparison_plots(self):
        class RunLog(dict):
            pass

        run_log = RunLog({
            'Hydro_v': ([0, 1], [.2, .3]),
            'M_Hydro': ([0, 1], [2., 3.]),
        })
        run_log._parameters = type('Parameters', (), {
            '__contains__': lambda self, name: name == 'Hydro_v_max',
            '__getitem__': lambda self, name: type('Parameter', (), {
                '__float__': lambda self: 1.,
                'is_constant': lambda self: True,
            })(),
        })()

        result = plot_checkpoint.CheckpointPlot(
            Checkpoint('ONE'), (None, None, run_log))
        names = plot_checkpoint.comparison_plot_names([result])
        figures = plot_checkpoint._comparison_reaction_rate_signatures([result])

        self.assertIn('reaction_rate_signature', names)
        self.assertIn('reaction_rate_signature_log', names)
        self.assertNotIn('Hydro_reaction_rate', names)
        self.assertEqual(figures['reaction_rate_signature'].axes[0].get_yscale(), 'linear')
        self.assertEqual(figures['reaction_rate_signature_log'].axes[0].get_yscale(), 'log')

    def test_pathway_plot_flag_removes_that_pathway_from_the_suite(self):
        class RunLog(dict):
            pass

        run_log = RunLog({
            'Hydro_v': ([0, 1], [.2, .3]),
            'M_Hydro': ([0, 1], [2., 3.]),
        })
        run_log._parameters = type('Parameters', (), {
            '__contains__': lambda self, name: name == 'Hydro_v_max',
            '__getitem__': lambda self, name: type('Parameter', (), {
                '__float__': lambda self: 1.,
                'is_constant': lambda self: True,
            })(),
        })()
        result = plot_checkpoint.CheckpointPlot(
            Checkpoint('ONE'), (None, None, run_log))

        with mock.patch.dict(
                plot_checkpoint.report.PATHWAY_PLOT_FLAGS, {'Hydro': False}):
            names = plot_checkpoint.comparison_plot_names([result])

        self.assertNotIn('Hydro_v', names)
        self.assertNotIn('reaction_rate_signature', names)
        self.assertNotIn('Hydro_limitation_factors', names)


if __name__ == '__main__':
    unittest.main()
