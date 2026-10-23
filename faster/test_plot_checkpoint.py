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
        self.assertIn('replica 1351-4', title)
        self.assertEqual(labels, ['cp-ONE', 'cp-TWO'])

    def test_only_differing_attributes_are_added_to_each_legend_label(self):
        title, labels = plot_checkpoint.comparison_labels([
            Checkpoint('ONE', model='model-A', replica='1351-4'),
            Checkpoint('TWO', model='model-B', replica='1351-5'),
        ])

        self.assertIn('fit split', title)
        self.assertNotIn('model ', title)
        self.assertIn('model model-A', labels[0])
        self.assertIn('replica 1351-5', labels[1])
        self.assertNotIn('fit split', labels[0])


class CommandLineTest(unittest.TestCase):
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


if __name__ == '__main__':
    unittest.main()
