import unittest

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt

import plot

# ``plot`` preserves the project's LaTex preference; use Matplotlib's built-in
# renderer here so this test also verifies that labels do not require LaTex.
matplotlib.rcParams['text.usetex'] = False


class Replica:
    sample = 'test sample'

    def CO2(self):
        return [0, 1], [1.0, 2.0]

    def CH4(self):
        return [0, 1], [0.1, 0.2]


class FigureSpecTest(unittest.TestCase):
    def tearDown(self):
        plt.close('all')

    def test_named_grid_axes_are_used_by_pool_plotters(self):
        figure, axes = plot.create_figure(plot.FigureSpec(
            ncols=2,
            axis_names=('CO2', 'CH4'),
            figsize=(8, 4),
            sharex=True,
        ))
        run_log = {
            'CO2': ([0, 1], [1.1, 2.1]),
            'CH4': ([0, 1], [0.11, 0.21]),
        }

        plot.plot_replica_pools(Replica(), axes)
        plot.plot_run_pools(run_log, axes)

        self.assertEqual(figure.get_size_inches().tolist(), [8.0, 4.0])
        self.assertEqual(list(axes), ['CO2', 'CH4'])
        self.assertEqual([len(axis.lines) for axis in axes.values()], [2, 2])
        figure.canvas.draw()

    def test_default_size_scales_with_the_number_of_axes(self):
        figure, _ = plot.create_figure(plot.FigureSpec(axis_names=('only',)))

        self.assertEqual(figure.get_size_inches().tolist(), [4.0, 4.0])

    def test_log_scale_requires_non_negative_data_and_sets_positive_limits(self):
        figure, axis = plt.subplots()
        axis.plot([0, 1], [0, 10])

        self.assertTrue(plot.format_ax(axis, log_scale=True))
        self.assertEqual(axis.get_yscale(), 'log')
        self.assertGreater(axis.get_ylim()[0], 0)

        _, negative_axis = plt.subplots()
        negative_axis.plot([0, 1], [-1, 10])
        self.assertFalse(plot.format_ax(negative_axis, log_scale=True))
        self.assertEqual(negative_axis.get_yscale(), 'linear')

    def test_mosaic_creates_axes_from_its_labels(self):
        _, axes = plot.create_figure(plot.FigureSpec(
            mosaic=[['CO2', 'CH4'], ['summary', 'summary']],
            figsize=(8, 6),
        ))
        self.assertEqual(list(axes), ['CO2', 'CH4', 'summary'])

    def test_model_variant_line_styles_are_centralized(self):
        self.assertEqual(plot.MODEL_VARIANT_LINESTYLES,
                         {'A': '-', 'B': '--', 'C': ':', 'D': '-.'})

    def test_figure_plan_reuses_named_secondary_axis(self):
        calls = []

        def curve(context, axis, *, value, label):
            calls.append((context.checkpoint, axis, value, label))
            axis.plot([0, 1], [value, value], label=label)

        plan = plot.FigurePlan(
            plot.FigureSpec(axis_names=('gas',), figsize=(4, 3)),
            calls=(
                plot.PlotCall(curve, 'gas', {'value': 1, 'label': 'primary'}),
                plot.PlotCall(curve, 'gas', {'value': 2, 'label': 'right one'},
                              secondary=True),
                plot.PlotCall(curve, 'gas', {'value': 3, 'label': 'right two'},
                              secondary=True),
            ),
            suptitle='planned figure',
        )

        figure, axes = plot.render_figure(
            plan, plot.PlotContext(checkpoint='cp-TEST'))

        self.assertEqual(list(axes), ['gas', 'gas_right'])
        self.assertEqual(len(axes['gas'].lines), 1)
        self.assertEqual(len(axes['gas_right'].lines), 2)
        self.assertTrue(all(checkpoint == 'cp-TEST'
                            for checkpoint, *_ in calls))
        self.assertEqual(figure._suptitle.get_text(), 'planned figure')

    def test_label_keeps_words_in_text_mode_and_units_in_math_mode(self):
        self.assertEqual(plot.label('time', 'd'), r'time [$\mathrm{d}$]')
        self.assertEqual(plot.math(r'CO_2'), r'$CO_2$')

    def test_plot_book_renders_each_named_figure_with_one_context(self):
        seen = []

        def record_context(context, axis):
            seen.append(context)
            axis.plot([0, 1], [0, 1])

        book = plot.PlotBook({
            'first': plot.FigurePlan(
                plot.FigureSpec(axis_names=('axis',)),
                (plot.PlotCall(record_context, 'axis'),)),
            'second': plot.FigurePlan(
                plot.FigureSpec(axis_names=('axis',)),
                (plot.PlotCall(record_context, 'axis'),)),
        })
        context = plot.PlotContext(checkpoint='cp-TEST')

        figures = plot.render_book(book, context)

        self.assertEqual(list(figures), ['first', 'second'])
        self.assertEqual(seen, [context, context])


if __name__ == '__main__':
    unittest.main()
