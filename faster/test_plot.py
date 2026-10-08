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

    def test_mosaic_creates_axes_from_its_labels(self):
        _, axes = plot.create_figure(plot.FigureSpec(
            mosaic=[['CO2', 'CH4'], ['summary', 'summary']],
            figsize=(8, 6),
        ))
        self.assertEqual(list(axes), ['CO2', 'CH4', 'summary'])


if __name__ == '__main__':
    unittest.main()
