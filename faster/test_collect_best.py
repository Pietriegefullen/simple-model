import unittest

import collect_best
import matplotlib.pyplot as plt
import system


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
        self.assertEqual(figures['Fe3_Aceto_deltaG_r'][1].lines[-1].get_linestyle(), '--')
        self.assertEqual(figures['Fe3_Aceto_MM'][1].get_ylim(), (0.0, 1.0))


if __name__ == '__main__':
    unittest.main()
