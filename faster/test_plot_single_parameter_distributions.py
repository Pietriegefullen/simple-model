import unittest
from unittest import mock

import matplotlib.pyplot as plt
import parameters
import plot_single_parameter_distributions as distributions


class Checkpoint:
    def __init__(self, model_id, fit_mode='single', value=1.0):
        self.model_id = model_id
        self.fit_mode = fit_mode
        self.loss_id = 'loss-test'
        self.replica = ('1351', '4')
        self.id = 'checkpoint-1'
        self.loss = 1.0
        self.parameters = parameters.ModelParameters({'Homo_v_max': value})
        self.run_config = {}


class SingleParameterDistributionsTest(unittest.TestCase):
    def tearDown(self):
        plt.close('all')

    def test_single_checkpoints_uses_retained_single_fits_only(self):
        single = Checkpoint('model-A')
        split = Checkpoint('model-B', fit_mode='split')
        with mock.patch.object(distributions.collect_best, 'collect_best',
                               return_value=[split, single]):
            selected = distributions.single_checkpoints('/ignored')
        self.assertEqual(selected, [single])

    def test_r2_filter_is_strict_and_ignores_nonfinite_scores(self):
        checkpoints = [Checkpoint('model-A'), Checkpoint('model-B'), Checkpoint('model-C')]
        with mock.patch.object(distributions, 'fit_r2', side_effect=[.8, .8001, float('nan')]):
            accepted = distributions.scored_checkpoints(checkpoints, .8)
        self.assertEqual(accepted, [(checkpoints[1], .8001)])

    def test_parameter_rows_keep_model_identity_and_fit_r2(self):
        checkpoint = Checkpoint('model-A', value=2.5)
        rows = distributions.parameter_rows([(checkpoint, .91)])
        self.assertEqual(rows, [{
            'model_id': 'model-A', 'loss_id': 'loss-test', 'sample': '1351',
            'replica': '4', 'checkpoint_id': 'checkpoint-1', 'loss': 1.0,
            'r2': .91, 'parameter': 'Homo_v_max', 'value': 2.5,
        }])

    def test_parameter_plot_labels_total_and_per_model_checkpoint_counts(self):
        figure, axis = plt.subplots()
        distributions._draw_parameter(
            axis, 'Homo_v_max', {'model-A': [1.0, 2.0], 'model-B': [3.0]})
        self.assertEqual(axis.get_title(), 'Homo_v_max (n=3)')
        self.assertEqual(
            [label.get_text() for label in axis.get_xticklabels()],
            ['model-A\n(n=2)', 'model-B\n(n=1)'])


if __name__ == '__main__':
    unittest.main()
