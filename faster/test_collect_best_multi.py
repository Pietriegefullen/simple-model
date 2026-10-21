import os
import tempfile
import unittest
from unittest import mock

import plot_best
import model
import parameters


def multi_checkpoint(loss, checkpoint_id, replicas):
    pathway_model = model.build_model(['Hydro'], normalized=True)
    run_config = {
        'model': pathway_model.get_config(only_structure=True),
        'chosen': {
            'sample': 0,
            'validation_replica': 0,
            't_start': None,
            't_end': None,
            'fit_mode': 'multi',
            'pathways': ['Hydro'],
            'parameter_override': {},
            'normalized_parameters': True,
            'algorithm': 'differential_evolution',
        },
        'objective': {
            'loss_weight': {'CO2': 1., 'CH4': 1.},
            'reduction': {'CO2': 'mse', 'CH4': 'mse'},
            'transform': {'CO2': ['normalize'], 'CH4': ['log', 'normalize']},
        },
        'algo': {},
        'range': {},
        'multi_replica': {'shared': ['Hydro_v_max'], 'replicas': replicas},
    }
    records = [
        {
            'sample': replica[:-1],
            'replica_number': replica[-1],
            'replica': replica,
            'loss': .1,
            'r2': .9,
            'checkpoint': f'/inner/{replica}',
            'parameters': {'M_Hydro': .01},
        }
        for replica in replicas
    ]
    return plot_best.Checkpoint(
        parameters.ModelParameters({'Hydro_v_max': 1.}), loss, run_config,
        checkpoint_id, f'/outer/{checkpoint_id}', {'replicas': records})


class MultiCheckpointCollectionTest(unittest.TestCase):
    def test_distinct_joint_configurations_are_not_collapsed_to_placeholder_replica(self):
        first = multi_checkpoint(3., 'cp-FIRST', ['13514', '13515'])
        replacement = multi_checkpoint(2., 'cp-BETTER', ['13514', '13515'])
        second = multi_checkpoint(1., 'cp-SECOND', ['13694', '13695'])

        with mock.patch.object(plot_best, 'collect_all',
                               return_value=[first, replacement, second]):
            selected = plot_best.collect_best()

        self.assertEqual({checkpoint.id for checkpoint in selected},
                         {'cp-BETTER', 'cp-SECOND'})

    def test_filters_match_replicas_inside_joint_fit(self):
        checkpoint = multi_checkpoint(1., 'cp-TEST', ['13514', '13515'])

        selected = plot_best.filter_checkpoints(
            [checkpoint], samples=['1351'], replicas=['5'], fit_modes=['multi'])

        self.assertEqual(selected, [checkpoint])

    def test_directory_and_manifest_identify_the_exact_joint_fit(self):
        checkpoint = multi_checkpoint(1., 'cp-TEST', ['13514', '13515'])
        directory = plot_best.multi_plot_directory(checkpoint)
        self.assertIn('_multi_', os.path.basename(directory))
        self.assertIn(checkpoint.run_id, directory)
        self.assertIn(checkpoint.id, directory)

        with tempfile.TemporaryDirectory() as temporary_directory:
            manifest = plot_best.write_multi_fit_manifest(
                temporary_directory, checkpoint)
            with open(manifest) as handle:
                contents = handle.read()

        self.assertIn('cp-TEST', contents)
        self.assertIn('/inner/13514', contents)
        self.assertIn('FULL FIT CONFIGURATION', contents)


if __name__ == '__main__':
    unittest.main()
