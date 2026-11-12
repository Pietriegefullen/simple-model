import unittest
import tempfile
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import checkpoint_report
import model_inspector


class CheckpointFilenameMetadataTest(unittest.TestCase):
    def test_extracts_catalogue_metadata_from_canonical_filename(self):
        metadata = checkpoint_report.checkpoint_filename_metadata(
            '/results/loss-ABCD_model-EFGH_1351_4_split_loss-1_run-IJK_cp-LMN')

        self.assertEqual(metadata, {
            'model_id': 'model-EFGH',
            'run_id': 'run-IJK',
            'checkpoint_id': 'cp-LMN',
            'sample': '1351',
            'replica': '4',
            'fit_mode': 'split',
        })

    def test_catalogue_reads_each_checkpoint_file_once(self):
        payload = {
            'parameters': {},
            'total_loss': .25,
            'run_config': {
                'model': {'version': '0.4'},
                'chosen': {'sample': '1351', 'validation_replica': '4',
                           'fit_mode': 'single'},
                'objective': {},
            },
        }
        with tempfile.TemporaryDirectory() as directory:
            Path(directory, 'loss-X_model-Y_1351_4_single_run-Z_cp-Q').touch()
            with patch.object(checkpoint_report.parameters, 'load_file',
                              return_value=payload) as load_file, \
                 patch.object(checkpoint_report.parameters, 'load_parameter_file',
                              side_effect=AssertionError('must not parse twice')):
                checkpoints = checkpoint_report.collect_all(directory, verbose=False)

        self.assertEqual(len(checkpoints), 1)
        load_file.assert_called_once()


class InspectorFilterTest(unittest.TestCase):
    class Checkpoint:
        def __init__(self, sample, loss=1.):
            self.loss = loss
            self.id = 'cp-from-file'
            self.run_config = {'model': {'version': '0.4'}}
            self.replica = (sample, '4')
            self.fit_mode = 'single'

        @property
        def model_id(self):
            raise AssertionError('Model ID must not be read for a sample filter.')

        @property
        def run_id(self):
            raise AssertionError('Run ID must not be read for a sample filter.')

    def test_sample_filter_does_not_evaluate_unrequested_ids(self):
        first = self.Checkpoint('1351', loss=.1)
        second = self.Checkpoint('1352', loss=.2)

        matched = model_inspector.filter_checkpoints(
            [first, second],
            {'sample': '1351', 'model_id': '', 'run_id': '',
             'replica': '', 'fit_mode': '', 'model_version': '',
             'checkpoint_id': ''})

        self.assertEqual(matched, [first])

    def test_new_matching_checkpoints_are_shown_before_older_results(self):
        older = self.Checkpoint('1351')
        older.source = 'older.json'
        newer = self.Checkpoint('1351')
        newer.source = 'newer.json'

        ordered = model_inspector.prioritise_new_checkpoints(
            [older, newer], {'newer.json'})

        self.assertEqual(ordered, [newer, older])


class DirectCheckpointTest(unittest.TestCase):
    def test_loads_requested_checkpoint_without_collecting_catalogue(self):
        source = Path('/results/example_cp-ABCD-EFGH')
        contents = {
            'total_loss': .25,
            'replicas': [],
        }
        run_config = {
            'model': {'version': '0.4'},
            'chosen': {'sample': '1351', 'validation_replica': '4',
                       'fit_mode': 'single'},
            'objective': {},
        }
        with patch.object(model_inspector, 'resolve_checkpoint', return_value=source) as resolve, \
             patch.object(model_inspector, 'replay_inputs',
                          return_value=(source, run_config, 'parameters', contents)) as replay, \
             patch.object(model_inspector.report, 'collect_all',
                          side_effect=AssertionError('must not load the catalogue')):
            checkpoint = model_inspector.load_checkpoint('ABCD')

        resolve.assert_called_once_with('ABCD')
        replay.assert_called_once_with(source)
        self.assertEqual(checkpoint.id, 'cp-ABCD-EFGH')
        self.assertEqual(checkpoint.parameters, 'parameters')
        self.assertEqual(checkpoint.loss, .25)

    def test_checkpoint_details_include_selected_checkpoint_provenance(self):
        checkpoint = SimpleNamespace(
            id='cp-ABCD', loss=.25,
            run_config={'model': {'version': '0.4'}}, model_id='model-XYZ',
            run_id='run-123', replica=('1351', '4'), fit_mode='single',
            source='/results/checkpoint')

        details = model_inspector.checkpoint_details(checkpoint)

        self.assertIn('Checkpoint ID: cp-ABCD', details)
        self.assertIn('Model ID: model-XYZ', details)
        self.assertIn('Source: /results/checkpoint', details)

    def test_legacy_converted_checkpoint_uses_v042_display_contract(self):
        source = Path('/converted/cp-ABCD-EFGH__liter__1351_4')
        contents = {
            'total_loss': .25,
            'unit_conversion': {
                'target_units': {'M_*': 'micromol/L'},
                'replica': {'sample': '1351', 'replica_number': '4'},
            },
        }
        run_config = {
            'model': {'version': '0.3'},
            'chosen': {'sample': '0', 'validation_replica': '0', 'fit_mode': 'multi',
                       'pathways': ['old'], 'normalized_parameters': True},
            'objective': {},
        }
        with patch.object(model_inspector, 'resolve_checkpoint', return_value=source), \
             patch.object(model_inspector, 'replay_inputs',
                          return_value=(source, run_config, 'parameters', contents)):
            checkpoint = model_inspector.load_checkpoint(str(source))

        self.assertEqual(model_inspector.replay_model_settings(checkpoint), (
            '0.4.2', ('Hydrolysis', 'Fermentation', 'Hydro', 'Aceto', 'Homo', 'Fe3'), False))
        self.assertEqual(checkpoint.run_config['chosen']['fit_mode'], 'single')
        self.assertEqual(checkpoint.replica, ('1351', '4'))

    def test_manual_converted_result_does_not_require_an_objective_call(self):
        parameters = object()
        pathway_model = SimpleNamespace(parameters=lambda: parameters)
        checkpoint = SimpleNamespace(
            loss=.25, run_config={'model': {}}, source='/converted/cp-ABCD',
            metadata={'converted_checkpoint': {'format': 'litre-normalized-v1'}})

        result = model_inspector.manual_converted_result(
            checkpoint, (pathway_model, object(), object()))

        self.assertEqual(result.comparison_name, 'manual values')
        self.assertEqual(result.checkpoint.id, 'manual')
        self.assertEqual(result.checkpoint.loss, .25)


if __name__ == '__main__':
    unittest.main()
