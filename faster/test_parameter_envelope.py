import unittest

import parameter_envelope as envelope
import parameter_bounds
import plot_envelope


class Parameters:
    def __init__(self, values):
        self.values = values

    def get_values(self):
        return self.values


class Checkpoint:
    def __init__(self, identifier, loss, run_id='run-A', restart=None,
                 sample='1351', replica='4', values=None, fit_mode='single'):
        self.id = identifier
        self.loss = loss
        self.run_id = run_id
        self.model_id = 'model-ABCD'
        self.fit_mode = fit_mode
        self.replica = (sample, replica)
        self.metadata = {} if restart is None else {'cma_restart': {'index': restart}}
        self.parameters = Parameters(values or {'rate': loss})
        self.source = f'/checkpoints/{identifier}'
        self.multi_replica_results = ()


class SelectionTest(unittest.TestCase):
    def test_combines_metadata_and_partial_id_exclusion(self):
        wanted = Checkpoint('cp-KEEP', 2.)
        excluded = Checkpoint('cp-DROP', 1.)
        selected = envelope.select_checkpoints(
            [wanted, excluded], model_ids=['abcd'], samples=['1351'], replicas=['4'],
            sites=['sam'], origins=['core'], excluded_ids=['drop'],
            metadata={('1351', '4'): {'site': 'Samoylov', 'origin': 'core'}})
        self.assertEqual(selected, [wanted])

    def test_fit_mode_is_combined_with_other_filters(self):
        single = Checkpoint('single', 1., fit_mode='single')
        split = Checkpoint('split', 1., fit_mode='split')

        selected = envelope.select_checkpoints(
            [single, split], model_ids=['abcd'], fit_modes=['split'])

        self.assertEqual(selected, [split])

    def test_sample_glob_matches_checkpoint_sample(self):
        selected = envelope.select_checkpoints(
            [Checkpoint('sample', 1., sample='1373')], samples=['137*'])

        self.assertEqual([checkpoint.id for checkpoint in selected], ['sample'])

    def test_loss_upper_discards_penalty_checkpoints(self):
        selected = envelope.select_checkpoints(
            [Checkpoint('kept', 1.), Checkpoint('penalty', 1e100)], loss_upper=1e99)

        self.assertEqual([checkpoint.id for checkpoint in selected], ['kept'])

    def test_best_restart_keeps_run_best_and_each_restart_best(self):
        retained = envelope.best_with_restarts([
            Checkpoint('ordinary', 5.), Checkpoint('run-best', 1.),
            Checkpoint('restart-0-old', 4., restart=0),
            Checkpoint('restart-0-best', 2., restart=0),
            Checkpoint('restart-1-best', 3., restart=1),
        ])
        self.assertEqual({checkpoint.id for checkpoint in retained}, {
            'run-best', 'restart-0-best', 'restart-1-best'})


class EnvelopeRowsTest(unittest.TestCase):
    def test_records_bounds_and_their_checkpoint_ids(self):
        rows = envelope.envelope_rows([
            Checkpoint('high', 1., values={'rate': 4., 'ignored': float('nan')}),
            Checkpoint('low', 2., values={'rate': 2., 'flag': True}),
        ])
        self.assertEqual(rows, [{
            'parameter': 'rate', 'low': 2., 'high': 4., 'value_count': 2,
            'low_checkpoint_id': 'low', 'high_checkpoint_id': 'high',
        }])

    def test_checkpoint_records_include_sample_replica_numbers(self):
        records = envelope.checkpoint_records([
            Checkpoint('first', 1., sample='1351', replica='4'),
            Checkpoint('second', 2., sample='1352', replica='6'),
        ])

        self.assertEqual(
            [record['sample_replicas'] for record in records],
            [['1351-4'], ['1352-6']])


class BoundsStoreTest(unittest.TestCase):
    def test_resolves_partial_id_and_applies_numeric_bounds(self):
        from pathlib import Path
        from tempfile import TemporaryDirectory
        import json
        import parameters

        with TemporaryDirectory() as directory:
            path = Path(directory) / 'bounds-ABCD-EF12.json'
            path.write_text(json.dumps({'bounds_id': 'ABCD-EF12', 'parameters': [{
                'parameter': 'rate', 'low': 2., 'high': 8.,
            }]}))
            resolved, document = parameter_bounds.load_bounds('abcd', directory)
            model = parameters.ModelParameters({
                'rate': {'name': 'rate', 'value': 1., 'range': (0.1, 10.),
                         'scale': 'log', 'normalize': True},
            })
            parameter_bounds.apply_bounds(model, document)

        self.assertEqual(resolved, path)
        self.assertEqual((model['rate'].low, model['rate'].high), (2., 8.))
        self.assertEqual(model['rate'].value, 4.)


class EnvelopePlotTest(unittest.TestCase):
    def test_parser_requires_exactly_one_parameter_catalogue(self):
        with self.assertRaises(SystemExit):
            plot_envelope.parse_args(['ABCD'])
        with self.assertRaises(SystemExit):
            plot_envelope.parse_args(['ABCD', '--default', '--hypothesis'])

        args = plot_envelope.parse_args(['ABCD', '--hypothesis'])

        self.assertTrue(args.hypothesis)

    def test_log_scale_requires_positive_log_reference_and_values(self):
        self.assertTrue(plot_envelope._is_logarithmic(
            {'low': 2., 'high': 8.},
            {'low': 1., 'high': 10., 'value': 2., 'scale': 'log'}))
        self.assertFalse(plot_envelope._is_logarithmic(
            {'low': 0., 'high': 8.},
            {'low': 1., 'high': 10., 'value': 2., 'scale': 'log'}))

    def test_loads_marker_values_from_checkpoint_sources(self):
        from pathlib import Path
        from tempfile import TemporaryDirectory
        import json

        with TemporaryDirectory() as directory:
            source = Path(directory) / 'checkpoint.json'
            source.write_text(json.dumps({'parameters': {'rate': 3.}, 'total_loss': 0.0}))
            values = plot_envelope.checkpoint_values({'checkpoints': [{
                'checkpoint_id': 'cp-TEST', 'source': str(source),
            }]})

        self.assertEqual(values, {'cp-TEST': {
            'parameters': {'rate': 3.}, 'loss': 0.0,
        }})

    def test_missing_checkpoint_is_skipped_with_a_warning(self):
        with self.assertWarnsRegex(UserWarning, 'no longer available'):
            values = plot_envelope.checkpoint_values({'checkpoints': [{
                'checkpoint_id': 'cp-MISSING', 'source': '/not/a/checkpoint',
            }]})

        self.assertEqual(values, {})

    def test_loss_normalizer_accepts_a_shared_loss(self):
        normalizer = plot_envelope.loss_normalizer({
            'cp-A': {'parameters': {}, 'loss': .2},
            'cp-B': {'parameters': {}, 'loss': .2},
        })

        self.assertLess(normalizer.vmin, normalizer.vmax)


if __name__ == '__main__':
    unittest.main()
