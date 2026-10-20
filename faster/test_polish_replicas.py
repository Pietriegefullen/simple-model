import unittest

import polish_replicas


def multi_checkpoint():
    return {
        'parameters': {'shared_rate': 0.4},
        'run_config': {
            'chosen': {
                'fit_mode': 'multi',
                'pathways': ['Hydro'],
                'normalized_parameters': True,
                'parameter_override': {'fixed_rate': 2.0},
            },
            'objective': {},
            'multi_replica': {'shared': ['shared_rate']},
        },
        'replicas': [
            {'replica': '13514', 'sample': '1351', 'replica_number': '4',
             'parameters': {'local_rate': 0.25}},
            {'replica': '13515', 'sample': '1351', 'replica_number': '5',
             'parameters': {'local_rate': 0.75}},
        ],
    }


class PolishSelectionTest(unittest.TestCase):
    def test_defaults_to_every_saved_replica(self):
        selected = polish_replicas.select_replica_records(multi_checkpoint()['replicas'])
        self.assertEqual([record['replica'] for record in selected], ['13514', '13515'])

    def test_selects_explicit_replicas_in_requested_order(self):
        selected = polish_replicas.select_replica_records(
            multi_checkpoint()['replicas'], ['1351/5', '1351:4', '13515'])
        self.assertEqual([record['replica'] for record in selected], ['13515', '13514'])

    def test_rejects_replica_not_in_source_fit(self):
        with self.assertRaisesRegex(ValueError, 'not present'):
            polish_replicas.select_replica_records(multi_checkpoint()['replicas'], ['13664'])

    def test_builds_single_fit_that_freezes_source_shared_values(self):
        chosen = polish_replicas.polish_run_configuration(
            multi_checkpoint(), multi_checkpoint()['replicas'][0])
        self.assertEqual(chosen['fit_mode'], 'single')
        self.assertEqual(chosen['sample'], 1351)
        self.assertEqual(chosen['validation_replica'], 4)
        self.assertEqual(chosen['algorithm'], 'cma_es')
        self.assertEqual(chosen['parameter_override'],
                         {'fixed_rate': 2.0, 'shared_rate': 0.4})

    def test_parser_has_a_small_polishing_sigma_by_default(self):
        args = polish_replicas.parse_args(['ABCD'])
        self.assertEqual(args.sigma, 0.10)
        self.assertEqual(args.max_evaluations, 3000)


if __name__ == '__main__':
    unittest.main()
