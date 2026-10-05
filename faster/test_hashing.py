import copy
import json
from pathlib import Path
import tempfile
import unittest

import hashing
import numpy as np
import reindex_checkpoints


def checkpoint_with(value_overrides=None):
    checkpoint = {
        'parameters': {
            'Hydro_thermodynamics': True,
            'rate': 1.0,
        },
        'total_loss': 0.25,
        'run_config': {
            'model': {
                'Hydro': {
                    'microbe': {
                        'name': 'M_Hydro', 'death_rate': 0.0,
                        'v_max': 'variable', 'Kmb': 0.0, 'CUE': 'variable',
                        'C_source': 'CO2',
                    },
                    'educts': {'CO2': {
                        'name': 'CO2', 'stoichiometry': 1.0,
                        'Km': 'variable', 'inhibition': float('inf'),
                    }},
                    'products': {},
                    'use_thermodynamics': True,
                },
                'version': '0.1',
            },
            'chosen': {
                'sample': '1351', 'validation_replica': '4',
                't_start': None, 't_end': None, 'fit_mode': 'split',
                'pathways': ['Hydro'], 'parameter_override': {
                    'Hydro_thermodynamics': True,
                },
                'normalized_parameters': True,
                'algorithm': 'differential_evolution',
            },
            'objective': {
                'loss_weight': {'CO2': 1.0},
                'reduction': {'CO2': 'mse'},
                'transform': {'CO2': ['normalize']},
            },
            'algo': {
                'strategy': 'rand1bin', 'updating': 'deferred', 'popsize': 10,
                'workers': -1, 'tol': 0.0001, 'init': 'sobol', 'polish': False,
                'recombination': 0.7, 'mutation': [0.5, 1.0],
            },
            'range': {
                'rate': {
                    'name': 'rate', 'range': [0.0, 1.0], 'scale': 'linear',
                    'normalize': False,
                },
            },
        },
    }
    for path, value in (value_overrides or {}).items():
        target = checkpoint
        *parents, key = path.split('.')
        for parent in parents:
            target = target[parent]
        target[key] = value
    return checkpoint


class HashingTest(unittest.TestCase):
    def test_equivalent_scalar_representations_have_the_same_checkpoint_id(self):
        canonical = checkpoint_with()
        alternate = checkpoint_with({
            'parameters.Hydro_thermodynamics': '1',
            'parameters.rate': '1.0',
            'total_loss': '0.25',
            'run_config.model.Hydro.use_thermodynamics': 'True',
            'run_config.model.Hydro.microbe.death_rate': '0',
            'run_config.chosen.sample': 1351,
            'run_config.chosen.validation_replica': 4.0,
            'run_config.chosen.normalized_parameters': 1,
            'run_config.chosen.parameter_override.Hydro_thermodynamics': 'true',
            'run_config.objective.loss_weight.CO2': '1',
            'run_config.algo.popsize': '10.0',
            'run_config.algo.polish': 'false',
            'run_config.range.rate.range': ['0', '1'],
            'run_config.range.rate.normalize': '0',
        })

        self.assertEqual(hashing.build_checkpoint_id(canonical),
                         hashing.build_checkpoint_id(alternate))
        self.assertEqual(hashing.build_run_id(canonical['run_config']),
                         hashing.build_run_id(alternate['run_config']))

    def test_boolean_conversion_is_explicit(self):
        self.assertTrue(hashing.as_bool('TRUE'))
        self.assertFalse(hashing.as_bool('0'))
        self.assertTrue(hashing.as_bool(np.bool_(True)))
        with self.assertRaises(ValueError):
            hashing.as_bool('not a bool')

    def test_model_hashing_does_not_mutate_the_caller_config(self):
        model = checkpoint_with()['run_config']['model']
        original = copy.deepcopy(model)
        hashing.build_model_id(model)
        self.assertEqual(model, original)

    def test_reindex_archives_original_and_writes_canonical_checkpoint(self):
        alternate = checkpoint_with({
            'parameters.Hydro_thermodynamics': 'false',
            'run_config.model.Hydro.use_thermodynamics': '0',
            'run_config.chosen.normalized_parameters': 'true',
        })
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / 'results'
            original = source / 'old-folder' / 'checkpoint.json'
            original.parent.mkdir(parents=True)
            with original.open('w') as handle:
                json.dump(alternate, handle)

            archive = root / 'archive'
            migrations, skipped, errors = reindex_checkpoints.plan_migration(source, archive)
            self.assertEqual(skipped, 0)
            self.assertEqual(errors, [])
            self.assertEqual(len(migrations), 1)

            removed_directories = reindex_checkpoints.execute_migration(migrations, source)
            self.assertTrue(archive.joinpath('old-folder', 'checkpoint.json').is_file())
            self.assertEqual(removed_directories, 1)
            self.assertFalse(source.joinpath('old-folder').exists())
            with migrations[0].destination.open() as handle:
                written = json.load(handle)
            self.assertIs(written['parameters']['Hydro_thermodynamics'], False)
            self.assertIs(written['run_config']['chosen']['normalized_parameters'], True)


if __name__ == '__main__':
    unittest.main()
