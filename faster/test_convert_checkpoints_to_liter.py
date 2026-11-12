import unittest
from unittest import mock

import CONSTANTS
import convert_checkpoints_to_liter as conversion


class ConversionTest(unittest.TestCase):
    def test_microbial_and_non_microbial_values_use_correct_factors(self):
        values, metadata = conversion.convert_values(
            {'M_Hydro': CONSTANTS.MOLAR_MASS_BIOMASS, 'Hydro_Km_H2': 2.,
             'Hydrolysis_Kmb': CONSTANTS.MOLAR_MASS_BIOMASS,
             'Hydro_v_max': 2., 'Hydro_CUE': .2, 'death_rate': .01,
             'Hydro_thermodynamics': True},
            dry_weight_g=2., water_content_ml=100.)
        self.assertAlmostEqual(values['M_Hydro'], 20.)
        self.assertAlmostEqual(values['Hydrolysis_Kmb'], 20.)
        self.assertAlmostEqual(values['Hydro_Km_H2'], 40.)
        self.assertAlmostEqual(values['Hydro_v_max'],
                               2. * CONSTANTS.MOLAR_MASS_BIOMASS)
        self.assertEqual(values['Hydro_CUE'], .2)
        self.assertEqual(values['death_rate'], .01)
        self.assertIs(values['Hydro_thermodynamics'], True)
        self.assertEqual(metadata['gdw_per_litre'], 20.)

    def test_joint_checkpoint_is_expanded_and_merges_shared_values(self):
        document = {
            'parameters': {'Hydro_Km_H2': 2.},
            'run_config': {'chosen': {'fit_mode': 'multi'}},
            'replicas': [
                {'sample': '1351', 'replica_number': '4',
                 'parameters': {'M_Hydro': CONSTANTS.MOLAR_MASS_BIOMASS}},
                {'sample': '1351', 'replica_number': '5',
                 'parameters': {'M_Hydro': CONSTANTS.MOLAR_MASS_BIOMASS * 2}},
            ],
        }
        expanded = list(conversion.expanded_solutions(document))
        self.assertEqual(len(expanded), 2)
        self.assertEqual(expanded[0][2], {
            'Hydro_Km_H2': 2., 'M_Hydro': CONSTANTS.MOLAR_MASS_BIOMASS})
        self.assertEqual(expanded[1][2]['M_Hydro'], CONSTANTS.MOLAR_MASS_BIOMASS * 2)

    def test_envelopes_are_compared_to_v04_default_bounds(self):
        rows = conversion.envelope_rows({'Hydro_Km_H2': [13., 200.], 'death_rate': [.1]})
        by_name = {row['parameter']: row for row in rows}
        self.assertEqual(by_name['Hydro_Km_H2']['comparison'], 'inside_default_bounds')
        self.assertEqual(by_name['death_rate']['comparison'], 'outside_default_bounds')

    def test_declared_numeric_ranges_are_converted(self):
        ranges = {'Hydro_Km_H2': {'range': [1., 2.]},
                  'M_Hydro': {'range': [CONSTANTS.MOLAR_MASS_BIOMASS,
                                        CONSTANTS.MOLAR_MASS_BIOMASS * 2]},
                  'Hydrolysis_Kmb': {'range': [CONSTANTS.MOLAR_MASS_BIOMASS,
                                                CONSTANTS.MOLAR_MASS_BIOMASS * 2]},
                  'Hydro_v_max': {'range': [1., 2.]},
                  'Hydro_CUE': {'range': [.1, .2]}}
        converted = conversion.convert_declared_ranges(ranges, 2., 100.)
        self.assertEqual(converted['Hydro_Km_H2']['range'], [20., 40.])
        self.assertEqual(converted['M_Hydro']['range'], [20., 40.])
        self.assertEqual(converted['Hydrolysis_Kmb']['range'], [20., 40.])
        self.assertEqual(converted['Hydro_v_max']['range'], [
            CONSTANTS.MOLAR_MASS_BIOMASS,
            2. * CONSTANTS.MOLAR_MASS_BIOMASS])
        self.assertEqual(converted['Hydro_CUE']['range'], [.1, .2])

    def test_converted_document_declares_its_v042_replay_contract(self):
        document = {
            'parameters': {'death_rate': .01},
            'run_config': {
                'model': {'version': '0.3'},
                'chosen': {'fit_mode': 'single', 'sample': '1351',
                           'validation_replica': '4'},
            },
        }
        replica = mock.Mock(dry_weight=2., water_content=100.)

        converted, _, _ = next(conversion.converted_documents(
            document, '/results/cp-OLD', {('1351', '4'): replica}))

        self.assertEqual(converted['run_config']['model']['version'], '0.3')
        self.assertEqual(converted['converted_checkpoint'], {
            'format': 'litre-normalized-v1',
            'source_model_version': '0.3',
            'replay_model_version': '0.4.2',
            'pathways': ['Hydrolysis', 'Fermentation', 'Hydro', 'Aceto', 'Homo', 'Fe3'],
            'normalized_parameters': False,
        })

    def test_only_lowest_loss_checkpoint_is_kept_per_run(self):
        documents = [
            ('z', {'run_config': {'run': 'one'}, 'total_loss': 2.}, 'model-X'),
            ('a', {'run_config': {'run': 'one'}, 'total_loss': 1.}, 'model-X'),
            ('b', {'run_config': {'run': 'two'}, 'total_loss': .5}, 'model-X'),
        ]
        with mock.patch('hashing.build_run_id', side_effect=lambda config: config['run']):
            selected = conversion.best_checkpoint_per_run(documents)
        self.assertEqual([str(path) for path, _, _ in selected], ['a', 'b'])


if __name__ == '__main__':
    unittest.main()
