import unittest

import matplotlib.pyplot as plt

import system
import variant_diagnostics


def pathway(name, coefficients):
    instance = type(name, (), {})()
    instance.pathway_vector = [0.0] * len(system.SYSTEM)
    for pool, value in coefficients.items():
        instance.pathway_vector[system.index(pool)] = value
    return instance


def run_log(include_thermodynamic_factor=True):
    log = {
        'H2': ([0, 1, 2], [10, 8, 6]),
        'Acetate': ([0, 1, 2], [8, 7, 6]),
        'CO2': ([0, 1, 2], [5, 6, 7]),
        'CH4': ([0, 1, 2], [1, 2, 3]),
    }
    for pathway_name, microbe in {
            'Homo': 'M_Homo', 'Hydro': 'M_Hydro',
            'Fe3': 'M_Fe3', 'Aceto': 'M_Ac'}.items():
        log[microbe] = ([0, 1, 2], [1, 1.5, 2])
        log[f'{pathway_name}_v'] = ([0, 1, 2], [.2, .3, .4])
        log[f'{pathway_name}_MM'] = ([0, 1, 2], [.2, .5, .8])
        log[f'{pathway_name}_deltaG_r'] = ([0, 1, 2], [-20, -21, -22])
        if include_thermodynamic_factor:
            log[f'{pathway_name}_thermodynamic_factor'] = ([0, 1, 2], [.4, .6, .8])
    return log


def model(include_hydro=True):
    pathways = [
        pathway('Homo', {'H2': -4, 'CO2': -2, 'Acetate': 1}),
        pathway('Fe3', {'Acetate': -1.5, 'CO2': 2, 'Fe3': -8}),
        pathway('Aceto', {'Acetate': -1.25, 'CO2': 1, 'CH4': 1}),
    ]
    if include_hydro:
        pathways.append(pathway('Hydro', {'H2': -4, 'CO2': -1, 'CH4': 1}))
    return type('Model', (), {'contributing_pathways': pathways})()


class VariantDiagnosticsTest(unittest.TestCase):
    def tearDown(self):
        plt.close('all')

    def test_pathway_pool_rate_includes_biomass_and_stoichiometry(self):
        time, rate = variant_diagnostics.pathway_pool_rate(
            model(), run_log(), 'Homo', 'H2')
        self.assertEqual(time.tolist(), [0.0, 1.0, 2.0])
        self.assertTrue((abs(rate - [-.8, -1.8, -3.2]) < 1e-12).all())

    def test_complete_abc_suite_has_all_six_mechanistic_figures(self):
        specification = lambda variant: type('Specification', (), {'variant': variant})()
        specifications = [specification(variant) for variant in ('A', 'B', 'C')]
        model_runs = {
            'A': (model(), None, run_log()),
            'B': (model(), None, run_log(include_thermodynamic_factor=False)),
            'C': (model(include_hydro=False), None, run_log()),
        }

        figures = dict(variant_diagnostics.diagnostic_figures(
            specifications, model_runs, 'test'))

        self.assertEqual(set(figures), {
            '05_flux_partition',
            '06_ab_thermodynamic_effect',
            '07_ac_pathway_removal_response',
            '08_substrate_flux_phase',
            '09_mass_balance_attribution',
            '10_ab_limitation_fingerprint',
        })
        self.assertEqual(len(figures['07_ac_pathway_removal_response'].axes), 4)
        self.assertEqual(len(figures['09_mass_balance_attribution'].axes), 4)

    def test_flux_partition_uses_pathway_color_and_variant_line_style(self):
        specification = lambda variant: type('Specification', (), {'variant': variant})()
        specifications = [specification(variant) for variant in ('A', 'B')]
        model_runs = {
            'A': (model(), None, run_log()),
            'B': (model(), None, run_log()),
        }

        figure = variant_diagnostics.plot_flux_partition(
            specifications, model_runs, 'test')
        homo_a, hydro_a, homo_b = figure.axes[0].lines[:3]

        self.assertEqual(homo_a.get_color(), 'tab:pink')
        self.assertEqual(hydro_a.get_color(), 'tab:cyan')
        self.assertEqual(homo_a.get_linestyle(), '-')
        self.assertEqual(homo_b.get_linestyle(), '--')


if __name__ == '__main__':
    unittest.main()
