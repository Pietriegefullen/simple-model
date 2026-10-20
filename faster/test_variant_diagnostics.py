import unittest

import matplotlib.pyplot as plt
import numpy as np

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
            '11_ac_homo_activity_delta_g',
        })
        self.assertEqual(len(figures['07_ac_pathway_removal_response'].axes), 4)
        self.assertEqual(len(figures['09_mass_balance_attribution'].axes), 4)

    def test_ac_homo_activity_and_delta_g_share_a_variant_comparison_plot(self):
        specification = lambda variant: type('Specification', (), {'variant': variant})()
        specifications = [specification(variant) for variant in ('A', 'C')]
        model_runs = {
            'A': (model(), None, run_log()),
            'C': (model(include_hydro=False), None, run_log()),
        }

        figure = variant_diagnostics.plot_ac_homo_activity_and_delta_g(
            specifications, model_runs, 'test')
        activity, delta_g = figure.axes

        self.assertEqual([line.get_linestyle() for line in activity.lines], ['-', ':'])
        self.assertEqual([line.get_linestyle() for line in delta_g.lines], ['-', ':'])
        np.testing.assert_allclose(activity.lines[0].get_ydata(), [.2, .45, .8])
        np.testing.assert_allclose(delta_g.lines[0].get_ydata(),
                                   [-.02, -.021, -.022])
        self.assertIn('M(t)v', activity.get_ylabel())
        self.assertIn(r'\Delta G_\mathrm{r}', delta_g.get_ylabel())

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

    def test_d_is_included_in_the_selected_diagnostics(self):
        specification = lambda variant: type('Specification', (), {'variant': variant})()
        specifications = [specification(variant) for variant in ('A', 'D')]
        model_runs = {
            'A': (model(), None, run_log()),
            'D': (model(), None, run_log(include_thermodynamic_factor=False)),
        }

        figures = dict(variant_diagnostics.diagnostic_figures(
            specifications, model_runs, 'test'))

        thermodynamic_effect = figures['06_ab_thermodynamic_effect']
        self.assertTrue(any(line.get_linestyle() == '-.'
                            for axis in thermodynamic_effect.axes
                            for line in axis.lines))
        counterfactual = figures['07_ac_pathway_removal_response']
        self.assertTrue(any(line.get_linestyle() == '-.'
                            for axis in counterfactual.axes
                            for line in axis.lines))
        fingerprint = figures['10_ab_limitation_fingerprint']
        self.assertEqual(fingerprint.axes[1].get_title(), r'$\mathrm{model\ D}$')

    def test_missing_logged_factor_is_marked_unavailable_in_diagnostics(self):
        specification = lambda variant: type('Specification', (), {'variant': variant})()
        specifications = [specification(variant) for variant in ('A', 'B')]
        model_runs = {
            'A': (model(), None, run_log()),
            'B': (model(), None, run_log(include_thermodynamic_factor=False)),
        }

        self.assertIsNone(variant_diagnostics._thermodynamic_factor(
            model_runs['B'][2], 'Homo'))

        effect = variant_diagnostics.plot_ab_thermodynamic_effect(
            specifications, model_runs, 'test')
        self.assertTrue(any('not logged' in text.get_text()
                            for text in effect.axes[0].texts))

        fingerprint = variant_diagnostics.plot_ab_limitation_fingerprint(
            specifications, model_runs, 'test')
        missing_factor_cells = fingerprint.axes[1].images[0].get_array()[0]
        self.assertTrue(np.ma.getmaskarray(missing_factor_cells).all())
        self.assertIn('not', fingerprint.axes[-1].get_ylabel())

    def test_fingerprint_preserves_a_constant_factor_of_one(self):
        log = run_log()
        for key, (time, _) in tuple(log.items()):
            if key.endswith('_thermodynamic_factor'):
                log[key] = time, [1., 1., 1.]
        specification = lambda variant: type('Specification', (), {'variant': variant})()
        specifications = [specification(variant) for variant in ('A', 'B')]
        model_runs = {
            'A': (model(), None, log),
            'B': (model(), None, log),
        }

        fingerprint = variant_diagnostics.plot_ab_limitation_fingerprint(
            specifications, model_runs, 'test')
        homo_thermodynamic_factor = fingerprint.axes[0].images[0].get_array()[0]
        np.testing.assert_allclose(homo_thermodynamic_factor, 1.)


if __name__ == '__main__':
    unittest.main()
