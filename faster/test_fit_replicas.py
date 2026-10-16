import unittest
from contextlib import redirect_stdout
import io

import fit_replicas
import parameters


class Replica:
    def __init__(self, sample, number):
        self.sample = sample
        self.replica_number = str(number)


class Sample:
    def __init__(self, name, site, origin, depth, replica_numbers):
        self.sample_name = str(name)
        self.site = site
        self.origin = origin
        self.depth = depth
        self.replicas = [Replica(self, number) for number in replica_numbers]


class Dataset:
    def __init__(self, samples):
        self.samples = samples


class SelectionTest(unittest.TestCase):
    def setUp(self):
        self.dataset = Dataset([
            Sample('1351', 'Samoylov', 'core', 1.5, (4, 5)),
            Sample('1352', 'Kurugnakh', 'cliff', 4.0, (4,)),
        ])

    def test_filters_site_origin_depth_and_full_replica_identifier(self):
        selected = fit_replicas.select_replicas(
            self.dataset, sites=['samoylov'], origins=['core'], depth=(1, 2),
            replicas=['13515'])

        self.assertEqual([fit_replicas.replica_identifier(replica)
                          for replica in selected], ['13515'])

    def test_split_selects_one_replica(self):
        selected = fit_replicas.select_replicas(self.dataset, split=('1352', '4'))

        self.assertEqual([fit_replicas.replica_identifier(replica)
                          for replica in selected], ['13524'])

    def test_random_selection_is_reproducible(self):
        first = fit_replicas.select_replicas(self.dataset, random_count=2, seed=7)
        second = fit_replicas.select_replicas(self.dataset, random_count=2, seed=7)

        self.assertEqual([fit_replicas.replica_identifier(replica) for replica in first],
                         [fit_replicas.replica_identifier(replica) for replica in second])


class SharedParameterParsingTest(unittest.TestCase):
    def test_expands_globs_without_duplicates(self):
        self.assertEqual(
            fit_replicas.expand_parameter_names(
                ['Hydro_*', 'Homo_v_max', 'Hydro_v_max'],
                ['Hydro_v_max', 'Hydro_CUE', 'Homo_v_max']),
            ['Hydro_v_max', 'Hydro_CUE', 'Homo_v_max'])

    def test_rejects_unknown_name(self):
        with self.assertRaisesRegex(ValueError, 'matches no model parameter'):
            fit_replicas.expand_parameter_names(['missing'], ['Hydro_v_max'])


class SetupSummaryTest(unittest.TestCase):
    def test_declared_default_ranges_remain_locally_optimizable(self):
        template = fit_replicas.build_parameter_template(fit_replicas.PATHWAYS)

        self.assertTrue(template.parameters()['M_Ferm'].is_variable())
        self.assertTrue(template.parameters()['Fe3'].is_variable())

    def test_reports_parameter_roles_and_nested_objective(self):
        sample = Sample('1351', 'Samoylov', 'core', 1.5, (4,))
        local = parameters.ModelParameters({
            'shared_rate': parameters.Parameter('shared_rate', 0.5, (0.1, 1.0)),
            'local_rate': parameters.Parameter('local_rate', 0.5, (0.1, 1.0)),
            'fixed_rate': parameters.Parameter('fixed_rate', 2.0),
            'TOC': parameters.Parameter('TOC', 0.0),
        })
        shared = parameters.ModelParameters({
            'shared_rate': parameters.Parameter('shared_rate', 0.5, (0.2, 0.8)),
        })
        output = io.StringIO()
        with redirect_stdout(output):
            fit_replicas.print_setup_summary(
                sample.replicas, ['Hydro'], ['Homo'], ['shared_rate'], shared,
                {'13514': local}, {'fixed_rate': 2.0}, {}, [], None)

        text = output.getvalue()
        self.assertIn('13514', text)
        self.assertIn('Omitted pathways: Homo', text)
        self.assertIn('global (outer DE)', text)
        self.assertIn('local (inner DE)', text)
        self.assertIn('constant', text)
        self.assertIn('OBJECTIVE STRUCTURE', text)
        self.assertIn('L(g) = Σ_r', text)


class LocalR2ThresholdTest(unittest.TestCase):
    def test_stops_only_when_the_threshold_is_reached(self):
        callback = fit_replicas.StopWhenLocalR2Reached(0.7)
        callback.set_objective(type('Objective', (), {'last_R2': lambda self: 0.7})())

        with self.assertRaises(fit_replicas.LocalR2ThresholdReached):
            callback()

    def test_parser_accepts_local_r2_threshold(self):
        args = fit_replicas.parse_args([
            '--shared', 'Hydro_v_max', '--local-r2-threshold', '0.7',
            '--local-workers', '-1'])

        self.assertEqual(args.local_r2_threshold, 0.7)
        self.assertEqual(args.local_workers, -1)


if __name__ == '__main__':
    unittest.main()
