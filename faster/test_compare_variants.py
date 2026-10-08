import unittest

import compare_variants


class Checkpoint:
    def __init__(self, model_id, loss):
        self.model_id = model_id
        self.loss = loss
        self.loss_id = 'loss-test'
        self.fit_mode = 'split'
        self.replica = ('1351', '4')


class ComparisonGroupsTest(unittest.TestCase):
    def test_groups_keep_the_lowest_loss_checkpoint_per_variant(self):
        first_a = Checkpoint('model-FA75-first', 2.0)
        best_a = Checkpoint('model-PIRL-best', 1.0)
        b = Checkpoint('model-EM33-b', 3.0)

        groups = compare_variants.comparison_groups((first_a, best_a, b))

        self.assertEqual(len(groups), 1)
        group = ('loss-test', 'split', ('1351', '4'))
        self.assertIs(groups[group]['A'], best_a)
        self.assertIs(groups[group]['B'], b)

    def test_target_is_separate_by_loss_mode_and_replica(self):
        group = ('loss-test', 'split', ('1351', '4'))
        target = compare_variants.comparison_target('/tmp/comparisons', group)
        self.assertEqual(
            str(target), '/tmp/comparisons/loss-test_split/1351-4')


if __name__ == '__main__':
    unittest.main()
