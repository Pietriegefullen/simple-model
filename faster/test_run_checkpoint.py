"""Tests for the shared checkpoint input helpers.

The former ``run_checkpoint`` command is deliberately only a compatibility
entry point; checkpoint input now belongs to ``checkpoint_io``.
"""

import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

import checkpoint_io


class CheckpointInputTest(unittest.TestCase):
    def test_resolves_a_unique_partial_checkpoint_id(self):
        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            checkpoint = root / 'fit_loss-1_cp-ABCD-EFGH'
            checkpoint.touch()

            resolved = checkpoint_io.resolve_checkpoint('abcd', root)

        self.assertEqual(resolved, checkpoint)

    def test_rejects_an_ambiguous_partial_checkpoint_id(self):
        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / 'fit_loss-1_cp-ABCD-EFGH').touch()
            (root / 'fit_loss-2_cp-ABCD-IJKL').touch()

            with self.assertRaisesRegex(ValueError, 'ambiguous'):
                checkpoint_io.resolve_checkpoint('abcd', root)

    def test_extracts_the_checkpoint_id(self):
        self.assertEqual(checkpoint_io.checkpoint_id(
            'fit_loss-1_cp-ABCD-EFGH'), 'ABCD-EFGH')


if __name__ == '__main__':
    unittest.main()
