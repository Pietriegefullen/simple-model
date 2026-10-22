"""Compatibility command for plotting the best saved checkpoints.

Use ``plot_checkpoint.py --best`` for the same interface directly.
"""

from __future__ import annotations

import sys

import plot_checkpoint


def main(argv=None):
    """Delegate all argument handling and rendering to ``plot_checkpoint``."""
    arguments = sys.argv[1:] if argv is None else list(argv)
    return plot_checkpoint.main(['--best', *arguments])


if __name__ == '__main__':
    main()
