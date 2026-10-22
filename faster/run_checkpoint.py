"""Deprecated compatibility entry point for ``plot_checkpoint.py``.

Checkpoint discovery and replay live in :mod:`checkpoint_io`; plotting is
handled by :mod:`plot_checkpoint`.  New commands should call
``python plot_checkpoint.py CHECKPOINT`` directly.
"""

from checkpoint_io import checkpoint_id, resolve_checkpoint
from plot_checkpoint import main


if __name__ == '__main__':
    main()
