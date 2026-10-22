"""Small, shared helpers for finding and replaying saved checkpoints."""

from __future__ import annotations

import copy
import re
from pathlib import Path

import parameters
from fit_sample import _requested_checkpoint
from USER_VARIABLES import RESULTS_DIRECTORY


CHECKPOINT_ID_PATTERN = re.compile(r'(?:^|[_-])cp-([a-z0-9-]+)$', re.IGNORECASE)


def checkpoint_id(path):
    """Return the short checkpoint ID embedded in a result filename."""
    match = CHECKPOINT_ID_PATTERN.search(Path(path).name)
    return None if match is None else match.group(1)


def resolve_checkpoint(identifier, results_directory=RESULTS_DIRECTORY):
    """Resolve a checkpoint path or a unique full/partial checkpoint ID."""
    direct_path = Path(identifier)
    if direct_path.is_file():
        return direct_path

    root = Path(results_directory)
    if not root.is_dir():
        raise ValueError(f'Results directory does not exist: {root}')
    query = str(identifier).casefold().removeprefix('cp-')
    matches = [
        (found_id, path)
        for path in root.rglob('*') if path.is_file()
        for found_id in [checkpoint_id(path)]
        if found_id is not None and query in found_id.casefold()
    ]
    if not matches:
        raise ValueError(f'No checkpoint ID matching {identifier!r} in {root}')
    if len(matches) > 1:
        suggestions = '\n'.join(
            f'  {found_id}: {path.relative_to(root)}'
            for found_id, path in matches[:12])
        remaining = '' if len(matches) <= 12 else f'\n  ... and {len(matches) - 12} more'
        raise ValueError(
            f'Checkpoint ID {identifier!r} is ambiguous; use more characters:\n'
            f'{suggestions}{remaining}')
    return matches[0][1]


def replay_inputs(checkpoint_path):
    """Load the saved run configuration and fitted parameter values."""
    source = resolve_checkpoint(checkpoint_path)
    contents = _requested_checkpoint(source)
    return source, copy.deepcopy(contents['run_config']), parameters.ModelParameters(
        contents['parameters']), contents
