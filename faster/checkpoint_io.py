"""Small, shared helpers for finding and replaying saved checkpoints."""

from __future__ import annotations

import copy
import re
from pathlib import Path

import parameters
from USER_VARIABLES import RESULTS_DIRECTORY


# A normal checkpoint ends with its ID; CMA restart minima append their
# restart label after it.  Keep the ID capture fixed-width so the label cannot
# be consumed as part of an otherwise hyphen-tolerant identifier.
CHECKPOINT_ID_PATTERN = re.compile(
    r'(?:^|[_-])cp-([a-z0-9]{4}-[a-z0-9]{4})(?=$|\.json$|_restart-)', re.IGNORECASE)


def checkpoint_id(path):
    """Return an ID from a normal or CMA-restart checkpoint filename."""
    match = CHECKPOINT_ID_PATTERN.search(Path(path).name)
    return None if match is None else match.group(1)


def resolve_checkpoint(identifier, results_directory=RESULTS_DIRECTORY, *, additional_roots=()):
    """Resolve a checkpoint path or a unique full/partial checkpoint ID."""
    direct_path = Path(identifier)
    if direct_path.is_file():
        return direct_path

    roots = [Path(root) for root in (results_directory, *additional_roots) if Path(root).is_dir()]
    if not roots:
        raise ValueError(f'Checkpoint directory does not exist: {results_directory}')
    query = str(identifier).casefold().removeprefix('cp-')
    if not query:
        raise ValueError('A checkpoint ID must contain at least one character.')
    matches = [
        (found_id, path)
        for root in roots
        for path in root.rglob('*')
        for found_id in [checkpoint_id(path)]
        if found_id is not None and query in found_id.casefold() and path.is_file()
    ]
    if not matches:
        raise ValueError(f'No checkpoint ID matching {identifier!r} in '
                         + ', '.join(map(str, roots)))
    if len(matches) > 1:
        suggestions = '\n'.join(
            f'  {found_id}: {path}'
            for found_id, path in matches[:12])
        remaining = '' if len(matches) <= 12 else f'\n  ... and {len(matches) - 12} more'
        raise ValueError(
            f'Checkpoint ID {identifier!r} is ambiguous; use more characters:\n'
            f'{suggestions}{remaining}')
    return matches[0][1]


def load_checkpoint_document(path):
    """Load an explicitly selected annotated checkpoint."""
    try:
        contents = parameters.load_file(path)
    except (OSError, ValueError) as error:
        raise ValueError(f'Could not read checkpoint {path}: {error}') from error
    required = {'parameters', 'total_loss', 'run_config'}
    present = set(contents) if isinstance(contents, dict) else set()
    if not required.issubset(present):
        raise ValueError(
            f'{path} is not an annotated checkpoint; missing '
            + ', '.join(sorted(required - present)))
    return contents


def replay_inputs(checkpoint_path):
    """Load the saved run configuration and fitted parameter values."""
    source = resolve_checkpoint(checkpoint_path)
    contents = load_checkpoint_document(source)
    return source, copy.deepcopy(contents['run_config']), parameters.ModelParameters(
        contents['parameters']), contents
