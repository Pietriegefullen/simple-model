"""Persist and apply reusable fitted-parameter bounds."""

from __future__ import annotations

import json
import math
import re
import uuid
from pathlib import Path

from USER_VARIABLES import PROJECT_DIRECTORY


BOUNDS_DIRECTORY = Path(PROJECT_DIRECTORY) / 'bounds'
_ID_PATTERN = re.compile(r'^bounds-([a-z0-9]{4}-[a-z0-9]{4})\.json$', re.IGNORECASE)


def new_bounds_id():
    token = uuid.uuid4().hex[:8].upper()
    return f'{token[:4]}-{token[4:]}'


def bounds_path(bounds_id, directory=BOUNDS_DIRECTORY):
    return Path(directory) / f'bounds-{bounds_id}.json'


def resolve_bounds(identifier, directory=BOUNDS_DIRECTORY):
    """Resolve a full or unique partial bounds ID to its stored JSON file."""
    directory = Path(directory)
    query = str(identifier).casefold().removeprefix('bounds-')
    matches = [path for path in directory.glob('bounds-*.json')
               if (match := _ID_PATTERN.match(path.name))
               and query in match.group(1).casefold()]
    if not matches:
        raise ValueError(f'No bounds ID matching {identifier!r} in {directory}')
    if len(matches) > 1:
        choices = ', '.join(path.stem.removeprefix('bounds-') for path in matches[:12])
        raise ValueError(f'Bounds ID {identifier!r} is ambiguous; use more characters: {choices}')
    return matches[0]


def load_bounds(identifier, directory=BOUNDS_DIRECTORY):
    path = resolve_bounds(identifier, directory)
    with path.open() as handle:
        document = json.load(handle)
    if not isinstance(document, dict) or not isinstance(document.get('parameters'), list):
        raise ValueError(f'{path} is not a parameter-bounds file.')
    return path, document


def apply_bounds(model_parameters, document):
    """Narrow compatible model parameters to the bounds in an envelope file."""
    for row in document['parameters']:
        name, low, high = row.get('parameter'), row.get('low'), row.get('high')
        if name not in model_parameters:
            raise ValueError(f'Bounds contain unknown model parameter {name!r}.')
        if (isinstance(low, bool) or isinstance(high, bool)
                or not isinstance(low, (int, float)) or not isinstance(high, (int, float))
                or not math.isfinite(low) or not math.isfinite(high) or low > high):
            raise ValueError(f'Bounds for {name!r} must be finite numeric low/high values.')
        parameter = model_parameters[name]
        if parameter.options is not None:
            raise ValueError(f'Bounds cannot be applied to categorical parameter {name!r}.')
        if low == high:
            parameter.constant(low)
        else:
            if parameter.scale == 'log':
                if low <= 0:
                    raise ValueError(f'Log-scaled parameter {name!r} requires positive bounds.')
                midpoint = math.sqrt(low * high)
            else:
                midpoint = (low + high) / 2
            parameter.variable(midpoint, (low, high))
            parameter.transformer = None
    return model_parameters
