"""Plot a stored parameter-bounds envelope against a selected model catalogue.

Run from the repository root, using a full or unique partial bounds ID::

    python faster/plot_envelope.py ABCD --hypothesis
"""

from __future__ import annotations

import argparse
import math
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize
import numpy as np

import checkpoint_io
from model_factory import (build_model, current_version,
                           envelope_converted_model_parameters,
                           hypothesis_model_parameters)
import parameter_bounds
import parameters
from USER_VARIABLES import PLOT_DIRECTORY


def converted_envelope_bounds():
    """Return the current converted-envelope range and scale for each parameter."""
    return {
        parameter.name: {
            'low': parameter.low,
            'high': parameter.high,
            'value': parameter.value,
            'scale': parameter.scale,
        }
        for parameter in envelope_converted_model_parameters(
            current_version(), normalized=True)
    }


def model_catalogue_bounds(catalogue):
    """Return range and scale metadata for a named model parameter catalogue."""
    version = current_version()
    if catalogue == 'default':
        catalogue_parameters = build_model(version, normalized=True).parameters()
    elif catalogue == 'converted-envelope':
        catalogue_parameters = envelope_converted_model_parameters(
            version, normalized=True)
    elif catalogue == 'hypothesis':
        catalogue_parameters = hypothesis_model_parameters(version, normalized=True)
    else:
        raise ValueError(f'Unknown parameter catalogue {catalogue!r}.')
    return {
        parameter.name: {
            'low': parameter.low,
            'high': parameter.high,
            'value': parameter.value,
            'scale': parameter.scale,
        }
        for parameter in catalogue_parameters
    }


def _is_logarithmic(bounds, reference):
    values = [value for value in (
        bounds['low'], bounds['high'],
        None if reference is None else reference['low'],
        None if reference is None else reference['high'])
        if value is not None]
    return (reference is not None and reference['scale'] == 'log'
            and all(value > 0 for value in values))


def checkpoint_values(document):
    """Load parameter values and losses from checkpoint files named by the artifact."""
    values = {}
    missing = []
    records = document.get('checkpoints', ())
    if not isinstance(records, list):
        return values
    for record in records:
        identifier = record.get('checkpoint_id')
        source = record.get('source')
        try:
            if source:
                path = Path(source)
                if not path.is_file():
                    raise FileNotFoundError(f'recorded source does not exist: {path}')
            else:
                path = checkpoint_io.resolve_checkpoint(identifier)
            contents = parameters.load_file(path)
            values[str(identifier)] = {
                'parameters': parameters.ModelParameters(contents['parameters']).get_values(),
                'loss': float(contents['total_loss']),
            }
        except (OSError, KeyError, TypeError, ValueError) as error:
            missing.append(f'{identifier} ({error})')
    if missing:
        warnings.warn(
            'Skipping checkpoint(s) no longer available for marker plotting: '
            + '; '.join(missing), stacklevel=2)
    return values


def draw_parameter(axis, bounds, reference, reference_label, values, colour_map, colour_norm):
    """Draw used checkpoint values and their selected-catalogue reference."""
    low, high = bounds['low'], bounds['high']
    log_scale = _is_logarithmic(bounds, reference)
    if reference is not None and reference['low'] is not None:
        axis.axvspan(reference['low'], reference['high'], color='tab:orange',
                     alpha=.24, label=f'{reference_label} bounds')
    elif reference is not None:
        axis.axvline(reference['value'], color='tab:orange', linestyle='--', linewidth=2,
                    label=f'{reference_label} fixed value')
    markers = [
        {'value': checkpoint['parameters'][bounds['parameter']], 'loss': checkpoint['loss']}
        for checkpoint in values.values()
        if bounds['parameter'] in checkpoint['parameters']]
    if not markers:
        raise ValueError(
            f"No checkpoint values found for parameter {bounds['parameter']!r}.")
    stored_values = [marker['value'] for marker in markers]
    offsets = np.linspace(.38, .62, len(markers))
    axis.scatter(stored_values, offsets, c=[marker['loss'] for marker in markers],
                 cmap=colour_map, norm=colour_norm, s=46, zorder=3,
                 label='used checkpoint value')
    if log_scale:
        axis.set_xscale('log')
    values = [low, high, *stored_values]
    if reference is not None:
        values.extend(value for value in (reference['low'], reference['high'], reference['value'])
                      if value is not None)
    minimum, maximum = min(values), max(values)
    if log_scale:
        axis.set_xlim(minimum / 1.7, maximum * 1.7)
    elif minimum == maximum:
        padding = max(abs(minimum) * .1, 1.)
        axis.set_xlim(minimum - padding, maximum + padding)
    else:
        padding = (maximum - minimum) * .12
        axis.set_xlim(minimum - padding, maximum + padding)
    axis.set_yticks([])
    axis.grid(axis='x', color='.9')
    return log_scale


def loss_normalizer(values):
    """Choose a readable loss scale for checkpoint-marker colours."""
    losses = np.asarray([checkpoint['loss'] for checkpoint in values.values()], dtype=float)
    low, high = float(losses.min()), float(losses.max())
    if low > 0 and high / low >= 100:
        return LogNorm(vmin=low, vmax=high)
    if low == high:
        high = low + max(abs(low) * .01, 1.)
    return Normalize(vmin=low, vmax=high)


def plot_envelope(document, output, catalogue='converted-envelope'):
    """Write a multi-panel PNG comparing stored and selected-catalogue bounds."""
    rows = document['parameters']
    if not rows:
        raise ValueError('Bounds file contains no parameters.')
    reference_bounds = model_catalogue_bounds(catalogue)
    reference_label = f'--{catalogue}'
    values = checkpoint_values(document)
    if not values:
        raise ValueError('Bounds file records no checkpoints to plot.')
    colour_map = plt.get_cmap('viridis_r')
    colour_norm = loss_normalizer(values)
    columns = 2
    rows_count = math.ceil(len(rows) / columns)
    figure, axes = plt.subplots(rows_count, columns, figsize=(13, 2.35 * rows_count),
                                squeeze=False)
    for axis, bounds in zip(axes.flat, rows):
        reference = reference_bounds.get(bounds['parameter'])
        draw_parameter(axis, bounds, reference, reference_label, values, colour_map, colour_norm)
        unavailable = '' if reference is not None else f' (not in {catalogue})'
        axis.set_title(bounds['parameter'] + unavailable, loc='left', fontsize=10)
    for axis in axes.flat[len(rows):]:
        axis.set_visible(False)
    handles, labels = axes.flat[0].get_legend_handles_labels()
    figure.legend(handles, labels, loc='upper center', ncol=2, frameon=False)
    colourbar = figure.colorbar(
        plt.cm.ScalarMappable(norm=colour_norm, cmap=colour_map), ax=list(axes.flat),
        pad=.02, fraction=.025)
    colourbar.set_label('Checkpoint loss')
    figure.suptitle(f"Stored bounds {document.get('bounds_id', '?')} vs {catalogue} bounds",
                    y=.995)
    figure.tight_layout(rect=(0, 0, 1, .94))
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=180)
    plt.close(figure)
    return output


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('bounds_id', help='Full or unique partial bounds ID.')
    catalogue = parser.add_mutually_exclusive_group(required=True)
    catalogue.add_argument('--default', action='store_true',
                           help='Compare against the default model parameters.')
    catalogue.add_argument('--converted-envelope', action='store_true',
                           help='Compare against the converted-envelope parameters.')
    catalogue.add_argument('--hypothesis', action='store_true',
                           help='Compare against the hypothesis model parameters.')
    parser.add_argument('--bounds-directory', type=Path,
                        default=parameter_bounds.BOUNDS_DIRECTORY,
                        help='Directory containing stored bounds files (default: bounds/).')
    parser.add_argument('--output', type=Path,
                        help='PNG destination (default: plots/bounds/bounds-<ID>.png).')
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    catalogue = ('default' if args.default else 'converted-envelope'
                 if args.converted_envelope else 'hypothesis')
    source, document = parameter_bounds.load_bounds(args.bounds_id, args.bounds_directory)
    bounds_id = document.get('bounds_id', source.stem.removeprefix('bounds-'))
    output = args.output or Path(PLOT_DIRECTORY) / 'bounds' / f'bounds-{bounds_id}.png'
    plot_envelope(document, output, catalogue)
    print(f'Wrote {output}')
    return output


if __name__ == '__main__':
    main()
