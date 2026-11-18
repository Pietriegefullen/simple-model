"""Map the conditional loss surface of a pathway's ``M(0)`` and ``v_max``.

The command replays a saved single- or split-fit checkpoint over a logarithmic
two-parameter grid.  It deliberately does *not* optimise at each grid point:
all parameters other than the selected pathway's initial biomass and specific
maximum rate retain their checkpoint values.  Consequently, the figure is a
direct local test of their conditional identifiability.

Examples
--------
    python faster/plot_m0_vmax_surface.py 6LND Hydro --grid-density 81
    python faster/plot_m0_vmax_surface.py IUPI Hydro --grid-density 61 --log-span 1.5

Results are written below ``results/m0_vmax_surfaces/`` unless ``--output-dir``
is supplied.  The PNG is accompanied by a CSV grid and a JSON provenance file.
"""

from __future__ import annotations

import argparse
import copy
import csv
import json
import math
import re
from datetime import datetime, timezone
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import data
from checkpoint_io import checkpoint_id, replay_inputs
from model_factory import build_model
import optimizer
import parameters
from USER_VARIABLES import PROJECT_DIRECTORY


# Hydrolysis and fermentation intentionally share M_Ferm.  They are excluded
# here: scaling it to compensate one rate would also alter the other pathway's
# initial capacity, so it is not the requested two-parameter experiment.
PATHWAY_PARAMETERS = {
    'aceto': ('Aceto', 'M_Ac', 'Ac_v_max'),
    'fe3': ('Fe3', 'M_Fe3', 'Fe3_v_max'),
    'homo': ('Homo', 'M_Homo', 'Homo_v_max'),
    'hydro': ('Hydro', 'M_Hydro', 'Hydro_v_max'),
}
UNSUPPORTED_PATHWAYS = {
    'fermentation': 'Fermentation shares M_Ferm with Hydrolysis. '
                    'A two-parameter surface would also change Hydrolysis.',
    'hydrolysis': 'Hydrolysis shares M_Ferm with Fermentation. '
                  'A two-parameter surface would also change Fermentation.',
}


def pathway_parameters(name):
    """Return canonical pathway, initial-biomass, and rate parameter names."""
    key = name.casefold().replace('_', '').replace('-', '').replace(' ', '')
    if key in PATHWAY_PARAMETERS:
        return PATHWAY_PARAMETERS[key]
    if key in UNSUPPORTED_PATHWAYS:
        raise ValueError(UNSUPPORTED_PATHWAYS[key] + ' Choose Aceto, Fe3, Homo, or Hydro.')
    choices = ', '.join(item[0] for item in PATHWAY_PARAMETERS.values())
    raise ValueError(f'Unknown pathway {name!r}. Choose one of: {choices}.')


def parameter_bounds(parameter):
    """Return finite positive lower/upper bounds for a log-scaled parameter."""
    lower, upper = parameter.lower(), parameter.upper()
    if not (math.isfinite(lower) and math.isfinite(upper) and 0 < lower < upper):
        raise ValueError(f'{parameter.name} needs finite, positive bounds; got {lower!r}, {upper!r}.')
    return float(lower), float(upper)


def logarithmic_grid(centre, bounds, span, density):
    """Return a centred log grid, clipped to both requested span and bounds."""
    if density < 3:
        raise ValueError('--grid-density must be at least 3.')
    if not math.isfinite(span) or span <= 0:
        raise ValueError('--log-span must be a positive finite number of decades.')
    if centre <= 0:
        raise ValueError('Checkpoint parameter values must be positive for a logarithmic grid.')
    lower, upper = bounds
    requested_low, requested_high = centre * 10 ** (-span), centre * 10 ** span
    low, high = max(lower, requested_low), min(upper, requested_high)
    if not low < high:
        raise ValueError(f'No usable grid range around {centre:g} within bounds {bounds}.')
    return np.geomspace(low, high, density)


def grid_loss(run_config, fitted_parameters, biomass_name, rate_name,
              biomass_grid, rate_grid):
    """Evaluate the saved fitting objective at every grid cell."""
    chosen = run_config['chosen']
    if chosen.get('fit_mode') == 'multi':
        raise ValueError('Joint multi-replica checkpoints require a selected local replica and are not yet '
                         'supported by this conditional surface command.')
    dataset = data.get_data_before_day()
    sample = dataset[chosen['sample']]
    split = sample.get_split(chosen['validation_replica'], chosen['fit_mode'])
    model = build_model(run_config['model']['version'], chosen['pathways'],
                        chosen['normalized_parameters'])
    model.parameters().set(copy.deepcopy(fitted_parameters))
    objective = optimizer.build_objective_function(
        model, split['fit'], run_config['objective'], t_start=chosen.get('t_start'),
        t_end=chosen.get('t_end'))

    losses = np.full((len(rate_grid), len(biomass_grid)), np.nan)
    cnt = 0
    for rate_index, rate in enumerate(rate_grid):
        for biomass_index, biomass in enumerate(biomass_grid):
            cnt += 1
            print(f'\r{cnt}/{losses.size} ({cnt/losses.size*100:.1f}%)', end = '')
            candidate = copy.deepcopy(fitted_parameters)
            candidate[biomass_name].set(float(biomass))
            candidate[rate_name].set(float(rate))
            try:
                losses[rate_index, biomass_index] = objective(candidate, transformed=False)
            except Exception as error:
                # Numerical failures are visible as blank cells while preserving
                # the rest of a potentially expensive grid evaluation.
                print(f'Grid cell M={biomass:.5g}, v={rate:.5g} failed: {error}')
    return losses


def metadata_fragment(run_config):
    """Build concise provenance included in both title and file name."""
    chosen = run_config['chosen']
    return (str(chosen['sample']), str(chosen['validation_replica']),
            str(chosen['fit_mode']), str(run_config['model']['version']))


def slug(value):
    return re.sub(r'[^A-Za-z0-9_.-]+', '-', str(value)).strip('-')


def write_heatmap(losses, biomass_grid, rate_grid, biomass_optimum, rate_optimum,
                  pathway, title, output):
    """Write the heatmap and constant-product reference line."""
    finite = np.isfinite(losses)
    if not np.any(finite):
        raise RuntimeError('Every grid evaluation failed; no heatmap can be drawn.')
    reference_loss = losses[np.unravel_index(np.nanargmin(losses), losses.shape)]
    delta_loss = losses - reference_loss
    figure, axis = plt.subplots(figsize=(9.2, 7.2), constrained_layout=True)
    mesh = axis.pcolormesh(biomass_grid, rate_grid, delta_loss, shading='auto', cmap='viridis')
    colorbar = figure.colorbar(mesh, ax=axis, pad=.02)
    colorbar.set_label(r'$\Delta$ objective (relative to best grid cell)')

    product = biomass_optimum * rate_optimum
    ridge_rate = product / biomass_grid
    usable = (ridge_rate >= rate_grid.min()) & (ridge_rate <= rate_grid.max())
    axis.plot(biomass_grid[usable], ridge_rate[usable], color='white', linewidth=1.7,
              linestyle='--', label=r'$M_0v_{\max}=M_{0,\mathrm{opt}}v_{\max,\mathrm{opt}}$')
    axis.scatter([biomass_optimum], [rate_optimum], color='white', edgecolor='black',
                 linewidth=.8, s=52, zorder=3, label='checkpoint')
    axis.set_xscale('log')
    axis.set_yscale('log')
    axis.set_xlabel(r'initial biomass $M_0$')
    axis.set_ylabel(r'specific maximum rate $v_{\max}$')
    axis.set_title(title)
    axis.legend(loc='best', framealpha=.9)
    figure.savefig(output, dpi=200)
    plt.close(figure)
    return float(reference_loss)


def write_grid_csv(path, biomass_grid, rate_grid, losses):
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=('M0', 'v_max', 'loss'))
        writer.writeheader()
        for rate_index, rate in enumerate(rate_grid):
            for biomass_index, biomass in enumerate(biomass_grid):
                writer.writerow({'M0': biomass, 'v_max': rate,
                                 'loss': losses[rate_index, biomass_index]})


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('checkpoint', help='full or uniquely resolvable partial checkpoint ID')
    parser.add_argument('pathway', help='Aceto, Fe3, Homo, or Hydro')
    parser.add_argument('--grid-density', type=int, default=61,
                        help='number of points along each logarithmic axis (default: 61)')
    parser.add_argument('--log-span', type=float, default=1.,
                        help='half-width of the grid in decades around the checkpoint (default: 1)')
    parser.add_argument('--output-dir', type=Path,
                        help='directory for this analysis; default is a dedicated results subdirectory')
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    pathway, biomass_name, rate_name = pathway_parameters(args.pathway)
    source, run_config, fitted_parameters, contents = replay_inputs(args.checkpoint)
    if biomass_name not in fitted_parameters or rate_name not in fitted_parameters:
        raise ValueError(f'Checkpoint does not contain {biomass_name} and {rate_name}.')
    biomass_optimum, rate_optimum = float(fitted_parameters[biomass_name]), float(fitted_parameters[rate_name])
    # Checkpoints retain fitted values but not the searchable bounds.  Recover
    # those bounds from the exact versioned model described by the checkpoint.
    reference_model = build_model(run_config['model']['version'], run_config['chosen']['pathways'],
                                  run_config['chosen']['normalized_parameters'])
    model_parameters = reference_model.parameters()
    biomass_grid = logarithmic_grid(biomass_optimum, parameter_bounds(model_parameters[biomass_name]),
                                    args.log_span, args.grid_density)
    rate_grid = logarithmic_grid(rate_optimum, parameter_bounds(model_parameters[rate_name]),
                                 args.log_span, args.grid_density)

    short_id = checkpoint_id(source) or source.name
    sample, replica, fit_mode, model_version = metadata_fragment(run_config)
    filename_stem = (f'sample-{slug(sample)}_replica-{slug(replica)}_fit-{slug(fit_mode)}_'
                     f'model-{slug(model_version)}_cp-{slug(short_id)}_pathway-{slug(pathway)}_'
                     f'grid-{args.grid_density}_span-{args.log_span:g}dec')
    output_dir = (args.output_dir or Path(PROJECT_DIRECTORY) / 'results' / 'm0_vmax_surfaces' /
                  f'sample-{slug(sample)}_replica-{slug(replica)}').resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f'Evaluating {args.grid_density ** 2} cells for {pathway} from checkpoint cp-{short_id}.')
    losses = grid_loss(run_config, fitted_parameters, biomass_name, rate_name, biomass_grid, rate_grid)
    title = (f'{pathway}: conditional $M_0$--$v_{{\\max}}$ loss surface\n'
             f'checkpoint cp-{short_id} | sample {sample}, replica {replica}\n'
             f'fit {fit_mode}, model {model_version}; all other parameters fixed at checkpoint')
    image_path = output_dir / f'{filename_stem}.png'
    best_grid_loss = write_heatmap(losses, biomass_grid, rate_grid, biomass_optimum, rate_optimum,
                                   pathway, title, image_path)
    csv_path = output_dir / f'{filename_stem}.csv'
    write_grid_csv(csv_path, biomass_grid, rate_grid, losses)
    metadata = {
        'created_at_utc': datetime.now(timezone.utc).isoformat(),
        'source_checkpoint': str(source.resolve()), 'checkpoint_id': short_id,
        'pathway': pathway, 'biomass_parameter': biomass_name, 'v_max_parameter': rate_name,
        'checkpoint_M0': biomass_optimum, 'checkpoint_v_max': rate_optimum,
        'checkpoint_product_M0_v_max': biomass_optimum * rate_optimum,
        'sample': sample, 'replica': replica, 'fit_mode': fit_mode, 'model_version': model_version,
        'grid_density': args.grid_density, 'log_span_decades': args.log_span,
        'M0_grid_bounds': [float(biomass_grid[0]), float(biomass_grid[-1])],
        'v_max_grid_bounds': [float(rate_grid[0]), float(rate_grid[-1])],
        'best_grid_loss': best_grid_loss,
        'method': 'conditional objective replay; all parameters except M0 and v_max fixed',
        'files': {'plot': image_path.name, 'grid': csv_path.name},
        'saved_checkpoint_loss': contents.get('total_loss'),
    }
    metadata_path = output_dir / f'{filename_stem}.json'
    with metadata_path.open('w') as handle:
        json.dump(metadata, handle, indent=2)
    print(f'Wrote heatmap: {image_path}')
    return image_path


if __name__ == '__main__':
    main()
