"""Collect the best checkpoint for each model, replica, and fit mode.

Running this module writes a data/model comparison plot for every selected
checkpoint below ``best_<loss-id>_<model-id>_<fit-mode>``.
"""

import argparse
import os
from dataclasses import dataclass
from functools import cached_property
from itertools import cycle
import traceback

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

import data
import hashing
import optimizer
import parameters
import plot
from fit_sample import run
from model_variants import model_variant_from_id
from USER_VARIABLES import PROJECT_DIRECTORY, RESULTS_DIRECTORY

# ``00`` is reserved for the data/model fit (CO2 and CH4).  The remaining
# system pools are grouped with the pathway diagnostics below, rather than
# duplicating the gases or a microbe pool in a separate group.
POOL_PLOT_ORDER = ('TOC', 'DOC', 'Acetate', 'Fe3', 'H2', 'H2O', 'Fe2')
PATHWAY_PLOT_ORDER = (
    'Hydrolysis', 'Fermentation', 'Fe3', 'Aceto', 'Hydro', 'Homo',
)
PATHWAY_MICROBE_POOL = {
    'Hydrolysis': 'M_Ferm',
    'Fermentation': 'M_Ferm',
    'Fe3': 'M_Fe3',
    'Aceto': 'M_Ac',
    'Hydro': 'M_Hydro',
    'Homo': 'M_Homo',
}
PATHWAY_V_MAX_PARAMETER = {
    'Hydrolysis': 'Hydrolysis_v_max',
    'Fermentation': 'Ferm_v_max',
    'Fe3': 'Fe3_v_max',
    'Aceto': 'Ac_v_max',
    'Hydro': 'Hydro_v_max',
    'Homo': 'Homo_v_max',
}
PATHWAY_QUANTITY_ORDER = (
    'thermodynamic_factor', 'deltaG_r', 'MM', 'microbe_pool',
)

@dataclass
class Checkpoint:
    """A loaded checkpoint and the metadata used to group it."""

    parameters: parameters.ModelParameters
    loss: float
    run_config: dict
    id: str
    source: str = None

    @cached_property
    def model_id(self):
        return hashing.build_model_id(self.run_config['model'])

    @cached_property
    def loss_id(self):
        return hashing.build_loss_id(self.run_config['objective'])

    @cached_property
    def replica(self):
        chosen = self.run_config['chosen']
        return str(chosen['sample']), str(chosen['validation_replica'])

    @cached_property
    def replica_name(self):
        return '-'.join(self.replica)

    @cached_property
    def fit_mode(self):
        return self.run_config['chosen']['fit_mode']

    @cached_property
    def key(self):
        return '_'.join((self.loss_id, self.model_id))
    
    def __str__(self):
        return self.key

def collect_all(results_directory=RESULTS_DIRECTORY):
    """Load every annotated checkpoint in ``results_directory``."""
    collected = []
    for root, _, files in os.walk(results_directory):
        for filename in files:
            if filename.startswith('.') or filename.endswith('.log') or filename.endswith('.json'):
                continue

            file_path = os.path.join(root, filename)
            try:
                loaded_parameters, loss, run_config = parameters.load_parameter_file(file_path)
            except Exception as ex:
                print(traceback.format_exc())
                input()
            checkpoint_ids = [part for part in filename.split('_')
                              if part.startswith('cp-')]
            if not checkpoint_ids:
                raise ValueError(f'Checkpoint ID missing from {file_path}')
            if not run_config:
                raise ValueError(f'Run configuration missing from {file_path}')

            collected.append(Checkpoint(loaded_parameters, loss, run_config,
                                        checkpoint_ids[0], file_path))
    return collected


def collect_best(results_directory=RESULTS_DIRECTORY):
    """Return the lowest-loss checkpoint per model/loss/replica/fit-mode."""
    best_by_group = {}
    for checkpoint in collect_all(results_directory):
        group = (checkpoint.key, checkpoint.fit_mode, checkpoint.replica)
        best = best_by_group.get(group)
        if best is None or checkpoint.loss < best.loss:
            best_by_group[group] = checkpoint
    return list(best_by_group.values())


def filter_checkpoints(checkpoints, samples=None, model_ids=None, fit_modes=None,
                       replicas=None, loss_ids=None):
    """Keep checkpoints matching every supplied command-line filter."""
    filters = {
        'sample': None if samples is None else set(map(str, samples)),
        'model_id': None if model_ids is None else set(model_ids),
        'fit_mode': None if fit_modes is None else set(fit_modes),
        'replica': None if replicas is None else set(map(str, replicas)),
        'loss_id': None if loss_ids is None else set(loss_ids),
    }
    selected = []
    for checkpoint in checkpoints:
        sample, replica = checkpoint.replica
        if filters['sample'] is not None and sample not in filters['sample']:
            continue
        if (filters['model_id'] is not None and
                checkpoint.model_id not in filters['model_id']):
            continue
        if (filters['fit_mode'] is not None and
                checkpoint.fit_mode not in filters['fit_mode']):
            continue
        if filters['replica'] is not None and replica not in filters['replica']:
            continue
        if filters['loss_id'] is not None and checkpoint.loss_id not in filters['loss_id']:
            continue
        selected.append(checkpoint)
    return selected


def plot_target(checkpoint, number=0, quantity=None):
    """Return the directory and basename for a numbered checkpoint plot."""
    directory = os.path.join(
        PROJECT_DIRECTORY,
        '_'.join(('best', checkpoint.key, checkpoint.fit_mode)),
        checkpoint.replica_name,
    )
    filename_parts = [
        f'{number:02d}',
    ]
    if quantity is not None:
        filename_parts.append(quantity)
    filename_parts.extend((
        checkpoint.replica_name, checkpoint.model_id,
        checkpoint.fit_mode, checkpoint.id,
    ))
    return directory, '_'.join(filename_parts)


def _run_log_keys(run_log):
    """Return run-log keys without depending on its private storage."""
    return set(run_log.keys())


def pathway_plot_quantities(run_log):
    """Yield ``(name, time, values)`` in the required plot-file order.

    Quantities added to a pathway log in the future are appended after the
    specified diagnostics for that pathway.  ``v_max`` is shown as a dashed
    reference line on the corresponding ``v`` plot, not as its own plot.
    """
    keys = _run_log_keys(run_log)

    for pool in POOL_PLOT_ORDER:
        if pool in keys:
            yield pool, *run_log[pool]

    for pathway in PATHWAY_PLOT_ORDER:
        prefix = f'{pathway}_'
        pathway_keys = {key for key in keys if key.startswith(prefix)}
        if not pathway_keys:
            continue

        for quantity in PATHWAY_QUANTITY_ORDER:
            if quantity == 'microbe_pool':
                microbe_pool = PATHWAY_MICROBE_POOL[pathway]
                if microbe_pool in keys:
                    yield microbe_pool, *run_log[microbe_pool]
            else:
                name = prefix + quantity
                if name in pathway_keys:
                    yield name, *run_log[name]

        handled = {
            prefix + quantity
            for quantity in PATHWAY_QUANTITY_ORDER
            if quantity != 'microbe_pool'
        }
        # Add every remaining pathway diagnostic (for example ``v``,
        # inhibition, and methane production) in a stable order.
        for name in sorted(pathway_keys - handled):
            yield name, *run_log[name]


def pathway_v_max_values(run_log):
    """Return the configured maximum rate for each active pathway."""
    parameters = getattr(run_log, '_parameters', None)
    if parameters is None:
        return {}

    keys = _run_log_keys(run_log)
    return {
        pathway: float(parameters[PATHWAY_V_MAX_PARAMETER[pathway]])
        for pathway in PATHWAY_PLOT_ORDER
        if any(key.startswith(f'{pathway}_') for key in keys)
    }


def quantity_ylabel(name):
    """Return a compact, quantity-appropriate y-axis label."""
    if name.endswith('_deltaG_r'):
        return r'$\Delta G_\mathrm{r}\ [J/mol]$'
    if name.endswith(('_thermodynamic_factor', '_MM', '_inhib')):
        return r'$\mathrm{factor}$'
    if name.endswith('_v') or name.endswith('_v_max'):
        return r'$\mathrm{rate}$'
    return r'$\mathrm{substance\ [\mu mol/g\ dry\ weight]}$'


def figure_suptitle(checkpoint, fit_replicas=(), curve_replica=None):
    """Return shared metadata, including the data behind a model curve."""
    fit_replica_numbers = ', '.join(
        str(replica.replica_number) for replica in fit_replicas)
    curve_replica_number = getattr(curve_replica, 'replica_number', None)
    first_line = ' | '.join((
        checkpoint.replica[0],
        checkpoint.fit_mode,
        f'model {model_variant(checkpoint)}',
    ))
    return first_line


def plot_quantity(name, time, values, checkpoint, v_max=None, suptitle=None):
    """Plot one modeled pool or pathway quantity in the fit-plot style."""
    figure, axis = plt.subplots(figsize=(4, 4))
    axis.plot(time, values, '-', color='k', label=name)
    if v_max is not None:
        axis.axhline(v_max, color='k', linestyle='--', label='v_max')
        axis.legend(fancybox=False, edgecolor='k', loc='best')
    axis.set_title(name)
    axis.set_xlabel(r'$\mathrm{time\ [d]}$')
    axis.set_ylabel(quantity_ylabel(name))
    if name.endswith(('_thermodynamic_factor', '_MM')):
        axis.set_ylim(0, 1)
    plot.format_ax(axis)
    figure.suptitle(figure_suptitle(checkpoint) if suptitle is None else suptitle)
    figure.tight_layout(rect=(0, 0, 1, 0.85))
    return figure, axis

def model_variant(cp):
    variant = model_variant_from_id(cp.model_id)
    if variant is not None:
        return variant
    raise Exception('Model variant not identified: ' + cp.model_id)


def sample_y_limits(sample, relative_padding=0.05):
    """Return fixed CO2 and CH4 limits from all measured sample replicas."""
    limits = {}
    for pool in ('CO2', 'CH4'):
        values = []
        for replica in sample.replicas:
            _, replica_values = getattr(replica, pool)()
            replica_values = np.asarray(replica_values, dtype=float)
            usable = np.isfinite(replica_values)
            if pool == 'CH4':
                usable &= replica_values > 0
            values.extend(replica_values[usable])

        if not values:
            continue
        values = np.asarray(values)
        if pool == 'CO2':
            low, high = np.min(values), np.max(values)
            padding = max((high - low) * relative_padding,
                          max(abs(low), abs(high), 1.) * relative_padding)
            limits[pool] = (low - padding, high + padding)
        else:
            log_values = np.log10(values)
            low, high = np.min(log_values), np.max(log_values)
            padding = max((high - low) * relative_padding, relative_padding)
            limits[pool] = (10**(low - padding), 10**(high + padding))
    return limits


def fit_replicas(checkpoint, dataset):
    """Return the sample, validation replica, and fitted replicas for a plot."""
    sample_number, validation_number = checkpoint.replica
    sample = dataset[sample_number]
    validation_replica = sample[validation_number]
    replicas = sample.get_split(validation_number, checkpoint.fit_mode)['fit']
    return sample, validation_replica, replicas


def create_fit_figure():
    """Create the standard two-pool figure used for fitted-model plots."""
    return plot.create_figure(plot.FigureSpec(
        ncols=2,
        axis_names=('CO2', 'CH4'),
        figsize=(10, 5),
        sharex=True,
    ))


def plot_fit_measurements(checkpoint, dataset, axes):
    """Draw the fit/validation measurements and return their legend handles."""
    _, validation_replica, fitted_replicas = fit_replicas(checkpoint, dataset)
    markersize = 3
    markers = cycle(('^', 's', 'D', 'v', 'P', 'X'))
    handles = []
    for fitted_replica in fitted_replicas:
        marker = next(markers)
        plot.plot_replica_pools(
            fitted_replica, axes, marker=marker, ms=markersize, mfc='none')
        handles.append(Line2D(
            [], [], color='k', marker=marker, linestyle='None',
            markerfacecolor='none', markersize=markersize,
            label=f'fit replica {fitted_replica.replica_number}',
        ))

    if checkpoint.fit_mode == 'split':
        plot.plot_replica_pools(
            validation_replica, axes, marker='o', ms=markersize)
        handles.append(Line2D(
            [], [], color='k', marker='o', linestyle='None', markersize=markersize,
            label=f'validation replica {validation_replica.replica_number}',
        ))
    return handles


def fit_score_label(checkpoint, objective, validation_objective, run_log,
                    model_label=None):
    """Return the R² label used for a model line in a fit comparison."""
    r2_values = [(objective.R2(run_log), 'fit')]
    if checkpoint.fit_mode != 'single':
        r2_values.append((validation_objective.R2(run_log), 'val'))
    r2_text = ', '.join(
        rf'$R^2_{{\mathrm{{{name}}}}} = {value:.2f}$'
        for value, name in r2_values
    )
    if model_label is None:
        model_label = rf'$\mathrm{{model\ {{{model_variant(checkpoint)}}}}}$'
    return ' '.join((model_label, r2_text))


def format_fit_axes(checkpoint, sample, axes):
    """Apply the standard fit-plot scales, limits, and axes style."""
    y_limits = sample_y_limits(sample)
    for pool, axis in axes.items():
        plot.format_ax(
            axis,
            log_scale='log' in checkpoint.run_config['objective']['transform'][pool],
        )
        if pool in y_limits:
            axis.set_ylim(y_limits[pool])


def finish_fit_figure(figure, suptitle):
    """Apply the fixed report layout shared by fit comparison figures."""
    figure.suptitle(suptitle)
    figure.subplots_adjust(
        left=0.12, right=0.97, bottom=0.14, top=0.78, wspace=0.30)



def plot_fit(checkpoint, dataset, model_run=None):
    """Plot fitted replicas, validation data, and a model run for a checkpoint."""
    if model_run is None:
        model_run = run(checkpoint.run_config, checkpoint.parameters)
    pathway_model, objective, run_log = model_run
    sample, validation_replica, fitted_replicas = fit_replicas(checkpoint, dataset)
    validation_objective = optimizer.build_objective_function(
        pathway_model,
        validation_replica,
        checkpoint.run_config['objective'],
        t_start=0,
        t_end=None,
    )
    figure, axes = create_fit_figure()
    handles = plot_fit_measurements(checkpoint, dataset, axes)
    plot.plot_run_pools(run_log, axes)

    handles.append(Line2D(
        [], [], color='k', linestyle='-',
        label=fit_score_label(checkpoint, objective, validation_objective, run_log),
    ))

    axes['CH4'].legend(
        handles=handles, fancybox=False, edgecolor='k', loc='lower right')
    format_fit_axes(checkpoint, sample, axes)
    finish_fit_figure(figure, figure_suptitle(
        checkpoint, fitted_replicas, getattr(run_log, 'replica', None)))

    return figure, axes


def plot_all(results_directory=RESULTS_DIRECTORY, **filters):
    """Generate plots for all best checkpoints and return their paths."""
    dataset = data.get_data_before_carex()
    paths = []
    checkpoints = filter_checkpoints(collect_best(results_directory), **filters)
    for checkpoint in checkpoints:
        model_run = run(checkpoint.run_config, checkpoint.parameters)
        figure, _ = plot_fit(checkpoint, dataset, model_run=model_run)
        directory, filename = plot_target(checkpoint)

        os.makedirs(directory, exist_ok=True)
        path = os.path.join(directory, filename)
        #print('saving to  ', os.path.split(directory)[-1], filename)
        figure.savefig(path, dpi=300)
        plt.close(figure)
        paths.append(path)

        _, _, run_log = model_run
        sample = dataset[checkpoint.replica[0]]
        fit_replicas = sample.get_split(
            checkpoint.replica[1], checkpoint.fit_mode)['fit']
        suptitle = figure_suptitle(
            checkpoint, fit_replicas, getattr(run_log, 'replica', None))
        v_max_values = pathway_v_max_values(run_log)
        for number, (quantity, time, values) in enumerate(
                pathway_plot_quantities(run_log), start=1):
            pathway = quantity[:-2] if quantity.endswith('_v') else None
            figure, _ = plot_quantity(
                quantity, time, values, checkpoint,
                v_max=v_max_values.get(pathway), suptitle=suptitle)
            _, filename = plot_target(checkpoint, number, quantity)
            path = os.path.join(directory, filename)
            # print('saving to  ', os.path.split(directory)[-1], filename)
            figure.savefig(path, dpi=300)
            plt.close(figure)
            paths.append(path)
    return paths


def parse_args():
    parser = argparse.ArgumentParser(
        description='Plot the best checkpoints, optionally limited to selected fits.')
    parser.add_argument('--results-directory', default=RESULTS_DIRECTORY,
                        help='Checkpoint root to scan.')
    parser.add_argument('--sample', dest='samples', nargs='+', metavar='SAMPLE',
                        help='Only plot these sample IDs.')
    parser.add_argument('--model-id', dest='model_ids', nargs='+', metavar='MODEL_ID',
                        help='Only plot these exact model IDs.')
    parser.add_argument('--fit-mode', dest='fit_modes', choices=('single', 'split'),
                        nargs='+', metavar='FIT_MODE',
                        help='Only plot these fit modes.')
    parser.add_argument('--replica', dest='replicas', nargs='+', metavar='REPLICA',
                        help='Only plot these validation replica numbers.')
    parser.add_argument('--loss-id', dest='loss_ids', nargs='+', metavar='LOSS_ID',
                        help='Only plot these exact loss IDs.')
    parser.add_argument('--dry-run', action='store_true',
                        help='Print selected checkpoints without creating plots.')
    return parser.parse_args()


def main(args=None):
    args = parse_args() if args is None else args
    filters = {
        'samples': args.samples,
        'model_ids': args.model_ids,
        'fit_modes': args.fit_modes,
        'replicas': args.replicas,
        'loss_ids': args.loss_ids,
    }
    if args.dry_run:
        checkpoints = filter_checkpoints(
            collect_best(args.results_directory), **filters)
        for checkpoint in checkpoints:
            print(checkpoint.model_id, checkpoint.loss_id,
                  checkpoint.replica_name, checkpoint.fit_mode)
        print(f'{len(checkpoints)} checkpoint(s) selected')
        return checkpoints
    return plot_all(args.results_directory, **filters)


if __name__ == '__main__':
    main()
