"""Plot one or more saved checkpoints, including aligned comparisons.

The command accepts the same full or partial checkpoint identifiers as
``run_checkpoint.py``.  A single checkpoint keeps the established
``plot_best`` report.  With several checkpoints, each result is added to
the corresponding fit and diagnostic figure instead of producing a separate
set of files.
"""

from __future__ import annotations

import argparse
import copy
from dataclasses import dataclass
from fnmatch import fnmatchcase
from pathlib import Path
from typing import Any, Iterable

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

import data
import parameters
import plot
from fit_sample import _requested_checkpoint, run
from run_checkpoint import (checkpoint_id, prepare_output_directory,
                            replay_configuration, resolve_checkpoint)
from USER_VARIABLES import RESULTS_DIRECTORY, PLOT_DIRECTORY
from user_input import parse


@dataclass(frozen=True)
class CheckpointPlot:
    """A checkpoint together with the one model evaluation used to plot it."""

    checkpoint: Any
    model_run: tuple[Any, Any, Any]
    comparison_name: str | None = None


@dataclass(frozen=True)
class ComparisonSpec:
    """A named, parameter-only counterfactual replay specification."""

    name: str
    overrides: tuple[tuple[str, Any], ...]


def add_comparison_argument(parser) -> None:
    """Add the repeatable command-line option shared by plot commands."""
    parser.add_argument(
        '--compare', action='append', nargs='+', metavar='SPEC', default=[],
        help=('add a named replay: NAME PARAMETER VALUE [PARAMETER VALUE ...]. '
              'Parameter names support shell-style wildcards, e.g. '
              '--compare D "*_thermodynamics" False'))


def parse_comparisons(specifications) -> list[ComparisonSpec]:
    """Parse ``--compare`` values while retaining wildcard expansion for replay."""
    comparisons = []
    for specification in specifications or []:
        if len(specification) < 3 or (len(specification) - 1) % 2:
            raise ValueError(
                '--compare requires NAME followed by one or more PARAMETER VALUE pairs.')
        name = specification[0]
        if not name.strip():
            raise ValueError('--compare NAME must not be empty.')
        overrides = tuple(
            (parameter, parse(value))
            for parameter, value in zip(specification[1::2], specification[2::2]))
        comparisons.append(ComparisonSpec(name, overrides))
    return comparisons


def _checkpoint_name(checkpoint) -> str:
    identifier = str(checkpoint.id)
    return identifier if identifier.casefold().startswith('cp-') else f'cp-{identifier}'


def _model_variant(checkpoint) -> str:
    """Return a readable variant while retaining unknown model identifiers."""
    import plot_best

    try:
        variant = str(plot_best.model_variant(checkpoint))
        return str(checkpoint.model_id) if variant == '?' else variant
    except Exception:
        return str(checkpoint.model_id)


def checkpoint_attributes(checkpoint, comparison_name=None) -> dict[str, str]:
    """Return human-readable checkpoint metadata for titles and legends."""
    return {
        'checkpoint': _checkpoint_name(checkpoint),
        'model': comparison_name or _model_variant(checkpoint),
        'replica': checkpoint.replica_name,
        'fit': checkpoint.fit_mode,
        'loss': checkpoint.loss_id,
    }


def comparison_labels(results: Iterable[Any]) -> tuple[str, list[str]]:
    """Return a common title fragment and minimal, unique curve labels.

    Checkpoint ID is always retained in the legend.  Model, replica, fitting
    mode, and objective appear there only when they differ; shared attributes
    are promoted to the title instead.  This prevents a comparison legend
    from repeating the same provenance on every line.
    """
    attributes = [
        checkpoint_attributes(
            result.checkpoint, result.comparison_name)
        if isinstance(result, CheckpointPlot) else checkpoint_attributes(result)
        for result in results
    ]
    if not attributes:
        return '', []

    names = ('model', 'replica', 'fit', 'loss')
    differing = {
        name for name in names
        if len({attribute[name] for attribute in attributes}) > 1
    }
    common = [
        f'{name} {attributes[0][name]}'
        for name in names if name not in differing
    ]
    labels = []
    for attribute in attributes:
        parts = [attribute['checkpoint']]
        parts.extend(
            f'{name} {attribute[name]}' for name in names if name in differing)
        labels.append(' | '.join(parts))
    return ' | '.join(common), labels


def load_checkpoint(identifier, results_directory=RESULTS_DIRECTORY) -> CheckpointPlot:
    """Resolve and replay one checkpoint without changing it on disk."""
    import plot_best

    source = resolve_checkpoint(identifier, results_directory)
    run_config, initial_parameters = replay_configuration(source)
    model_run = run(run_config, initial_parameters)
    pathway_model, objective, _ = model_run
    contents = _requested_checkpoint(source)
    effective_config = copy.deepcopy(run_config)
    effective_config['model'] = pathway_model.get_config(only_structure=True)
    metadata = {'replicas': contents.get('replicas', [])} \
        if isinstance(contents, dict) else {}
    checkpoint = plot_best.Checkpoint(
        parameters=pathway_model.parameters(),
        loss=float(objective.last_call()[-1]),
        run_config=effective_config,
        id=checkpoint_id(source) or source.name,
        source=str(source),
        metadata=metadata,
    )
    return CheckpointPlot(checkpoint, model_run)


def load_checkpoints(identifiers, results_directory=RESULTS_DIRECTORY):
    """Resolve and replay checkpoint identifiers in the requested order."""
    return [load_checkpoint(identifier, results_directory) for identifier in identifiers]


def _expanded_overrides(checkpoint, specification: ComparisonSpec) -> dict[str, Any]:
    """Expand override patterns against parameters actually present in a checkpoint."""
    available = set(checkpoint.parameters.get_values())
    overrides = {}
    for pattern, value in specification.overrides:
        matches = sorted(name for name in available if fnmatchcase(name, pattern))
        if not matches:
            raise ValueError(
                f'--compare {specification.name!r}: parameter pattern {pattern!r} '
                f'matches no parameters in {_checkpoint_name(checkpoint)}.')
        overrides.update({name: value for name in matches})
    return overrides


def replay_comparison(result: CheckpointPlot,
                      specification: ComparisonSpec) -> CheckpointPlot:
    """Replay one checkpoint with named parameter overrides.

    Thermodynamics switches also have a structural model setting.  Keeping it
    in sync with a ``*_thermodynamics=False`` override makes model D-style
    counterfactuals faithfully disable the multiplier rather than merely
    changing stored parameter metadata.
    """
    import plot_best

    source = result.checkpoint
    overrides = _expanded_overrides(source, specification)
    run_config = copy.deepcopy(source.run_config)
    parameter_values = parameters.ModelParameters(source.parameters.get_config())
    for name, value in overrides.items():
        parameter_values[name].constant(value)
        if name.endswith('_thermodynamics') and value is False:
            pathway = name.removesuffix('_thermodynamics')
            if pathway in run_config['model']:
                run_config['model'][pathway]['use_thermodynamics'] = False
    run_config['chosen'].setdefault('parameter_override', {}).update(overrides)

    model_run = run(run_config, parameter_values)
    pathway_model, objective, _ = model_run
    effective_config = copy.deepcopy(run_config)
    effective_config['model'] = pathway_model.get_config(only_structure=True)
    checkpoint = plot_best.Checkpoint(
        parameters=pathway_model.parameters(),
        loss=float(objective.last_call()[-1]),
        run_config=effective_config,
        id=source.id,
        source=source.source,
        metadata=source.metadata,
    )
    return CheckpointPlot(checkpoint, model_run, specification.name)


def add_comparisons(results: Iterable[CheckpointPlot],
                    specifications: Iterable[ComparisonSpec]) -> list[CheckpointPlot]:
    """Return original results followed by every requested replay of each one."""
    originals = list(results)
    return originals + [
        replay_comparison(result, specification)
        for specification in specifications for result in originals
    ]


def _comparison_suptitle(results, subject: str) -> str:
    common, _ = comparison_labels(results)
    return subject if not common else f'{subject} | {common}'


def _legend(axis, labels) -> None:
    """Add the result legend only when there is more than one curve."""
    if len(labels) > 1:
        handles, existing_labels = axis.get_legend_handles_labels()
        label_set = set(labels)
        result_handles = [
            handle for handle, existing_label in zip(handles, existing_labels)
            if existing_label in label_set
        ]
        result_labels = [
            existing_label for existing_label in existing_labels
            if existing_label in label_set
        ]
        axis.legend(result_handles, result_labels,
                    fancybox=False, edgecolor='k', loc='best')


def comparison_fit(results: list[CheckpointPlot], dataset):
    """Draw one shared CO2/CH4 fit figure for all checkpoint results."""
    import plot_best

    _, labels = comparison_labels(results)
    figure, axes = plot_best.create_fit_figure()
    colours = plt.rcParams['axes.prop_cycle'].by_key()['color']
    measurements_drawn = set()
    samples = []

    for index, (result, label) in enumerate(zip(results, labels)):
        checkpoint = result.checkpoint
        pathway_model, _, run_log = result.model_run
        del pathway_model
        sample, _, _ = plot_best.fit_replicas(checkpoint, dataset)
        samples.append(sample)
        measurement_key = (checkpoint.replica, checkpoint.fit_mode)
        if measurement_key not in measurements_drawn:
            plot_best.plot_fit_measurements(checkpoint, dataset, axes)
            measurements_drawn.add(measurement_key)
        colour = colours[index % len(colours)]
        for pool, axis in axes.items():
            time, values = run_log[pool]
            axis.plot(time, values, color=colour, linewidth=1.5, label=label)

    # Use the broadest data limits when checkpoints use several samples.
    for pool, axis in axes.items():
        limits = [plot_best.sample_y_limits(sample).get(pool) for sample in samples]
        limits = [limit for limit in limits if limit is not None]
        if limits:
            axis.set_ylim(min(limit[0] for limit in limits),
                          max(limit[1] for limit in limits))
        if any('log' in result.checkpoint.run_config['objective']['transform'][pool]
               for result in results):
            axis.set_yscale('log')
        plot.format_ax(axis, log_scale=axis.get_yscale() == 'log')
    _legend(axes['CH4'], labels)
    plot_best.finish_fit_figure(
        figure, _comparison_suptitle(results, 'model fit'))
    return figure, axes


def _comparison_quantities(results: list[CheckpointPlot]):
    """Return figures for every ordinary diagnostic logged by any result."""
    import plot_best

    _, labels = comparison_labels(results)
    quantities = {}
    for index, result in enumerate(results):
        pathway_model, _, run_log = result.model_run
        vmax = plot_best.pathway_v_max_values(run_log)
        for name, time, values, pathway in plot_best.pathway_plot_quantities(
                run_log, include_pathway=True, pathway_model=pathway_model):
            capacity = (pathway if name.startswith('M_') and pathway is not None
                        and plot_best.pathway_uses_capacity_parameter(
                            run_log, pathway) else None)
            key = f'{capacity}_{name}' if capacity is not None else name
            quantities.setdefault(key, []).append(
                (index, name, time, values, vmax.get(pathway), capacity))

    figures = {}
    colours = plt.rcParams['axes.prop_cycle'].by_key()['color']
    for key, curves in quantities.items():
        _, first_name, _, _, _, capacity = curves[0]
        figure, axes = plot.create_figure(plot.FigureSpec(
            axis_names=('diagnostic',), figsize=(4, 4)))
        axis = axes['diagnostic']
        for index, name, time, values, vmax, _ in curves:
            values = np.asarray(values, dtype=float)
            if name.endswith('_deltaG_r'):
                values = values * plot_best.GIBBS_ENERGY_SCALE
            colour = colours[index % len(colours)]
            axis.plot(time, values, color=colour, linewidth=1.5,
                      label=labels[index])
            if vmax is not None and name.endswith('_v'):
                axis.axhline(vmax, color=colour, linestyle='--', linewidth=1)
        if first_name.endswith('_deltaG_r'):
            axis.axhline(plot_best.GIBBS_MINIMUM * plot_best.GIBBS_ENERGY_SCALE,
                         color='k', linestyle=':', linewidth=1,
                         label=r'$\mathrm{Gibbs\ minimum}$')
        axis.set_title(plot_best.quantity_title(first_name, capacity))
        axis.set_xlabel(plot.label('time', 'd'))
        axis.set_ylabel(plot_best.quantity_ylabel(first_name, capacity))
        if first_name.endswith(('_thermodynamic_factor', '_MM', '_inhib')):
            axis.set_ylim(0, 1)
        plot.format_ax(axis)
        _legend(axis, labels)
        figure.suptitle(_comparison_suptitle(results, key))
        figure.tight_layout(rect=(0, 0, 1, .88))
        figures[key] = figure
    return figures


def _comparison_reaction_rates(results: list[CheckpointPlot]):
    """Return one comparison figure per biomass-scaled pathway rate."""
    import plot_best

    _, labels = comparison_labels(results)
    rates = {}
    for index, result in enumerate(results):
        _, _, run_log = result.model_run
        vmax = plot_best.pathway_v_max_values(run_log)
        for pathway, time, actual, maximum in plot_best.pathway_reaction_rate_curves(
                run_log, vmax):
            rates.setdefault(pathway, []).append((index, time, actual, maximum))

    figures = {}
    colours = plt.rcParams['axes.prop_cycle'].by_key()['color']
    for pathway, curves in rates.items():
        figure, axes = plot.create_figure(plot.FigureSpec(
            axis_names=('diagnostic',), figsize=(4, 4)))
        axis = axes['diagnostic']
        for index, time, actual, maximum in curves:
            colour = colours[index % len(colours)]
            axis.plot(time, actual, color=colour, label=labels[index])
            axis.plot(time, maximum, color=colour, linestyle='--', linewidth=1)
        axis.set_title(f'{pathway} reaction rate')
        axis.set_xlabel(plot.label('time', 'd'))
        axis.set_ylabel(plot.label('reaction rate', r'\mu mol/g dry weight/d'))
        plot.format_ax(axis)
        _legend(axis, labels)
        figure.suptitle(_comparison_suptitle(results, f'{pathway} reaction rate'))
        figure.tight_layout(rect=(0, 0, 1, .88))
        figures[f'{pathway}_reaction_rate'] = figure
    return figures


def _comparison_limitation_factors(results: list[CheckpointPlot]):
    """Return one all-factor comparison for every pathway logged by a result."""
    import plot_best

    _, labels = comparison_labels(results)
    series = {}
    for index, result in enumerate(results):
        _, _, run_log = result.model_run
        keys = set(run_log.keys())
        for pathway in plot_best.PATHWAY_PLOT_ORDER:
            if f'{pathway}_v' not in keys:
                continue
            factors = []
            for suffix, style in (('_MM', '-'), ('_thermodynamic_factor', '--'),
                                  ('_inhib', '-.')):
                key = f'{pathway}{suffix}'
                if key in keys:
                    factors.append((suffix, style, run_log[key]))
            if factors:
                series.setdefault(pathway, []).append((index, factors))

    figures = {}
    colours = plt.rcParams['axes.prop_cycle'].by_key()['color']
    factor_handles = [
        Line2D([], [], color='k', linestyle='-', label='Michaelis-Menten'),
        Line2D([], [], color='k', linestyle='--', label='thermodynamic'),
        Line2D([], [], color='k', linestyle='-.', label='inhibition'),
        Line2D([], [], color='k', linewidth=2, label='all factors'),
    ]
    for pathway, result_factors in series.items():
        figure, axes = plot.create_figure(plot.FigureSpec(
            axis_names=('diagnostic',), figsize=(4, 4)))
        axis = axes['diagnostic']
        checkpoint_handles = []
        for index, factors in result_factors:
            colour = colours[index % len(colours)]
            active = []
            for _, style, values in factors:
                time, factor = values
                axis.plot(time, factor, color=colour, linestyle=style, linewidth=1.2)
                active.append(values)
            time, total = plot_best._factor_product(active)
            axis.plot(time, total, color=colour, linewidth=2)
            checkpoint_handles.append(Line2D([], [], color=colour, linewidth=2,
                                             label=labels[index]))
        axis.axhline(1, color='0.65', linestyle=':', linewidth=1)
        axis.set_title(f'{pathway} limitation factors')
        axis.set_xlabel(plot.label('time', 'd'))
        axis.set_ylabel('factor')
        axis.set_ylim(0, 1.05)
        plot.format_ax(axis)
        handles = factor_handles if len(labels) == 1 else checkpoint_handles + factor_handles
        axis.legend(handles=handles, fancybox=False, edgecolor='k', loc='best')
        figure.suptitle(_comparison_suptitle(results, f'{pathway} limitation'))
        figure.tight_layout(rect=(0, 0, 1, .88))
        figures[f'{pathway}_limitation_factors'] = figure
    return figures


def comparison_figures(results: list[CheckpointPlot], dataset):
    """Build the shared figures that can be compared curve-for-curve."""
    if not results:
        raise ValueError('At least one checkpoint is required.')
    figures = {'fit': comparison_fit(results, dataset)}
    figures.update({name: (figure, None)
                    for name, figure in _comparison_quantities(results).items()})
    figures.update({name: (figure, None)
                    for name, figure in _comparison_reaction_rates(results).items()})
    figures.update({name: (figure, None)
                    for name, figure in _comparison_limitation_factors(results).items()})
    return figures


def save_comparison(results: list[CheckpointPlot], output_directory, dataset=None):
    """Save the aligned fit and diagnostic figures for several checkpoints."""
    output_directory = prepare_output_directory(output_directory)
    dataset = data.get_data_before_carex() if dataset is None else dataset
    figures = comparison_figures(results, dataset)
    paths = []
    for number, (name, (figure, _)) in enumerate(figures.items()):
        path = Path(output_directory) / f'{number:02d}_{name}.svg'
        figure.savefig(path, dpi=300)
        plt.close(figure)
        paths.append(path)
    return paths


def save_single(result: CheckpointPlot, output_directory, dataset=None, *,
                target=None, prepare=True):
    """Save the complete established ``plot_best`` report for one result."""
    import plot_best

    output_directory = (prepare_output_directory(output_directory)
                        if prepare else Path(output_directory))
    dataset = data.get_data_before_carex() if dataset is None else dataset
    checkpoint, model_run = result.checkpoint, result.model_run
    figure, _ = plot_best.plot_fit(checkpoint, dataset, model_run=model_run)
    if target is None:
        def target(number, quantity=None):
            stem = 'fit' if quantity is None else quantity
            return Path(output_directory) / f'{number:02d}_{stem}.svg'

    fit_path = Path(target(0, None))
    figure.savefig(fit_path, dpi=300)
    plt.close(figure)

    sample = dataset[checkpoint.replica[0]]
    fitted_replicas = sample.get_split(
        checkpoint.replica[1], checkpoint.fit_mode)['fit']
    suptitle = plot_best.figure_suptitle(
        checkpoint, fitted_replicas, getattr(model_run[2], 'replica', None))

    return [fit_path, *plot_best._save_diagnostic_plots(
        checkpoint, model_run, str(output_directory), target, suptitle)]


def save_plots(results: list[CheckpointPlot], output_directory, dataset=None):
    """Save all standard plots, combined when several checkpoints are given."""
    if len(results) == 1:
        return save_single(results[0], output_directory, dataset)
    return save_comparison(results, output_directory, dataset)


def default_output_directory(checkpoints, comparisons=()) -> Path:
    """Return a readable, non-project-root default below ``PLOT_DIRECTORY``."""
    def component(identifier):
        name = Path(str(identifier)).name
        return ''.join(
            character if character.isalnum() or character in '-._' else '-'
            for character in name).strip('-') or 'checkpoint'

    names = [component(checkpoint) for checkpoint in checkpoints]
    names.extend(component(comparison.name) for comparison in comparisons)
    category = 'checkpoint' if len(names) == 1 else 'checkpoint-comparison'
    return Path(PLOT_DIRECTORY) / category / '__'.join(names)


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description='Plot one or more saved checkpoints; comparisons share figures.')
    parser.add_argument('checkpoints', nargs='+', metavar='CHECKPOINT',
                        help='checkpoint path or unique full/partial cp-ID')
    parser.add_argument('--output-dir', metavar='DIRECTORY',
                        help='new or empty directory for SVG figures; defaults '
                             'below USER_VARIABLES.PLOT_DIRECTORY')
    parser.add_argument('--results-directory', default=RESULTS_DIRECTORY,
                        help='checkpoint root used for partial IDs')
    add_comparison_argument(parser)
    args = parser.parse_args(argv)
    args.comparisons = parse_comparisons(args.compare)
    if args.output_dir is None:
        args.output_dir = default_output_directory(args.checkpoints, args.comparisons)
    return args


def main(argv=None):
    args = parse_args(argv)
    results = load_checkpoints(args.checkpoints, args.results_directory)
    results = add_comparisons(results, args.comparisons)
    paths = save_plots(results, args.output_dir)
    print(f'Wrote {len(paths)} plot(s) to {args.output_dir}')
    return paths


if __name__ == '__main__':
    main()
