"""Plot one or more saved checkpoints, including aligned comparisons.

The command accepts full or partial checkpoint identifiers.  A single
checkpoint keeps the established report.  With several checkpoints, each result is added to
the corresponding fit and diagnostic figure instead of producing a separate
set of files.
"""

from __future__ import annotations

import argparse
import copy
import math
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
import checkpoint_report as report
import optimizer
from checkpoint_io import checkpoint_id, resolve_checkpoint, replay_inputs
from fit_sample import run
from USER_VARIABLES import PROJECT_DIRECTORY, RESULTS_DIRECTORY, PLOT_DIRECTORY
from user_input import parse

from CONSTANTS import GIBBS_MINIMUM


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


def _checkpoint_directory_name(checkpoint) -> str:
    """Return a checkpoint output-folder name that makes fit quality visible."""
    try:
        loss = float(checkpoint.loss)
    except (AttributeError, TypeError, ValueError):
        # Retain a useful output name for lightweight callers that provide
        # checkpoint identity but not the numerical objective value.
        return _checkpoint_name(checkpoint)
    return f'{_checkpoint_name(checkpoint)}_loss-{loss:.8g}'


def _model_variant(checkpoint) -> str:
    """Return a readable variant while retaining unknown model identifiers."""
    try:
        variant = str(report.model_variant(checkpoint))
        return str(checkpoint.model_id) if variant == '?' else variant
    except Exception:
        return str(checkpoint.model_id)


def checkpoint_attributes(checkpoint, comparison_name=None) -> dict[str, str]:
    """Return human-readable checkpoint metadata for titles and legends."""
    sample, separator, replica = checkpoint.replica_name.partition('-')
    if not separator:
        replica = checkpoint.replica_name
    return {
        'checkpoint': _checkpoint_name(checkpoint),
        'model': comparison_name or _model_variant(checkpoint),
        'sample': sample,
        'replica': replica,
        'fit': checkpoint.fit_mode,
        'loss': checkpoint.loss_id,
    }


def plotted_checkpoint(identifier, source, metadata, run_config, model_run,
                       comparison_name=None) -> CheckpointPlot:
    """Pair one completed model run with the metadata needed by report figures."""
    pathway_model, objective, _ = model_run
    effective_config = copy.deepcopy(run_config)
    effective_config['model'] = pathway_model.get_config(only_structure=True)
    checkpoint = report.Checkpoint(
        parameters=pathway_model.parameters(),
        loss=float(objective.last_call()[-1]),
        run_config=effective_config,
        id=identifier,
        source=str(source),
        metadata=metadata,
    )
    return CheckpointPlot(checkpoint, model_run, comparison_name)


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

    names = ('model', 'sample', 'replica', 'fit', 'loss')
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
    source = resolve_checkpoint(identifier, results_directory)
    _, run_config, initial_parameters, contents = replay_inputs(source)
    model_run = run(run_config, initial_parameters, log_prediction=True)
    metadata = {'replicas': contents.get('replicas', [])} \
        if isinstance(contents, dict) else {}
    return plotted_checkpoint(
        checkpoint_id(source) or source.name, source, metadata, run_config, model_run)


def load_checkpoints(identifiers, results_directory=RESULTS_DIRECTORY):
    """Resolve and replay checkpoint identifiers in the requested order."""
    return [load_checkpoint(identifier, results_directory) for identifier in identifiers]


def replay_saved_checkpoint(checkpoint) -> CheckpointPlot:
    """Evaluate one checkpoint already selected from the checkpoint catalogue."""
    return CheckpointPlot(checkpoint, run(
        checkpoint.run_config, checkpoint.parameters, log_prediction=True))


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

    model_run = run(run_config, parameter_values, log_prediction=True)
    return plotted_checkpoint(
        source.id, source.source, source.metadata, run_config, model_run,
        specification.name)


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


def _legend(axis, labels, handles=()) -> None:
    """Add a legend for measurements and every displayed result curve."""
    existing_handles, existing_labels = axis.get_legend_handles_labels()
    label_set = set(labels)
    result_handles = [
        handle for handle, existing_label in zip(existing_handles, existing_labels)
        if existing_label in label_set
    ]
    result_labels = [
        existing_label for existing_label in existing_labels
        if existing_label in label_set
    ]
    handles = [*handles, *result_handles]
    if handles:
        axis.legend(handles, [handle.get_label() for handle in handles],
                    fancybox=False, edgecolor='k', loc='best')


def comparison_fit(results: list[CheckpointPlot], dataset):
    """Draw one shared CO2/CH4 fit figure for all checkpoint results."""
    _, labels = comparison_labels(results)
    figure, axes = report.create_fit_figure()
    colours = plt.rcParams['axes.prop_cycle'].by_key()['color']
    measurements_drawn = set()
    measurement_handles = []
    samples = []

    for index, (result, label) in enumerate(zip(results, labels)):
        checkpoint = result.checkpoint
        pathway_model, _, run_log = result.model_run
        del pathway_model
        sample, _, _ = report.fit_replicas(checkpoint, dataset)
        samples.append(sample)
        measurement_key = (checkpoint.replica, checkpoint.fit_mode)
        if measurement_key not in measurements_drawn:
            measurement_handles.extend(report.plot_fit_measurements(checkpoint, dataset, axes))
            measurements_drawn.add(measurement_key)
        colour = colours[index % len(colours)]
        for pool, axis in axes.items():
            time, values = run_log[pool]
            axis.plot(time, values, color=colour, linewidth=1.5, label=label)

    # Use the broadest measured-data limits when checkpoints use several samples.
    y_limits = {}
    for pool, axis in axes.items():
        limits = [report.sample_y_limits(sample).get(pool) for sample in samples]
        limits = [limit for limit in limits if limit is not None]
        if limits:
            y_limits[pool] = (min(limit[0] for limit in limits),
                              max(limit[1] for limit in limits))
    report.format_fit_axes(results[0].checkpoint, samples[0], axes,
                           y_limits=y_limits)
    for axis in axes.values():
        _legend(axis, labels, measurement_handles)
    report.finish_fit_figure(
        figure, _comparison_suptitle(results, 'model fit'))
    return figure, axes


def comparison_fit_pool(results: list[CheckpointPlot], dataset, pool: str):
    """Draw one gas axis from the standard fit figure as its own figure."""
    if pool not in ('CO2', 'CH4'):
        raise KeyError(f'Unknown fit pool: {pool}')
    _, labels = comparison_labels(results)
    figure, axes = plot.create_figure(plot.FigureSpec(
        axis_names=(pool,), figsize=(5, 5)))
    axis = axes[pool]
    colours = plt.rcParams['axes.prop_cycle'].by_key()['color']
    measurements_drawn = set()
    measurement_handles = []
    samples = []

    for index, (result, label) in enumerate(zip(results, labels)):
        checkpoint = result.checkpoint
        _, _, run_log = result.model_run
        sample, _, _ = report.fit_replicas(checkpoint, dataset)
        samples.append(sample)
        measurement_key = (checkpoint.replica, checkpoint.fit_mode)
        if measurement_key not in measurements_drawn:
            measurement_handles.extend(report.plot_fit_measurements(
                checkpoint, dataset, axes, pools=(pool,)))
            measurements_drawn.add(measurement_key)
        time, values = run_log[pool]
        axis.plot(time, values, color=colours[index % len(colours)],
                  linewidth=1.5, label=label)

    limits = [report.sample_y_limits(sample).get(pool) for sample in samples]
    limits = [limit for limit in limits if limit is not None]
    y_limits = ({pool: (min(limit[0] for limit in limits),
                        max(limit[1] for limit in limits))} if limits else {})
    report.format_fit_axes(results[0].checkpoint, samples[0], axes, y_limits=y_limits)
    _legend(axis, labels, measurement_handles)
    report.finish_fit_figure(
        figure, _comparison_suptitle(results, f'model fit {pool}'))
    return figure, axes


def _comparison_quantities(results: list[CheckpointPlot], only=None):
    """Return figures for every ordinary diagnostic logged by any result."""
    if not report.plot_enabled('quantities'):
        return {}
    _, labels = comparison_labels(results)
    quantities = {}
    for index, result in enumerate(results):
        pathway_model, _, run_log = result.model_run
        vmax = report.pathway_v_max_values(run_log)
        for name, time, values, pathway in report.pathway_plot_quantities(
                run_log, include_pathway=True, pathway_model=pathway_model):
            if not report.plot_enabled('quantities', pathway):
                continue
            if name in report.MODELLED_GAS_POOLS and not report.plot_enabled('gas_pools'):
                continue
            capacity = (pathway if name.startswith('M_') and pathway is not None
                        and report.pathway_uses_capacity_parameter(
                            run_log, pathway) else None)
            if (capacity is not None and
                    not report.plot_enabled('pathway_capacity', capacity)):
                continue
            key = f'{capacity}_{name}' if capacity is not None else name
            quantities.setdefault(key, []).append(
                (index, name, time, values, vmax.get(pathway), capacity))

    figures = {}
    colours = plt.rcParams['axes.prop_cycle'].by_key()['color']
    for key, curves in quantities.items():
        if only is not None and key != only:
            continue
        _, first_name, _, _, _, capacity = curves[0]
        figure, axes = plot.create_figure(plot.FigureSpec(
            axis_names=('diagnostic',), figsize=(4, 4)))
        axis = axes['diagnostic']
        for index, name, time, values, vmax, _ in curves:
            values = np.asarray(values, dtype=float)
            if name.endswith('_deltaG_r'):
                values = values * report.GIBBS_ENERGY_SCALE
            colour = colours[index % len(colours)]
            axis.plot(time, values, color=colour, linewidth=1.5,
                      label=labels[index])
            if vmax is not None and name.endswith('_v'):
                axis.axhline(vmax, color=colour, linestyle='--', linewidth=1)
        if first_name.endswith('_deltaG_r'):
            axis.axhline(GIBBS_MINIMUM * report.GIBBS_ENERGY_SCALE,
                         color='k', linestyle=':', linewidth=1,
                         label=r'$\mathrm{Gibbs\ minimum}$')
        axis.set_title(report.quantity_title(first_name, capacity))
        axis.set_xlabel(plot.label('time', 'd'))
        axis.set_ylabel(report.quantity_ylabel(first_name, capacity))
        if first_name.endswith(('_thermodynamic_factor', '_MM', '_inhib')):
            axis.set_ylim(0, 1)
        plot.format_ax(axis)
        _legend(axis, labels)
        figure.suptitle(_comparison_suptitle(results, key))
        figure.tight_layout(rect=(0, 0, 1, .88))
        figures[key] = figure
    return figures


def _comparison_reaction_rate_signatures(results: list[CheckpointPlot], only=None):
    """Return linear and logarithmic FPCA-style pathway-rate signatures."""
    if not report.plot_enabled('reaction_rate_signature'):
        return {}
    _, labels = comparison_labels(results)
    curves = []
    for index, result in enumerate(results):
        _, _, run_log = result.model_run
        vmax = report.pathway_v_max_values(run_log)
        curves.extend((index, pathway, time, actual)
                      for pathway, time, actual, _ in
                      report.pathway_reaction_rate_curves(run_log, vmax)
                      if report.plot_enabled('reaction_rate_signature', pathway))
    if not curves:
        return {}

    styles = ('-', '--', ':', '-.')
    figures = {}
    for log_scale in (False, True):
        key = 'reaction_rate_signature' + ('_log' if log_scale else '')
        if only is not None and key != only:
            continue
        figure, axes = plot.create_figure(plot.FigureSpec(
            axis_names=('diagnostic',), figsize=(8, 4.5)))
        axis = axes['diagnostic']
        for index, pathway, time, actual in curves:
            axis.plot(time, actual, color=plot.pathway_color(pathway),
                      linestyle=styles[index % len(styles)], linewidth=1.5)
        axis.set_title(r'$M(t)\,v(t)$ pathway-activity signature')
        axis.set_xlabel(plot.label('time', 'd'))
        axis.set_ylabel(plot.label('reaction rate', plot.REACTION_RATE_UNIT))
        plot.format_ax(axis, log_scale=log_scale)
        pathway_handles = [
            Line2D([], [], color=plot.pathway_color(pathway), label=report.pathway_display_name(pathway))
            for pathway in report.PATHWAY_PLOT_ORDER
            if any(curve_pathway == pathway for _, curve_pathway, _, _ in curves)
        ]
        handles = pathway_handles
        if len(labels) > 1:
            handles += [
                Line2D([], [], color='k', linestyle=styles[index % len(styles)], label=label)
                for index, label in enumerate(labels)
            ]
        axis.legend(handles=handles, fancybox=False, edgecolor='k', loc='best',
                    fontsize='small')
        figure.suptitle(_comparison_suptitle(results, 'pathway-activity signature'))
        figure.tight_layout(rect=(0, 0, 1, .88))
        figures[key] = figure
    return figures


def _comparison_limitation_factors(results: list[CheckpointPlot], only=None):
    """Return one all-factor comparison for every pathway logged by a result."""
    if not report.plot_enabled('limitation_factors'):
        return {}
    _, labels = comparison_labels(results)
    series = {}
    for index, result in enumerate(results):
        _, _, run_log = result.model_run
        keys = set(run_log.keys())
        for pathway in report.PATHWAY_PLOT_ORDER:
            if (not report.plot_enabled('limitation_factors', pathway) or
                    f'{pathway}_v' not in keys):
                continue
            factors = []
            for suffix, style in (('_MM', '--'), ('_thermodynamic_factor', '-.'),
                                  ('_inhib', ':')):
                key = f'{pathway}{suffix}'
                if key in keys:
                    factors.append((suffix, style, run_log[key]))
            if factors:
                series.setdefault(pathway, []).append((index, factors))

    figures = {}
    colours = plt.rcParams['axes.prop_cycle'].by_key()['color']
    factor_handles = [
        Line2D([], [], color='k', linestyle='--', label=r'$f_{\mathrm{MM}}$'),
        Line2D([], [], color='k', linestyle='-.', label=r'$f_{\mathrm{G}}$'),
        Line2D([], [], color='k', linestyle=':', label=r'$f_{\mathrm{I}}$'),
        Line2D([], [], color='k', linestyle='-', linewidth=2, label='all factors'),
    ]
    for pathway, result_factors in series.items():
        key = f'{pathway}_limitation_factors'
        if only is not None and key != only:
            continue
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
            time, total = report._factor_product(active)
            axis.plot(time, total, color='k', linestyle='-', linewidth=2)
            checkpoint_handles.append(Line2D([], [], color='k', linestyle='-', linewidth=2,
                                             label=labels[index]))
        axis.axhline(1, color='0.65', linestyle=':', linewidth=1)
        axis.set_title(rf'${plot.pathway_name(pathway)}\ \mathrm{{limitation\ factors}}$')
        axis.set_xlabel(plot.label('time', 'd'))
        axis.set_ylabel('factor')
        axis.set_ylim(0, 1.05)
        plot.format_ax(axis)
        handles = factor_handles if len(labels) == 1 else checkpoint_handles + factor_handles
        axis.legend(handles=handles, fancybox=False, edgecolor='k', loc='best')
        figure.suptitle(_comparison_suptitle(results, f'{pathway} limitation'))
        figure.tight_layout(rect=(0, 0, 1, .88))
        figures[key] = figure
    return figures


def comparison_figures(results: list[CheckpointPlot], dataset):
    """Build the shared figures that can be compared curve-for-curve."""
    if not results:
        raise ValueError('At least one checkpoint is required.')
    figures = ({'fit': comparison_fit(results, dataset)}
               if report.plot_enabled('fit') else {})
    figures.update({name: (figure, None)
                    for name, figure in _comparison_quantities(results).items()})
    figures.update({name: (figure, None)
                    for name, figure in _comparison_reaction_rate_signatures(results).items()})
    figures.update({name: (figure, None)
                    for name, figure in _comparison_limitation_factors(results).items()})
    samples = [report.fit_replicas(result.checkpoint, dataset)[0]
               for result in results]
    x_limits = plot.time_limits([
        time for sample in samples for replica in sample.replicas
        for time, _ in (replica.CO2(),)
    ])
    for figure, _ in figures.values():
        for axis in figure.axes:
            plot.format_time_axis(axis, x_limits)
    return figures


def comparison_plot_names(results: list[CheckpointPlot], *, split_fit=False) -> list[str]:
    """Return standard comparison plot keys without creating any figures."""
    if not results:
        raise ValueError('At least one checkpoint is required.')
    names = (['fit_CO2', 'fit_CH4'] if split_fit else ['fit']) if report.plot_enabled('fit') else []
    quantities = set()
    has_reaction_rate_signature = False
    factors = set()
    for result in results:
        pathway_model, _, run_log = result.model_run
        for quantity, _, _, pathway in report.pathway_plot_quantities(
                run_log, include_pathway=True, pathway_model=pathway_model):
            if not report.plot_enabled('quantities', pathway):
                continue
            if (quantity in report.MODELLED_GAS_POOLS and
                    not report.plot_enabled('gas_pools')):
                continue
            capacity = (pathway if quantity.startswith('M_') and pathway is not None
                        and report.pathway_uses_capacity_parameter(run_log, pathway)
                        else None)
            if (capacity is not None and
                    not report.plot_enabled('pathway_capacity', capacity)):
                continue
            quantities.add(f'{capacity}_{quantity}' if capacity is not None else quantity)
        vmax = report.pathway_v_max_values(run_log)
        has_reaction_rate_signature |= (
            report.plot_enabled('reaction_rate_signature') and any(
                report.plot_enabled('reaction_rate_signature', pathway)
                for pathway, *_ in report.pathway_reaction_rate_curves(run_log, vmax)))
        for pathway in report.PATHWAY_PLOT_ORDER:
            if (report.plot_enabled('limitation_factors', pathway) and
                    f'{pathway}_v' in run_log.keys()):
                factors.add(f'{pathway}_limitation_factors')
    names.extend(sorted(quantities))
    if has_reaction_rate_signature:
        names.extend(('reaction_rate_signature', 'reaction_rate_signature_log'))
    names.extend(sorted(factors))
    return names


def comparison_figure(results: list[CheckpointPlot], dataset, name: str):
    """Render exactly one named comparison figure.

    Model runs are supplied by the caller, so switching plots does not replay
    a checkpoint.  Diagnostic generators also accept a requested key, avoiding
    the former cost of constructing every report figure just to show one.
    """
    if name == 'fit' and report.plot_enabled('fit'):
        return comparison_fit(results, dataset)
    if name.startswith('fit_') and report.plot_enabled('fit'):
        return comparison_fit_pool(results, dataset, name.removeprefix('fit_'))
    for builder in (_comparison_quantities, _comparison_reaction_rate_signatures,
                    _comparison_limitation_factors):
        figures = builder(results, only=name)
        if name in figures:
            return figures[name], None
    raise KeyError(f'Unknown comparison plot: {name}')


def prepare_output_directory(directory, *, allow_existing=False) -> Path:
    """Create a report directory, optionally retaining its existing figures."""
    target = Path(directory)
    if not target.is_absolute():
        target = Path(PROJECT_DIRECTORY) / target
    if target.exists() and not target.is_dir():
        raise ValueError(f'Output path is not a directory: {target}')
    if target.exists() and any(target.iterdir()) and not allow_existing:
        raise ValueError(f'Output directory is not empty: {target}')
    target.mkdir(parents=True, exist_ok=True)
    return target


def save_comparison(results: list[CheckpointPlot], output_directory, dataset=None,
                    *, overwrite=False):
    """Save the aligned fit and diagnostic figures for several checkpoints."""
    output_directory = prepare_output_directory(output_directory, allow_existing=True)
    dataset = data.get_data_before_day() if dataset is None else dataset
    figures = comparison_figures(results, dataset)
    paths = []
    for number, (name, (figure, _)) in enumerate(figures.items()):
        path = Path(output_directory) / f'{number:02d}_{name}.png'
        if overwrite or not path.exists():
            figure.savefig(path, dpi=300)
            paths.append(path)
        plt.close(figure)
    return paths


def save_single(result: CheckpointPlot, output_directory, dataset=None, *,
                target=None, prepare=True, overwrite=False):
    """Save the complete standard report for one result."""

    output_directory = (prepare_output_directory(output_directory, allow_existing=True)
                        if prepare else Path(output_directory))
    dataset = data.get_data_before_day() if dataset is None else dataset
    checkpoint, model_run = result.checkpoint, result.model_run
    if target is None:
        def target(number, quantity=None):
            return Path(output_directory) / single_plot_filename(
                number, checkpoint, quantity)

    paths = []
    if report.plot_enabled('fit'):
        figure, _ = report.plot_fit(checkpoint, dataset, model_run=model_run)
        fit_path = Path(target(0, None))
        if overwrite or not fit_path.exists():
            figure.savefig(fit_path, dpi=300)
            paths.append(fit_path)
        plt.close(figure)

    sample = dataset[checkpoint.replica[0]]
    fitted_replicas = sample.get_split(
        checkpoint.replica[1], checkpoint.fit_mode)['fit']
    suptitle = report.figure_suptitle(
        checkpoint, fitted_replicas, getattr(model_run[2], 'replica', None))

    return [*paths, *report._save_diagnostic_plots(
        checkpoint, model_run, str(output_directory), target, suptitle,
        x_limits=report.sample_time_limits(sample), overwrite=overwrite)]


def single_plot_filename(number, checkpoint, quantity=None) -> str:
    """Return a sort-friendly report filename with its sample and replica."""
    stem = 'fit' if quantity is None else quantity
    replica = getattr(checkpoint, 'replica', None)
    replica_name = '-'.join(replica) if replica is not None else checkpoint.replica_name
    return f'{number:02d}_{replica_name}_{stem}.png'


def save_plots(results: list[CheckpointPlot], output_directory, dataset=None,
               *, overwrite=False):
    """Save all standard plots, combined when several checkpoints are given."""
    if len(results) == 1:
        return save_single(results[0], output_directory, dataset, overwrite=overwrite)
    return save_comparison(results, output_directory, dataset, overwrite=overwrite)


def default_output_directory(checkpoints, comparisons=()) -> Path:
    """Return a readable, non-project-root default below ``PLOT_DIRECTORY``."""
    def component(identifier):
        name = Path(str(identifier)).name
        return ''.join(
            character if character.isalnum() or character in '-._' else '-'
            for character in name).strip('-') or 'checkpoint'

    checkpoint_objects = [
        item.checkpoint if isinstance(item, CheckpointPlot) else item
        for item in checkpoints
    ]
    names = []
    for checkpoint in checkpoint_objects:
        if hasattr(checkpoint, 'id'):
            name = component(_checkpoint_directory_name(checkpoint))
        else:
            name = component(checkpoint)
        if name not in names:
            names.append(name)
    names.extend(component(comparison.name) for comparison in comparisons)
    category = 'checkpoint' if len(names) == 1 else 'checkpoint-comparison'
    directory = Path(PLOT_DIRECTORY) / category
    if checkpoint_objects and all(
            hasattr(checkpoint, attribute)
            for checkpoint in checkpoint_objects
            for attribute in ('fit_mode', 'replica_name')):
        fit_modes = {checkpoint.fit_mode for checkpoint in checkpoint_objects}
        replicas = {checkpoint.replica_name for checkpoint in checkpoint_objects}
        directory /= next(iter(fit_modes)) if len(fit_modes) == 1 else 'mixed-fit-mode'
        directory /= next(iter(replicas)) if len(replicas) == 1 else 'mixed-replica'
    return directory / '__'.join(names)


def _add_best_filters(parser) -> None:
    """Add the selectors used only with ``--best``."""
    parser.add_argument('--sample', dest='samples', nargs='+', metavar='SAMPLE')
    parser.add_argument('--model-id', dest='model_ids', nargs='+', metavar='MODEL_ID')
    parser.add_argument('--fit-mode', dest='fit_modes', nargs='+',
                        choices=('single', 'split', 'multi'), metavar='FIT_MODE')
    parser.add_argument('--replica', dest='replicas', nargs='+', metavar='REPLICA')
    parser.add_argument('--loss-id', dest='loss_ids', nargs='+', metavar='LOSS_ID')
    parser.add_argument('--dry-run', action='store_true',
                        help='list selected best checkpoints without plotting')


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description='Plot one or more saved checkpoints; comparisons share figures.')
    parser.add_argument('checkpoints', nargs='*', metavar='CHECKPOINT',
                        help='checkpoint path or unique full/partial cp-ID')
    parser.add_argument('--best', action='store_true',
                        help='plot every best checkpoint selected by the filters below')
    parser.add_argument('--best-restart', action='store_true',
                        help=('plot the best checkpoint per CMA-ES restart; where no '
                              'restarts were saved, plot the ordinary best checkpoint. '
                              'Invalid-objective penalty checkpoints are excluded'))
    parser.add_argument('--output-dir', metavar='DIRECTORY',
                        help='new or empty directory for SVG figures; defaults '
                             'below USER_VARIABLES.PLOT_DIRECTORY')
    parser.add_argument('--overwrite', action='store_true',
                        help='replace figures that already exist in the target directory')
    parser.add_argument('--results-directory', default=RESULTS_DIRECTORY,
                        help='checkpoint root used for partial IDs')
    _add_best_filters(parser)
    add_comparison_argument(parser)
    args = parser.parse_args(argv)
    if args.best and args.best_restart:
        parser.error('--best and --best-restart cannot be combined')
    if (args.best or args.best_restart) == bool(args.checkpoints):
        parser.error('give checkpoint IDs, or use --best or --best-restart, but not both')
    args.comparisons = parse_comparisons(args.compare)
    return args


def best_checkpoints(args):
    """Return the best checkpoints selected by the common command interface."""
    selected = report.filter_checkpoints(
        report.collect_all(args.results_directory),
        samples=args.samples, model_ids=args.model_ids, fit_modes=args.fit_modes,
        replicas=args.replicas, loss_ids=args.loss_ids)
    return report.choose_best(selected)


def _best_group(checkpoint):
    """Return the standard ``--best`` grouping key for a checkpoint."""
    return (checkpoint.model_id, checkpoint.fit_mode,
            checkpoint.run_id if checkpoint.is_multi_replica else checkpoint.replica)


def _restart_index(checkpoint):
    """Return a checkpoint's CMA-ES restart index, when it has one."""
    restart = (getattr(checkpoint, 'metadata', None) or {}).get('cma_restart')
    if not isinstance(restart, dict) or restart.get('index') is None:
        return None
    return str(restart['index'])


def _usable_checkpoint(checkpoint):
    """Whether a checkpoint has a replayable, non-penalty objective value."""
    try:
        loss = float(checkpoint.loss)
    except (AttributeError, TypeError, ValueError):
        return False
    return math.isfinite(loss) and loss < optimizer.INVALID_OBJECTIVE_PENALTY


def best_restart_checkpoints(checkpoints):
    """Select restart minima, falling back to ordinary best checkpoints.

    CMA-ES records a minimum for every restart, including restart zero.  For a
    standard ``--best`` group with any such records, retain the best usable
    result from every restart (and keep independent runs separate).  A group
    with no restart records retains its ordinary best result instead.
    """
    checkpoints = list(checkpoints)
    restart_groups = {
        _best_group(checkpoint) for checkpoint in checkpoints
        if _restart_index(checkpoint) is not None
    }
    usable = [checkpoint for checkpoint in checkpoints if _usable_checkpoint(checkpoint)]

    selected = report.choose_best([
        checkpoint for checkpoint in usable
        if _best_group(checkpoint) not in restart_groups
    ])
    best_by_restart = {}
    for checkpoint in usable:
        index = _restart_index(checkpoint)
        if index is None:
            continue
        # The same restart numbers recur in independently configured runs.
        key = (_best_group(checkpoint), checkpoint.run_id, index)
        prior = best_by_restart.get(key)
        if prior is None or checkpoint.loss < prior.loss:
            best_by_restart[key] = checkpoint
    return [*selected, *best_by_restart.values()]


def selected_best_checkpoints(args):
    """Return checkpoints chosen by either best-selection command mode."""
    selected = report.filter_checkpoints(
        report.collect_all(args.results_directory),
        samples=args.samples, model_ids=args.model_ids, fit_modes=args.fit_modes,
        replicas=args.replicas, loss_ids=args.loss_ids)
    return (best_restart_checkpoints(selected) if args.best_restart
            else report.choose_best(selected))


def save_best(checkpoints, args):
    """Write one report per best checkpoint, with optional named replays."""
    checkpoints = list(checkpoints)
    paths = []
    for checkpoint in checkpoints:
        print(checkpoint)
        if checkpoint.is_multi_replica:
            if args.comparisons:
                raise ValueError(
                    '--compare does not support joint multi-replica checkpoints; '
                    'select --fit-mode single or split.')
            paths.extend(report.plot_multi_checkpoint(
                checkpoint, data.get_data_before_day()))
            continue
        results = add_comparisons([replay_saved_checkpoint(checkpoint)],
                                  args.comparisons)
        default_directory = default_output_directory([checkpoint], args.comparisons)
        directory = default_directory if args.output_dir is None else args.output_dir
        if args.output_dir is not None and len(checkpoints) > 1:
            directory = Path(args.output_dir) / default_directory.name
        paths.extend(save_plots(results, directory, overwrite=args.overwrite))
    return paths


def main(argv=None):
    args = parse_args(argv)
    if args.best or args.best_restart:
        checkpoints = selected_best_checkpoints(args)
        if args.dry_run:
            for checkpoint in checkpoints:
                print(checkpoint.id, checkpoint.model_id, checkpoint.replica_name,
                      checkpoint.fit_mode)
            print(f'{len(checkpoints)} checkpoint(s) selected')
            return checkpoints
        paths = save_best(checkpoints, args)
        print(f'Wrote {len(paths)} plot(s) below {PLOT_DIRECTORY}')
        return paths
    results = add_comparisons(
        load_checkpoints(args.checkpoints, args.results_directory), args.comparisons)
    directory = args.output_dir or default_output_directory(results, args.comparisons)
    paths = save_plots(results, directory, overwrite=args.overwrite)
    print(f'Wrote {len(paths)} plot(s) to {directory}')
    return paths



if __name__ == '__main__':
    main()
