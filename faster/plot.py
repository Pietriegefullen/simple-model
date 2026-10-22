"""Reusable, declarative plotting helpers.

The plotting interface has three small, composable parts:

* :class:`FigureSpec` describes a layout and names its primary axes.
* :class:`PlotCall` assigns a plotting routine to one of those axes (or its
  optional right-hand secondary axis).
* :func:`render_figure` creates the figure and executes the calls with one
  shared :class:`PlotContext`.  The context can be built directly from a
  saved checkpoint with :func:`replay_context`.

For example::

    context = replay_context("ABCD")
    plan = FigurePlan(
        FigureSpec(ncols=2, axis_names=("CO2", "CH4")),
        calls=(
            PlotCall(measured_pool, "CO2", {"pool": "CO2"}),
            PlotCall(modelled_pool, "CO2", {"pool": "CO2"}),
            PlotCall(modelled_pool, "CO2", {"pool": "CH4"}, secondary=True),
            PlotCall(measured_pool, "CH4", {"pool": "CH4"}),
            PlotCall(modelled_pool, "CH4", {"pool": "CH4"}),
        ),
    )
    figure, axes = render_figure(plan, context)

A routine used by ``PlotCall`` has the explicit signature
``routine(context, axis, **options)``.  This keeps layout, data selection,
and drawing separate, while allowing a panel to contain any number of curves.
``subplot_mosaic`` remains available through ``FigureSpec(mosaic=...)``.
``plot_fit`` remains as a compatibility wrapper for the one legacy caller.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Callable, Mapping, Sequence

import numpy as np

import matplotlib.pyplot as plt

CO2_COLOR = 'tab:blue'
CH4_COLOR = 'tab:orange'
POOL_NAMES = ('CO2', 'CH4')
GIBBS_ENERGY_SCALE = 1e-3
GIBBS_ENERGY_UNIT = r'\mathrm{kJ}\,\mathrm{mol}^{-1}'
# Model identity is always encoded by line style in variant-comparison plots.
# Pathway identity is encoded by ``pathways.pathway_color``.
MODEL_VARIANT_LINESTYLES = {
    'A': '-',
    'B': '--',
    'C': ':',
    'D': '-.',
}


def configure_latex() -> None:
    """Use LaTeX for all Matplotlib text without resetting user styling.

    Plain strings are rendered by LaTeX in text mode.  Use :func:`math` only
    for the mathematical portions of a label, for example
    ``label('time', 'd')`` gives ``time [$\\mathrm{d}$]``.
    """
    plt.rcParams['text.usetex'] = True


def math(expression: str) -> str:
    """Return a LaTeX math fragment suitable for a Matplotlib label."""
    return rf'${expression}$'


def label(text: str, unit: str | None = None) -> str:
    """Build a text-mode label with an optional LaTeX-rendered unit."""
    if unit is None:
        return text
    unit_text = rf'\mathrm{{{unit}}}'
    return f'{text} [{math(unit_text)}]'


configure_latex()


@dataclass(frozen=True)
class FigureSpec:
    """Configuration for a figure and its named axes.

    ``figsize`` is in inches.  ``box_aspect`` is the height/width ratio of
    every axes box (use it for visual shape, rather than ``set_aspect``, which
    changes the data-coordinate scaling).  Give either a rectangular grid via
    ``nrows``/``ncols`` or a Matplotlib ``subplot_mosaic`` via ``mosaic``.
    """

    nrows: int = 1
    ncols: int = 1
    axis_names: tuple[str, ...] | None = None
    mosaic: str | Sequence[Sequence[str]] | None = None
    figsize: tuple[float, float] | None = None
    sharex: bool | str = False
    sharey: bool | str = False
    box_aspect: float | None = None
    # ``None`` keeps the builder neutral; report-specific code can choose
    # ``tight_layout`` after it has added titles, labels, and legends.
    layout: str | None = None
    subplot_kw: Mapping[str, Any] = field(default_factory=dict)
    gridspec_kw: Mapping[str, Any] = field(default_factory=dict)


@dataclass(frozen=True)
class PlotContext:
    """The data shared by all routines in a rendered figure.

    ``checkpoint`` deliberately remains opaque: it may be a checkpoint path,
    its short identifier, or a loaded checkpoint object.  Plot routines only
    rely on the data they need, normally ``replica`` and ``run_log``.
    """

    checkpoint: Any = None
    replica: Any = None
    run_log: Any = None
    model_run: tuple[Any, ...] | None = None


@dataclass(frozen=True)
class PlotCall:
    """Add one plotting routine to a named primary or secondary axis."""

    routine: Callable[..., Any]
    axis: str
    options: Mapping[str, Any] = field(default_factory=dict)
    secondary: bool = False


@dataclass(frozen=True)
class FigurePlan:
    """A figure layout together with the routines that populate it."""

    spec: FigureSpec
    calls: tuple[PlotCall, ...] = ()
    suptitle: str | None = None


@dataclass(frozen=True)
class PlotBook:
    """Named figures that share one plotting context and checkpoint replay."""

    figures: Mapping[str, FigurePlan]


def create_figure(spec: FigureSpec) -> tuple[plt.Figure, dict[str, plt.Axes]]:
    """Create a configured figure and return its axes by stable names.

    Naming axes avoids positional coupling at call sites.  A two-panel figure
    can therefore be changed from a row to a column without changing plotting
    functions or the code that selects ``axes['CO2']``.
    """
    if spec.figsize is None:
        if spec.mosaic is None:
            rows, columns = spec.nrows, spec.ncols
        elif isinstance(spec.mosaic, str):
            mosaic_rows = spec.mosaic.split(';') if ';' in spec.mosaic \
                else spec.mosaic.splitlines()
            rows = len(mosaic_rows)
            columns = max(map(len, mosaic_rows))
        else:
            rows, columns = np.asarray(spec.mosaic, dtype=object).shape
        figsize = (4 * columns, 4 * rows)
    else:
        figsize = spec.figsize

    if spec.mosaic is not None:
        if spec.axis_names is not None:
            raise ValueError('axis_names cannot be used with a mosaic; use its labels.')
        figure = plt.figure(figsize=figsize, layout=spec.layout)
        axes = figure.subplot_mosaic(
            spec.mosaic,
            sharex=spec.sharex,
            sharey=spec.sharey,
            subplot_kw=dict(spec.subplot_kw),
            gridspec_kw=dict(spec.gridspec_kw),
        )
    else:
        if spec.nrows < 1 or spec.ncols < 1:
            raise ValueError('nrows and ncols must be positive.')
        count = spec.nrows * spec.ncols
        names = spec.axis_names or tuple(f'ax{i}' for i in range(count))
        if len(names) != count:
            raise ValueError(
                f'Expected {count} axis names for a {spec.nrows}x{spec.ncols} grid, '
                f'got {len(names)}.')
        if len(set(names)) != len(names):
            raise ValueError('axis_names must be unique.')
        figure, grid = plt.subplots(
            spec.nrows,
            spec.ncols,
            figsize=figsize,
            sharex=spec.sharex,
            sharey=spec.sharey,
            squeeze=False,
            layout=spec.layout,
            subplot_kw=dict(spec.subplot_kw),
            gridspec_kw=dict(spec.gridspec_kw),
        )
        axes = dict(zip(names, grid.flat))

    if spec.box_aspect is not None:
        for axis in axes.values():
            axis.set_box_aspect(spec.box_aspect)
    return figure, dict(axes)


def render_figure(plan: FigurePlan, context: PlotContext | None = None,
                  *, secondary_axis_offset: float | None = None
                  ) -> tuple[plt.Figure, dict[str, plt.Axes]]:
    """Render a declarative figure plan and return all named axes.

    Primary axes retain their names from ``plan.spec``.  A secondary axis is
    stored as ``"<name>_right"``.  Multiple calls marked ``secondary=True``
    share that same axis, so secondary data and curves naturally use one
    scale and one legend.

    ``secondary_axis_offset`` can move every right spine outwards, which is
    useful when a narrow panel has a long secondary label.
    """
    figure, axes = create_figure(plan.spec)
    context = PlotContext() if context is None else context
    secondary_axes: dict[str, plt.Axes] = {}

    for call in plan.calls:
        if call.axis not in axes:
            available = ', '.join(axes)
            raise KeyError(
                f'PlotCall refers to unknown axis {call.axis!r}; '
                f'available axes are: {available}.')
        target = axes[call.axis]
        if call.secondary:
            target = secondary_axes.get(call.axis)
            if target is None:
                secondary_name = f'{call.axis}_right'
                if secondary_name in axes:
                    raise ValueError(
                        f'Cannot create the secondary axis {secondary_name!r}; '
                        'that name is already used by the figure layout.')
                target = axes[call.axis].twinx()
                if secondary_axis_offset is not None:
                    target.spines['right'].set_position(
                        ('outward', secondary_axis_offset))
                secondary_axes[call.axis] = target
                axes[secondary_name] = target
        call.routine(context, target, **dict(call.options))

    if plan.suptitle is not None:
        figure.suptitle(plan.suptitle)
    return figure, axes


def replay_context(checkpoint: Any, *, parameter_overrides=None,
                   t_start=None, t_end=None) -> PlotContext:
    """Replay a checkpoint once and return a context for one or more figures.

    The expensive model evaluation happens here, not once per ``PlotCall``.
    Passing the returned context to several :func:`render_figure` calls lets a
    report reuse the exact same model trajectory.  Imports are intentionally
    lazy so normal plotting does not depend on checkpoint tooling.
    """
    from checkpoint_io import replay_inputs
    from fit_sample import run

    _, run_config, parameter_values, _ = replay_inputs(checkpoint)
    overrides = parameter_overrides or {}
    run_config['chosen'].setdefault('parameter_override', {}).update(overrides)
    if t_start is not None:
        run_config['chosen']['t_start'] = t_start
    if t_end is not None:
        run_config['chosen']['t_end'] = t_end
    model_run = run(run_config, parameter_values)
    return PlotContext(
        checkpoint=checkpoint,
        replica=getattr(model_run[-1], 'replica', None),
        run_log=model_run[-1],
        model_run=model_run,
    )


def render_book(book: PlotBook, context: PlotContext | None = None
                ) -> dict[str, tuple[plt.Figure, dict[str, plt.Axes]]]:
    """Render every named figure in a book with the same data context."""
    return {
        name: render_figure(plan, context)
        for name, plan in book.figures.items()
    }


def render_checkpoint(checkpoint: Any, book: PlotBook, *,
                      parameter_overrides=None, t_start=None, t_end=None
                      ) -> dict[str, tuple[plt.Figure, dict[str, plt.Axes]]]:
    """Replay ``checkpoint`` once, then render every figure in ``book``.

    This is the high-level entry point for a report.  It deliberately returns
    figures rather than saving them: output formats, filenames, and ownership
    stay with the caller.
    """
    context = replay_context(
        checkpoint,
        parameter_overrides=parameter_overrides,
        t_start=t_start,
        t_end=t_end,
    )
    return render_book(book, context)


def _pool_axes(ax: object, separate: bool) -> tuple[plt.Figure, dict[str, plt.Axes]]:
    """Normalise legacy axes inputs to the named-axes interface."""
    if ax is None:
        names = POOL_NAMES if separate else ('pool',)
        figure, axes = create_figure(FigureSpec(
            ncols=2 if separate else 1,
            axis_names=names,
            figsize=(8, 4),
        ))
        if separate:
            return figure, axes
        return figure, {pool: axes['pool'] for pool in POOL_NAMES}
    if isinstance(ax, Mapping):
        missing = set(POOL_NAMES) - set(ax)
        if missing:
            raise ValueError(f'Missing pool axes: {sorted(missing)}.')
        return next(iter(ax.values())).figure, dict(ax)
    if isinstance(ax, (list, tuple, np.ndarray)):
        flat_axes = np.asarray(ax, dtype=object).flat
        axes = list(flat_axes)
        if len(axes) != 2:
            raise ValueError('Expected exactly two axes for CO2 and CH4.')
        return axes[0].figure, dict(zip(POOL_NAMES, axes))
    if isinstance(ax, plt.Axes):
        return ax.figure, {pool: ax for pool in POOL_NAMES}
    raise TypeError('ax must be an Axes, a two-axis sequence, or a pool-axis mapping.')

def pool_color(pool):
    if pool == 'CO2':
        c = CO2_COLOR
    elif pool == 'CH4': 
        c = CH4_COLOR
    else:
        raise NotImplementedError()
    return c

def format_ax(ax, log_scale=False):
    """Apply shared axes styling and an optional safe logarithmic y scale.

    A logarithmic scale is meaningful only when every drawn value is
    non-negative and at least one is positive.  Its limits are set from the
    plotted data, replacing any zero lower bound with a small positive one.
    """
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    if not log_scale:
        return False
    lower, upper = ax.dataLim.intervaly
    if not (np.isfinite(lower) and np.isfinite(upper) and lower >= 0 and upper > 0):
        return False
    positive_values = [
        np.asarray(line.get_ydata(), dtype=float)
        for line in ax.lines
    ]
    positive_values = [
        values[np.isfinite(values) & (values > 0)]
        for values in positive_values
    ]
    positive_values = [values for values in positive_values if values.size]
    positive_lower = (min(np.min(values) for values in positive_values)
                      if positive_values else upper * 1e-6)
    ax.set_yscale('log')
    ax.set_ylim(positive_lower / 1.05, upper * 1.05)
    return True

def plot_pool_data(replica, pool, ax, **kwargs):
    """Add one replica's measured pool to ``ax`` and return that same axes."""
    if pool not in POOL_NAMES:
        raise ValueError(f'Unknown measured pool: {pool}.')
    time, values = getattr(replica, pool)()
    options = {'marker': 'x', 'linestyle': 'None'}
    options.update(kwargs)
    ax.plot(time, values, color=pool_color(pool), label='incubation data',
            clip_on=False, **options)
    ax.set_title(rf'$\mathrm{{{pool[:-1]}_{pool[-1]}}}$')
    ax.set_ylabel(r'$\mathrm{substance\ [\mu mol/g\ dry\ weight]}$')
    ax.set_xlabel(r'$\mathrm{time\ [d]}$')
    return ax


def plot_replica_pools(replica, axes: Mapping[str, plt.Axes], **kwargs):
    """Add CO2 and CH4 measurements to explicitly supplied named axes."""
    for pool in POOL_NAMES:
        plot_pool_data(replica, pool, axes[pool], **kwargs)
    return axes


def plot_run_pool(run_log, pool, ax, **kwargs):
    """Add one modeled pool from a run log to ``ax`` and return it."""
    if pool not in POOL_NAMES:
        raise ValueError(f'Unknown modeled pool: {pool}.')
    time, values = run_log[pool]
    ax.plot(time, values, '-', color=pool_color(pool), **kwargs)
    return ax


def plot_run_pools(run_log, axes: Mapping[str, plt.Axes], **kwargs):
    """Add modeled CO2 and CH4 curves to explicitly supplied named axes."""
    for pool in POOL_NAMES:
        plot_run_pool(run_log, pool, axes[pool], **kwargs)
    return axes


def measured_pool(context: PlotContext, axis: plt.Axes, *, pool: str,
                  replica=None, **kwargs):
    """A :class:`PlotCall` routine for one measured CO2 or CH4 pool.

    ``replica`` is optional to support a figure that intentionally combines
    observations from several replicas while keeping one shared context.
    """
    source = context.replica if replica is None else replica
    if source is None:
        raise ValueError('measured_pool needs context.replica or replica=...')
    return plot_pool_data(source, pool, axis, **kwargs)


def modelled_pool(context: PlotContext, axis: plt.Axes, *, pool: str,
                  run_log=None, **kwargs):
    """A :class:`PlotCall` routine for one modelled CO2 or CH4 pool."""
    source = context.run_log if run_log is None else run_log
    if source is None:
        raise ValueError('modelled_pool needs context.run_log or run_log=...')
    return plot_run_pool(source, pool, axis, **kwargs)


def log_quantity(context: PlotContext, axis: plt.Axes, *, quantity: str,
                 run_log=None, **kwargs):
    """Plot any ``(time, values)`` quantity logged by a model run.

    This deliberately performs no semantic formatting.  It is the compact
    building block for diagnostic figures which put several arbitrary logged
    quantities on the same primary or secondary axis.
    """
    source = context.run_log if run_log is None else run_log
    if source is None:
        raise ValueError('log_quantity needs context.run_log or run_log=...')
    time, values = source[quantity]
    return axis.plot(time, values, **kwargs)


def plot_fit(run_log, ax=None, separate=True):
    """Compatibility wrapper around :func:`plot_run_pools`."""
    figure, axes = _pool_axes(ax, separate)
    plot_run_pools(run_log, axes)
    return figure, axes
