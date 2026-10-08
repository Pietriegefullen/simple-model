"""Reusable, axes-oriented plotting helpers.

The public boundary in this module is deliberately small:

1. describe a figure once with :class:`FigureSpec`;
2. create it with :func:`create_figure`; and
3. pass the named axes to functions that add data or model curves.

For example::

    spec = FigureSpec(ncols=2, axis_names=("CO2", "CH4"), figsize=(8, 4))
    figure, axes = create_figure(spec)
    plot_replica_pools(replica, axes, marker="o")
    plot_run_pools(run_log, axes)

``subplot_mosaic`` is available through ``FigureSpec(mosaic=...)`` when a
rectangular grid is not expressive enough.  The older ``plot_data`` and
``plot_fit`` functions remain as compatibility wrappers.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Mapping, Sequence

import numpy as np

import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams.update(mpl.rcParamsDefault)
plt.rcParams['text.usetex'] = True
import matplotlib.ticker as ticker

from chemistry import GIBBS_MINIMUM as DGmin
from pathways import pathway_color

CO2_COLOR = 'tab:blue'
CH4_COLOR = 'tab:orange'
POOL_NAMES = ('CO2', 'CH4')


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
    figsize: tuple[float, float] = (8, 4)
    sharex: bool | str = False
    sharey: bool | str = False
    box_aspect: float | None = None
    # ``None`` keeps the builder neutral; report-specific code can choose
    # ``tight_layout`` after it has added titles, labels, and legends.
    layout: str | None = None
    subplot_kw: Mapping[str, Any] = field(default_factory=dict)
    gridspec_kw: Mapping[str, Any] = field(default_factory=dict)


def create_figure(spec: FigureSpec) -> tuple[plt.Figure, dict[str, plt.Axes]]:
    """Create a configured figure and return its axes by stable names.

    Naming axes avoids positional coupling at call sites.  A two-panel figure
    can therefore be changed from a row to a column without changing plotting
    functions or the code that selects ``axes['CO2']``.
    """
    if spec.mosaic is not None:
        if spec.axis_names is not None:
            raise ValueError('axis_names cannot be used with a mosaic; use its labels.')
        figure = plt.figure(figsize=spec.figsize, layout=spec.layout)
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
            figsize=spec.figsize,
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

def get_axes(ax, separate = True):
    """Compatibility wrapper; new code should use ``create_figure``."""
    return _pool_axes(ax, separate)

def pool_color(pool):
    if pool == 'CO2':
        c = CO2_COLOR
    elif pool == 'CH4': 
        c = CH4_COLOR
    else:
        raise NotImplementedError()
    return c

def format_ax(ax, log_scale = False):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    if log_scale:
        ax.set_yscale('log')

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


def plot_data(replica, ax=None, separate=True, **kwargs):
    """Compatibility wrapper around :func:`plot_replica_pools`.

    It returns a named axes dictionary (rather than a positional tuple), which
    is accepted unchanged by later calls to ``plot_data`` and ``plot_fit``.
    """
    figure, axes = _pool_axes(ax, separate)
    plot_replica_pools(replica, axes, **kwargs)
    sample_label = str(replica.sample).replace(' ', r'\ ')
    figure.suptitle(rf'$\mathrm{{{sample_label}}}$')
    return figure, axes


def plot_fit(run_log, ax=None, separate=True):
    """Compatibility wrapper around :func:`plot_run_pools`."""
    figure, axes = _pool_axes(ax, separate)
    plot_run_pools(run_log, axes)
    return figure, axes


def design(ax):
    
    is_right = ax.yaxis.get_ticks_position() == 'right'
    if is_right:
        ax.spines['left'].set_visible(False)
        ax.spines['bottom'].set_visible(False)

    else:
        ax.spines['right'].set_visible(False)

    ax.spines['top'].set_visible(False)
    
    ax.spines['left'].set_position(('data', 0))
    
    ax.tick_params(axis="y", direction='out')
    if ax.get_yscale() == 'linear':
        ax.ticklabel_format(axis='y', style='plain')
        ax.spines['bottom'].set_position(('data', 0))
        ax.yaxis.set_major_formatter(ticker.StrMethodFormatter('{x:,.1f}'))
        
    return ax

def xaxis_time(ax):
    ax.set_xlabel('t [d]',  loc = 'right')
    ax.xaxis.set_major_formatter(ticker.StrMethodFormatter('{x:,.0f}'))
    
    maxy = None
    for line in ax.get_lines():
        y = line.get_ydata()
        if maxy is None or max(y) > maxy:
            maxy = max(y)
    if maxy <= 0:
        ax.tick_params(axis="x", direction='in', labeltop = True, labelbottom = False)
        ax.set_xlabel('t [d]',  loc = 'right', labelpad = -20)
        ax.xaxis.set_label_coords(1.1, 1.06)
        
def title(ax, log, pathway = ''):
    r = log.replica
    sample = ' '.join([r.sample.sample_name, r.sample.site, f'({r.sample.origin})'])
    pwy = ''
    if not pathway == '':
        pwy = f': {pathway} pathway'
    ax.set_title(f'{sample}, replica {r.replica_number}{pwy}')
    return ax

def plot_Gibbs(log, pathway, ax = None):
    if ax is None:
       fig, ax = plt.subplots()
    
    t, DGr = log[f'{pathway}_deltaG_r']
    
    ax.plot(t, DGr)
    
    ax.plot([min(t), max(t)], [DGmin, DGmin], 'k--')
    ax.annotate('Gibbs minimum', xy = (max(t), DGmin), ha = 'right', va = 'top',
                xytext = (0, -6), textcoords = 'offset points')
        
    ax = design(ax)
    xaxis_time(ax)
    
    ax.set_ylabel('ΔG [J/mol]', rotation = 0, loc = 'top', labelpad = -30)
    
    ax = title(ax, log, pathway)

    return ax

def legend(ax):
    handles, labels = [], []
    for axs in plt.gcf().axes:
        handles += axs.get_legend_handles_labels()[0]
        labels += axs.get_legend_handles_labels()[1]
        
    plt.gcf().axes[0].legend(handles, labels, 
                             loc = 'best', 
                             fancybox = False, 
                             edgecolor = 'k')


def plot_data2(x, y, ax = None, **kwargs):
    if ax is None:
        fig, ax = plt.subplots()
    
    ax.plot(x, y, 'x', **kwargs)
    xaxis_time(ax)
    return ax

def plot_model(x,y, ax = None, **kwargs):
    if ax is None:
        fig, ax = plt.subplots()
    
    ax.plot(x,y, '-', **kwargs)
    xaxis_time(ax)
    return ax

def plot_fit2(log, plot_CO2 = True, plot_CH4 = True, ax = None, log_co2 = False, log_ch4 = True):
    if ax is None:
        fig, ax = plt.subplots()
    
    axs = []
    if plot_CO2:
        t, data = log.replica.CO2()
        plot_data(t, data, ax = ax, color = CO2_COLOR)
        ylim = ax.get_ylim()
        r2 = log.R2('CO2', log_fit = log_co2)
        plot_model(*log['CO2'], ax = ax, color = CO2_COLOR, label = f'R² = {r2:.2f}')
        ax.set_ylim(ylim)
        design(ax)
        axs.append(ax)
        
        ax.set_ylabel(f'CO2 [μmol]')

    if plot_CH4:
        if plot_CO2:
            ax = ax.twinx()
        t, data = log.replica.CH4()
        ax = plot_data(t, data, ax = ax, color = CH4_COLOR)
        ylim = ax.get_ylim()
        r2 = log.R2('CH4', log_fit = log_ch4)
        plot_model(*log['CH4'], ax = ax, color = CH4_COLOR, label = f'R² = {r2:.2f}')
        ax.set_ylim(ylim)

        ax.set_ylabel(f'CH4 [μmol]')

        ax.set_yscale('log')
        design(ax)
        axs.append(ax)
        
        if plot_CO2:
            axs[0].tick_params(axis='y', labelcolor=CO2_COLOR)
            axs[0].yaxis.label.set_color(CO2_COLOR)
            axs[1].tick_params(axis='y', labelcolor=CH4_COLOR)
            axs[1].yaxis.label.set_color(CH4_COLOR)

    design(ax)
    title(ax, log)
    legend(ax)


def plot_pathways(log, ax = None):
    if ax is None:
        fig, ax = plt.subplots()
        
    # TODO: 
    # for each pathway get biomass, 
    # get v (not v_max)
    # compute pathway 'activity' as biomass * v
    # plot as stacked?
    stack = []
    pathways = set([name.split('_')[0] for name in log._log.keys() if '_' in name])
    pathways.remove('CO2')
    pathways.remove('CH4')
    pathways.remove('M')
    labels = []
    for pathway in pathways:
        microbe = 'M_' + pathway
        v = pathway + '_v'
        
        if microbe == 'M_Aceto':
            microbe = 'M_Ac'
        elif microbe == 'M_Fermentation' or  microbe == 'M_Hydrolysis':
            microbe = 'M_Ferm'
        
        if not v in log._log:
            continue
    
        t, biomass = log[microbe]
        t_v, _v = log[v]
        _v = np.interp(t, t_v, _v)
        stack.append(biomass*_v)
        labels.append(pathway)
        
    ax.stackplot(t, *stack, labels = labels, 
                 colors = [pathway_color(p) for p in labels])
    #ax.set_yscale('log')
    #ax.set_ylim([1e-8, 1e0])
    design(ax)
    ax.legend()
    
    return ax


def plot_thermodynamics(log, pathway, ax = None):
    if ax is None:
        fig, ax = plt.subplots()
    
    t, MM = log[f'{pathway}_MM']
    t, f = log[f'{pathway}_thermodynamic_factor']
    
    plot_model(t,MM, ax = ax)
    plot_model(t,f, ax = ax)
    design(ax)
    
    return ax
