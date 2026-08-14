
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams.update(mpl.rcParamsDefault)
#plt.rcParams['text.usetex'] = True
import matplotlib.ticker as ticker

from chemistry import GIBBS_MINIMUM as DGmin


CO2_COLOR = 'tab:blue'
CH4_COLOR = 'tab:orange'

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


def plot_data(x, y, ax = None, **kwargs):
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

def plot_fit(log, plot_CO2 = True, plot_CH4 = True, ax = None, log_co2 = False, log_ch4 = True):
    if ax is None:
        fig, ax = plt.subplots()
    
    axs = []
    if plot_CO2:
        t, data = log.replica.CO2()
        plot_data(t, data, ax = ax, color = CO2_COLOR)
        r2 = log.R2('CO2', log_fit = log_co2)
        plot_model(*log['CO2'], ax = ax, color = CO2_COLOR, label = f'R² = {r2:.2f}')
        design(ax)
        axs.append(ax)
        
        ax.set_ylabel(f'CO2 [μmol]')

    if plot_CH4:
        if plot_CO2:
            ax = ax.twinx()
        t, data = log.replica.CH4()
        ax = plot_data(t, data, ax = ax, color = CH4_COLOR)
        r2 = log.R2('CH4', log_fit = log_ch4)
        plot_model(*log['CH4'], ax = ax, color = CH4_COLOR, label = f'R² = {r2:.2f}')

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
    
    
def plot_thermodynamics(log, pathway, ax = None):
    if ax is None:
        fig, ax = plt.subplots()
    
    t, MM = log[f'{pathway}_MM']
    t, f = log[f'{pathway}_thermodynamic_factor']
    
    plot_model(t,MM, ax = ax)
    plot_model(t,f, ax = ax)
    design(ax)
    
    return ax