from cycler import cycler
import matplotlib
import numpy as np
import seaborn as sns


linewidth=0.4 #0.4
majorticksize=3.5
minorticksize=2.0
labelsize=9 #6
fontsize=9 #7
figsize=(3.5, 2.8) #(7, 5.6) (3.5, 2.8)

#colors = sns.color_palette(palette='RdYlGn')
#colors=sns.color_palette(palette='tab20')
#colors = sns.color_palette(palette='Spectral', n_colors=15)
# colors=[
#     '#515151', #grey
#     '#1A6FDF', #blue
#     '#D32F2F', #dark red
#     '#303F9F', #indigo
#     '#008176', #teal
#     '#eecc16', #yellow
#     "#7847a2", #purple
#     '#FF4081', #pink
#     '#FB6501', #orange
#     #'#F0746E', #orange
#     #'#455A64', #blue grey
#     '#b3b3b3', #grey
#     '#045275', #dark blue
#     '#1A6FDF', #blue
#     #'#089099', #dark cyan
#     '#7CCBA2', #light green
#     '#F6EDAA', #light yellow
#     #'#FCDE9C', #light orange
#     '#BdBeDC', #light purple
#     '#EAA5C2', #light pink
#     #'#DC3977', #dark pink
#     #'#7C1D6F', #dark purple
# ]

colors = [ 
    # origin colors
    '#515151', #grey
    '#F14040', #red
    '#1A6FDF', #blue
    #'#303F9F', #indigo
    #'#D32F2F', #dark red
    '#37AD6B', #green
    '#B177DE', #purple
    '#FEC211', #yellow
    '#999999', #grey
    '#FF4081', #hot pink
    '#FF6666', #baby pink
    '#FB6501', #orange
    '#6699CC', #light blue
    '#97D7F2', #light blue
    '#EAA5C2', #light pink
    '#B3D49D', #light green
    '#BdBeDC', #light purple
    '#F6EDAA', #light yellow
    '#CC9900', #yellow
    '#00CBCC', #cyan
    '#7D4E4E', #brown
    '#8E8E00', #olive
    '#6FB802', #light green
    '#07AEE3', #light blue
]


preamble  = r"""
              \usepackage{color}
              \usepackage[tx]{sfmath}
              \usepackage{helvet}
              \usepackage{sansmath}
           """

options={
    'font.size' : fontsize,

    'axes.labelsize': labelsize, #'medium'
    'axes.linewidth': linewidth, #0.8 ; 1.2
    'axes.prop_cycle': cycler('color',colors ),
    #'default': cycler('color', ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']),
    #'nature':['#E64B35', '#4DBBD5', '#00A087', '#3C5488', '#F39B7F', '#8491B4', '#91D1C2FF', '#DC0000', '#7E6148', '#B09C85']
    'axes.titlesize': fontsize, #large

    'boxplot.boxprops.color': '9ae1f9', #'black'
    'boxplot.boxprops.linewidth': 0, #1.0
    'boxplot.capprops.color': 'C0', #'black'
    'boxplot.flierprops.color': 'C0', #'black'
    'boxplot.flierprops.markeredgecolor': 'white', #'black'
    'boxplot.flierprops.markerfacecolor': 'auto', #'none'
    'boxplot.flierprops.markersize': 7, #6.0
    'boxplot.meanprops.color': 'C1', #'C2'
    'boxplot.meanprops.markeredgecolor': 'C1', #'C2'
    'boxplot.meanprops.markerfacecolor': 'C1', #'C2'
    'boxplot.meanprops.markersize': 7, #6.0
    'boxplot.medianprops.color': '9ae1f9',
    'boxplot.patchartist': True, #False
    'boxplot.whiskerprops.color': 'C0', #'black'

    'savefig.bbox':'tight',
    'savefig.pad_inches':0.01,

    'figure.figsize' : figsize,
    'figure.dpi': 300, #250
    'figure.edgecolor': 'white', #(1, 1, 1, 0)
    'figure.facecolor': 'white', #(1, 1, 1, 0)
    'figure.frameon': True, #True

    'font.family': ['sans-serif'], #['sans-serif']
    'font.sans-serif': ['Helvetica','Arial','Avenir','Computer Modern Sans Serif','DejaVu Sans','Bitstream Vera Sans','Lucida Grande','Verdana','Geneva','Lucid','Avant Garde','sans-serif'], 
    #'Helvetica, Computer Modern Sans Serif, DejaVu Sans, 'Bitstream Vera Sans, Lucida Grande, Verdana, Geneva, Lucid, Arial, Avant Garde, sans-serif', 
    
    #'text.usetex': True,
    'text.latex.preamble': preamble,

    'grid.alpha': 0.2, #1.0
    'grid.color': '93939c', # #b0b0b0
    'grid.linestyle': '--', #'-'
    'grid.linewidth':linewidth,

    'legend.edgecolor': '0.8', 
    'legend.fancybox': False, #True
    'legend.frameon': False, #True
    'legend.fontsize':labelsize,

    'lines.markeredgecolor': 'black', #'auto'
    'lines.markersize':3,
    'lines.linewidth':linewidth*2,

    'patch.linewidth':linewidth,
    'patch.facecolor': '000000', #'C0'
    'patch.force_edgecolor': True, #False
    'scatter.edgecolors': '000000', #'face'

    'xtick.labelsize': labelsize, #'medium'
    'xtick.direction': 'out',
    'xtick.minor.visible': False,
    'xtick.major.size': majorticksize, #3.5
    'xtick.minor.size': minorticksize,
    'xtick.major.top': False, #True
    'xtick.major.width': linewidth, #0.8
    'xtick.minor.width': linewidth,
    'xtick.minor.bottom': True, #True
    'xtick.minor.top': True, #True
    'xtick.alignment': 'center',

    'ytick.labelsize': labelsize, #'medium'
    'ytick.minor.visible' : False,
    'ytick.major.right': True, #True
    'ytick.major.size': majorticksize, #3.5  distance to the minor tick label in points
    'ytick.minor.size': minorticksize,
    'ytick.major.width': linewidth, #0.8
    'ytick.minor.width':linewidth,
    'ytick.minor.pad':3.4, #3.4
    'ytick.minor.right': True #True
} 

matplotlib.rcParams.update(options)

from matplotlib.ticker import AutoMinorLocator
from matplotlib import rcParams


def pretty_plot(width=None, height=None, plt=None, dpi=None):
    if plt is None:
        plt = matplotlib.pyplot
        if width is None:
            width = matplotlib.rcParams["figure.figsize"][0]
        if height is None:
            height = matplotlib.rcParams["figure.figsize"][1]
        if dpi is not None:
            matplotlib.rcParams["figure.dpi"] = dpi

        fig = plt.figure(figsize=(width, height))
        ax = fig.add_subplot(1, 1, 1)
        matplotlib.rcParams["xtick.minor.visible"]=True
        matplotlib.rcParams["ytick.minor.visible"]=True
        ax.yaxis.set_minor_locator(AutoMinorLocator(2))
        ax.xaxis.set_minor_locator(AutoMinorLocator(2)) 
    return plt

def pretty_subplot(
    nrows=None,
    ncols=None,
    width=None,
    height=None,
    sharex=False,
    sharey=True,
    dpi=None,
    plt=None,
    gridspec_kw=None,
    total_num=None,
    ):
    if nrows is None and ncols is None:
        if total_num == 1:
            nrows, ncols = 1, 1
        elif total_num > 1:
            nrows, ncols = 1, total_num
        elif total_num == 4:
            nrows, ncols = 2, 2
        elif total_num > 4:
            nrows, ncols = 2, total_num // 2
        elif total_num == 9:
            nrows, ncols = 3, 3
        elif total_num > 9:
            nrows, ncols = total_num // 3, 3
        else:
            raise ValueError("total_num must be a positive integer")
    elif nrows is None:
        nrows = total_num // ncols
    elif ncols is None:
        ncols = total_num // nrows
    
    if width is None:
        width = rcParams["figure.figsize"][0]*ncols
    if height is None:
        height = rcParams["figure.figsize"][1]*nrows


    plt = matplotlib.pyplot
    fig, axes = plt.subplots(
        nrows,
        ncols,
        sharex=sharex,
        sharey=sharey,
        dpi=dpi,
        figsize=(width, height),
        facecolor="w",
        gridspec_kw=gridspec_kw,
    )
    # if nrows is None and ncols is None:
    #     if total_num == 1:
    #         axes = [axes]
    
    # for ax in axes:
    #     ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    #     ax.xaxis.set_minor_locator(AutoMinorLocator(2))

    axs = axs.flatten() if isinstance(axes, np.ndarray) else [axes]

    return fig, axes

def draw_themed_line(y, ax, orientation="horizontal"):
    """Draw a horizontal line using the theme settings
    Args:
        y (float): Position of line in data coordinates
        ax (Axes): Matplotlib Axes on which line is drawn
    """

    # Note to future developers: feel free to add plenty more optional
    # arguments to this to mess with linestyle, zorder etc.
    # Just .update() the options dict

    themed_line_options = dict(
        color=rcParams["grid.color"],
        linestyle="--",
        dashes=(5, 2),
        zorder=0,
        linewidth=rcParams["ytick.major.width"],
    )

    if orientation == "horizontal":
        ax.axhline(y, **themed_line_options)
    elif orientation == "vertical":
        ax.axvline(y, **themed_line_options)
    else:
        raise ValueError(f'Line orientation "{orientation}" not supported')