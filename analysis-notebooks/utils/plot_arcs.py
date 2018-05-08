import pandas as pd

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import gridspec
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.path as mpath
Path = mpath.Path
import matplotlib.patches as mpatches
from matplotlib.ticker import FormatStrFormatter
from matplotlib.legend_handler import HandlerPatch
matplotlib.style.use('classic')

import numpy as np
import scipy.stats as stats
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster, fclusterdata
from scipy import interp
from sklearn.preprocessing import LabelEncoder, OneHotEncoder, StandardScaler, MinMaxScaler

from utils.calc_clv_sc import get_strand, get_num_sc, get_property


COLOR_PC = 'black'
COLOR_NMD = (0.34, 0.86, 0.5792)
# COLOR_BOTH = 'red'
# COLOR_BOTH = 'orange'
# COLOR_BOTH = '#be29ec'          # purple-like
COLOR_BOTH = '#FF7F00'          # dark orange
# COLOR_BOTH = '#F5BE41'          # lighter orange
# COLOR_BOTH = '#FF9F00'          # dark orange
COLOR_INSIGNIFICANT = '#4d4d4d'  # this appears to be a good level of gray for clv arrows


MAIN_PLOT_GD_PAIRS_SIMPLE = [
    ('FGF2', 'LUAD'),
    ('CCNE1', 'LUAD'),
    ('RNF43', 'KIRC'),
    ('RNF43', 'UCEC'),
]

MAIN_PLOT_GD_PAIRS_COMPLEX = [
    ('CDKN2A', 'KIRC'),
    ('CDKN2A', 'HNSC'),
    ('EZH2', 'LUAD'),
    ('PTCH1', 'BRCA')
]

MAIN_PLOT_GD_PAIRS = MAIN_PLOT_GD_PAIRS_SIMPLE + MAIN_PLOT_GD_PAIRS_COMPLEX


# ordered in the way to arrange panels
# MAIN_PLOT_GD_PAIRS = [
#     ('FGF2', 'LUAD'),
#     # LUSC has fewer sample numbers than LUAD, just go with LUAD
#     # ('FGF2', 'LUSC'),
#     ('RNF43', 'KIRC'),
#     ('CCNE1', 'LUAD'),
#     # ('CCNE1', 'LUSC'),
#     ('RNF43', 'UCEC'),
#     ('CDKN2A', 'KIRC'),
#     ('EZH2', 'LUAD'),
#     ('CDKN2A', 'HNSC'),
#     ('PTCH1', 'BRCA')
# ]



def process_arrow_df(idf, xlim=None, arrow_head_length=0.07, y_interval=0):
    odf = idf.copy()

    # for those diff is 0, keep it 0, otherwise, scaled to a minimum diff
    abs_diffs = odf.N2T_ratio_diff.abs().values
    # mm_scaler = MinMaxScaler(
    #     [arrow_head_length, 0.5]
    # ).fit(
    #     abs_diffs[abs_diffs > 0].reshape(-1, 1)
    # )
    # odf['abs_N2T_ratio_diff_mmscaled'] = odf.N2T_ratio_diff.abs().apply(
    #     lambda x: 0 if x == 0 else mm_scaler.transform([[x]])[0][0])
    odf['abs_N2T_ratio_diff'] = abs_diffs

    odf['color'] = odf.apply(calc_arrow_color, axis=1)
    odf['alpha'] = odf.diff_is_significant.apply(lambda v: 1 if v else 0.7)
    odf['linestyle'] = odf.diff_is_significant.apply(lambda v: '-')
    odf['head_length'] = arrow_head_length
    odf['y_interval'] = y_interval
    if xlim is not None:
        xmin, xmax = xlim
        odf['head_width'] = (xmax - xmin) * 0.018
    return odf


def calc_xlim(arc_df, arw_df):
    xmin1 = arw_df.mclv_t.min().min()
    xmax1 = arw_df.mclv_t.max().max()
    xmin2 = arc_df[['aclv_t', 'sc_t']].min().min()
    xmax2 = arc_df[['aclv_t', 'sc_t']].max().max()

    xmin = min(xmin1, xmin2)
    xmax = max(xmax1, xmax2)
    # added buffer so that arc ends won't be at the edge
    bf = (xmax - xmin) * .06
    if bf == 0:
        bf += 10                # to avoid error
    xmin -= bf
    xmax += bf
    return xmin, xmax


def calc_arc_color(source):
    if source == 'protein_coding':
        return COLOR_PC
    elif source == 'nonsense_mediated_decay':
        return COLOR_NMD
    else:
        raise


def calc_arc_linewidth(source):
    if source == 'protein_coding':
        return 1
    elif source == 'nonsense_mediated_decay':
        return 1.4
    else:
        raise


def process_arc_df(idf, ymax, y_interval=0):
    """idf should be gene-specific"""
    odf = idf.copy()
    odf['y_interval'] = y_interval
    odf['mid_point'] = (odf['sc_t'] + odf['aclv_t']) / 2

    order = odf['alen'].argsort().argsort()
    odf['order'] = order
    # take log to squash the difference in between
    # arc_height = np.log(order + 1) + 0.2
    arc_height = order + 1
    # for inverted y axis
    arc_height /= max(arc_height)
    # so that the height utilize the whole axes height
    odf['arc_height'] = arc_height * 1.97 * (ymax - y_interval)

    if odf.source.unique().shape[0] == 1:
        odf['color'] = 'black'
    else:
        odf['color'] = odf.source.apply(calc_arc_color)
    odf['linewidth'] = odf.source.apply(calc_arc_linewidth)
    return odf


def plot_arcs(ax, arc_df):
    for k, row in arc_df.iterrows():
        path = Path([
            (row.sc_t, row.y_interval),
            (row.mid_point, row.arc_height),
            (row.aclv_t, row.y_interval)
        ], [Path.MOVETO, Path.CURVE3, Path.CURVE3])

        path_patch = mpatches.PathPatch(
            path, fc='none', color=row.color, lw=row.linewidth,
            transform=ax.transData)

        ax.add_patch(path_patch)


def plot_stop_codons(ax, arc_df):
    sc_df = arc_df[['sc_t', 'source', 'y_interval']]\
            .drop_duplicates()\
            .sort_values('sc_t')

    symbols, labels = [], []
    for k, row in sc_df.sort_values('source').iterrows():
        if row.source == 'protein_coding':
            fc = COLOR_PC
            size = 20
            mk = 'o'
        elif row.source == 'nonsense_mediated_decay':
            fc = COLOR_NMD
            size = 100
            mk = 'o'
        # elif row.source == 'retained_intron':
        #     fc = 'orange'
        #     size = 20
        #     mk = '^'
        # elif row.source == 'processed_transcript':
        #     fc = 'purple'
        #     size = 20
        #     mk = 'v'
        # elif row.source == 'non_stop_decay':
        #     fc = 'pink'
        #     size = 40
        #     mk = '*'
        else:
            raise
        lbl = row.source.replace('_', ' ') + ' stop codon'
        symb = ax.scatter([row.sc_t], [row.y_interval],
                          alpha=0.7,
                          color=fc,
                          edgecolors='black',
                          linewidth=0.5,
                          marker=mk,
                          s=size,
                          zorder=10,
                          clip_on=False)
        symbols.append(symb)
        labels.append(lbl)
    return symbols, labels


def calc_arrow_color(row):
    if not row.diff_is_significant:
        return COLOR_INSIGNIFICANT
    if row.num_sc == 1:
        if 'protein_coding' in row.src_list:
            # return COLOR_PC
            if row.N2T_ratio_diff > 0:
                return 'red'
            elif row.N2T_ratio_diff < 0:
                return 'blue'
            else:
                return 'gray'
        elif 'nonsense_mediated_decay' in row.src_list:
            return COLOR_NMD
    elif row.num_sc > 1:
        return COLOR_BOTH
    else:
        raise ValueError('cannot decide arrow color')


def plot_strand(ax, xlim, strand, y_interval=0, head_width=0.08):
    xmin, xmax = xlim
    head_length = (xmax - xmin) * 0.022

    arrow_params = dict(
        head_width=head_width,
        head_length=head_length,
        color='black', lw=0.5,
    )
    if strand == '+':
        _ = ax.arrow(xmin, y_interval,
                     xmax - xmin - head_length, 0,
                     **arrow_params)
    elif strand == '-':
        _ = ax.arrow(xmax, y_interval,
                     xmin - xmax + head_length, 0,
                     **arrow_params)
    else:
        raise


def clean_arc_ax(ax):
    for _ in ['left', 'top', 'right', 'bottom']:
        ax.spines[_].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_visible(False)
    ax.spines['bottom'].set_position(('data', 0))
    # so the tick label won't be hided by the second subplot
    ax.tick_params(direction='out')


def plot_connector(ax, arc_df):
    """connect annotated clv to the below ax where predicted clvs are"""
    # gene = get_property(arc_df, 'gene_name')
    clvs = arc_df['aclv_t'].values.tolist()
    scs = arc_df['sc_t'].values.tolist()

    # lowest = -0.1 if (gene, dise) in MAIN_PLOT_GD_PAIRS_COMPLEX else -0.8
    lowest = -1
    # for i in clvs + scs:
    for i in clvs:
        # -0.25 is fragile, but seems to be good for now. It's affect by hspace
        # when specifying the grid and ymin, but the exact relationship is unclear, not sure
        # what unit does hspace use

        # or if zorders of different axes are pro, then it's fine.
        ax.plot([i, i], [lowest, 0], ':', linewidth=0.5, color='#333333', clip_on=False)


def plot_arcs_stop_codons_and_strand(ax, gd_arc_df, xlim, ylim=None):
    plot_arcs(ax, gd_arc_df)
    plot_connector(ax, gd_arc_df)
    symbols, labels = plot_stop_codons(ax, gd_arc_df)
    xmin, xmax = xlim
    strand = get_strand(gd_arc_df)
    plot_strand(ax, xlim, strand)
    clean_arc_ax(ax)

    ax.set_xlim(xlim)
    if ylim is None:
        ylim = [-0.2, 1]
    ax.set_ylim(ylim)
    ax.ticklabel_format(style='sci', axis='x', scilimits=(0, 2))
    # modifying y seems to have no effect
    # http://stackoverflow.com/questions/39620700/positioning-the-exponent-of-tick-labels-when-using-scientific-notation-in-matplo
    # ax.get_xaxis().get_offset_text().set_position((1, 0))
    ax.xaxis.offsetText.set_visible(False)


def make_legend_arrow_up(
        legend, orig_handle, xdescent, ydescent, width, height, fontsize):
    # Trying to put arrows in the legend, too. Get pretty complicated, may try
    # later, not as important, learned from
    # http://stackoverflow.com/questions/22348229/matplotlib-legend-for-an-arrow

    # print('aa', legend.get_label(), legend.get_texts(), orig_handle.get_label())
    # print(type(legend))
    # print(orig_handle.get_label())
    # print(legend, orig_handle, xdescent, ydescent)
    # print(width, height)
    p = mpatches.FancyArrow(
        0.5 * width, 0, 0, height,
        length_includes_head=True,
        head_width=0.75 * height,
        head_length=0.3 * height)
    return p
