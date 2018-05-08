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
from utils.plot_arcs import MAIN_PLOT_GD_PAIRS_SIMPLE, MAIN_PLOT_GD_PAIRS_COMPLEX


def adjust_clv_text_coord(gene, dise, ki, x, y):
    """ki: 0-based index of the clv"""

    # Be specific, use both gene and dise for conditional match

    if gene == 'FGF2' and dise == 'LUAD':
        if ki == 0:
            x += 120

    if gene == 'CDKN2A':
        if dise == "HNSC":
            if ki == 0:
                x -= 90
            if ki == 2:
                x += 90
            if ki in [3, 4]:
                y += 0.15       # to avoid clash with genomic coordinates
        elif dise == 'KIRC':
            if ki == 0:
                x -= 70
            if ki == 2:
                x += 120

    if gene == 'EZH2' and dise == "LUAD":
        if ki == 0:
            x -= 500
        elif ki == 1:
            x += 500

    if gene == 'PTCH1' and dise == "BRCA":
        if ki == 0:
            x -= 500
        if ki == 1:
            y += 0.15
        if ki == 2:
            x += 500
        # if ki == 4:
        #     x -= 200
        if ki == 5:
            x -= 100
        if ki == 6:
            x += 600
    return x, y


def label_clvs(ax, arrow_df, text_arrow_interval, supp=False):
    """text_arrow_interval: to reduce overlap between text and arrow"""
    gene = get_property(arrow_df, 'gene_name')
    dise = get_property(arrow_df, 'disease')
    # clv_prefix = 'D' if gene == 'CDKN2A' else gene[0]
    clv_prefix = gene[0]
    for ki, (_, row) in enumerate(arrow_df.iterrows()):
        delta = row.abs_N2T_ratio_diff
        rat_chg = row.N2T_ratio_change
        hdl = row.head_length

        x = row.mclv_t
        y, dy = calc_arrow_y_coord(rat_chg, delta, hdl)

        text = '{0}{1}'.format(clv_prefix, ki + 1)
        if supp:
            fontsize = 7
        else:
            fontsize = 12
        params = dict(horizontalalignment='center',
                      color=row.color,
                      alpha=row.alpha,
                      fontsize=fontsize
        )


        if rat_chg == 'up':
            y += text_arrow_interval
        elif rat_chg == 'down':
            y = y + dy + hdl + text_arrow_interval
        else:
            raise

        x, y = adjust_clv_text_coord(gene, dise, ki, x, y)
        ax.text(x, y, text, **params)


def plot_clv_arrows(ax, arrow_df, supp=False):
    """
    mmscaler: MinMaxScaler to scale arrow heights
    """
    arrow_patches, arrow_labels, sort_keys = [], [], []
    gene = get_property(arrow_df, 'gene_name')
    dise = get_property(arrow_df, 'disease')

    if (gene, dise) in MAIN_PLOT_GD_PAIRS_SIMPLE:
        text_arrow_interval = 0.22
    else:
        text_arrow_interval = 0.14

    label_clvs(ax, arrow_df, text_arrow_interval, supp)

    for ki, (_, row) in enumerate(arrow_df.iterrows()):
        delta = row.abs_N2T_ratio_diff
        begx, dx = row.mclv_t, 0
        begy, dy = calc_arrow_y_coord(
            row.N2T_ratio_change, delta, row.head_length)

        arrow_pat = ax.arrow(
            # it's actually x, y, dx, dy, refer to doc for details
            begx, begy, dx, dy,
            head_length=row.head_length,
            head_width=row.head_width,
            color=row.color,
            # ec='black',
            alpha=row.alpha,
            ls='-',
            lw=1)

        # if row.color == COLOR_INSIGNIFICANT:
        #     arrow_labels.append('insignificant frequency change')
        # elif row.color == COLOR_PC:
        #     arrow_labels.append('significant frequency change of PC')
        # elif row.color == COLOR_NMD:
        #     arrow_labels.append('significant frequency change of NMD')
        # elif row.color == COLOR_BOTH:
        #     arrow_labels.append('significant frequency change of both')
        # else:
        #     raise

        arrow_patches.append(arrow_pat)
        sort_keys.append(row.color)

    # take one of each in order
    idxes = []
    # for i in ['red', 'blue', 'green', 'yellow', 'black', 'gray']:
    for i in arrow_df.color.unique().tolist():
        try:
            idx = sort_keys.index(i)
            idxes.append(idx)
        except:
            pass
            # print('not found {0} arrows'.format(i), end=',')

    arrow_patches = [_ for k, _ in enumerate(arrow_patches) if k in idxes]
    arrow_labels = [_ for k, _ in enumerate(arrow_labels) if k in idxes]
    # if not (arrow_df.disease.unique()[0] == 'THCA' and arrow_df.gene_name.unique()[0] == 'FGF2'):
    # print(arrow_patches)
    assert len(arrow_patches) >= 1
    return arrow_patches, arrow_labels


def calc_arrow_y_coord(ratio_change, delta, arrow_head_length):
    rat_chg = ratio_change
    ahl = arrow_head_length

    # use a small epsilon to make the direction right, otherwise, it could
    # result in wrong direction due to numerical error
    epsilon = 1e-6
    if rat_chg == 'up':
        y = ahl if delta < ahl else delta
        dy = min(-epsilon, - (delta - ahl))
    elif rat_chg == 'down':
        y = 0
        dy = max(epsilon,   (delta - ahl))
    return y, dy


def clean_arrow_ax(ax):
    for _ in ['top', 'right', 'bottom']:
        ax.spines[_].set_visible(False)
    ax.xaxis.set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.yaxis.grid(True, linestyle=':', color='gray', linewidth=0.350)
    ax.set_axisbelow(True)      # to put grid below other elements


def plot_clvs(ax, gd_arw_df, xlim, supp):
    xmin, xmax = xlim
    # ax.plot([xmin, xmax], [0, 0], ':', lw=0.5, color='black')
    mclv_arrow_patches, mclv_arrow_labels = plot_clv_arrows(ax, gd_arw_df, supp)
    clean_arrow_ax(ax)

    ax.set_xlim(xmin, xmax)
    # a little negative buffer, otherwise, up arrows look truncated
    ax.set_ylim(-0.005, 0.7)
    ax.yaxis.set_ticks([0, 0.2, 0.4, 0.6])
    ax.tick_params(axis='y', labelsize=9)
    ax.invert_yaxis()
    ax.set_ylabel('$\Delta$', fontsize=13, labelpad=7, rotation=0, va='center')
    # ax.ticklabel_format(style='sci', axis='x', scilimits=(0, 2))


def assertEqual(a, b):
    try:
        assert a == b
    except AssertionError:
        print('{0} != {1}'.format(a, b))
        raise


def test_calc_arrow_y_coord():
    arw_len = 1
    assertEqual(calc_arrow_y_coord('up', 4, arw_len), (4, -3))
    assertEqual(calc_arrow_y_coord('down', 4, arw_len), (0, 3))

    # tricky cases when delta < arw_len
    arw_len = 2
    assertEqual(calc_arrow_y_coord('up', 5, arw_len), (5, -3))
    assertEqual(calc_arrow_y_coord('down', 5, arw_len), (0, 3))

    assertEqual(calc_arrow_y_coord('up', 1, arw_len)[0], 2)
    assertEqual(calc_arrow_y_coord('down', 1, arw_len)[0], 0)


test_calc_arrow_y_coord()
