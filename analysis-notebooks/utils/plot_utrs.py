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


# the import doesn't work for some reason
from utils.plot_arcs import COLOR_PC, COLOR_NMD




def plot_utrs(ax, arc_df, xlim):
    """connect annotated clv to the below ax where predicted clvs are"""
    gene = get_property(arc_df, 'gene_name')
    arc_df = arc_df.sort_values(['sc_t', 'alen'], ascending=[False, True])

    num_utrs = arc_df.shape[0]
    # height per transcript: stop codon height, scaled to between [0, 1]
    box_h = 1 / (num_utrs)
    # last exon height, the rest are used to separate neighbouring utrs
    if gene == 'CDKN2A':
        sc_h = box_h * 0.6
    else:
        sc_h = box_h * 0.7
    # utr height
    utr_h = sc_h / 4
    # stop codon/last exon width
    sc_w = (xlim[-1] - xlim[0]) / 50
    strand = get_property(arc_df, 'strand')

    for k, (key, row) in enumerate(arc_df.iterrows()):
        sc = row.sc_t
        clv = row.aclv_t

        if strand == '+':
            sc_x = sc - sc_w
            utr_x = sc
        else:
            sc_x = sc
            utr_x = clv

        # this is equivalent to one side padding for simplicity
        sc_y = box_h * k
        utr_y = (sc_h - utr_h) / 2 + sc_y
        utr_w = abs(clv - sc)

        if row.source == 'nonsense_mediated_decay':
            color = COLOR_NMD
        else:
            color = COLOR_PC

        # plot last exon
        # print(sc_y, utr_y)
        add_rect(ax, sc_x, sc_y, sc_w, sc_h, color)
        # plot utr
        add_rect(ax, utr_x, utr_y, utr_w, utr_h, color)
    clean_ax(ax, xlim)


def add_rect(ax, x, y, width, height, color):
    # ref on rect plotting http://matthiaseisen.com/pp/patterns/p0203/
    rect = mpatches.Rectangle((x, y), width, height, color=color)
    # rect = mpatches.Rectangle((x, y), width, height, facecolor=color, edgecolor='gray', linewidth=0.3)
    ax.add_patch(rect)


def clean_ax(ax, xlim):
    ax.set_xlim(xlim)
    ax.set_ylim([-0.05, 1])
    ax.axis('off')


    # for i in clvs:
    #     # -0.25 is fragile, but seems to be good for now. It's affect by hspace
    #     # when specifying the grid and ymin, but the exact relationship is unclear, not sure
    #     # what unit does hspace use
    #     ax.plot([i, i], [-0.25, 0], ':', color='gray', clip_on=False)

