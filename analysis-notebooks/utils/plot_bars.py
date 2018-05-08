import re

import matplotlib
import numpy as np
from utils.calc_clv_sc import get_property


COUNTS_TEXT_X_LOC_DD = {
    # in terms of percentage of the xlim
    'FGF2': 0.08,
    'CCNE1': 0.08,
    'RNF43': 0.65,
    'EZH2': 0.6,
    'CDKN2A': 0.6,
    'PTCH1': 0.6,
    'DRAM1': 0.1,
    'GNAS': 0.1,
    'TERT': 0.6,
    'REL': 0.4,
    'WT1': 0.6,
    'CHURC1': 0.6,
    'HNF1A': 0.2,
    'MET': 0.4,
    'MITF': 0.52,
    'MDM2': 0.4,
    'KIT': 0.1,
    'CDKN2C': 0.1,
}

# NEEDS_ROTATE_GB_PAIRS = [
#     ('GNAS', 'KICH'),
#     ('GNAS', 'BRCA'),
# ]

NEEDS_ROTATE_GENES = [
    'GNAS'
]

BS_BLUE = (np.array([2, 117, 216]) / 255).tolist()
BS_GREEN = (np.array([92, 184, 89]) / 255).tolist()
BS_RED  = (np.array([217, 83, 79]) / 255).tolist()
BS_ORANGE = (np.array([240, 173, 78]) / 255).tolist()


def init_ax2(ax):
    ax2 = ax.twiny()
    ax2.xaxis.tick_bottom()
    ax2.xaxis.set_label_position('bottom')
    # http://stackoverflow.com/questions/22942984/matplotlib-bug-data-is-shown-above-spine-with-twinx
    # don't need these two lines any more as bars are set to be transparent

    # ax.set_zorder(ax2.get_zorder()+1) # put ax in front of ax2
    # ax.patch.set_visible(False) # hide the 'canvas'
    ax2.grid(False)
    return ax2


def filter_expr_df(df_expr, df_clv, gene_name, disease):
    aids = df_clv\
           .query('disease == "{0}"'.format(disease))\
           .query('gene_name == "{0}"'.format(gene_name))\
           .analysis_id.unique().tolist()

    res = df_expr.query('disease == "{0}"'.format(disease))\
                 .query('analysis_id in {0}'.format(aids))

    assert len(aids) == res.shape[0]
    return res[['analysis_id', 'sstype', gene_name]]


def plot_expr(ax, expr_df, gene_name):
    # xleft = expr_df[gene_name].min()
    xleft = 0
    xright = expr_df[gene_name].quantile(0.98)

    bins = np.linspace(xleft, xright, 25)
    for st, col, ls in zip(
            ['normal', 'tumour'],
            # [BS_BLUE, BS_RED],
            ['blue', 'red'],
            ['-', '-']
    ):
        vals = expr_df.query('sstype == "{0}"'.format(st))[gene_name].values

        hist, binedges = np.histogram(vals, bins=bins, normed=False)
        ys = hist / vals.shape[0] # convert to percentage
        xs = (binedges[:-1] + binedges[1:]) / 2

        # simulate an edge color
        lne, = ax.plot(xs, ys, color='white', linestyle=ls, lw=3.3)
        # real plot
        lne, = ax.plot(xs, ys, color=col, linestyle=ls, lw=1.5)


def rotate_xtick_labels(ax):
    for tick in ax.get_xticklabels():
        tick.set_rotation(30)


def calc_counts_box_width(count_text):
    width_dd = {
        # 3: 0.32,
        # 4: 0.36,

        # after using latex, the font width changed
        3: 0.27,
        4: 0.31,
    }
    ss = re.findall(r'\d+', count_text)
    max_num = max(len(_) for _ in ss)
    return width_dd.get(max_num, 0.28)


def gen_counts_text(df_expr):
    """returns e.g. N: 56\nT: 498"""
    sstype_counts = df_expr[['analysis_id', 'sstype']].sstype.value_counts()
    # text = 'N: {0}\nT: {1}'.format(
    #     sstype_counts.loc['normal'],
    #     sstype_counts.loc['tumour'])

    # for bold latex rendered string
    text = '\\textbf{{N: {0}}}\n\\textbf{{T: {1}}}'.format(
        sstype_counts.loc['normal'],
        sstype_counts.loc['tumour'])
    # print(text)
    return text


def add_counts_txt(ax, df_expr, gene_name=None,
                   x=0.5, y=0.735, bg_color=None):
    if bg_color is None:
        bg_color = 'white'

    if gene_name is not None:
        x = COUNTS_TEXT_X_LOC_DD.get(gene_name, x)

    text = gen_counts_text(df_expr)

    # if bg_color == BS_BLUE:
    #     text = 'Short.\n' + text  # shortening
    # elif bg_color == BS_RED:
    #     text = 'Len.\n' + text  # lengthening
    # else:
    #     text = 'Cmplx.\n' + text  # complex

    height = 0.22
    # height = 0.28; y = 0.68     # three lines, e.g. with trend
    rw, rh = calc_counts_box_width(text), height  # width & height
    # bg_color = 'white'          # overwrite bg_color with text
    rect_params = dict(
        color=bg_color,
        ec='black' if bg_color == 'white' else 'none',
        # ec='none',
        transform=ax.transAxes
    )

    rect = matplotlib.patches.Rectangle((x, y), rw, rh, **rect_params)
    ax.add_patch(rect)

    rx, ry = rect.get_xy()
    cx = rx
    cy = ry + rect.get_height()/2.0

    # + 0.035, it's like a padding-left, while also control the box width to
    # make the text look center
    text_color = 'black' if bg_color == 'white' else 'white'
    # print(cx, cy, text_color)
    ax.text(cx + 0.035, cy, text, color=text_color, weight='bold',
            ha='left', va='center', transform=ax.transAxes)

    # # transform=ax.transAxes: seem to make x as a percentage
    # # http://matplotlib.org/1.5.1/api/pyplot_api.html#matplotlib.pyplot.text
    # ax.text(x, y, text, ha='left', va='center', fontsize=15, transform=ax.transAxes)


def calc_edgecolors(adf, hcolor, pval_cutoff=0.01):
    # cc = lambda v: hcolor if v < pval_cutoff else 'grey'
    return adf['fisher_exact_p'].apply(
        lambda v: 'black' if v < pval_cutoff else 'grey')


def calc_facecolors(bar_df, hcolor, pval_cutoff=0.01):
    res =  bar_df['fisher_exact_p'].apply(
        lambda v: ['blue', 'red'] if v < pval_cutoff else ['grey', 'grey'])
    res = np.array([np.array(_) for _ in res]).T
    # e.g. np.array([['grey', 'blue', 'blue'], ['grey', 'red', 'red']])
    return res


def apply_hatches(ax, num_grp):
    """
    # make bars prettier

    grp_size is 2 (normal & tumor)
    num_grp is the same as the number of cleavage sites per gene
    """
    bars = ax.patches
    # first normal and then tumor
    # hats = '\\' * num_grp + '/' * num_grp
    # for bar, hatch in zip(bars, hats):
    #     bar.set_hatch(hatch)

    # cols = ['blue'] * num_grp + ['red'] * num_grp
    # for bar, col in zip(bars, cols):
    #     print(bar.get_edgecolor()
    #     if bar.get_edgecolor() == 'gray':
    #         bar.set_color('gray')
    #     else:
    #         bar.set_color(col)


def plot_freq(ax, bar_df, highlight_color='black', do_plot_expr=False):
    for i in ['gene_name', 'disease', 'strand']:
        assert bar_df[i].unique().shape[0] == 1

    bar_df.set_index(bar_df['mclv'], inplace=True)
    bar_df.sort_index(inplace=True)

    # highlight cleavage sites with significant shifts
    # ecs = calc_edgecolors(bar_df, highlight_color)
    facecolors = calc_facecolors(bar_df, highlight_color)
    bar_df[['N_on_ratio', 'T_on_ratio']].plot.bar(
        ax=ax,
        edgecolor='white',
        color=facecolors,
        lw=1.6,
        alpha=1,
        rot=0, legend=False)

    apply_hatches(ax, bar_df.shape[0])
    assert bar_df.gene_name.unique().shape[0] == 1
    gene_name = bar_df.gene_name.unique()[0]
    dise = bar_df.disease.unique()[0]
    if do_plot_expr:
        label_bars_at_top(ax, gene_name, dise)
    else:
        label_bars_at_bottom(ax, gene_name, dise)


def label_bars_at_bottom(ax, gene_name, dise):
    clv_prefix = gene_name[0]
    bars = ax.patches
    num_clvs = len(bars) // 2
    txts = []
    for k in range(num_clvs):
        txt = '{0}{1}'.format(clv_prefix, k + 1)
        txts.append(txt)
    print(txts)
    ax.set_xticklabels(txts)


def label_bars_at_top(ax, gene_name, dise):
    clv_prefix = gene_name[0]
    bars = ax.patches
    num_clvs = len(bars) // 2
    for k in range(num_clvs):
        # each group is two bar, normal and tumour are partitioned, so do the
        # following slicing
        grp = [bars[k], bars[k + num_clvs]]
        mx = np.mean([(_.get_x() + _.get_width() / 2) for _ in grp])
        my = np.max([_.get_height() for _ in grp]) + 0.05

        # if gene_name == 'FGF2' and dise == 'LUAD' and k == 0:
        #     my += 0.07          # to aovid hiding by expr lines
        # if gene_name == 'CDKN2A' and dise in ['KIRC', 'HNSC'] and k == 0:
        #     mx += 0.05

        txt = '{0}{1}'.format(clv_prefix, k + 1)
        gray = (0.5019607843137255, 0.5019607843137255, 0.5019607843137255, 1)
        if grp[0].get_facecolor() == gray:
            color = 'gray'
        else:
            color = 'black'
        ax.text(mx, my, txt, color=color, va='center', ha='center',
                zorder=20,
                fontsize=7
        )


def set_title(ax, bar_df, has_diffpref=False, fontsize=15):
    gene_name = bar_df.gene_name.unique()[0]
    disease = bar_df.disease.unique()[0]
    strand = bar_df.strand.values[0]
    # # unfortunately, seqname was not left in the df
    # seqname = bar_df.mkid.values[0].split('|')[0]
    # title = '{gene_name} ({seqname}, {strand}), {disease}'.format(**locals())
    title = r'\textit{{{gene_name}}}, {disease}'.format(**locals())
    # print(title)
    if has_diffpref:
        ax.set_title(title, backgroundcolor='orange', fontsize=15)
    elif (gene_name, disease) == ('FGF2', 'THCA'):
        # reference
        ax.set_title(title, color='white', fontweight='bold', backgroundcolor='black')
    else:
        # ax.set_title(title)
        ax.set_title(title, fontsize=15)


def clean_ax_bar(ax):
    ax.set_xlabel('')
    ax.tick_params(top=False, bottom=False)
    ax.yaxis.grid(True, linestyle=':', color='gray', linewidth=0.35)
    ax.set_axisbelow(True)


def plot_bar_and_expr(ax, gd_bar_df, gd_expr_df,
                      set_expr_xlabel=True, set_expr_ylabel=True,
                      bg_color_dd=None,
                      ax_bar_ylim=[0, 1],
                      do_plot_expr=False):
    assert gd_bar_df.gene_name.unique().shape[0] == 1
    gene = gd_bar_df.gene_name.unique()[0]
    assert gd_bar_df.disease.unique().shape[0] == 1
    dise = gd_bar_df.disease.unique()[0]
    if bg_color_dd is None:
        bg_color_dd = {}

    ax_bar = ax
    plot_freq(ax_bar, gd_bar_df, highlight_color='blue',
              do_plot_expr=do_plot_expr)
    clean_ax_bar(ax_bar)

    if do_plot_expr:
        ax.set_xticklabels([])
        ax_expr = init_ax2(ax_bar)
        plot_expr(ax_expr, gd_expr_df, gene)

        if set_expr_xlabel:
            ax_expr.set_xlabel('RPKMS')
        # if set_expr_ylabel:
        #     ax_expr.set_ylabel('Fraction of samples')

        if gene in NEEDS_ROTATE_GENES:
            rotate_xtick_labels(ax_expr)

        bcol = bg_color_dd.get((gene, dise))
        add_counts_txt(ax_expr, gd_expr_df, gene, bg_color=bcol)
        set_title(ax_expr, gd_bar_df)

        for _ in [ax_bar, ax_expr]:
            _.set_ylim(ax_bar_ylim)
            # if ticklabel_fontsize is not None:
            #     _.tick_params(axis='both', labelsize=ticklabel_fontsize)
        label_delta(ax)
    else:
        if gene in NEEDS_ROTATE_GENES:
            rotate_xtick_labels(ax_expr)

        bcol = bg_color_dd.get((gene, dise))
        add_counts_txt(ax_bar, gd_expr_df, gene, bg_color=bcol)
        set_title(ax_bar, gd_bar_df)

        ax_bar.set_ylim(ax_bar_ylim)

        # if gene == 'FGF2' and dise == 'LUAD':
        # if (
        #         (gene == 'FGF2' and dise == 'LUAD') or
        #         (gene == 'CCNE1' and dise == 'LUAD') or
        #         (gene == 'RNF43' and dise in ['KIRC', 'UCEC'])
        # ):
            # label_delta(ax)

        label_delta(ax)

    # not sure why previous line doesn't work: 
    # ax_expr.set_ylabel('Fraction of samples')
    ax_bar.set_ylabel('Fraction of samples')


def label_delta(ax):
    """add delta label above the lower bar"""
    bars = ax.patches
    # there is an additional patch corresponding to the count box, should be
    # ignored
    num_grp = len(bars) // 2
    if num_grp < 6:
        fontsize = 12
        delta_height = 0.1     # the height of delta label
    else:
        fontsize = 8
        delta_height = 0.072

    for i in range(num_grp):
        bar_grp = [bars[i], bars[i + num_grp]]
        heights = [_.get_height() for _ in bar_grp]
        diff = np.abs(np.diff(heights))
        if diff < delta_height * 1:
            continue
        idx = np.argmin([_.get_height() for _ in bar_grp])

        bar = bar_grp[idx]
        width = bar.get_width()
        # * 0.97: take edge line width into effect, and try to make the line
        # * look at the center as much as possible
        xloc = bar.get_x() + (width) / 2

        y1 = bar.get_height() - 0.0015
        y2 = y1 + (diff - delta_height) / 2
        y3 = y2 + delta_height
        y4 = y1 + diff - 0.006

        # y1 = bar.get_height() - 0.001
        # y2 = y1 + diff * 0.33
        # y3 = y1 + diff * 0.65
        # y4 = y1 + diff * 0.98
        ax.plot([xloc - width / 4, xloc + width / 4], [y1, y1], lw=0.6, color='k')
        ax.plot([xloc, xloc], [y1, y2], linewidth=0.6, color='k' )
        # * 1.007: otherwise, the text doesn't seem to align to the center
        if num_grp < 6:
            xloc *= 1.007
        else:
            xloc *= 1.002
        ax.text(xloc, y2 + 0.004, r'$\Delta$', ha='center', va='bottom', fontsize=fontsize)
        ax.plot([xloc, xloc], [y3, y4], linewidth=0.6, color='k', )
        ax.plot([xloc - width / 4, xloc + width / 4], [y4, y4], lw=0.6, color='k')
        # print(idx)
