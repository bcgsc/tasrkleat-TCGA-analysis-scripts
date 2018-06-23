import os
import sys

import pandas as pd
import numpy as np

sys.path.insert(0, '../apa_manuscript/')


def load_polya_df(polya_csv):
    return pd.read_csv(polya_csv)


def load_cm2_df(cm2_bed):
    df = pd.read_csv(cm2_bed, header=None, usecols=[0, 1, 5], sep='\t')
    df.columns = ['seqname', 'clv', 'strand']
    return df


def load_kleat_df(kleat_csv):
    cols = ['strand', 'seqname', 'mclv']
    df = pd.read_csv(kleat_csv, usecols=cols)
    df = df.rename(columns={'mclv': 'clv'})
    df = df.drop_duplicates().copy()
    return df


def load_raw_kleat_df(kleat_csv):
    cols = ['transcript_strand', 'chromosome', 'cleavage_site']
    df = pd.read_csv(kleat_csv, usecols=cols, sep='\t')
    df = df.rename(columns={
        'chromosome': 'seqname',
        'transcript_strand': 'strand',
        'cleavage_site': 'clv',
    })
    df = df.drop_duplicates().copy()
    return df


def map_clvs(df_pred, df_ref):
    ref_dd = df_ref.groupby(['seqname', 'strand']).apply(
        lambda g: g.clv.values).to_dict()

    _dfs = []
    for k, g in df_pred.groupby(['seqname', 'strand']):
        _df = g.copy()

        ref_clvs = ref_dd.get(k)
        if ref_clvs is None:
            _df['mapped_ref_clv'] = np.nan
        else:
            _df['mapped_ref_clv'] = g.clv.apply(
                lambda v: ref_clvs[np.argmin(np.abs(v - ref_clvs))])
        _dfs.append(_df)

    df_mapped = pd.concat(_dfs)

    df_mapped['dist'] = df_mapped.mapped_ref_clv - df_mapped.clv
    df_mapped['abs_dist'] = np.abs(df_mapped['dist'])
    return df_mapped


def compare(df_pred, df_ref, dist_cutoff=50):
    df_mapped = map_clvs(df_pred, df_ref)
    df_mapped['is_tp'] = df_mapped.abs_dist < dist_cutoff

    sensitivity = df_mapped.query('is_tp').shape[0] / df_ref.shape[0]
    precision = df_mapped.query('is_tp').shape[0] / df_mapped.shape[0]
    f1 = (2 * sensitivity * precision) / (sensitivity + precision)

    return sensitivity, precision, f1


def replace_seqname(df):
    df['seqname'] = df.seqname.str.replace('chr', '').replace('M', 'MT')


if __name__ == "__main__":
    pred_file, pred_type = sys.argv[1:]

    truth_sample = 'UHRC1'
    truth_csv = './UHR/C1/polyA-Seq/polyA-Seq-truth-114-genes.csv'

    df_ref = load_polya_df(truth_csv)

    pred_file_loader = eval(f'load_{pred_type}_df')
    df_pred = pred_file_loader(pred_file)

    replace_seqname(df_ref)
    replace_seqname(df_pred)

    se, pr, f1 = compare(df_pred, df_ref)
    print(f'compared to\tSensitivity\tPrecision\tF1')
    print(f'{truth_sample}\t{se}\t{pr}\t{f1}')

    output = f'{os.path.dirname(pred_file)}/{os.path.basename(pred_file)}.vs-polyA-Seq.csv'
    with open(output, 'wt') as opf:
        opf.write('{0} {1} {2} {3}\n'.format(pred_file, se, pr, f1))
    print('saved to {0}'.format(output))
