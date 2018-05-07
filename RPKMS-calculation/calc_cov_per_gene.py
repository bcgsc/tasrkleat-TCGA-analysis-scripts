import os
import sys

import pandas as pd


def read_genomecov(infile):
    df = pd.read_csv(infile, sep='\t', header=None, names=['chromosome', 'start', 'end', 'coverage'])
    df = df.query('coverage > 0').copy()
    df['len'] = abs(df.end - df.start) + 1
    df['num_bases'] = df.len * df.coverage
    return df


def calc_norm_cov(grp, genomecov_df):
    _df = genomecov_df[
        (genomecov_df.chromosome == grp.chromosome.values[0]) 
        & (genomecov_df.start >= grp.start.values[0]) 
        & (genomecov_df.end <= grp.end.values[0])
    ]
    if _df.len.sum() == 0: # also means _df.shape[0] is 0
        return 0
    else:
        return _df.num_bases.sum() / _df.len.sum()


if __name__ == "__main__":
    df_gene_coords = pd.read_csv('/dev/shm/df_gene_coords.csv')

    # the genomecov file produced by bedtools genomecov
    input_file = sys.argv[1]
    output_file = os.path.join(os.path.dirname(input_file), 'genomecov-by-gene.csv')

    df = read_genomecov(input_file)
    res = df_gene_coords.groupby('gene_name').apply(
        lambda grp: calc_norm_cov(grp, df)).rename('cov_average_num_bases').reset_index()
    print('saving to {0} ...'.format(output_file))
    res.to_csv(output_file, index=False)
