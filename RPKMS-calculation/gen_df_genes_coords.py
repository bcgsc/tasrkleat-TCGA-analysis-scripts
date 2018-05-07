from collections import Counter

import pandas as pd


def read_target_genes(target_genes_file):
    res = []
    with open(target_genes_file) as inf:
        for line in inf:
            res.append(line.strip())
    return res


def read_gtf(input_gtf, target_genes):
    df_gtf = pd.read_csv(input_gtf, compression='gzip', sep='\t', header=None,
                         names=['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes'])
    gene_names = df_gtf.attributes.str.extract(
        r'gene_id\ \"(?P<gene_name>.+?)\";', expand=False)
    df_gtf['gene_name'] = gene_names
    df_gtf.drop('attributes', axis=1, inplace=True)
    df_gtf = df_gtf[df_gtf.gene_name.isin(target_genes)]
    return df_gtf


def gen_gene_coords(input_gtf, target_genes):
    df_gtf = read_gtf(input_gtf, target_genes)

    df_gene_coords = df_gtf.groupby('gene_name').apply(
        lambda g: pd.Series([
            g.seqid.unique()[0],
            g.strand.unique()[0],
            g.start.min(),
            g.end.max()]))
    df_gene_coords.columns = ['chromosome', 'strand', 'start', 'end']
    df_gene_coords.reset_index(inplace=True)
    df_gene_coords['len'] = abs(df_gene_coords.end - df_gene_coords.start)

    # sanity check
    assert (df_gene_coords.query('gene_name == "BRCA1"').strand == '-').all()
    assert (df_gene_coords.query('gene_name == "BRCA1"').start == 41196312).all()
    assert (df_gene_coords.query('gene_name == "BRCA1"').end == 41277500).all()

    return df_gene_coords


if __name__ == "__main__":
    target_genes_file = '/projects/btl2/zxue/tasrkleat/test-results-gmap/utrtargets/targets-for-tasrkleat/target_genes.txt'
    target_genes = read_target_genes(target_genes_file)

    gtf = '/projects/btl2/zxue/tasrkleat/test-results-gmap/reference/KLEAT-2.5.0/ensembl.fixed.sorted.gtf.gz'
    df_gene_coords = gen_gene_coords(gtf, target_genes)

    # df_gene_coords.to_csv('/dev/shm/df_gene_coords.csv', index=False)
    df_gene_coords.to_csv('./df_gene_coords.csv', index=False)
