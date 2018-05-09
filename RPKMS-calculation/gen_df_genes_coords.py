import pandas as pd

import logging
logging.basicConfig(
    level=logging.DEBUG, format='%(asctime)s|%(levelname)s|%(message)s')


def read_target_genes(target_genes_file):
    logging.info('reading {0}'.format(target_genes_file))
    res = []
    with open(target_genes_file) as inf:
        for line in inf:
            res.append(line.strip())
    return res


def read_gtf(input_gtf, target_genes):
    logging.info('reading {0}'.format(input_gtf))
    # http://uswest.ensembl.org/info/website/upload/gff.html
    cols = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand',
            'frame', 'attribute']
    df_gtf = pd.read_csv(input_gtf, header=None, sep='\t', comment='#',
                         low_memory=False, names=cols)
    df_gtf = df_gtf[-df_gtf.seqname.str.contains('PATCH')]
    gene_names = df_gtf.attribute.str.extract(
        r'gene_name\ \"(?P<gene_name>.+?)\";', expand=False)
    df_gtf['gene_name'] = gene_names
    df_gtf.drop('attribute', axis=1, inplace=True)
    df_gtf = df_gtf[df_gtf.gene_name.isin(target_genes)]
    return df_gtf


def gen_gene_coords(input_gtf, target_genes):
    df_gtf = read_gtf(input_gtf, target_genes)

    logging.info('generating gene coordinates')
    df_gene_coords = df_gtf.groupby('gene_name').apply(
        lambda g: pd.Series([
            g.seqname.unique()[0],
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
    target_genes_file = './reference_data/target_genes.txt'
    target_genes = read_target_genes(target_genes_file)

    gtf = './reference_data/Homo_sapiens.GRCh37.75.gtf'
    df_gene_coords = gen_gene_coords(gtf, target_genes)

    df_gene_coords.sort_values('gene_name').to_csv(
        'results/df_gene_coords.csv', index=False)
