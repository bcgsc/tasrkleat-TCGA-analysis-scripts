import os
import sys
sys.path.insert(0, '/projects/btl/zxue/tasrkleat-TCGA-results/apa_manuscript/')

import pandas as pd
import numpy as np
import matplotlib
import pysam

from hexamer_search import search_hexamer
from utils.cluster import cluster_clv_sites


TARGET_GENE_TSV = './reference_data/target_genes_with_type.tsv'
REF_FA = './reference_data/Homo_sapiens.GRCh37.75.dna_sm.primary_assembly.fa'

ANNOTATED_CLV_SC_MAPPING = './reference_data/annotated-clv-sc-mapping.csv.gz'
df_clv_sc = pd.read_csv(ANNOTATED_CLV_SC_MAPPING)
ANNOT_GENE2CLV_DD = df_clv_sc.groupby(['gene_name']).apply(lambda g: g.aclv.unique())


def get_target_genes():
    return pd.read_csv(TARGET_GENE_TSV, sep='\t').Gene.values.tolist()


def filter_by_target_genes(adf):
    wanted_cols = [
        'gene',
        'transcript_strand',
        'chromosome',
        'cleavage_site',
        'length_of_tail_in_contig',
        'number_of_bridge_reads',
        'max_bridge_read_tail_length',
        'tail+bridge_reads',
        'flag'
    ]

    bdf = adf[wanted_cols].copy()
    cdf = bdf.rename(columns={
        'gene': 'gene_name',
        'cleavage_site': 'clv',
        'transcript': 'transcript_id',
        'transcript_strand': 'strand',
        'chromosome': 'seqname'})

    target_genes = get_target_genes()
    ddf = cdf[cdf.gene_name.isin(target_genes)]
    return ddf


def search_hexamer_wrapper(refseq, chrm, clv, strand, window=50):
    chrm = chrm.replace('chr', '')
    res = search_hexamer.search(refseq, chrm, clv, strand, window)
    if res is None:
        res = ['NA', -1, -1]
    return pd.Series(res, index=['hexamer', 'hexamer_id', 'hexamer_loc0'])


def research_hexamer(adf):
    adf['clv0'] = adf['clv'] - 1  # convert 1-based to 0-based

    refseq = pysam.FastaFile(REF_FA)

    _cols = ['seqname', 'gene_name', 'clv0', 'strand']
    clv0s_df = adf[_cols].drop_duplicates()

    hexm_df = clv0s_df.apply(
        lambda row: search_hexamer_wrapper(
            refseq, row.seqname, row.clv0, row.strand), axis=1)

    hexm_df['hexamer_loc'] = hexm_df.hexamer_loc0.apply(
        lambda v: v + 1 if v != -1 else v)  # convert to 1-based

    hexm_ndf = pd.concat([clv0s_df, hexm_df], axis=1)

    bdf = adf.merge(
        hexm_ndf, on=['seqname', 'gene_name', 'clv0', 'strand'], how='left')

    assert bdf.shape[0] == adf.shape[0]

    cdf = bdf.drop(['clv0', 'hexamer_loc0'], axis=1)
    return cdf


def map_kclv2aclv(row):
    """map KLEAT predicted clv to a closest annotated clv"""
    gene = row.gene_name
    kclv = row['clv']
    poss_aclvs = ANNOT_GENE2CLV_DD[gene]
    return poss_aclvs[np.argmin(np.abs(poss_aclvs - kclv))]


def map_to_annotated_cs(adf):
    _df = adf[['gene_name', 'clv']].drop_duplicates()

    _df['aclv'] = _df.apply(map_kclv2aclv, axis=1)

    bdf = adf.merge(_df, on=['gene_name', 'clv'], how='left')

    assert bdf.shape[0] == adf.shape[0]

    bdf['signed_NDA'] = bdf['clv'] - bdf['aclv']
    bdf['NDA'] = bdf['signed_NDA'].abs()

    return bdf


HEXMAERS = ['AATAAA', 'ATTAAA', 'AGTAAA', 'TATAAA', 'CATAAA', 'GATAAA',
            'AATATA', 'AATACA', 'AATAGA', 'AAAAAG', 'ACTAAA', 'AAGAAA',
            'AATGAA', 'TTTAAA', 'AAAACA', 'GGGGCT']


def cluster_to_stable(adf):
    # hierarchical clustering help reveal local structure, so clustering over
    # all relevant data first. In other words, cluster-then-filter instead of
    # filter-then-cluster
    _df = adf.copy()
    for i in range(2):
        print('{0} clustering...'.format(i))
        # cluster twice to final results more stable, see the experiment below
        # NOTE: %time magic inside function may mess up with the returned value
        _df = _df.groupby('gene_name').apply(cluster_clv_sites, 20).reset_index(drop=True)
    return _df


def filter_by_confidence(adf, nda, ltc=4, nbr=2, maxbrtl=4, hexamers=HEXMAERS):
    cdf = adf[
        (
            # (adf.NDA <= 25) & (adf.flag == 1)
            adf.NDA <= nda
        )
        |
        (
            (
                (adf.length_of_tail_in_contig >= ltc) |
                (adf.number_of_bridge_reads >= nbr) |
                (adf.max_bridge_read_tail_length >= maxbrtl)
            )
            &
            (
                # adf.hexamer.isin(["AATAAA", "ATTAAA"])
                adf.hexamer.isin(hexamers)
            )
        )
    ].copy()
    return cdf


def filter_by_tbr(adf, tbr):
    cdf = adf[adf['tail+bridge_reads'] >= tbr].copy()
    return cdf


if __name__ == "__main__":
    input_kleat_output = sys.argv[1]
    filter_style = sys.argv[2]

    df = pd.read_csv(input_kleat_output, sep='\t')
    odf = filter_by_target_genes(df)
    odf = research_hexamer(odf)
    odf = map_to_annotated_cs(odf)

    if filter_style == 'A':
        # my way of filtering
        outdir = os.path.dirname(input_kleat_output) + '/postproc-styleA-polyA-confidence'
        nda, ltc, nbr, maxbrtl = [int(_) for _ in sys.argv[3:7]]
        hxm = sys.argv[7:]
        out_name = 'nda{nda}-ltc{ltc}-nbr{nbr}-maxbrtl{maxbrtl}-hxm{0}.csv'.format(len(hxm), **locals())
        filtered_df = filter_by_confidence(odf, nda, ltc, nbr, maxbrtl, hxm)
    elif filter_style == 'B':
        # CM2 way of filtering by #tail+bridge_reads
        outdir = os.path.dirname(input_kleat_output) + '/postproce-styleB-tbr-tuning'
        tbr = int(sys.argv[3])
        out_name = 'tbr_gt_{0}.csv'.format(tbr)
        filtered_df = filter_by_tbr(odf, tbr)
    else:
        raise ValueError('unknown filter style: {0}'.format(filter_style))

    filtered_df = cluster_to_stable(filtered_df)
    out_csv = os.path.join(outdir, out_name)
    filtered_df[['seqname', 'strand', 'mclv']]\
        .drop_duplicates()\
        .rename(columns={'mclv': 'clv'})\
        .to_csv(out_csv, index=False)
    print('saved to {0}'.format(out_csv))
