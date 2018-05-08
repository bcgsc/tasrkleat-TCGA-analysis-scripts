import os

import pandas as pd
import numpy as np
import unittest


# This is a specialized module for extracting (if it's there in the gtf) or
# calculating annotated cleavage sites and stop codons from annotation GTF file

def get_strand(grp):
    return get_property(grp, 'strand')


def get_source(grp):
    return get_property(grp, 'source')


def get_property(grp, prop_name):
    """prop_name: e.g. strand, num_sc, source, gene_source, transcript_source, is_cds_end_NF"""
    try:
        vals = grp[prop_name].unique()
        assert vals.shape[0] == 1
        return vals[0]
    except:
        print('multiple {0} values found for {1}'.format(prop_name, grp))
        raise


def get_num_sc(sub_clv_sc_df):
    try:
        assert sub_clv_sc_df.gene_name.unique().shape[0] == 1
        return sub_clv_sc_df.sc.unique().shape[0]
    except:
        print(sub_clv_sc_df)
        raise


def extract_sc_per_transcript(scs):
    # 5.4% of the annotated transcripts don't have UTR
    if scs.shape[0] == 0:
        return np.nan

    strand = get_strand(scs)
    if strand == '+':
        return scs.end.max()
    else:
        return scs.start.min()


def extract_clv_per_transcript(utrs):
    # 5.4% of the annotated transcripts don't have UTR
    if utrs.shape[0] == 0:
        return np.nan

    strand = get_strand(utrs)
    if strand == '+':
        return utrs.end.max()
    elif strand == '-':
        return utrs.start.min()
    else:
        raise ValueError('unknown strand: {0}'.format(strand))


def calc_clv_per_transcript(grp):
    """infer 3'UTR end/CS when no 3'UTR is annotated"""
    tran = grp.query('feature == "transcript"')
    assert tran.shape[0] == 1
    strand = get_strand(grp)
    # since the last bp is supposed to represents a sc, just +/- 1
    if strand == '+':
        return tran.end.max() + 1
    elif strand == '-':
        return tran.start.min() - 1
    else:
        raise


# ref: https://github.com/bcgsc/utrtargets/blob/master/extract_targets.py
def get_cds_utr_iterator(group, strand):
    res = group[group.feature.isin(['CDS', 'UTR'])].sort_values('start')
    if strand == '+':
        return res.iterrows()
    elif strand == '-':
        return res.ix[::-1].iterrows()


def calc_sc_per_transcript(grp):
    """infer stop codon from 3'UTR when no sc is annotated"""
    tran = grp.query('feature == "transcript"')
    assert tran.shape[0] == 1

    state = '5UTR'              # starting from 5' UTR
    strand = get_strand(grp)
    for idx, row in get_cds_utr_iterator(grp, strand):
        # print(row.feature)
        if row.feature == 'CDS':
            state = 'CDS'
        else:
            if state == '5UTR':
                pass
            elif state == 'CDS' or '3UTR':
                state = '3UTR'
                return row.start - 1 if strand == '+' else row.end + 1
    # otherwise, return the tail base anyway as the stop codon tail
    return tran.end.max() if strand == '+' else tran.start.min()


def gen_sc_clv_per_transcript_lean(grp):
    has_sc, has_utr = False, False
    strand = get_strand(grp)

    sc, clv = np.nan, np.nan

    scs = grp.query('feature == "stop_codon"')
    if scs.shape[0] > 0:
        has_sc = True
        sc = extract_sc_per_transcript(scs)

    # this could be 5' UTR, needs a way to exclude it!!
    utrs = grp.query('feature == "UTR"')
    if utrs.shape[0] > 0:
        clv = extract_clv_per_transcript(utrs)
        # compare sc and clv to make sure it's indeed 3' UTR
        if has_sc:
            if strand == '+':
                if clv > sc:
                    has_utr = True
            elif strand == '-':
                if clv < sc:
                    has_utr = True

    if not has_utr:
        # clv can always be easily calculated
        clv = calc_clv_per_transcript(grp)
        has_utr = True

    if not has_sc:
        sc = calc_sc_per_transcript(grp)

    source = get_source(grp)
    strand = get_strand(grp)
    return source, sc, clv, strand


def gen_sc_clv_per_transcript(grp):
    """
    refactored out gen_sc_clv_per_transcript_lean mainly to keep api consistent
    while also test cases valid
    """
    res_lean = gen_sc_clv_per_transcript_lean(grp)
    return pd.Series(
        list(res_lean) + [
            get_property(grp, 'gene_source'),
            get_property(grp, 'transcript_source'),
            get_property(grp, 'is_cds_end_NF'),
            get_property(grp, 'is_cds_start_NF'),
        ], index=[
            'source', 'sc', 'clv', 'strand',
            'gene_source',
            'transcript_source',
            'is_cds_end_NF',
            'is_cds_start_NF',
        ])


class TestCases(unittest.TestCase):
    def test_chr9_minus_strand_CDKN2A_no_3UTR(self):
        tdf = pd.read_csv(os.path.join(os.path.dirname(__file__), './tests/CDKN2A-ENST00000446177.csv'))
        source, sc, clv, strand = gen_sc_clv_per_transcript_lean(tdf)
        self.assertEqual(clv, sc - 1)
        self.assertEqual(clv, 21968723)

    def test_chr13_plus_strand_BRCA2_no_3UTR(self):
        # A case where there is no 3' UTR or stop codon. It stops at AAA (chr9: 32907428)
        tdf = pd.read_csv(os.path.join(os.path.dirname(__file__), './tests/BRCA2-ENST00000530893.csv'))
        source, sc, clv, strand = gen_sc_clv_per_transcript_lean(tdf)
        self.assertEqual(clv, sc + 1)
        self.assertEqual(clv, 32907429)

    def test_chr13_plus_strand_BRCA2_with_3UTR(self):
        # A typical case where stop codon is there but not annotated.
        # The 3 bases right before 32970230 (27-29) are TAG, a classic stop codon.
        tdf = pd.read_csv(os.path.join(os.path.dirname(__file__), './tests/BRCA2-ENST00000470094.csv'))
        source, sc, clv, strand = gen_sc_clv_per_transcript_lean(tdf)
        self.assertEqual(sc, 32970229)
        self.assertEqual(clv, 32972410)

    def test_chr13_plus_strand_BRCA2_with_3UTR_case2(self):
        # no stop_codon, but with 3'UTR
        tdf = pd.read_csv(os.path.join(os.path.dirname(__file__), './tests/BRCA2-ENST00000528762.csv'))
        source, sc, clv, strand = gen_sc_clv_per_transcript_lean(tdf)
        self.assertEqual(sc, 32950807)
        self.assertEqual(clv, 32953633)

    def test_chr9_minus_strand_GNAQ_with_3UTR(self):
        tdf = pd.read_csv(os.path.join(os.path.dirname(__file__), './tests/GNAQ-ENST00000397476.csv'))
        source, sc, clv, strand = gen_sc_clv_per_transcript_lean(tdf)
        self.assertEqual(sc, 80336239)
        self.assertEqual(clv, 80335834)


if __name__ == "__main__":
    unittest.main()
