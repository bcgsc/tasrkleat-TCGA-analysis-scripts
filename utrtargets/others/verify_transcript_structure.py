"""
This is a piece of UNCLEANED code, which leads to discovery of a bug in gtf
file, see this post for details

https://www.biostars.org/p/206362/#208593


Usage:

python verify_transcript_structure.py | tee output.txt

Then inspect the lines starting with ERROR

"""


import sys
import re

import pandas as pd
pd.set_option('display.max_columns', 250)
# Don't cut off long string
# http://stackoverflow.com/questions/26277757/pandas-to-html-truncates-string-contents
pd.set_option('display.max_colwidth', -1)


# adf: annotation df
col_names = ['seqname', 'source', 'feature', 'start', 'end', 'score',
             'strand', 'frame', 'attribute']
adf = pd.read_csv(
    # e.g. '../reference/Homo_sapiens.GRCh37.75.gtf'
    sys.argv[1],
    header=None, sep='\t', comment='#', low_memory=False,
    names=col_names)

re_gene_id = re.compile(r'gene_id "(?P<gene_id>[-$\.\w]+)"')
re_transcript_id = re.compile('transcript_id "(?P<transcript_id>[\w]+)"')
re_gene_name = re.compile(r'gene_name "(?P<gene_name>[-$\.\w]+)"')


def extract_transcript_id(attribute):
    res = re_transcript_id.search(attribute)
    # if feature is gene, the transcript_id will be empty
    if res is None:
        return ''
    else:
        return res.group('transcript_id')


def extract_gene_id(attribute):
    res = re_gene_id.search(attribute)
    return res.group('gene_id')


def extract_gene_name(attribute):
    res = re_gene_name.search(attribute)
    return res.group('gene_name')

adf['len'] = adf.end - adf.start + 1


adf['transcript_id'] = adf.attribute.apply(extract_transcript_id)
adf['gene_id'] = adf.attribute.apply(extract_gene_id)
adf['gene_name'] = adf.attribute.apply(extract_gene_name)
adf.drop('attribute', axis=1, inplace=True)


def assertEqual(a, b):
    try:
        assert a == b
    except AssertionError as err:
        print(('{0} != {1}'.format(a, b)))


def verify_transcript_structure(transcript_df):
    sum_len = transcript_df.groupby('feature').sum().len
    a = 0
    for f in ['CDS', 'UTR', 'stop_codon']:
        if f in sum_len.index:
            a += sum_len.ix[f]
    return a == sum_len.ix['exon']


adf = adf[(adf.source == 'protein_coding') & (adf.feature != 'gene') & (adf.feature != 'Selenocysteine')]
groups = adf.groupby('transcript_id')
for tid, tdf in groups:
    if not verify_transcript_structure(tdf):
        print('ERROR: {0}'.format(tid))
    else:
        print('CHECK: {0}'.format(tid))
