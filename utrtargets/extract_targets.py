#!/usr/bin/env python


import sys
import re

import pandas as pd

import pysam


# tested on a few test cases, made sure the regex string passed all gtf entries
# gene_name "MIR1302-10"; gene_id "ENST00000473358";
# gene_name "DDX11L1"; gene_id "ENST00000456328";
# gene_name "RP11-34P13.7"; gene_id "ENST00000466430";
# gene_name "HLA-DRB3$0301"; gene_id "ENST00000426847";

RE_GENE_ID = re.compile(r'gene_id "(?P<gene_id>[-$\.\w]+)"')
RE_TRANSCRIPT_ID = re.compile('transcript_id "(?P<transcript_id>[\w]+)"')
RE_GENE_NAME = re.compile(r'gene_name "(?P<gene_name>[-$\.\w]+)"')


# def get_coords_selective(group, strand):
#     """get the coordinates of exons that contain 3' UTRs and one extra neighboring
#     exon, much more selective than get_coords. The later is preferred, but this
#     function is kept for later reference
#     """
#     state = '5UTR'              # starting from 5' UTR
#     three_UTRs = []
#     for idx, row in get_cds_utr_iterator(group, strand):
#         if row.feature == 'CDS':
#             state = 'CDS'
#         else:
#             if state == '5UTR':
#                 pass
#             elif state == 'CDS' or '3UTR':
#                 state == '3UTR'
#                 three_UTRs.append(row)

#     if strand == '+':
#         tip = 'end'
#     else:
#         tip = 'start'

#     # found tips of all 3' UTRs, then based on the tips look for the
#     # corresponding exons
#     UTR_tips = [_[tip] for _ in three_UTRs]

#     coords, remains = [], []
#     if not UTR_tips:
#         return coords

#     for idx, row in get_exon_iterator(group, strand):
#         if row[tip] in UTR_tips:
#             coords.append((row.start, row.end))
#         else:
#             remains.append((row.start, row.end))
#     # append the last exon preceeding the one that contains UTR
#     if remains:
#         coords.append(remains[-1])

#     # sort to be left to right so that the corresponding sequences can be
#     # extracted correctly from fasta
#     coords.sort(key=lambda x: x[0])
#     return coords


def get_coords(group, buffer_size=300):
    """get the coordinates of all exons plus buffer_size to ends of 5' UTR and 3'
    UTRs

    group is a dataframe of all subgene elements for a particular transcript_id
    """
    coords = group[group.feature == 'exon'][['start', 'end']].sort_values('start').values.tolist()
    coords[0][0] -= buffer_size
    coords[-1][1] += buffer_size
    coords = [tuple(_) for _ in coords]
    return coords


def get_cds_utr_iterator(group, strand):
    res = group[group.feature.isin(['CDS', 'UTR'])].sort_values('start')
    if strand == '+':
        return res.iterrows()
    elif strand == '-':
        return res.ix[::-1].iterrows()


def get_exon_iterator(group, strand):
    res = group[group.feature.isin(['exon'])].sort_values('start')
    if strand == '+':
        return res.iterrows()
    elif strand == '-':
        return res.ix[::-1].iterrows()


def extract_transcript_id(attribute):
    res = RE_TRANSCRIPT_ID.search(attribute)
    # if feature is gene, the transcript_id will be empty
    if res is None:
        return ''
    else:
        return res.group('transcript_id')


def extract_gene_id(attribute):
    res = RE_GENE_ID.search(attribute)
    return res.group('gene_id')


def extract_gene_name(attribute):
    res = RE_GENE_NAME.search(attribute)
    return res.group('gene_name')


def extract_info(adf):
    """extract interested information"""
    print('extracting length...')
    adf['len'] = adf.end - adf.start + 1

    print('extracting transcript id...')
    adf['transcript_id'] = adf.attribute.apply(extract_transcript_id)

    print('extracting gene id...')
    adf['gene_id'] = adf.attribute.apply(extract_gene_id)

    print('extracting gene name...')
    adf['gene_name'] = adf.attribute.apply(extract_gene_name)

    adf.drop('attribute', axis=1, inplace=True)


def load_gtf(gtf, target_genes=None):
    """
    target_genes: a list of interested gene names
    """
    # http://uswest.ensembl.org/info/website/upload/gff.html
    names = ['seqname', 'source', 'feature', 'start', 'end',
             'score', 'strand', 'frame', 'attribute']
    # adf: annotation df
    print('reading {0}...'.format(gtf))
    adf = pd.read_csv(gtf, header=None, sep='\t', comment='#',
                      low_memory=False, names=names)

    # only interested in a subset of the gtf
    adf = adf[adf.source == 'protein_coding']
    adf = adf[-adf.feature.isin(['transcript', 'gene'])]
    extract_info(adf)
    if target_genes is not None:
        adf = adf[adf.gene_name.isin(target_genes)]
    return adf


def sanity_check(group):
    # assert all seq sections have the same strand
    assert group.strand.unique().shape[0] == 1
    assert group.strand.unique()[0] in ['+', '-']
    # assert all seq sections are in the same chromosome
    assert group.seqname.unique().shape[0] == 1
    assert group.gene_name.unique().shape[0] == 1
    assert group.gene_id.unique().shape[0] == 1


def fetch_seq(fasta, chrom, coords):
    chunks, header_parts = [], []
    for k, (start, end) in enumerate(coords):
        # start & end from GTF are inclusive and 1-based (ensembl)
        # https://www.biostars.org/p/84686/ while pysam.FastaFile.fetch is
        # -0-based and right-exclusive
        try:
            chunks.append(fasta.fetch(chrom, start - 1, end))
        except KeyError as err:
            print(err)
        header_parts.append('chunk{0}:{1}-{2}'.format(k + 1, start, end))
    seq = ''.join(chunks)
    return seq, header_parts


def main(gtf, fasta, target_genes, output_fa='targets.fa'):
    adf = load_gtf(gtf, target_genes)
    fasta = pysam.FastaFile(fasta)

    gs = adf.groupby('transcript_id')

    with open(output_fa, 'wt') as opf:
        for transcript_id, group in gs:
            chrom = group.seqname.unique()[0]
            gene_name = group.gene_name.unique()[0]
            gene_id = group.gene_id.unique()[0]

            # coords = get_coords_selective(group, strand)
            coords = get_coords(group)
            if not coords:      # meaning no 3' UTR chunks found
                continue

            header = [gene_name, gene_id, transcript_id, chrom]
            seq, header_parts = fetch_seq(fasta, chrom, coords)
            header.extend(header_parts)
            header.append(len(seq))
            header = '|'.join(map(str, header))

            if not seq:
                # sometimes seq can be '' because it's in a patch chromosome.
                # e.g. HG1257_PATCH, HG122_PATCH
                continue

            opf.write('>{0}\n{1}\n'.format(header, seq))


if __name__ == "__main__":
    if len(sys.argv[1:]) != 4:
        print('''
Usage: python extract_targets.py <gtf-file> <fasta-file> <target-genes-list> <output-fa>

<target-genes-file> is just a list of target genes, see share/target_genes.template.txt for format
''')
        sys.exit(1)

    gtf, fa, target_genes_file, output_fa = sys.argv[1:]

    with open(target_genes_file) as inf:
        target_genes = [_.strip() for _ in inf.readlines()]
        print('{0} genes read from {1}'.format(
            len(target_genes), target_genes_file))

    main(gtf, fa, target_genes, output_fa)
