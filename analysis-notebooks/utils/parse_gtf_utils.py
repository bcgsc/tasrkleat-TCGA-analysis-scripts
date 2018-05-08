import re


# tested on a few test cases, made sure the regex string passed all gtf entries
# gene_name "MIR1302-10"; gene_id "ENST00000473358";
# gene_name "DDX11L1"; gene_id "ENST00000456328";
# gene_name "RP11-34P13.7"; gene_id "ENST00000466430";
# gene_name "HLA-DRB3$0301"; gene_id "ENST00000426847";

RE_GENE_ID = re.compile(r'gene_id "(?P<gene_id>[-$\.\w]+)"')
RE_TRANSCRIPT_ID = re.compile('transcript_id "(?P<transcript_id>[\w]+)"')
RE_GENE_NAME = re.compile(r'gene_name "(?P<gene_name>[-$\.\w]+)"')
RE_GENE_SOURCE = re.compile(r'gene_source "(?P<gene_source>[-$\.\w]+)"')
RE_TRANSCRIPT_SOURCE = re.compile(r'transcript_source "(?P<transcript_source>[-$\.\w]+)"')
# the coding region end could not be confirmed.
RE_IS_CDS_END_NF = re.compile(r'tag "(?P<is_cds_end_nf>cds_end_NF)";')
RE_IS_CDS_START_NF = re.compile(r'tag "(?P<is_cds_end_nf>cds_start_NF)";')


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


def extract_gene_source(attribute):
    res = RE_GENE_SOURCE.search(attribute)
    return res.group('gene_source')


def extract_transcript_source(attribute):
    res = RE_TRANSCRIPT_SOURCE.search(attribute)
    if res is None:
        return ''
    else:
        return res.group('transcript_source')


def extract_cds_end_NF(attribute):
    res = RE_IS_CDS_END_NF.search(attribute)
    return True if res else False


def extract_cds_start_NF(attribute):
    res = RE_IS_CDS_START_NF.search(attribute)
    return True if res else False


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

    print('extracting gene source...')
    adf['gene_source'] = adf.attribute.apply(extract_gene_source)

    print('extracting transcript source...')
    adf['transcript_source'] = adf.attribute.apply(extract_transcript_source)

    print('extracting tag cds_end_NF...')
    adf['is_cds_end_NF'] = adf.attribute.apply(extract_cds_end_NF)

    print('extracting tag cds_start_NF...')
    adf['is_cds_start_NF'] = adf.attribute.apply(extract_cds_start_NF)

    adf.drop('attribute', axis=1, inplace=True)


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
