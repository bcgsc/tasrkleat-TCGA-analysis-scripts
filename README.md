## Raw data inside `gtfs`

* `ucsc.gtf.gz` was downloaded from http://genome.ucsc.edu/cgi-bin/hgTables, and
then sorted and indexed. See `how-to-download-ucsc-gtf-gz.png` for the
parameters selected when it's downloaded.

* `Homo_sapiens.GRCh37.75.gtf.gz` was downloaded from
http://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz,
and then sorted and indexed.

* `ensembl.fixed.sorted.gz` and `ensembl.fixed.sorted.gz.tbi` were downloaded from
https://github.com/bcgsc/KLEAT/. The former is basically the same as
`ucsc.gtf.gz` but with its `gene_id` attribute modified.

Commands used for sorting and indexing:
```bash
# ref: http://www.htslib.org/doc/tabix.html
for GTF_PREFIX in ucsc.gtf.gz Homo_sapiens.GRCh37.75.gtf; do 
    unpigz -c ${GTF_PREFIX}.gz | \grep -v ^"#" | sort -k1,1 -k4,4n | bgzip > ${GTF_PREFIX}.sorted.gz
    tabix Homo_sapiens.GRCh37.75.gtf.sorted.gz
done
```

## To reproduce the results as in `diff`, run 

The comparison is limited to Chromosome 1.
```
bash compare.sh
```

Then, verify they are indeed the same. This should output all 0 return codes.

```bash
for i in  \
    klt_CDS.tsv \
    klt_exon.tsv \
    klt_start_codon.tsv \
    klt_stop_codon.tsv \
    ucsc_CDS.tsv \
    ucsc_exon.tsv \
    ucsc_start_codon.tsv \
    ucsc_stop_codon.tsv \
    ens_CDS.tsv \
    ens_exon.tsv \
    ens_start_codon.tsv \
    ens_stop_codon.tsv \
    sdiff_CDS.log \
    sdiff_exon.log \
    sdiff_start_codon.log \
    sdiff_stop_codon.log \
    ; do
    diff __diff/$i diff/$i;
    echo $?
done
```

### Findings after comparison:

1. All exon coordinates in `ens_exon.tsv` and `klt_exon.tsv` are the same.

2. Not all CDS, start\_codon, stop\_codon coordinates in `ens.gtf.gz` and
`klt_gtf.gz` are the same. In general, some of the CDS coordinates in
`ens.gtf.gz` have 3 more bases than corresponding ones in `klt.gtf.gz`.

3. Problem: `klt_stop_codons.tsv` is a *superset* of `ens_stop_codons.tsv`.
However, for some of the annotations in `klt_stop_codons.tsv`, it doesn't really
look like a stop codon when inspected closely. e.g. chr1:738532-738534, minus
strand:

```
$:zcat klt.gtf.gz | grep 'stop_codon.*738532.*ENST00000599533'
chr1    hg19_ensGene    stop_codon      738532  738534  0.000000        -       .       gene_id "AL669831.1"; transcript_id "ENST00000599533";
$:zcat ens.gtf.gz | grep 'stop_codon.*738532.*ENST00000599533'
# not output
```

The corresponding bases are CTG with reverse complementary as CAG, which is a
codon for Glutamine (Q), and does not correspond to any of the
[stop codons](https://en.wikipedia.org/wiki/Stop_codon): **TAG**, **TAA**, and
**TGA** in DNA, or **UAG**, **UAA**, and **UGA** in RNA.

<!-- This is very likely due to frame/phase -->
<!-- 4. Problem: this is common to both `klt_stop_codons.tsv` and -->
<!-- `ens_stop_codons.tsv`: some of the stop codons are only two bp (e.g. -->
<!-- chr1:1203242-1203243). Not sure why. -->


### Conclusion

The two, `ensembl.fixed.sorted.gz` and `Homo_sapiens.GRCh37.75.gtf.gz` are not
the same, and shouldn't be used interchangably. The annotations for MT and
unassembled contigs could be more different based on observation.
