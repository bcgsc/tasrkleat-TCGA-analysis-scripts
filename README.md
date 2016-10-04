`ucsc.gtf.gz` is downloaded from http://genome.ucsc.edu/cgi-bin/hgTables. See
`how-to-download-ucsc-gtf-gz.png` for the parameters selected when it's
downloaded.

`ensembl.fixed.sorted.gz` is downloaded from
https://github.com/bcgsc/KLEAT/raw/master/ensembl.fixed.sorted.gz. It's
basically the same as `ucsc.gtf.gz` but with gene_id attribute modified.


`Homo_sapiens.GRCh37.75.gtf.gz` is downloaded from
http://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz.


# To reproduce the results as in `diff`

## `ensembl.fixed.sorted.gz` vs `ucsc.gtf.gz`

To verify `ensembl.fixed.sorted.gz` and `ucsc.gtf.gz` are largely the same.

```
ln -sf ensembl.fixed.sorted.gz dan.gtf.gz

for c in exon CDS start_codon stop_codon; do
    zcat dan.gtf.gz   | \grep -P "^chr1\t.*\t${c}\t" | sed 's/^chr//g' | awk '{print $1,$3,$4,$5,$7}' | tr ' ' '\t' | sort -k 1n -k 3n -k 4n  > dan_${c}.tsv &
done
wait

for c in exon CDS start_codon stop_codon; do
    zcat ucsc.gtf.gz  | \grep -P "^chr1\t.*\t${c}\t" | sed 's/^chr//g' | awk '{print $1,$3,$4,$5,$7}' | tr ' ' '\t' | sort -k 1n -k 3n -k 4n  > ucsc_${c}.tsv
done
wait

# Should print all 0
for c in exon CDS start_codon stop_codon; do
    diff ucsc_${c}.tsv dan_${c}.tsv
    echo $?
done
```

## `ensembl.fixed.sorted.gz` vs `Homo_sapiens.GRCh37.75.gtf.gz`

To show `ensembl.fixed.sorted.gz` and `Homo_sapiens.GRCh37.75.gtf.gz` downloaded
from Ensembl directly have different coordinates.

The following analysis are limited to Chromosome 1.

```
ln -sf Homo_sapiens.GRCh37.75.gtf.gz ens.gtf.gz

for c in exon CDS start_codon stop_codon; do
    zcat ens.gtf.gz  | \grep -P "^1\t.*\t${c}\t" | awk '{print $1,$3,$4,$5,$7}' | tr ' ' '\t' | sort -k 1n -k 3n -k 4n > ens_${c}.tsv &
done
wait

for c in exon CDS start_codon stop_codon; do
    sdiff ens_${c}.tsv dan_${c}.tsv > sdiff_${c}.txt
    echo $?
done
```

### Finding:

1. All exon coordinates in them are the same.

2. Not all CDS, start\_codon, stop\_codon coordinates in them are the same. In
general, some of the CDS coordinates in ucsc.gtf.gz have 3 more bases than
corresponding ones in ensembl.fixed.sorted.gz. 

3. `dan_stop_codons.tsv` is a *subset* of `ens_stop_codons.tsv`. For some of the
annotations in `dan_stop_codons.tsv`, but they don't really look like so when
inspected closely. e.g. chr1:738532-738534, minus strand, the bases are GTC,
which doesn't correspond to any of the stop codons (UAA, UAG, UGA)

4. This is common to both `dan_stop_codons.tsv` and `ens_stop_codons.tsv`: some
of the stop codons are only two bp. Not sure why. e.g. chr1:1203242-1203243

### Conclusion

The two, `ensembl.fixed.sorted.gz` and `Homo_sapiens.GRCh37.75.gtf.gz` are not
the same, and shouldn't be used interchangably. The annotations for MT and
unassembled contigs are expected to be more different based on observation
