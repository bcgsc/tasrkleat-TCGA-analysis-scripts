# UTRTargets

This repo includes the relevant scripts and files used for generating the
targets sequences for the TASRKleat pipeline. The basic process:

1. Collect a list of gene names
1. Download and preprocess the raw reference FASTA file and GTF files
1. Obtain genomic coordniates of all exons of all transcripts by filtering the
   annotation records in the GTF file against the list of genes. Relationships
   among exon, UTR, CDS, start codon & stop codon are shown below (https://www.biostars.org/p/206362/#208593):
   
   1. UTR is part of exon
   1. CDS is part of exon
   1. start codon is part of CDS, hence part of exon, too
   1. stop codon is neither part of UTR or part of CDS, but it's still part of exon.

1. Add buffer region to the first and last exons, respectively, i.e. minus the starting
   coordinate of the first exon and plus the ending coordinate by a specified
   buffer size (e.g. 300 bps). Pay attention to strand (+/-).
1. Extract the sequences from the FASTA file based on the coordinates collected
   for all transcripts of all genes. Note: the collected sequences can be highly
   repetitive since one can have multiple transcripts via alternative splicing,
   which is inconvenient for alignment.

You can then do whatever your want with the collected sequences, e.g. build a
bloomfilter based on them.

## To reproduce `targets-for-tasrkleat/targets.fa`

```
git clone https://github.com/bcgsc/utrtargets

cd utrtargets/reference && ./wget.sh

cd ..


# create an virtual environment with Anaconda
conda create -p venv -y --file conda-package-list.txt
source activate venv/

# pysam-0.9.1 is not available on the default conda chanel, and putting it in
# conda-package-list.txt doesn't work. conda would complain when creating a new
# environment
conda install -c bioconda pysam=0.9.1


# this creates the interested sequences in FASTA format
python extract_targets.py \
    reference/Homo_sapiens.GRCh37.75.gtf \
    reference/Homo_sapiens.GRCh37.75.dna_sm.primary_assembly.fa \
    targets-for-tasrkleat/target_genes.txt \
    output.fa
    
# make sure it's the same
diff output.fa targets-for-tasrkleat/targets.fa
```


<!-- ### Source fasta files -->

<!-- 1. `targetUTRcell2009.fa`: from Ewan 6 targets -->
<!-- 1. `oncogene_targets.160531.fa`: from Ewan, 99 targets -->
<!-- 1. `rchiu_targets.fa`: from Readman, 10 targets -->

<!-- ### Combine -->

<!-- ``` -->
<!-- cat src/*.fa > combined.fa -->
<!-- ``` -->

<!-- ### biobloommaker -->

<!-- version: 2.0.12 -->

<!-- ``` -->
<!-- biobloommaker -p bf/combined combined.fa -->
<!-- ``` -->
