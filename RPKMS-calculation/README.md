# Steps

### Run bedtools `genome_cov` to calculate coverage per base

1. use `calc_genome_cov.sh` to generate commands for `genome_cov` calculation in
   parallel.
1. then execute the generated commands to run `genome_cov`.

### Use custom script to calculate coverage per gene.

1. run `gen_df_genes_coords.py` to obtain genomic coordinates of candidate
   genes, and save them into `df_gene_coords.csv`, which will be used by
   `calc_expression_level.py`.
1. run `calc_cov_per_gene.py`, and calculate coverage per gene.

### Run `concat_all_genomecov_by_gene.ipynb` to concatenate all `genome-by-gene.csv`
files, and save to `all-genome-cov-by-genes.csv`.


# About kmers

Just put here for future reference. kmer analysis isn't part of the manuscript.

kmers files in tcga

e.g.

```
hpce705[0]:l ${PWD}/tcga/BRCA/RNA/RNA-Seq/UNC-LCCC/ILLUMINA/UNCID_2170988.76f171dd-1980-4504-8dc2-d8c8264458c9.130513_UNC14-SN744_0340_BD250CACXX_7_TAGCTT.tar.gz/tasrkleat-results/align_reads2genome/
total 99M
 31M -rw-r--r-- 1 zxue users  28M Oct 27  2016 genomecov.bed.csv
 26K -rw-r--r-- 1 zxue users 2.7K Oct 28  2016 genomecov-by-gene.csv
 30M -rw-r--r-- 1 zxue users  26M Nov 10  2016 mer_counts.jf
 25M -rw-r--r-- 1 zxue users  22M Nov 10  2016 kmer_count.csv.gz
 16M -rw-r--r-- 1 zxue users  14M Nov 17  2016 kmer_count.sorted.csv.gz
 2.0K drwxr-sr-x 2 zxue users  182 Nov 18  2016 ./
 2.0K drwxr-sr-x 4 zxue users   74 Nov 19  2016 ../
```
