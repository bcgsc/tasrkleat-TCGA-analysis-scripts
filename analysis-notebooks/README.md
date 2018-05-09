To reproduce the results from the paper, please follow the instructions in the
notebooks:

Please download needed data from
http://bcgsc.ca/downloads/tasrkleat-static/off-cloud/

**Note**: *the location of output directories may be outdated in the notebooks.
Please modify the path to those directories or add symblinks accordingly. In
general, it looks for the three directories, `metadata`, `reference_data`,
`results_data` for data
[already previously prepared](http://bcgsc.ca/downloads/tasrkleat-static/off-cloud/),
and outputs analysis results in a separate directory `results`*.

# Key notebooks and their functions:

## Core

Key notebooks for generating the results reported in the manuscript, including
supplemental materials

##### `diff-APA-analysis.ipynb`

Analysis of tumor-specific APA cases, their recurrence, and trends (multi or
not). Also, other stats include:
  * 40% (31/77) of the cases report here have more than 2 CSs that undergo significant
    change cases. Write out `aca_pval_with_sc_info.csv`.
  * 11% (33/297) of the predicted CSs undergoing siginificant frequency shift
    from normal to tumor in 33 genes are mapped to more than one stop codons
  * ~48% (16/33) of the reported APA genes have stop codons involved in NMD

##### `resolve-APA-trend.ipynb`

To resolve the trend in terms of shortening/lengthening for the identified
tumor-specific APA cases. Write out `aca_trends.csv`.

##### `map-clv2sc-based-on-annotation-for-114-genes.ipynb`

Extract the mapping between cleavage sites and stop codons based on GTF
annotation for 114 genes.

##### `plot-arc-plots-for-supp-mat-for-114-genes-and-14-diseases-without-percent-cutoff-with-expr.ipynb`

In contrast to the figures in main, there is cutoff imposed on the percent
change, so the plots might be noiser, but more comprehensive in terms of
information. Also included is expression distribution.


In general, the notebook names are indicative of their functions.

## CS-postprocessing

##### `concate-all-cba.KLEAT.ipynb`

Simply concatenate all cba.KLEAT files from the online, id each file with
kleat_fname column, which has one-to-one correspondence to TCGA barcode.

##### `cleanup-all-cba.KLEAT.ipynb`

Filter CSs by target genes, map kleat_fname to analysis_id, remove unused
columns, and rename remaining columns accordingly.

##### `hexamer-research.ipynb`

Redo hexamer search as that in the pipeline doesn't consider lower case
sequences.

##### `filter-based-on-APA-evidence-and-annotation.ipynb` 

Filter CSs based on APA evidence and annotation.

##### `cluster-filtered-cleavage-sites.ipynb`

Cluster cleavage sites using single-linkage.


## Misc

Miscellaneous results that may not be reported in the manuscript, but still
may be useful for reproduction or insights.

##### `GTF-analysis-general-properties.ipynb`

* Count stop codons per gene: about 65% of genes have multiple stop codons
* % of protein coding genes over whole genome: ~2%, not including introns and
  union overlapping exons from multiple transcripts of the same gene.
* Even chr1 has just a bit over 2k (2076) genes.
* Over 73% of protein coding genes have more than one transcripts

##### `GTF-analysis-general-properties-for-target-genes.ipynb`

* No overlap among target genes
* Total length of target genes is ~93.5M bp, including introns

##### `GTF-analysis-of-overlapped-genes.ipynb`

Overall, 18.6% of neighbouring genes overlap across all chromosomes.
