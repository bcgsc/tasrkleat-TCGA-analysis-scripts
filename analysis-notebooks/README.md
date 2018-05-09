To reproduce the results from the paper, please follow the instructions in the
notebooks:

Please download needed data from
http://bcgsc.ca/downloads/tasrkleat-static/off-cloud/

**Note**: the location of output directories may be outdated in the notebooks.
Please add symblinks accordingly. In general, it looks for the three
directories, `metadata`, `reference_data`, `results_data` for data
[already previously prepared](http://bcgsc.ca/downloads/tasrkleat-static/off-cloud/),
and outputs analysis results in a separate directory `results`.

# Key notebooks and their functions:

## Core

Key notebooks for generating the results reported in the manuscript, including
supplemental materials

1. `diff-APA-analysis`: analysis of tumor-specific APA cases, their
  recurrence, and trends (multi or not). Also, other stats include:
  * 40% (31/77) of the cases report here have more than 2 CSs that undergo significant
    change cases. Write out `aca_pval_with_sc_info.csv`.
  * 11% (33/297) of the predicted CSs undergoing siginificant frequency shift
    from normal to tumor in 33 genes are mapped to more than one stop codons
  * ~48% (16/33) of the reported APA genes have stop codons involved in NMD
1. `resolve-APA-trend.ipynb`: trying to resolve the trend in terms of
   shortening/lengthening for the identified tumor-specific APA cases. Write out
   `aca_trends.csv`.
1 `map-clv2sc-based-on-annotation-for-114-genes.ipynb`: extract the mapping
  between cleavage sites and stop codons based on GTF annotation for 114 genes.

In general, the notebook names are indicative of their functions.

## Misc

Miscellaneous results that may not be reported in the manuscript, but still
may be useful for reproduction or insights.


## CS-postprocessing

1. `concate-all-cba.KLEAT.ipynb`: Simply concatenate all cba.KLEAT
  files from the online, id each file with kleat_fname column, which has
  one-to-one correspondence to TCGA barcode
1. `cleanup-all-cba.KLEAT.ipynb`: filter CSs by target genes, map kleat_fname to
  analysis_id, remove unused columns, and rename remaining columns accordingly
1. `hexamer-research.ipynb`: redo hexamer search as that in the pipeline doesn't
  consider lower case sequences.
1. `filter-based-on-APA-evidence-and-annotation.ipynb`: do as the file name
  indicated
1. `cluster-filtered-cleavage-sites.ipynb`: cluster cleavage sites with single-linkage


# CS postprocess order:

The result of postprocess is prepared for you here
http://bcgsc.ca/downloads/tasrkleat-static/off-cloud/results_data/all_cba.KLEAT.on-target-cleaned.filtered.clustered.csv.gz,

so you don't need to reproduce it yourself.

concat => clenaup => filter => cluster
