* `concate-all-cba.KLEAT.ipynb`: Simply concatenate all cba.KLEAT
  files from the online, id each file with kleat_fname column, which has
  one-to-one correspondence to TCGA barcode
* `cleanup-all-cba.KLEAT.ipynb`: filter CSs by target genes, map kleat_fname to
  analysis_id, remove unused columns, and rename remaining columns accordingly
* `hexamer-research.ipynb`: redo hexamer search as that in the pipeline doesn't
  consider lower case sequences.
* `filter-based-on-APA-evidence-and-annotation.ipynb`: do as the file name
  indicated
* `cluster-filtered-cleavage-sites.ipynb`: cluster cleavage sites with single-linkage
