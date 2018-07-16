This repo records the notebooks that are used in the off-cloud analysis of the
targeted APA study of TCGA samples.

Notebooks that related to result reproduction are mostly in
`analysis-notebooks`. For details, please see the `README.md` inside each
folder.


### Citation

<a href=https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-018-4903-7
, target=_blank> Xue Z, Warren RL, Gibb EA, MacMillan D, Wong J, Chiu R, et al.
<cite>Recurrent tumor-specific regulation of alternative polyadenylation of
cancer-related genes</cite>. BMC Genomics. 2018;19:536</a>


### Run codes

Most of the codes are run in jupyter notebooks format. Most of the packages used
are listed in `eda_import.py`. Some extra include pysam, etc.

To get started, you may not need all the packages. I used Miniconda to manage
python environments. Try getting start with

```
conda create -p venv -y python=3 pandas jupyter matplotlib scikit-learn
source activate venv
```

To run jupyter notebook server

```
jupyter notebook --no-browser --ip 0.0.0.0
```

### Data

The data used as input to the off-cloud analysis are produced by on-cloud run
with [tasrkleat](https://github.com/bcgsc/tasrkleat). All (i.e. 10K-TCGA-sample
run) are generated with tasrkleat-v0.1 except the benchmark part which used
tasrkleat-0.1.1 with adjustment in the reads2genome alignment step. The
reads2genome alignment step doesn't affect the results of KLEAT, but may have an
effect for the other tools.


### Q&A

If you have any question, please ask by [opening a new
issue](https://github.com/bcgsc/tasrkleat-TCGA-analysis-scripts/issues/new?template=-your-question-.md).
