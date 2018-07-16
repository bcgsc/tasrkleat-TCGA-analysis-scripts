### Citation

<a href=https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-018-4903-7
, target=_blank> Xue Z, Warren RL, Gibb EA, MacMillan D, Wong J, Chiu R, et al.
<cite>Recurrent tumor-specific regulation of alternative polyadenylation of
cancer-related genes</cite>. BMC Genomics. 2018;19:536</a>

Notebooks related to result reproduction are mostly in `analysis-notebooks`. For
details, please see the `README.md` inside each folder.


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

### Q&A

If you have any question, please ask by
[opening a new issue](https://github.com/bcgsc/tasrkleat-TCGA-analysis-scripts/issues/new?template=-your-question-.md).
