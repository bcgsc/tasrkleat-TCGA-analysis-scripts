### Citation

Xue, Z. et al. 2018. *Recurrent tumor-specific regulation of alternative
polyadenylation of cancer-related genes*. Under review. Preprint is available at
https://doi.org/10.1101/160960.

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
