# CGHub MANIFEST Summary

The notebooks in this project summarize interested statistics from from
[Cancer Genomics Hub](https://cghub.ucsc.edu/) (CGHub) based on the
[LATEST_MANIFEST.tsv](https://cghub.ucsc.edu/reports/SUMMARY_STATS/LATEST_MANIFEST.tsv)
file using [pandas](http://pandas.pydata.org/).

### Development

Install anacdona, then

    conda create -p venv -y pandas jupyter matplotlib
	source activate venv

or use virtualenv + pip

    virtualenv venv
	. venv/bin/activate
	pip install pandas jupyter matplotlib
	# install more packages as needed

### Run jupyter

    jupyter notebook 
	
Go to localhost:8888, and start playing with the notebooks.

After all, commit, push, and you may want to create a pull request.

