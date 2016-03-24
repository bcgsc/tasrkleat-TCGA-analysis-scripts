# Set up development environment

Install anacdona, then

    conda create -p venv --file requirements.conda
	source activate venv

or use virtualenv + pip

    virtualenv venv
	. venv/bin/activate
	pip install pandas jupyter matplotlib
	# install more packages as needed

# Run jupyter

    jupyter notebook 
	
Go to localhost:8888, and start playing with the notebooks.

After all, commit, push, and you may want to create a pull request.

