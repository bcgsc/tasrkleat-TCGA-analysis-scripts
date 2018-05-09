import os
import sys
import glob
import gzip
import re
import datetime
import tempfile
import itertools
import functools
import pickle
import multiprocessing
from collections import Counter

import numpy as np
import pandas as pd
import scipy.stats as stats
from scipy import interp
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster, fclusterdata

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.patheffects as PathEffects
from matplotlib import gridspec
from mpl_toolkits.mplot3d import Axes3D
from IPython.display import display, Image, HTML
import seaborn as sns

# related to preprocessing 
from sklearn.preprocessing import LabelEncoder, OneHotEncoder, StandardScaler, label_binarize
from sklearn.feature_selection import SelectKBest, chi2
from sklearn.model_selection import GridSearchCV

# related to unsupervised learning 
from sklearn.decomposition import TruncatedSVD, PCA
from sklearn.manifold import TSNE
from sklearn.neighbors import KNeighborsClassifier
from sklearn.cluster import KMeans

# related to learning
from sklearn.linear_model import LogisticRegression, LogisticRegressionCV
from sklearn import svm
from sklearn.tree import DecisionTreeClassifier, export_graphviz
from sklearn.ensemble import  GradientBoostingClassifier, AdaBoostClassifier, RandomForestClassifier
from sklearn.multiclass import OneVsRestClassifier, OneVsOneClassifier

# related to evaluation
from sklearn.model_selection import StratifiedKFold, train_test_split, cross_val_score
from sklearn import metrics
kappa_scorer = metrics.make_scorer(metrics.cohen_kappa_score)

pd.set_option('display.max_columns', 250)
# Don't cut off long string
# http://stackoverflow.com/questions/26277757/pandas-to-html-truncates-string-contents
pd.set_option('display.max_colwidth', -1)

mpl.style.use('classic')
mpl.rcParams['figure.figsize'] = (8, 4.5)
get_ipython().magic('matplotlib inline')

# for auto-reloading external modules
# see http://stackoverflow.com/questions/1907993/autoreload-of-modules-in-ipython
get_ipython().magic('load_ext autoreload')
get_ipython().magic('autoreload 2')
