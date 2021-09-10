# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.11.3
#   kernelspec:
#     display_name: conda-ml
#     language: python
#     name: conda-ml
# ---

# %% [markdown]
# # Clinical data and EDA for GEMDeCan

# %% [markdown]
# ## Imports

# %%
# Convenient IPython console
# comment it if you execute Jupyterlab on a remote server
# %qtconsole --style monokai
# and use instead
# # %connect_info
# then jupyter qtconsole --style monokai --existing ./remote_kernel-1 --ssh alexis@CRCT2112

# %% imports
import numpy as np
import pandas as pd
import random
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import os
from pathlib import Path
from itertools import combinations

from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import PowerTransformer
from sklearn.decomposition import PCA
from sklearn.decomposition import FastICA
from sklearn.decomposition import KernelPCA
import umap

from tableone import TableOne
from scipy import stats

pd.options.display.max_rows = 999
pd.options.display.max_columns = 999
# %% Load deconvolution and response data
path_interm = "../../data/intermediate/"

list_deconv = [
    path_interm + "/Riaz_Cell_2017/deconv-Riaz_Cell_2017.csv",
    path_interm + "/Hugo_Cell_2016/deconv-Hugo_Cell_2016.csv",
    path_interm + "/Gide_Cancer-Cell_2019/deconv-Gide_Cancer-Cell_2019.csv",
    path_interm + "/Maha/deconv-Maha.csv",
]
list_response = [
    path_interm + "/Riaz_Cell_2017/response-Riaz_Cell_2017.csv",
    path_interm + "/Hugo_Cell_2016/response-Hugo_Cell_2016.csv",
    path_interm + "/Gide_Cancer-Cell_2019/response-Gide_Cancer-Cell_2019.csv",
    path_interm + "/Maha/response-Maha.csv",
]
list_clinic = [
    path_interm + "/Riaz_Cell_2017/clinic-Riaz_Cell_2017.csv",
    path_interm + "/Hugo_Cell_2016/clinic-Hugo_Cell_2016.csv",
    path_interm + "/Gide_Cancer-Cell_2019/clinic-Gide_Cancer-Cell_2019.csv",
#     path_interm + "/Maha/.csv",
]
list_dataset = [
    'Riaz',
    'Hugo',
    'Gide',
    'Maha',
]
list_cancer = [
    'melanoma',
    'melanoma',
    'melanoma',
    'lung',
]

nb_samples = len(list_deconv)
# %% [markdown]
# ### third level

# %%
# comment code
1+1

# %% [markdown]
# #### fourth level

# %%
