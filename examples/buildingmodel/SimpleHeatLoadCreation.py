# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.4.2
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %% [markdown]
# ### Validation of the heat model with the Tabula Database

# %% [markdown]
# Import all relevant modules

# %%
import os
import pandas as pd
# %matplotlib inline
# %load_ext autoreload
# %autoreload 2
import tsib
import tsib
import tsib.buildingmanager as bm
import tsib
import tsib.thermal.utils as utils

# %%
import matplotlib.pyplot as plt

# %%
# %matplotlib inline
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})

# %%
import numpy as np

# %%
EXPORT_PATH = os.path.join('plots')

# %% [markdown]
# ### Define the building

# %% [markdown]
# Init building 

# %%
bdgcfg = tsib.BuildingConfiguration({ 'refurbishment': False, 'nightReduction': True, 'buildingYear': 2005,
                                       'occControl': False, 'capControl': True,'n_persons': 2,
                                       'comfortT_lb': 20.0, 'comfortT_ub': 26.0, 'roofOrientation': 0.0,
                                       'longitude': 7., 'latitude': 54. })

# %%
bdgObj = tsib.Building(configurator = bdgcfg)

# %% [markdown]
# ## Get the load

# %%
bdgObj.get_load()

# %% [markdown]
# Plot

# %%
bdgObj.detailed_results['Heating Load'].plot()

# %%
bdgObj.detailed_results['Heating Load']['20100101':'20100107'].plot()

# %%
