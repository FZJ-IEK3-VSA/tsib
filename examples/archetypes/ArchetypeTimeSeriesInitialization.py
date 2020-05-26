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
# # Validation of the tsib.buildingconfig
# Author: Leander Kotzur
#
# Date: 21.09.2019

# %%
# %matplotlib inline
# %load_ext autoreload
# %autoreload 2

# %%
import tsib
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import copy
import os

# %% [markdown]
# ### Read building configurations data

# %% [markdown]
# Building configurations for 5 archetypes

# %%
bdgs_raw = pd.read_csv(os.path.join('5_archetypes_config.csv'), index_col = 0)

# %%
bdgs_raw

# %%
bdgs_dict = bdgs_raw.T.to_dict()

# %% [markdown]
# ### Loop over the buildings, define them and get their profile

# %%
bdg_cfgs = {}
bdg_cfgs_print = {}

# %%
for bdg_ix in bdgs_dict:
    bdg_cfgs[bdg_ix] = tsib.BuildingConfiguration(bdgs_dict[bdg_ix])
    bdg_cfgs_print[bdg_ix] = bdg_cfgs[bdg_ix].getBdgCfg()

# %%
bdg_cfgs_print[0]

# %%
import tsib

# %%
bdgs = {}

# %%
for bdg_ix in bdgs_dict:
    bdgs[bdg_ix] = tsib.Building(configurator = bdg_cfgs[bdg_ix])
    bdgs[bdg_ix].getLoad()

# %% [markdown]
# ### Plot example results and export to CSV

# %%
bdgs[0].timeseries['Heating Load'].plot()

# %%
bdgs[0].timeseries['Electricity Load'].plot()

# %%
bdgs[0].toCSV(filename="example_result")
