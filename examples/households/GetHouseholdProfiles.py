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
# ## Create household profiles in parallel
#
# Author: Leander Kotzur

# %%
# %matplotlib inline
# %load_ext autoreload
# %autoreload 2

# %% [markdown]
# Import the profile generation module

# %%
import tsib

# %% [markdown]
# Create profiles for three, three person households

# %%
res = tsib.simHouseholdsParallel(3,2010, 3, cores = 4, get_hot_water = True)

# %%
res

# %%
