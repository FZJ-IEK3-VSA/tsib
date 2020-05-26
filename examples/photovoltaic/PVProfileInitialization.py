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
# # Test weather data and pvlib

# %%

import matplotlib.pyplot as plt
# %matplotlib inline
# %load_ext autoreload
# %autoreload 2
import tsib

# %% [markdown]
# ## 1. Simulate PV with TRY data

# %%
try_data, loc = tsib.readTRY(5)

# %%
tmy_data = tsib.TRY2TMY(try_data)

# %%
pv_try, space_cov = tsib.simPhotovoltaic(tmy_data, latitude = loc['latitude'], longitude = loc['longitude'],
                                    losses = 0.1,integrateInverter = True)

# %%
pv_try.sum()

# %%
pv_try.plot()

# %% [markdown]
# ## 2. Simulate with COSMO6 data

# %%
try_data, identifier = tsib.readCosmo(os.path.join(os.environ['DATA_SHARE'],
                                              'weather','cosmo','rea6','processed'),
                                loc['longitude'], loc['latitude'], 2010 )

# %%
tmy_data = tsib.TRY2TMY(try_data)

# %%
pv_try, space_cov = tsib.simPhotovoltaic(tmy_data, latitude = loc['latitude'], longitude = loc['longitude'],
                                    losses = 0.1)

# %%
pv_try.sum()

# %%
