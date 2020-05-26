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
# # tsib - Example
# Example usage of the Time Series Initialization of Buildings (tsib) library
# Date: 18.11.2019
#
# Author: Leander Kotzur

# %% [markdown]
# Import required packages

# %%
# %load_ext autoreload
# %autoreload 2
import copy
import tsib
import matplotlib.pyplot as plt
import tsam.timeseriesaggregation as tsam
# %matplotlib inline

# %% [markdown]
# Get a plotting function

# %%
def plotTS(data, periodlength, label):
    fig, axes = plt.subplots(figsize = [6, 2], dpi = 100, nrows = 1, ncols = 1)
    stacked, timeindex = tsam.unstackToPeriods(copy.deepcopy(data), periodlength)
    cax = axes.imshow(stacked.values.T, interpolation = 'nearest', 
                      vmin = stacked.min().min(), vmax = stacked.max().max())
    axes.set_aspect('auto')  
    axes.set_ylabel('Hour')
    plt.xlabel('Day')

    fig.subplots_adjust(right = 1.2)
    cbar=plt.colorbar(cax)    
    cbar.set_label(label)


# %% [markdown]
# ### Define the building

# %% [markdown]
# Get a building configuration

# %%
cfg = tsib.BuildingConfiguration({
    "country": 'DE', 
    "buildingYear": 1990,
    "latitude": 50.0,
    "longitude": 8.0,
    "n_persons": 2,
    "a_ref": 150.,
    "n_apartments":1,
    "surrounding":"Detached",
    "mean_load": True,
    "occControl": False,
    
})

# %% [markdown]
# Parameterize the building model itself with the configuration

# %%
bdg = tsib.Building(configurator = cfg)

# %% [markdown]
# ### 2. Show the weather data as basis

# %% [markdown]
# Get current time series data. Until here, only weather data is provided.

# %%
bdg.timeseries.head()

# %% [markdown]
# Plot the weather data

# %%
for name in bdg.timeseries:
    plotTS(bdg.timeseries[name], 24, 
           str(name) +  ' [' + bdg.units[name] + ']')

# %% [markdown]
# ### 3. Get occupancy data

# %% [markdown]
# Get the data

# %%
occData = bdg.getOccupancy()

# %% [markdown]
# Show occupancy related data

# %%
for name in occData:
    plotTS(bdg.timeseries[name], 24, 
           str(name) +  ' [' + bdg.units[name] + ']')

# %% [markdown]
# ### 4. Get heat load data

# %% [markdown]
# Simulate and get the heat load data

# %%
heatLoad = bdg.getHeatLoad()

# %% [markdown]
# Plot the heat load data

# %%
for name in heatLoad:
    if heatLoad[name].max() > 1e-4:
        plotTS(bdg.timeseries[name], 24, 
               str(name) +  ' [' + bdg.units[name] + ']')

# %% [markdown]
# ### 5. Get renewable potential data

# %% [markdown]
# Simulate and get the data

# %%
renewables = bdg.getRenewables()

# %% [markdown]
# Plot the renewable potential

# %%
for name in renewables:
    plotTS(bdg.timeseries[name], 24, 
           str(name) +  ' [' + bdg.units[name] + ']')
