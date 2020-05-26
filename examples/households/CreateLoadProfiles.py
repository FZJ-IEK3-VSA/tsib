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
# # Creates all maximal required load profiles with tsorb and tsib

# %%
import os
import time
import datetime

# %%
import pandas as pd

# %%
import tsib
import tsib

# %%
import multiprocessing as mp

# %% [markdown]
# Number of household size

# %%
HH_SIZE = [1,2,3,4,5]

# %% [markdown]
# Get time series for Potsdam

# %%
weather, loc = tsib.readTRY(try_num=4)

# %% [markdown]
# Get a time dataframe for performance checking

# %%
timeFrame = pd.DataFrame(0, columns = ['Start time','End time','Duration [s]'], index = HH_SIZE)

# %% [markdown]
# Generate for all seeds the profiles, if not existant

# %%
SEEDS = [x for x in range(tsib.buildingmodel.TOTAL_PROFILE_NUM+1)]

# %%
for n_persons in HH_SIZE:
    # get starttime
    startTime = time.time()
    timeFrame.loc[n_persons,'Start time'] = datetime.datetime.now()
    
    tsib.getHouseholdProfiles(n_persons,
                               weather, '4', seeds=SEEDS,
                               ignore_weather=True,
                               mean_load=True,
                              cores = mp.cpu_count() - 2)
    
    # get endtime
    endTime = time.time()
    timeFrame.loc[n_persons,'End time'] = datetime.datetime.now()
    timeFrame.loc[n_persons,'Duration [s]'] = endTime-startTime
    
    print('FINISHED - Duration for '+ str(n_persons) + ' persons: ' + str(endTime-startTime))

# %%
timeFrame.to_csv(os.path.join(os.getcwd(),'timeDuration.csv'))
