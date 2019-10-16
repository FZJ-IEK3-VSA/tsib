# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.4'
#       jupytext_version: 1.2.4
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# # Validation of the tsib.timeseriesmanager.get_era5 module
# Author: Leander Kotzur
#
# Date: 21.09.2019

# %matplotlib inline
# %load_ext autoreload
# %autoreload 2

import tsib.timeseriesmanager as tsm
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import copy
import os

from netCDF4 import Dataset,num2date

import cdsapi


# Get weather data

try_data,location = tsm.readTRY()

try_data.head()

location

# ## Get single profiles

filename = 'download2.nc'

# +

c = cdsapi.Client()

c.retrieve(
    'reanalysis-era5-land',
    {
        'format':'netcdf',
        'variable':[
            '2m_temperature',
            'surface_net_solar_radiation',
            'total_sky_direct_solar_radiation_at_surface'
        ],
        'area':str(location['latitude']) + 
                '/' + str(location['longitude']) +
                '/' + str(location['latitude']+0.1) + 
                '/' + str(location['longitude']+0.1),
        'year':'2014',
        'month':[
            '01','02','03',
            '04','05','06',
            '07','08','09',
            '10','11','12'
        ],
        'day':[
            '01','02','03',
            '04','05','06',
            '07','08','09',
            '10','11','12',
            '13','14','15',
            '16','17','18',
            '19','20','21',
            '22','23','24',
            '25','26','27',
            '28','29','30',
            '31'
        ],
        'time':[
            '00:00','01:00','02:00',
            '03:00','04:00','05:00',
            '06:00','07:00','08:00',
            '09:00','10:00','11:00',
            '12:00','13:00','14:00',
            '15:00','16:00','17:00',
            '18:00','19:00','20:00',
            '21:00','22:00','23:00'
        ]
    },
    filename,
)
# -

nc_pvFile = Dataset(filename,'r')
lats = nc_pvFile.variables['latitude'][:]  # extract/copy the data
lons = nc_pvFile.variables['longitude'][:]
lats = lats[:].squeeze()
lons = lons[:].squeeze()

# initialize a pandas.DataFrame to store all the data
time_index = pd.date_range(
    start=str(year) + "-01-01 00:30:00",
    end=str(year) + "-12-31 23:30:00",
    freq="H",
    tz="Europe/Berlin",
)



def get_era5(lon, lat, year):
    """
    Get data from the cosmo data set.

    Parameters
    ----------
    lon: float
        Longitude in degree.
    lat: float
        Latitude in degree.
    year: int
        Year from 2010 to 2015 as integer.

    Returns
    -------
    df: pd.DataFrame with the data.
    datapath: Identifier for the weather data.
    """
    try:
        import res
        import netCDF4 as nc
    except ImportError:
        raise ImportError('Internal packages are required.')

    # required variables and a translation to the terms used in tsib
    VARIABLES2TRY = {
        "SWDIFDS_RAD": "DHI",
        "SWDIRS_RAD": "DirectHI",
        "2t": "T",
        "windspeed_10": "WS",
        "sp": "p",
    }

    # get the data allocation
    data = res.weather.CosmoSource(
        [os.path.join(datapath, str(year), "*.nc")], bounds=(lon, lat), verbose=False
    )

    # get the index for the netCDF4 file
    ix = res.weather.CosmoSource.loc2Index(None, (lon, lat))

    # get a filepath for the resulting weather data set, ix[1]
    identifier = "COSMO_Year_" + str(year) + "_ix_" + str(ix[0]) + "_" + str(ix[1])
    filepath = os.path.join(DATA_PATH, "weatherdata", "COSMO", identifier)

    # get file if existing
    if os.path.isfile(filepath + ".csv"):
        df = pd.read_csv(filepath + ".csv", index_col=0, parse_dates=True)
        df.index = df.index.tz_localize("UTC").tz_convert("Europe/Berlin")

    # get data
    else:
        # initialize a pandas.DataFrame to store all the data
        time_index = pd.date_range(
            start=str(year) + "-01-01 00:30:00",
            end=str(year) + "-12-31 23:30:00",
            freq="H",
            tz="Europe/Berlin",
        )
        df = pd.DataFrame(index=time_index)

        # loop over required variables and load the data into the dataframe
        for var in VARIABLES2TRY.keys():
            # load dataset via netCDF4
            ds = nc.Dataset(data.variables.loc[var, "path"])
            df[VARIABLES2TRY[var]] = ds.variables[var][:, ix[0], ix[1]]

        # convert temperature from K to °C
        df["T"] = df["T"] - 273.15

        # get GHI
        df["GHI"] = df["DirectHI"] + df["DHI"]

        # calclulate DNI
        df["DNI"] = calculateDNI(df["DirectHI"], lon, lat)

        # remove leap day
        sel = np.logical_and((time_index.day == 29), (time_index.month == 2))
        df = df.loc[~sel]

        # save as .csv
        df.to_csv(filepath + ".csv")
    return df, identifier


nc_pvFile.variables

t2m_raw = nc_pvFile.variables['t2m'][:]

t2m = t2m_raw[:,0,0]


len(t2m)

# from netCDF4 import Dataset,num2date
# import cartopy.crs as ccrs
# from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
# from cartopy.util import add_cyclic_point
#
#
# preciPlot = nc_pvFile.variables['tp'][:]
#
# pp = preciPlot[0,:,:]
#
#
# ax1 = plt.axes(projection=ccrs.PlateCarree(central_longitude=180))
# clevs = np.arange(min(pp.flatten()),max(pp.flatten())*1000,1)
#
#
# shear_fill = ax1.contourf(lons,lats,pp*1000,clevs,
#                       transform=ccrs.PlateCarree(),cmap=plt.get_cmap('hsv'),
#                                                                       linewidth=(10,),levels=100,extend='both')
#
# ax1.coastlines()
# ax1.gridlines()
# ax1.set_xticks([0, 10,20,30,40,50,60,70,80,90,100], crs=ccrs.PlateCarree())
# ax1.set_yticks([0, 10,20,30,40,50,60], crs=ccrs.PlateCarree())
# lon_formatter = LongitudeFormatter(zero_direction_label=True,
#                                number_format='.0f')
# lat_formatter = LatitudeFormatter()
# ax1.xaxis.set_major_formatter(lon_formatter)
# ax1.yaxis.set_major_formatter(lat_formatter)
# cbar = plt.colorbar(shear_fill, orientation='horizontal')
# plt.title('Total precipitation', fontsize=16)
# plt.savefig('precip_era.png')
