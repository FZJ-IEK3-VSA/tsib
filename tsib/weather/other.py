
"""
Created on Thu May 12 10:22:47 2016

Read in weather data

@author: l.kotzur
"""
import os
import pytz
import warnings
import pvlib
import math

import numpy as np
import pandas as pd

import tsib.data

def readLindenberg(filepath, year=2000):
    """
    Reads the minutely resolved solar irradiation data from the Lindenberg
    Weather station.
    
    Parameters
    ---------
    filpath: str, required
        Path to the Lindenberg datenbasis folder which has two subfolders:
            "Irradiance" and "Weather"
    year: int, optional (default = 2000)
        Year to take the data. Only year 2000 is checked for missing data.
        
    Returns
    -------
    data, pd.DataFrame
        The whole weather data as time indexed DataFrame.
    location, dict
        Dictionary with geographical informations about the side.
    """
    # read in all 12 month
    raw1 = pd.DataFrame()
    raw2 = pd.DataFrame()
    for month in range(1, 13):
        read1 = pd.read_csv(
            os.path.join(
                filepath,
                "Irradiance",
                "LIN_radiation_" + str(year) + "-" + "{:02d}".format(month) + ".tab",
            ), sep="\t",
            header=31,
            index_col=0,
        )
        raw1 = raw1.append(read1)
        read2 = pd.read_csv(
            os.path.join(
                filepath,
                "Weather",
                "LIN_SYNOP_" + str(year) + "-" + "{:02d}".format(month) + ".tab",
            ), sep="\t",
            header=31,
            index_col=0,
        )
        raw2 = raw2.append(read2)

    # clean indices
    raw1.index = pd.to_datetime(raw1.index)
    raw1 = raw1[~raw1.index.duplicated(keep="first")]
    raw2.index = pd.to_datetime(raw2.index)
    raw2 = raw2[~raw2.index.duplicated(keep="first")]

    # merge weather and irradiance data
    data = pd.concat([raw1, raw2], join="outer", axis=1)

    # resample to one minute and interpolate missing values in between
    data = data.resample("1min").mean().interpolate()

    # drop redundant columns
    data = data.loc[:, ~data.columns.duplicated(keep="first")]

    # localize data to lindenberg
    location = {"name": "Lindenberg", "latitude": 52.21, "longitude": 14.12}

    # get correct time index
    data.index = (
        data.index.tz_localize(pytz.utc).tz_convert("Europe/Berlin").shift(-1, freq="H")
    )

    # rename relevant columns
    data = data.rename(
        columns={
            "DIR [W/m**2]": "DNI",
            "DIF [W/m**2]": "DHI",
            "SWD [W/m**2]": "GHI",
            "dd [deg]": "Wdir",
            "ff [m/s]": "Wspd",
            "TTT [°C]": "DryBulb",
        }
    )

    # manage missing data
    if year == 2000:
        data["2000-08-19":"2000-08-22"] = data["2000-08-15":"2000-08-18"].values
    else:
        warnings.warn(
            "Other years than 2010 are not validated in terms " + "of missing data"
        )

    return data, location




def readCosmo(datapath, lon, lat, year):
    """
    Get data from the cosmo data set.

    Parameters
    ----------
    datapath: str
        Path to the netCDF4 files.
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
        raise ImportError('TSA-Internal packages are required.')

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
    filepath = os.path.join(tsib.data.PATH, "weatherdata", "COSMO", identifier)

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
