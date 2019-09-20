# -*- coding: utf-8 -*-
"""
Created on Thu May 12 10:22:47 2016

Radiation Converter

@author: l.kotzur
"""
import pvlib
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import copy
import pytz
import warnings
import pvlib

from datetime import date

# ignore pv lib warnings
np.seterr(divide="ignore")
np.seterr(invalid="ignore")

DATA_PATH = os.path.join(os.path.dirname(os.path.realpath(__file__)), "data")


def calcGHI(timeSeries, longitude, latitude):
    """
        Calculates the global horizontal irradiation (GHI) by time, diffuse
        horizontal irradiation (DHI) and direct normal irradiation (DNI).
        """
    solarPos = pvlib.solarposition.get_solarposition(
        timeSeries.index, latitude, longitude
    )
    if "DNI" in timeSeries:
        timeSeries["GHI"] = timeSeries["DHI"] + timeSeries["DNI"] * math.cos(
            math.radians(solarPos["apparent_zenith"])
        )
    elif "B" in timeSeries:  # Direct horizontal irradiation TRY
        timeSeries["GHI"] = timeSeries["DHI"] + timeSeries["B"]
    return timeSeries


def readTMY(filepath=os.path.join("TMY", "Germany DEU Koln (INTL).csv")):
    """
    Reads a typical meteorological year file and gets the GHI, DHI and DNI from it.
    """
    # get data
    data = pd.read_table(
        os.path.join(DATA_PATH, "weatherdata", filepath),
        skiprows=([0, 1]),
        delimiter=",",
    )
    data.index = pd.date_range(
        "2010-01-01 00:30:00", periods=8760, freq="H", tz="Europe/Berlin"
    )
    data = data.rename(
        columns={"Beam": "DNI", "Diffuse": "DHI", "Tdry": "T", "Wspd": "WS"}
    )
    location_data = pd.read_table(
        os.path.join(DATA_PATH, "profiles", filepath), nrows=1, delimiter=","
    )

    location = {
        "name": location_data["City"].values[0],
        "latitude": location_data["Latitude"].values[0],
        "longitude": location_data["Longitude"].values[0],
    }
    return data, location


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
        read1 = pd.read_table(
            os.path.join(
                filepath,
                "Irradiance",
                "LIN_radiation_" + str(year) + "-" + "{:02d}".format(month) + ".tab",
            ),
            header=31,
            index_col=0,
        )
        raw1 = raw1.append(read1)
        read2 = pd.read_table(
            os.path.join(
                filepath,
                "Weather",
                "LIN_SYNOP_" + str(year) + "-" + "{:02d}".format(month) + ".tab",
            ),
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


def readTRY(try_num=4, year=2010):
    """
    Reads a test refence year file and gets the GHI, DHI and DNI from it.
    
    Parameters
    -------
    try_num: int (default: 4)
        The region number of the test reference year.
    year: int (default: 2010)
        The year. Only data for 2010 and 2030 available
    """
    # get the correct file path
    filepath = os.path.join(
        DATA_PATH,
        "weatherdata",
        "TRY",
        "TRY" + str(year) + "_" + str(try_num).zfill(2) + "_Jahr",
    )

    # get the geoposition
    with open(filepath + ".dat", encoding="utf-8") as fp:
        lines = fp.readlines()
        location_name = lines[1][9:-18].encode("utf-8").rstrip()
        lat = float(lines[2][6:8]) + float(lines[2][9:11]) / 60.0
        lon = float(lines[2][21:23]) + float(lines[2][24:26]) / 60.0
    location = {"name": location_name, "latitude": lat, "longitude": lon}

    # check if time series data already exists as .csv with DNI
    if os.path.isfile(filepath + ".csv"):
        data = pd.read_csv(filepath + ".csv", index_col=0, parse_dates=True)
        data.index = pd.to_datetime(data.index, utc=True).tz_convert("Europe/Berlin")
    # else read from .dat and calculate DNI etc.
    else:
        # get data
        data = pd.read_table(
            filepath + ".dat", sep="\s+", skiprows=([i for i in range(0, 36)] + [37])
        )
        data.index = pd.date_range(
            "2010-01-01 00:30:00", periods=8760, freq="H", tz="Europe/Berlin"
        )
        data["GHI"] = data["D"] + data["B"]
        data = data.rename(columns={"D": "DHI", "t": "T", "WG": "WS"})

        # calculate direct normal
        data["DNI"] = calculateDNI(data["B"], lon, lat)

        # save as .csv
        data.to_csv(filepath + ".csv")
    return data, location


def readtCosmoNetCDF4(datapath, lon, lat, year):
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


def calculateDNI(directHI, lon, lat, zenith_tol=87.0):
    """
    Calculates the direct NORMAL irradiance from the direct horizontal irradiance with the help of the PV lib.

    Parameters
    ----------
    directHI: pd.Series with time index
        Direct horizontal irradiance
    lon: float
        Longitude of the loaction
    lat: float
        Latitude of the location
    zenith_tol: float, optional
        Avoid cosinus of values above a certain zenith angle of in order to avoid division by zero.

    Returns
    -------
    DNI: pd.Series
    """
    solarPos = pvlib.solarposition.get_solarposition(directHI.index, lat, lon)
    solarPos["apparent_zenith"][solarPos.apparent_zenith > zenith_tol] = zenith_tol
    DNI = directHI.div(solarPos["apparent_zenith"].apply(math.radians).apply(math.cos))
    if DNI.isnull().values.any():
        raise ValueError("Something went wrong...")
    return DNI


def TRY2TM2(trydata):
    """
    Takes a pd.DataFrame from the readTRY function and translate the column
    names to the TM2 format.
    """
    trydata["Year"] = 2010
    trydata = trydata.rename(
        columns={
            "MM": "Month",
            "DD": "Day",
            "HH": "Hour",
            "T": "Tdry",
            "p": "Pres",
            "WR": "Wdir",
            "WS": "Wspd",
        }
    )
    trydata["Tdew"] = trydata["Tdry"]

    return trydata[
        [
            "Year",
            "Month",
            "Day",
            "Hour",
            "GHI",
            "DHI",
            "Tdry",
            "Tdew",
            "Pres",
            "Wdir",
            "Wspd",
        ]
    ]


def TRY2TMY(trydata):
    """
    Takes a pd.DataFrame from the readTRY function and translate the column
    names to the TM2 format.
    """
    return trydata.rename(
        columns={
            "MM": "Month",
            "DD": "Day",
            "HH": "Hour",
            "T": "DryBulb",
            "p": "Pressure",
            "WR": "Wdir",
            "WS": "Wspd",
        }
    )


def simPV_PV_Lib(
    tmy_data,
    surface_tilt=30,
    surface_azimuth=180,
    albedo=0.2,
    latitude=55,
    longitude=7,
    losses=0.1,
    load_module_data=False,
    module_name="Hanwha_HSL60P6_PA_4_250T__2013_",
    integrateInverter=True,
    inverter_name="ABB__MICRO_0_25_I_OUTD_US_208_208V__CEC_2014_",
):
    """
    Simulates a defined PV array with the Sandia PV Array Performance Model.
    The implementation is done in accordance with following tutorial:
    https://github.com/pvlib/pvlib-python/blob/master/docs/tutorials/tmy_to_power.ipynb
    
    Parameters
    ----------
    tmy_data: pandas.DataFrame(), required
        Weatherfile in the format of a tmy file.
    surface_tilt: int or float, optional (default:30)
        Tilt angle of of the array in degree.
    surface_azimuth: int or float, optional (default:180)
        Azimuth angle of of the array in degree. 180 degree means south,
        90 degree east and 270 west.
    albedo: float, optional (default: 0.2)
        Reflection coefficient of the sorounding area.    
    latitude: float, optional (default: 55)
        Latitude of the position in degree.
    longitude: float, optional (default: 7)
        Longitude of the position in degree.
    losses: float, optional (default: 0.1)
        Losses due to soiling, mismatch, diode connections, dc wiring etc.
    load_module_data: Boolean, optional (default: False)
        If True the module data base is loaded from the Sandia Website.
        Otherwise it is loaded from this relative path
            '\profiles\PV-Modules\sandia_modules.csv'.
    module_name: str, optional (default:'Hanwha_HSL60P6_PA_4_250T__2013_')
        Module name. The string must be existens in Sandia Module database.
    integrateInverter: bool, optional (default: True)
        If an inverter shall be added to the simulation, providing the photovoltaic output after the inverter.
    inverter_name: str, optional (default: 'ABB__MICRO_0_25_I_OUTD_US_208_208V__CEC_2014_')
        Type of inverter.
        
    Returns
    --------
    specific_load: pandas.DataFrame
        Timeseries of the specific load of such a PV-module per installed kWp. [kW/kWp]
    space_coverage: float
        Space required for the an installed capacity of 1 kWp [m^2/kWp]
    """

    # calculate the solar position for all times in the TMY file
    solpos = pvlib.solarposition.get_solarposition(tmy_data.index, latitude, longitude)

    # calculate extra terrestrial radiation- n eeded for perez array diffuse irradiance models
    dni_extra = pd.Series(
        pvlib.irradiance.get_extra_radiation(tmy_data.index), index=tmy_data.index
    )  # automatic pd time series in future pvlib version
    # calculate airmass
    airmass = pvlib.atmosphere.get_relative_airmass(solpos["apparent_zenith"])
    # use pereze model to calculate the plane of array diffuse sky radiation
    poa_sky_diffuse = pvlib.irradiance.perez(
        surface_tilt,
        surface_azimuth,
        tmy_data["DHI"],
        tmy_data["DNI"],
        dni_extra,
        solpos["apparent_zenith"],
        solpos["azimuth"],
        airmass,
    )
    # calculate ground diffuse with specified albedo
    poa_ground_diffuse = pvlib.irradiance.get_ground_diffuse(
        surface_tilt, tmy_data["GHI"], albedo=albedo
    )
    # calculate angle of incidence
    aoi = pvlib.irradiance.aoi(
        surface_tilt, surface_azimuth, solpos["apparent_zenith"], solpos["azimuth"]
    )
    # calculate plane of array irradiance
    poa_irrad = pvlib.irradiance.poa_components(
        aoi, tmy_data["DNI"], poa_sky_diffuse, poa_ground_diffuse
    )
    # calculate pv cell and module temperature
    pvtemps = pvlib.pvsystem.sapm_celltemp(
        poa_irrad["poa_global"], tmy_data["Wspd"], tmy_data["DryBulb"]
    )

    # load the sandia data
    if load_module_data:
        # load module data online
        modules = pvlib.pvsystem.retrieve_sam(name="SandiaMod")
        module = modules[module_name]
        # get inverter data
        inverters = pvlib.pvsystem.retrieve_sam("cecinverter")
        inverter = inverters[inverter_name]
    else:
        # load module and inverter data from csv
        modules = pd.read_csv(
            os.path.join(
                DATA_PATH, "additionaldata", "pvmodules", "sandia_modules.csv"
            ),
            index_col=0,
        )
        module = modules[module_name]
        module = pd.to_numeric(module, errors="coerce")

        inverters = pd.read_csv(
            os.path.join(
                DATA_PATH, "additionaldata", "inverters", "sandia_modules.csv"
            ),
            index_col=0,
        )
        inverter = inverters[inverter_name]
        inverter = pd.to_numeric(inverter, errors="coerce")

    # calculate effective irradiance on pv module
    sapm_irr = pvlib.pvsystem.sapm_effective_irradiance(
        module=module,
        poa_direct=poa_irrad['poa_direct'],
        poa_diffuse=poa_irrad['poa_diffuse'],
        airmass_absolute=airmass,
        aoi=aoi,
    )
    # calculate pv performance
    sapm_out = pvlib.pvsystem.sapm(
        sapm_irr, module=module, temp_cell=pvtemps["temp_cell"]
    )

    # calculate peak load of single module [W]
    peak_load = module.ix["Impo"] * module.ix["Vmpo"]

    if integrateInverter:
        # calculate load after inverter
        inv_load = pvlib.pvsystem.snlinverter(
            inverter=inverter, v_dc=sapm_out["v_mp"], p_dc=sapm_out["p_mp"]
        )
        # load in [kW/kWp]
        specific_load = inv_load / peak_load
    else:
        # load in [kW/kWp]
        specific_load = sapm_out["p_mp"] / peak_load

    # check for NaN entries and fill with zero
    if specific_load.isnull().any():
        # warnings.warn(str(sapm_out['p_mp'].isnull().sum()) + ' NaNs get replaced with zero.')
        specific_load = specific_load.fillna(0)

    # replace negative values
    specific_load[specific_load < 0] = 0

    # round small values to avoid numerical trouble in optimization
    specific_load = specific_load.round(decimals=5)

    # add losses (soiling, mismatch, diodes and connections, dc wiriing)
    specific_load = specific_load * (1 - losses)

    # space coverage in [m^2/kW]
    space_coverage = module.ix["Area"] / peak_load * 1000.0
    return specific_load, space_coverage


def simSolarThermal(
    tmy_data,
    surface_tilt=30,
    surface_azimuth=180,
    albedo=0.2,
    latitude=55,
    longitude=7,
    module_temp=45.0,
    c_0=0.791,
    c_1=4.47,
    c_2=0.0069,
):
    """
    Simulates a defined solar thermal panel. The irradiance calculation
    is based on the PV lib.
    
    Parameters
    ----------
    tmy_data: pandas.DataFrame(), required
        Weatherfile in the format of a tmy file.
    surface_tilt: int or float, optional (default:30)
        Tilt angle of of the array in degree.
    surface_azimuth: int or float, optional (default:180)
        Azimuth angle of of the array in degree. 180 degree means south,
        90 degree east and 270 west.
    albedo: float, optional (default: 0.2)
        Reflection coefficient of the sorounding area.    
    latitude: float, optional (default: 55)
        Latitude of the position in degree.
    longitude: float, optional (default: 7)
        Longitude of the position in degree.
    module_temp: float, optional (default: 30.)
        Temperature of the solar thermal module or collector
    c_0, c_1, c_2: floats, optional (default: 0.791, 4.47, 0.0069)
        Performance constants of the solar thermal panel.
        Check http://www.spf.ch/index.php?id=111 for performance coefficients.
        Defailt is SPF-Nr.: C1734
        http://www.spf.ch/index.php?id=111
        
    Returns
    --------
    spec_load: pandas.Series
        Timeseries of the specific load of such a solar thermal-module per 
        installed m^2. [kW/m^2]
    """
    import pvlib

    # calculate the solar position for all times in the TMY file
    solpos = pvlib.solarposition.get_solarposition(tmy_data.index, latitude, longitude)

    # calculate extra terrestrial radiation- n eeded for perez array diffuse irradiance models
    dni_extra = pd.Series(
        pvlib.irradiance.get_extra_radiation(tmy_data.index), index=tmy_data.index
    )  # automatic pd time series in future pvlib version
    # calculate airmass
    airmass = pvlib.atmosphere.get_relative_airmass(solpos["apparent_zenith"])
    # use pereze model to calculate the plane of array diffuse sky radiation
    poa_sky_diffuse = pvlib.irradiance.perez(
        surface_tilt,
        surface_azimuth,
        tmy_data["DHI"],
        tmy_data["DNI"],
        dni_extra,
        solpos["apparent_zenith"],
        solpos["azimuth"],
        airmass,
    )
    # calculate ground diffuse with specified albedo
    poa_ground_diffuse = pvlib.irradiance.get_ground_diffuse(
        surface_tilt, tmy_data["GHI"], albedo=albedo
    )
    # calculate angle of incidence
    aoi = pvlib.irradiance.aoi(
        surface_tilt, surface_azimuth, solpos["apparent_zenith"], solpos["azimuth"]
    )

    # calculate plane of array irradiance
    poa_irrad = pvlib.irradiance.poa_components(
        aoi, tmy_data["DNI"], poa_sky_diffuse, poa_ground_diffuse
    )
    # get global poa irradiance in W
    irr = poa_irrad["poa_global"]

    # get delta T between collector and ambient
    d_T = module_temp - tmy_data["DryBulb"]

    # calculate pv cell and module temperature
    eta = c_0 - (c_1 * d_T / irr) - (c_2 * (d_T ** 2) / irr)

    # replace negative efficiencies with zero
    eta[eta < 0.0] = 0.0

    # load in [kW/m^2]
    spec_load = eta * irr / 1e3

    # replace nans with zeros
    return spec_load.fillna(0.0)


def createWoodFireProfile(
    temperature,
    occ_act,
    n_ovens=1,
    T_oven_on=5,
    t_cool=5.0,
    fullloadSteps=450,
    seed=None,
):
    """
    Creates the profile of the heating of wood ovens based on the outside
    temperature and the activity of the occupants. The profile is generated
    based on stochastics.
    
    Parameters
    ----------
    temperature: pandas.Series(), required
        Outside temperature profile of the location.
    occ_act: pandas.Series(), required
        Series of values between 0 and 1 desribing the share of overall active
        occpants at every time step.
    n_ovens: int, optional (default:1)
        Number of ovens in the building.
    T_oven_on: int or float, optional (default:5.)
        Threeshold outside temperature [°C] when the oven is turned on.
    t_cool: int, optional (default:5)
        Number of timesteps till when the oven is cooled down again.
    fullloadSteps: int or float, optional (default:450)
        Resulting number of full load timesteps. Attention: This value is not
        exact since it is a stochastic profile generation.
    seed: int, optional (default:None)
        Seed required for reproduceability. 
        If none, it is completely random.
        
    Returns
    ----------
    overAllLoad: pandas.Series()
        Relative load profile of the ovens in kW/kWp.
    """

    # Oven is only turned under a temperature threeshold
    tempBool = temperature < T_oven_on

    # Increase probability that the oven is turned on in case it is colder outside
    relCold = (T_oven_on - temperature) / (T_oven_on - temperature.min())

    # Caclulate fire activation probability
    prob = occ_act.values * tempBool.values * relCold.values

    # avoid rounding errors
    prob[prob < 0] = 0

    overAllLoad = pd.Series(0, index=temperature.index)
    for n in range(int(n_ovens)):

        # Modifier to reduce probability in order to fit the full load hours
        p_mod = fullloadSteps / (prob.sum() * t_cool / 2)

        overallProb = prob * p_mod
        overallProb[overallProb > 1.0] = 1.0

        # Binary decision if an oven can be activated
        initLogArr = pd.Series(np.random.RandomState(seed).binomial(1, overallProb))

        # create the profile
        heatLoad = []
        loadBefore = 0
        for initLog in initLogArr:
            if initLog:
                load = 1.0
            else:
                if loadBefore > 0.001:
                    load = loadBefore - 1.0 / t_cool
                else:
                    load = 0
            heatLoad.append(load)
            loadBefore = load
        overAllLoad += pd.Series(heatLoad, index=temperature.index)
    profile = overAllLoad / n_ovens
    if abs(profile.sum() - fullloadSteps) > 100:
        warnings.warn(
            "Fullload hour deviation is higher than 100. "
            + "Input parameters make it difficult or impossible "
            + "to generate the expected profile"
        )
    return overAllLoad / n_ovens


def groupToPeriods(timeSeries, periodStepNumber):
    """
    Extend the timeseries to an integer multiple of the period length and
    groups the time series to the periods.

    Parameters
    -----------
    timeSeries, required
        pandas.DataFrame()
    periodStepNumber: integer, required
        The number of discrete timesteps which describe one period.
        
    Returns
    -------
    timeseries
        pandas.DataFrame() which is stacked such that each row represents a
        candidate period
    """
    group_index = []
    column_index = []
    if len(timeSeries) % periodStepNumber == 0:
        attached_timesteps = 0
    else:
        #            raise ValueError()
        attached_timesteps = periodStepNumber - len(timeSeries) % periodStepNumber
        rep_data = timeSeries.head(attached_timesteps)
        timeSeries = timeSeries.append(rep_data, ignore_index=True)

    for ii in range(0, len(timeSeries)):
        group_index.append(int(ii / periodStepNumber))
        column_index.append(ii - int(ii / periodStepNumber) * periodStepNumber)

    timeSeries.index = pd.MultiIndex.from_arrays(
        [column_index, group_index], names=["TimeStep", "PeriodNum"]
    )
    timeSeries = timeSeries.unstack(level="TimeStep")
    return timeSeries


def plotTimeSeriesFFT(timeSeries, limiter=None):
    """
    Does a Fast Fourier Transformation and returns the absolute values
    over the resulting frequencies.
    
    Parameters
    ------------
    timeSeries: pandas.DataFrame, required
        A pandas.DataFrame with timestamps as index.
    
    """
    # Number of sample points
    N = len(timeSeries)
    # sample spacing
    T = timeSeries.index[1] - timeSeries.index[0]

    x = timeSeries.index
    # normalized
    #    y = (timeSeries.values.astype(float) - timeSeries.min()) / (timeSeries.max() - timeSeries.min())
    y = timeSeries.values.astype(float)
    yf = np.fft.fft(y)
    xf = np.linspace(0.0, 1.0 / (2.0 * T.total_seconds()), N / 2.0)
    spd = 2.0 / N * np.abs(yf[0 : N / 2.0])
    # xf = np.fft.fftfreq(N,T)
    # plot results
    color = np.array([0.0, 83.0, 130.0]) / 255
    ac_color = np.array([192.0, 0.0, 0.0, 255.0]) / 255
    #    plt.figure()
    plt.scatter(xf[1:] * 3600, spd[1:], color=color, marker="o", s=1)
    import math

    if max(spd[1:]) < 0.1:
        ylim = math.ceil(max(spd[1:]) * 100) / 100
    else:
        ylim = math.ceil(max(spd[1:]) * 10) / 10
    ax = plt.gca()  #    ax.set_xscale('log')
    limiter = sorted(copy.deepcopy(spd))[-4]
    for ii, xf_val in enumerate(xf):
        if spd[ii] > limiter and ii > 0:  #
            #            if xf_val>1./(8750.*T.total_seconds()):
            #            raise ValueError()
            plt.scatter(
                xf_val * 3600,
                spd[ii],
                marker="o",
                s=20,
                edgecolors=ac_color,
                facecolors="None",
                linewidth=1,
            )
            ax.text(
                xf_val * 3600,
                ylim / 40 + spd[ii],
                str(int(round(1 / (xf_val * T.total_seconds())))) + " h",
                ha="left",
                va="bottom",
                color=ac_color,
            )
    #            else:
    #                plt.scatter(xf_val*3600,spd[ii],marker = 'o', edgecolors='r',facecolors='None')
    #                ax.text( xf_val*3600 , ylim/40 +spd[ii],
    #                        str(round(T.total_seconds()/3600*len(timeSeries),1))+'h',
    #                        ha='left', va='bottom',color='r')

    #    ax.set_xlabel('Frequency [1/h]')
    ax.set_xlim(
        [1.0 / (T.total_seconds() / 3600 * (N + 1)) / 2, 3600 / T.total_seconds() * 0.5]
    )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ylim_max = 10 ** (np.floor(np.log10(spd.max())) + 1)
    ax.set_ylim([ylim_max * 1e-5, ylim_max])
    #    ax.set_ylim([1e-5,1])
    #    ax.set_ylabel('Amplitude [-]')
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    return


if __name__ == "__main__":
    #    tmy_data, loc = readTMY()

    #    plt.style.use(os.path.join(DATA_PATH, 'utils', 'matplotlibrc.mplstyle'))
    try_data, loc = readTRY()
    #    htw_data = readHTWProfiles()
    raw = pd.DataFrame()
    raw["GHI"] = try_data["GHI"]
    tmy_data = TRY2TMY(try_data)
    specific_load_se, space_coverage_se = simPV_PV_Lib(
        tmy_data,
        surface_tilt=30,
        surface_azimuth=135,
        albedo=0.2,
        latitude=loc["latitude"],
        longitude=loc["longitude"],
        losses=0.1,
        load_module_data=False,
        module_name="Hanwha_HSL60P6_PA_4_250T__2013_",
    )
    spec_load_st = simSolarThermal(tmy_data)
    plotTimeSeriesFFT(raw)
#    raw['T'] =  try_data['T']
#    raw['E-Load'] =  htw_data['Profile 1']+htw_data['Profile 2']+htw_data['Profile 3']+htw_data['Profile 4']
#    raw['E-Load'] =  htw_data.sum(axis=1)

#    tmy_data = TRY2TMY(try_data)
#    specific_load, space_coverage = simPV_PV_Lib(tmy_data,
#                 surface_tilt = 45, surface_azimuth = 90, albedo = 0.2,
#                 latitude = loc['latitude'], longitude = loc['longitude'],
#                 load_module_data = False,
#                 module_name = 'Canadian_Solar_CS5P_220M___2009_')
#    plt.figure()
#    specific_load.plot()
#    plt.figure
#    specific_load.sum()
#    # plot FFT
#    plt.figure(figsize=[9,3])
#    plotTimeSeriesFFT( raw['GHI'], limiter = 0.05 )
#    plt.figure(figsize=[9,3])
#    plotTimeSeriesFFT( raw['T'], limiter = 0.05 )
#    plt.figure(figsize=[9,3])
#    plotTimeSeriesFFT( raw['E-Load'], limiter = 0.01 )
