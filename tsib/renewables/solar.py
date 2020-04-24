"""

Created on Sat Dec 10 12:40:17 2017

@author: Leander Kotzur
"""

import os
import pvlib

import pandas as pd

import tsib.data

def simPhotovoltaic(
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
            '\\profiles\\PV-Modules\\sandia_modules.csv'.
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
    temp_model = pvlib.temperature.TEMPERATURE_MODEL_PARAMETERS["sapm"]["open_rack_glass_glass"]
    pvtemps = pvlib.temperature.sapm_cell(
        poa_irrad["poa_global"], tmy_data["DryBulb"], tmy_data["Wspd"], **temp_model
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
                tsib.data.PATH, "sandia", "pvmodules", "sandia_modules.csv"
            ),
            index_col=0,
        )
        module = modules[module_name]
        module = pd.to_numeric(module, errors="coerce")

        inverters = pd.read_csv(
            os.path.join(
                tsib.data.PATH, "sandia", "inverters", "sandia_modules.csv"
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
        sapm_irr, module=module, temp_cell=pvtemps,
    )

    # calculate peak load of single module [W]
    peak_load = module.loc["Impo"] * module.loc["Vmpo"]

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
    space_coverage = module.loc["Area"] / peak_load * 1000.0
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
