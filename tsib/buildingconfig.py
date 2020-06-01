# -*- coding: utf-8 -*-
"""
Created on Wed Oct 11 21:17:37 2017

@author: Leander Kotzur
"""

import os
import copy
import warnings
import logging

import pandas as pd
import numpy as np

import tsib
import tsib.data


HEAT_TECHS = [
    "Oil boiler",
    "Gas boiler",
    "Heat pump",
    "Pellet boiler",
    "Electric heater",
    "District heating",
]

KWARG_TYPES = {
    "country": ['AT', 'BE', 'BG', 'CY', 'CZ', 'DE', 'DK', 'ES', 'FR', 'GB', 'GR',
       'HU', 'IE', 'IT', 'NL', 'NO', 'PL', 'RS', 'SE', 'SI', 'XX'], # country code used for the choice of the building parameters
    "n_storey": "NOT_IMPLEMENTED",  # number of storeys in the building
    "a_ref_app": float,  # living area of a single flat
    "a_ref": float,  # living area of the whole building
    "n_apartments": int,  # number of flats
    "surrounding": [
        "Detached",
        "Semi",
        "Terraced",
    ],  # if the house is standing alone or is surrounded
    "ID": str,  # direct identification of an IWU building by its ID
    "buildingYear": int,  # construction year
    "buildingType": [
        "AB",
        "SFH",
        "MFH",
        "TH",
    ],  # type of building (can be inherited by n_apartments)
    "buildingClassification": ["Gen", "HR", "Lightframe"],
    "eastOrOverall": ["N", "East"],
    "roofOrientation": float,  # roof azimuth in degree with 180 as south
    "roofTilt": float,  # rooftile angle in degree
    "refurbished": bool,  # if the building is already refurbished
    "buildnew": bool,  # if the building gets completely new constructed
    "onlyEnergyInvest": bool,  # if the cost of refurbishing the walls and the roof are only energy related
    "thermalClass": "NOT_IMPLEMENTED",  # ['very light', 'light', 'medium', 'heavy', 'very heavy'],
    "refurbishment": bool,  # if refurbishment options (changing insulation or ventilation) shall be considered
    "force_refurbishment": bool,  # if refurbishment must be selected
    "hotWaterElec": bool,  # if hot water is electrically provided
    "existingHeatSupply": [
        "Oil boiler",
        "Gas boiler",
        "Heat pump",
        "Pellet boiler",
        "Electric heater",
        "CHP",
        "District heating",
    ],
    "replaceHeatSupply": bool,  # if heat supply is at the end if the life time
    "T_sup": float,  # design supply temperature of the building
    "floorHeating": bool,  # if a floor heating is available --> set the supply temperature
    "ownership": bool,  # if the occupant is also the owner
    "WACC": float,  # interest rate - otherwise inherited from ownership
    "weatherData": pd.DataFrame,  # time series with the weather
    "weatherID": str,  # identifier of the chosen weather data
    "getNetCDF4": bool,  # if chosen
    "year": int,  # year for which it shall get optimized
    "longitude": float,  # longitude in degree
    "latitude": float,  # latitude in degree
    "nightReduction": bool,  # night reduction of comfort temperature
    "occControl": bool,  # if the heating is adapted to occupancy activity
    "capControl": bool,  # if the heat storage capacity of the bdg can be used
    "comfortT_lb": float,  # lower bound of the comfortable temperature if active
    "comfortT_ub": float,  # upper bound of the comfortable temperature if active
    "n_persons": int,  # number of persons living in a single flat
    "elecLoad": pd.Series,  # electricity load profile of a single flat with correct time index
    "elecLoadID": str,  # identifier of the load
    "hasFirePlace": bool,  # if the building has a fire place
    "hasSolarThermal": bool,  # if the building has solar thermal to provide hot water
    "hasPhotovoltaic": bool,  # if the building has a photovoltaic panel
    "varyoccupancy": int,  # for how many occupancy profiles the building shall be optimized
    "mean_load": bool,  # if the fluctuative profile or the mean hourly profile should be taken
    "a_roof": "NOT_IMPLEMENTED",  # the total roof area
    "windows_refurbished": "NOT_IMPLEMENTED",  # if the windows have allready been replaced
    "walls_refurbished": "NOT_IMPLEMENTED",  # if the walls have already gotton an additional insulation
    "roof_refurbished": "NOT_IMPLEMENTED",  # if the roof area has already gotton an additional insulation
    "costdata": str,  # file identifier with the related cost data
    "ventControl": bool, # if the ventilation system can be smart controlled
}

KWARG_DEFAULTS = {
    "country": 'DE', 
    "roofOrientation": 135.0,  # roof azimuth with 180 as south
    "refurbished": False,  # if the building is already refurbished
    "buildnew": False,  # if the building is newly constructed
    "onlyEnergyInvest": False,  # if the cost of refurbishing the walls and the roof are only energy rlated
    "thermalClass": "medium",
    "refurbishment": False,  # if refurbishment options (changing insulation) shall be considered
    "force_refurbishment": False,  # if refurbishment must be selected
    "hotWaterElec": False,  # if hot water is electrically provided
    "existingHeatSupply": "Oil boiler",
    "buildingYear": 1990,  # construction year
    "replaceHeatSupply": True,  # if heat supply is at the end if the life time
    "hasPhotovoltaic": False,  # if it exists already a photovoltaic panel
    "floorHeating": False,  # if a floor heating is available --> set the supply temperature
    "ownership": True,  # if the occupant is also the owner
    "getNetCDF4": False,  # if chosen
    "year": 2010,  # year for which it shall get optimized
    "longitude": 8.0,  # longitude in degree
    "latitude": 50.0,  # latitude in degree
    "nightReduction": True,  # night reduction of comfort temperature
    "capControl": True,  # if the heating capacity of the building can be used
    "occControl": False,  # if the heating is adapted to occupancy activity
    "comfortT_lb": 21.0,  # lower bound of the comfortable temperature if active
    "comfortT_ub": 24.0,  # upper bound of the comfortable temperature if active
    "n_persons": 2,  # number of persons living in a single flat
    "varyoccupancy": 1,  # for how many occupancy profiles the building shall be optimized
    "mean_load": False,  # if the fluctuative profile or the mean hourly profile should be taken
    "costdata": "default_2016",
    "ventControl": False, # if the ventilation system can be intelligently operated
}


class BuildingConfiguration(object):
    """
    Class which configures building types and decides based on the provision
    of parameters dynamically how to set the building ID and which parameters
    from which database need to get added.


    """

    def __init__(self, kwargs, database=None, ignore_profiles=False):
        """
        Initialize a unique building configuration which provides
        database identifier and configuration parameters based on the user
        input.
        
        Parameters
        ----------
        ignore_profiles: bool, optional (default = False)
            If the profiles should or should not be loaded during
            initialization. Basically this is only required for database
            access.
        """
        if not isinstance(kwargs, dict):
            raise ValueError("kwargs needs to be dictionary with all the parameters")

        # validate if the input parameters are correct
        for kwarg in list(kwargs):
            if not isinstance(kwarg, str):
                raise ValueError("Keyword identifiers need to be strings")
            # drop not defined keywords from dictionary
            elif kwargs[kwarg] is None:
                kwargs.pop(kwarg)
            elif kwarg in KWARG_TYPES:
                # check if implemented
                if KWARG_TYPES[kwarg] == "NOT_IMPLEMENTED":
                    raise NotImplementedError(
                        kwarg + " is not yet considered as keyword"
                    )
                # check for keywords which can only be a limited set of values
                elif type(KWARG_TYPES[kwarg]) is list:
                    if not kwargs[kwarg] in KWARG_TYPES[kwarg]:
                        raise ValueError(
                            "'"
                            + kwarg
                            + "' needs to be one of "
                            + "the following: "
                            + "{}".format(KWARG_TYPES[kwarg])
                        )
                    else:
                        pass  # is valid
                # check the data type of the keyword argument
                elif KWARG_TYPES[kwarg] is float:
                    if not isinstance(
                        kwargs[kwarg], (np.float16, np.float32, np.float64, float)
                    ):
                        raise ValueError(
                            "'" + kwarg + "' needs to be of " + str(KWARG_TYPES[kwarg])
                        )
                    else:
                        # format from numpy to python float
                        if not isinstance(kwargs[kwarg], float):
                            kwargs[kwarg] = kwargs[kwarg].item()  
                        else:
                            pass # is valid
                elif KWARG_TYPES[kwarg] is int:
                    if not isinstance(
                        kwargs[kwarg], (np.int8, np.int16, np.int32, np.int64, int)
                    ):
                        raise ValueError(
                            "'" + kwarg + "' needs to be of " + str(KWARG_TYPES[kwarg])
                        )
                    else:
                        # format from numpy to python int
                        if not isinstance(kwargs[kwarg], int):
                            kwargs[kwarg] = kwargs[kwarg].item()  
                        else:
                            pass # is valid
                elif KWARG_TYPES[kwarg] is bool:
                    if not isinstance(kwargs[kwarg], (bool, np.bool, np.bool_)):
                        raise ValueError(
                            "'" + kwarg + "' needs to be of " + str(KWARG_TYPES[kwarg])
                        )
                    else:
                        # format from numpy to python bool
                        if not isinstance(kwargs[kwarg], bool):
                            kwargs[kwarg] = kwargs[kwarg].item()  
                        else:
                            pass # is valid
                elif not isinstance(kwargs[kwarg], KWARG_TYPES[kwarg]):
                    raise ValueError(
                        "'" + kwarg + "' needs to be of " + str(KWARG_TYPES[kwarg])
                    )
                else:
                    pass  # is valid
            else:
                raise ValueError(kwarg + " is not a valid keyword argument")

        # fill some of the other kwargs with default values
        self.inputKwargs = copy.deepcopy(kwargs)
        for def_kwarg in KWARG_DEFAULTS:
            if not def_kwarg in self.inputKwargs:
                self.inputKwargs[def_kwarg] = KWARG_DEFAULTS[def_kwarg]

        self.IDentries = {}
        # init building configurator
        self.cfg = {}

        self.ignore_profiles = ignore_profiles

        if not database is None:
            raise NotImplementedError()

        return


    def getBdgCfg(self, includeSupply=True):
        """
        Returns a dictionary which is used, either to parameterize
        the buildingmodel.Building or the buildingopt.BuildingOpt

        Parameters
        ----------
        includeSupply: bool, optional (default: True)
            Decides wether also parameters for the supply optimization
            should be included.
        """
        if self.cfg:
            return self.cfg
        else:
            # get the iwu database
            raw = pd.read_csv(
                os.path.join(tsib.data.PATH, "episcope", "episcope.csv"), index_col=1,
            )
            # reduce to the country list of buildings
            self.iwu_bdg = raw[raw["Code_Country"] == self.inputKwargs.pop("country")]

            # call all functions which populate the building configurator
            cfg = self.cfg
            cfg = self._get_form(cfg, self.inputKwargs)
            cfg = self._get_fabric(cfg, self.inputKwargs)
            cfg = self._get_operation(cfg, self.inputKwargs)
            if includeSupply:
                cfg = self._get_equipment(cfg, self.inputKwargs)
            cfg = self._get_finance(cfg, self.inputKwargs)

            # check if cost data file exists
            if not os.path.exists(
                os.path.join(
                    tsib.data.PATH, "costdata", self.inputKwargs["costdata"] + ".xlsx"
                )
            ):
                raise ValueError(
                    "'costdata' file with name '"
                    + self.inputKwargs["costdata"]
                    + "' does not exists"
                )
            # add cost data
            cfg["costdata"] = self.inputKwargs.pop("costdata")
            cfg["costdatapath"] = os.path.join(tsib.data.PATH, "costdata", self.cfg["costdata"] + ".xlsx")
            self.IDentries["costdata"] = cfg["costdata"]

            # check if unused kwargs are left
            for remaining_kwg in self.inputKwargs:
                if remaining_kwg in KWARG_DEFAULTS:
                    if not self.inputKwargs[remaining_kwg] == KWARG_DEFAULTS[remaining_kwg]:
                        warnings.warn('Keyword ' + str(remaining_kwg) + ' is not used for the building parameterization.\n')
                    else:
                        logging.info('Keyword ' + str(remaining_kwg) + ' is not used. Nevertheless, it just holds the default value.\n')
                else:
                        warnings.warn('Keyword ' + str(remaining_kwg) + ' is not used for the building parameterization.\n')
            self._has_cfg = True

            return cfg

    def _get_operation(self, cfg, kwgs):
        """
        Gets all parameters which are relevant for the operation of the supply and heating system.
        """
        cfg["latitude"] = kwgs.pop("latitude")
        cfg["longitude"] = kwgs.pop("longitude")

        # required weatherdata
        weather_units = {"DHI": 'W/m^2', "T": '°C', "DNI": 'W/m^2'}
        cfg["weatherUnits"] = weather_units
        
        if "weatherData" in kwgs:
            if not "weatherID" in kwgs:
                raise ValueError(
                    "If weatherData is defined, also " + "weatherID has to be defined"
                )
            else:
                cfg["weather"] = kwgs.pop("weatherData")
                
                # check if it is the correct weather data
                for key in weather_units:
                    if not key in cfg["weather"].columns:
                        raise ValueError('Column "' + key + '" required in weatherData')

                cfg["design_T_min"] = cfg["weather"].min()["T"]
                cfg["weatherID"] = kwgs.pop("weatherID")
                if (
                    cfg["longitude"] == KWARG_DEFAULTS["longitude"]
                    or cfg["latitude"] == KWARG_DEFAULTS["latitude"]
                ):
                    warnings.warn(
                        "longitude and latitude are set to "
                        + "default values. It can cause an error in "
                        "the solar irration simulation"
                    )

        else:
            # get TRY weather and ISO
            cfg["weather"], cfg["design_T_min"], cfg[
                "weatherID"
            ] = tsib.getISO12831weather(
                cfg["longitude"],
                cfg["latitude"],
                year=kwgs.pop("year"),
                cosmo=kwgs.pop("getNetCDF4"),
            )

        # save relevant ID entries
        self.IDentries["T_min"] = cfg["design_T_min"]
        self.IDentries["weather"] = cfg["weatherID"]

        # get controller booleans
        for control in ["nightReduction", "occControl", "capControl","ventControl"]:
            cfg[control] = kwgs.pop(control)
            self.IDentries[control] = cfg[control]

        # get comfort zone
        cfg["comfortT_lb"] = kwgs.pop("comfortT_lb")
        cfg["comfortT_ub"] = kwgs.pop("comfortT_ub")
        self.IDentries["ComfortZone"] = (
            str(cfg["comfortT_lb"]) + "-" + str(cfg["comfortT_ub"])
        )
        if cfg["comfortT_lb"] >= cfg["comfortT_ub"] - 0.5:
            warnings.warn(
                "If the gap is too small, the solver can run into"
                + " numerical troubles",
                UserWarning,
            )

        # get occupancy
        cfg["n_persons"] = kwgs.pop("n_persons")

        if "n_apartments" in kwgs:
            cfg["n_apartments"] = kwgs.pop("n_apartments")
        else:
            cfg["n_apartments"] = self.iwu_bdg.loc[
                self.IDentries["Shape"], "n_Apartment"
            ]
            logging.info('number of app. "n_apartments" is inherited from IWU')
        self.IDentries["n_persons"] = cfg["n_persons"]
        self.IDentries["n_apartments"] = cfg["n_apartments"]

        if "elecLoad" in kwgs:
            if not "elecLoadID" in kwgs:
                raise ValueError(
                    "If elecLoad is defined, also " + " elecLoadID has to be defined"
                )
            else:
                cfg["elecLoad"] = kwgs.pop("elecLoad")
                self.IDentries["elecLoad"] = kwgs.pop("elecLoadID")
                cfg["tsorb_device_load"] = False
        else:
            self.IDentries["elecLoad"] = (
                "CREST_" + str(cfg["n_persons"]) + "x" + str(cfg["n_apartments"])
            )
            cfg["tsorb_device_load"] = True

        # define if an fire place is in the building
        if "hasFirePlace" in kwgs:
            cfg["hasFirePlace"] = kwgs.pop("hasFirePlace")
        else:
            # define oven only for building where the appartments are bigger
            # than 100 m^2 -> self chosen value
            if float(cfg["A_ref"]) / cfg["n_apartments"] >= 100.0:
                cfg["hasFirePlace"] = True

                # define size of the oven - atm just a fix value of 10 kW per flat
                cfg["fireplaceSize"] = 10.0 * cfg["n_apartments"]
            else:
                cfg["hasFirePlace"] = False
        self.IDentries["hasFirePlace"] = cfg["hasFirePlace"]

        # if a varying occupancy profile should get integrated
        if "varyoccupancy" in kwgs:
            cfg["varyoccupancy"] = kwgs.pop("varyoccupancy")
            if cfg["varyoccupancy"] < 1:
                raise ValueError('"varyoccupancy" needs to be at least 1.')
        else:
            cfg["varyoccupancy"] = 1
        self.IDentries["varyoccupancy"] = cfg["varyoccupancy"]

        # if the mean profile or the fluctuation pad profile should be taken
        cfg["mean_load"] = kwgs.pop("mean_load")
        self.IDentries["mean_load"] = cfg["mean_load"]

        # create seed for every building
        state_seed = (
            str(int(cfg["n_persons"]))
            + str(int(cfg["longitude"] * 100))[2:]
            + str(int(cfg["A_ref"]))
        )
        if len(state_seed) > 8:
            state_seed = state_seed[:8]
        state_seed = int(state_seed)
        cfg['state_seed'] = state_seed

        return cfg

    def _get_form(self, cfg, kwgs):
        '''
        Derives the form of the building, in terms of the size of the exterior walls etc.
        '''

        ### Either use predefined building types
        if "ID" in kwgs:
            iwu_bdg = self.iwu_bdg.loc[kwgs["ID"],:].to_dict()
            cfg['a_ref'] = iwu_bdg["A_C_Ref"]
            cfg["n_apartments"] = iwu_bdg["n_Apartment"]
            iwu_idx = kwgs["ID"]
        

        elif "buildingYear" and "buildingType" in kwgs:

            if "buildingClassification" in kwgs:
                bdg_class = kwgs.pop("buildingClassification")
            else:
                bdg_class = "Gen"
            if "eastOrOverall" in kwgs:
                bdg_eastwest = kwgs.pop("eastOrOverall")
            else:
                bdg_eastwest = "N"
            cfg['buildingYear'] = kwgs['buildingYear']
            iwu_bdg = self.iwu_bdg[
                (self.iwu_bdg["Year1_Building"] <= cfg["buildingYear"])
                & (cfg["buildingYear"] <= self.iwu_bdg["Year2_Building"])
                & (self.iwu_bdg["Code_BuildingSizeClass"] == kwgs.pop("buildingType"))
            ]
            # ugly but pycharm complains otherwise
            is_class = [
                bdg_class in btype.split(".") for btype in iwu_bdg["Code_BuildingType"]
            ]
            iwu_bdg = iwu_bdg[is_class]
            is_eastwest = [
                bdg_eastwest in btype.split(".")
                for btype in iwu_bdg["Code_BuildingType"]
            ]
            iwu_bdg = iwu_bdg[is_eastwest]
            # test if a building was chosen
            if not iwu_bdg.index.any():
                raise ValueError(
                    "Building with this definition does not exist in IWU-Database"
                )
            if len(iwu_bdg.index) > 1:
                raise ValueError("Internal Error: To many building data to unpack")

            iwu_idx = iwu_bdg.index[0]

            iwu_bdg = iwu_bdg.loc[iwu_idx, :].to_dict()

            cfg['a_ref'] = iwu_bdg["A_C_Ref"]
            cfg["n_apartments"] = iwu_bdg["n_Apartment"]
        
        ### Or derive from the shape
        else:

            # save the number of apps
            cfg["n_apartments"] = kwgs["n_apartments"]

            # get living area
            if ("a_ref_app" in kwgs and "n_apartments" in kwgs):
                if "a_ref" in kwgs:
                    cfg['a_ref'] = kwgs.pop("a_ref")
                else:
                    # calculate full reference area based on appartments
                    cfg['a_ref'] = kwgs.pop("a_ref_app") * cfg["n_apartments"]
            elif "a_ref" in kwgs:
                cfg['a_ref'] = kwgs.pop("a_ref")
            else:
                warnings.warn(
                    'No suffucient sufficient keyword arguments provided for the living area "a_ref". It is set to 150 m²'
                )
                cfg['a_ref'] = 150

            # get the surrounding
            if "surrounding" in kwgs:
                cfg["surrounding"] = kwgs.pop("surrounding")
            else:
                warnings.warn(
                    'No "surrounding" provided. Falling back to "detached".'
                )
                cfg['surrounding'] = "Detached"


            # reduce to buildings with equivalent surrounding
            surDict = {"B_Alone": "Detached", "B_N1": "Semi", "B_N2": "Terraced"}
            self.iwu_bdg.replace({"Code_AttachedNeighbours": surDict}, inplace=True)
            iwu_sur = self.iwu_bdg[
                self.iwu_bdg["Code_AttachedNeighbours"] == cfg['surrounding'] 
            ]
            
            # get the most similar building
            diff_area = abs(iwu_sur["A_C_Ref"] - cfg['a_ref'])

            iwu_idx = diff_area.idxmin()

            iwu_bdg = self.iwu_bdg.xs(iwu_idx).to_dict()

        # adaot the shape values from the chosen iwu bdg
        cfg = get_shape(cfg, iwu_bdg, cfg['a_ref'])

        self.IDentries["Shape"] = iwu_idx
        self.IDentries["A_ref"] = cfg["A_ref"]

        # integrate roof tilt
        if "roofTilt" in kwgs:
            cfg["roofTilt"] = kwgs.pop("roofTilt")
        else:
            if iwu_bdg["Code_RoofType"] in ["FR"]:
                cfg["roofTilt"] = 0.0
            else:
                cfg["roofTilt"] = 45.0
        self.IDentries["roofTilt"] = cfg["roofTilt"]

        # integrate roof orientation
        cfg["roofOrientation"] = kwgs.pop("roofOrientation")
        self.IDentries["roofOrientation"] = cfg["roofOrientation"]

        return cfg

    def _get_fabric(self, cfg, kwgs):
        """
        Get the fabric of the surrounding walls and derives the thermal conductivity
        """
        if "ID" in kwgs:
            self.IDentries["Fabric"] = kwgs["ID"]
            cfg = get_fabric(cfg, self.iwu_bdg.xs(kwgs["ID"]).to_dict())
            cfg["buildingYear"] = self.iwu_bdg.xs(kwgs["ID"]).to_dict()[
                "Year2_Building"
            ]
            # drop ID kwg after usage
            kwgs.pop('ID')
        else:
            # get the buildingyear
            year = None            
            if kwgs.pop("buildnew"):
                year = 2020
            elif "buildingYear" in kwgs:
                year = kwgs.pop("buildingYear")
                if kwgs.pop("refurbished"):
                    year_before = copy.deepcopy(year)
                    year = min(max(year + 40, 1995), 2020)
                    warnings.warn(
                        '"refurbished" just overwrites the buildingyear from '
                        + str(year_before)
                        + " to "
                        + str(year) 
                        + " for the chose type of fabric"
                    )
            else:
                raise ValueError('"buildingYear" is required as argument')

            cfg["buildingYear"] = year

            # get all buildings with this year
            iwu_year = self.iwu_bdg[
                (self.iwu_bdg["Year1_Building"] <= year)
                & (year <= self.iwu_bdg["Year2_Building"])
            ]

            # first try to get all with similar surrounding
            surDict = {"B_Alone": "Detached", "B_N1": "Semi", "B_N2": "Terraced"}
            self.iwu_bdg.replace({"Code_AttachedNeighbours": surDict}, inplace=True)
            iwu_sur = iwu_year[
                iwu_year["Code_AttachedNeighbours"] == cfg['surrounding'] 
            ]
            # just use those if such a building exist
            if len(iwu_sur) > 1:
                iwu_year = iwu_sur


            # get the one with the most similar reference area
            diff_area = abs(
                iwu_year["A_C_Ref"] * iwu_year["n_Apartment"]
                - cfg['a_ref']
            )
            iwu_idx = diff_area.idxmin()
            iwu_bdg = iwu_year.xs(iwu_idx).to_dict()

            cfg = get_fabric(cfg, iwu_bdg)

            # set chjosen buildingType as fabric entry
            self.IDentries["Fabric"] = iwu_idx

        # define thermal class
        cfg["thermalClass"] = kwgs.pop("thermalClass")
        self.IDentries["thermalClass"] = cfg["thermalClass"]

        # if refurbishment is an option for the optimization
        cfg["refurbishment"] = kwgs.pop("refurbishment")
        self.IDentries["refurbishment"] = cfg["refurbishment"]

        # if refurbishment is a forced into the optimization
        cfg["force_refurbishment"] = kwgs.pop("force_refurbishment")
        if cfg["force_refurbishment"] and not cfg["refurbishment"]:
            raise ValueError(
                'If "force_refurbishment" is activated, "refurbishment" must be activated as well.'
            )
        # force_refurbishment to ID entries
        self.IDentries['force_refurbishment'] = cfg['force_refurbishment']

        return cfg

    def _get_equipment(self, cfg, kwgs):
        """
        Gives the configuration for existing energy supply technologies
        in the building.
        """
        # if the hot water is supplied by an electricity boiler or by the
        # heating system
        cfg["hotWaterElec"] = kwgs.pop("hotWaterElec")
        self.IDentries["hotWaterElec"] = cfg["hotWaterElec"]

        # if hot water is generated electrically, correct the hot water demand (BDEW table)
        if not self.ignore_profiles:
            if cfg["hotWaterElec"]:
                cfg["hotWaterLoad"] = cfg["hotWaterLoad"] * 0.6

        # get existing heat supply
        cfg["existingHeatSupply"] = kwgs.pop("existingHeatSupply")
        self.IDentries["existingHeatSupply"] = cfg["existingHeatSupply"]

        # TODO: replace heat supply with heat equipment age
        cfg["replaceHeatSupply"] = kwgs.pop("replaceHeatSupply")
        self.IDentries["replaceHeatSupply"] = cfg["replaceHeatSupply"]

        # define if it has already solar thermal
        if "hasSolarThermal" in kwgs:
            cfg["hasSolarThermal"] = kwgs.pop("hasSolarThermal")
        else:
            # define solar thermal for all post enev 2009 gas boiler buildings
            if cfg["existingHeatSupply"] == "Gas boiler" and (
                cfg["buildingYear"] >= 2009
            ):
                cfg["hasSolarThermal"] = True

            else:
                cfg["hasSolarThermal"] = False
        self.IDentries["hasSolarThermal"] = cfg["hasSolarThermal"]

        # define if it has already photovoltaic
        cfg["hasPhotovoltaic"] = kwgs.pop("hasPhotovoltaic")
        self.IDentries["hasPhotovoltaic"] = cfg["hasPhotovoltaic"]

        # determine the design supply temperature depending on the size of the
        # building and the age
        if "T_sup" in kwgs:
            T_sup = kwgs.pop("T_sup")
        else:
            # TODO: find data for this
            T_sup = 70
            if cfg["n_apartments"] > 6:
                T_sup += 5
            if cfg["buildingYear"] > 1990:
                T_sup -= 10
            if cfg["buildingYear"] > 2000:
                T_sup -= 10
            if cfg["buildingYear"] > 2010:
                T_sup -= 10
            if "floorHeating" in kwgs:
                if kwgs.pop("floorHeating"):
                    T_sup = 40.0
        cfg["T_sup"] = T_sup
        cfg["T_ret"] = T_sup - 20.
        self.IDentries["T_sup"] = cfg["T_sup"]

        return cfg

    def _get_finance(self, cfg, kwgs):
        """
        Get the interest rate and the ownership structure of the building
        """
        cfg["ownership"] = kwgs.pop("ownership")

        cfg["onlyEnergyInvest"] = kwgs.pop("onlyEnergyInvest")

        if "WACC" in kwgs:
            cfg["WACC"] = kwgs.pop("WACC")
        else:
            if cfg["ownership"]:
                cfg["WACC"] = 0.03
            else:
                cfg["WACC"] = 0.06

        self.IDentries["WACC"] = cfg["WACC"]
        self.IDentries["ownership"] = cfg["ownership"]
        return cfg


def get_fabric(bdg, iwu_bdg):
    # merge two window types to one window type
    bdg["U_Window"] = (
        iwu_bdg["A_Window_1"] * iwu_bdg["U_Window_1"]
        + iwu_bdg["A_Window_2"] * iwu_bdg["U_Window_2"]
    ) / bdg["A_Window"]
    bdg["g_gl_n_Window"] = (
        iwu_bdg["A_Window_1"] * iwu_bdg["g_gl_n_Window_1"]
        + iwu_bdg["A_Window_2"] * iwu_bdg["g_gl_n_Window_1"]
    ) / bdg["A_Window"]

    # get u values of the walls etc.
    for wall in ["Wall_1", "Wall_2", "Wall_3"]:
        bdg["U_" + wall] = iwu_bdg["U_Actual_" + wall]
        bdg["b_Transmission_" + wall] = iwu_bdg["b_Transmission_" + wall]
    for roof in ["Roof_1", "Roof_2"]:
        bdg["U_" + roof] = iwu_bdg["U_Actual_" + roof]
        bdg["b_Transmission_" + roof] = iwu_bdg["b_Transmission_" + roof]
    for door in ["Door_1"]:
        bdg["U_" + door] = iwu_bdg["U_Actual_" + door]
    for floor in ["Floor_1", "Floor_2"]:
        bdg["U_" + floor] = iwu_bdg["U_Actual_" + floor]
        bdg["b_Transmission_" + floor] = iwu_bdg["b_Transmission_" + floor]

    # TODO: make them as own argument and move them to operation
    bdg["n_air_infiltration"] = iwu_bdg["n_air_infiltration"]
    bdg["n_air_use"] = iwu_bdg["n_air_use"]

    # include specific space heat demand by episcope
    bdg["q_h_nd"] = iwu_bdg["q_h_nd"]

    return bdg


def get_shape(bdg, iwu_bdg, a_ref):
    """
    Takes the IWU reference building and scales its shape such
     that it fits the given reference area

    :param bdg: Dictionary with the relevant parameters
    :param iwu_bdg: IWU-Reference building as dictionary
    :param a_ref: float
    :return: bdg
    """
    # get new reference area
    bdg["A_ref"] = a_ref

    # get ratio
    ratio = a_ref / iwu_bdg["A_C_Ref"]

    # keep the number of storeys and the height of the room
    bdg["n_Storey"] = iwu_bdg["n_Storey"]
    bdg["h_room"] = iwu_bdg["h_room"]

    # get specific and an aggregated windows area
    for di in ["North", "East", "South", "West"]:
        bdg["A_Window_" + di] = iwu_bdg["A_Window_" + di] * (ratio ** 0.5)
    bdg["A_Window_Horizontal"] = (
        iwu_bdg["A_Window_Horizontal"] * ratio / iwu_bdg["n_Storey"]
    )
    bdg["A_Window"] = sum(
        bdg["A_Window_" + di] for di in ["North", "East", "South", "West", "Horizontal"]
    )

    # get shading relevance
    bdg["F_sh_vert"] = iwu_bdg["F_sh_vert"]
    bdg["F_sh_hor"] = iwu_bdg["F_sh_hor"]
    bdg["F_f"] = iwu_bdg["F_f"]
    bdg["F_w"] = iwu_bdg["F_w"]

    # get new shape of wall, roof and floor
    for wall in ["Wall_1", "Wall_2", "Wall_3"]:
        bdg["A_" + wall] = iwu_bdg["A_" + wall] * (ratio ** 0.5)
    for roof in ["Roof_1", "Roof_2"]:
        bdg["A_" + roof] = iwu_bdg["A_" + roof] * (
            1 + (ratio - 1) / iwu_bdg["n_Storey"]
        )
    for floor in ["Floor_1", "Floor_2"]:
        bdg["A_" + floor] = iwu_bdg["A_" + floor] * (
            1 + (ratio - 1) / iwu_bdg["n_Storey"]
        )
    for door in ["Door_1"]:
        bdg["A_" + door] = iwu_bdg["A_" + door] * (ratio ** 0.5)
    return bdg


