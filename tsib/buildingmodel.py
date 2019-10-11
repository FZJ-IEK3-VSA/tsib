# -*- coding: utf-8 -*-
"""
Created on Sat Dec 10 12:40:17 2016

@author: Leander Kotzur
"""

import os
import warnings

import pandas as pd
import numpy as np

import tsib.household.profiles as household
import tsib.thermal.model5R1C as model5R1C
import tsib.timeseriesmanager as tsm
import tsib.buildingconfig as config
import tsib.data


TOTAL_PROFILE_NUM = 20



class Building(object):

    def __init__(
        self,
        buildingYear=1976,
        buildingType="SFH",
        buildingClassification="Gen",
        weatherData=None,
        weatherID=None,
        eastOrOverall="N",
        n_persons=None,
        n_apartments=None,
        latitude=50.0,
        longitude=8.0,
        roofOrientation=135.0,
        roofTilt=45.0,
        a_ref=None,
        thermalClass=None,
        nightReduction=True,
        occControl=False,
        capControl=True,
        ventControl=False,
        refurbishment=True,
        isRefurbished=False,
        comfortT_lb=21.0,
        comfortT_ub=24.0,
        costdata="default_2016",
        elecLoad=None,
        elecLoadID=None,
        configurator=None,
    ):
        """
        A building model which uses the IWU-Buildingtopology for parameterizing
        the physical building of a model. This can then be used to first run an
        occupancy simulation based on the CREST model (tsorb) and then a thermal
        simulation based on a 5R1C model to predict the heat loads. 
        
        Parameters
        ----------
        buildingYear: int, optional (default: 1976)
            Gets the buildingclasses for these years -> used for IWU
        buildingType: str, optional (default: SFH)
            Size category of the buildings: 'AB', 'SFH', 'SFH2', 'MFH', 'MFH2', 
            'TH' -> used for IWU
        buildingClassification: str, optional (default: Gen)
            If it is a generic or a lightframe house. -> used for IWU
        weatherData: pd.DataFrame, optional (default: None)
            If None, the weather data is derived from the longitude and 
            latitude and the related TRY region.
        weatherID: pd.DataFrame, optional (default: TRY-Zone)
            If weather data is externally set, please define an identifier.
            Otherwise it is automatically set to the TRY-Zone.
        eastOrOverall: str, optional (default: N)
            Generic 'N' or buildings from east Germany 'East' -> used for IWU
        n_persons: int, optional (default: None)
            Number of persons living in a single flat in the house. If not set
            it is randomly chosen.
        n_apartments: int, optional (default: None)
            Number of apartments in the house in the house. If not set
            it is chosen from the IWU-Data.
        latitude: float, optional (default: 50.)
            Location latitude of the building.
        longitude: float, optional (default: 8.)
            Location longitude of the building.
        roofOrientation: float, optional (default: 135.)
            Orientation of the gable in degree (only relevant for tilted roof)
            from North = 0
        roofTilt: float, optional (default: 45.)
            Tilt of the roof in degree. 0 is flat. The parameter gets 
            overwritten for buildings with flat roofs.
        thermalClass: str, optional (default: 'medium')
            The themal class is used for the heat simulation and determines
            the heat capacitiy of the building itself:
                'very light','light', 'medium', 'heavy', 'very heavy'
        nightReduction: boolean, optional (default: True)
            If activated, the lower bound of the temperature tolerance between 22:00 and 6:59
            is set to 18°C.
        occControl: boolean, optional (default: False)
            Occupance controller which increases the temperature tolerance to 10-30°C in case
            nobody is at home.
        capControl: boolean, optional (default: True)
            Controller which uses the thermal capacity of the building.
        ventControl: bool, optional (default: False)
            Allows a control of the ventilation system.
            Attention: The correct equation is nonlinear, why here a linearized
            version is used, which results in a small error.
        refurbishment: bool, optional (default: True)
            Activates the degree of freedom to refurbish to true.
        isRefurbished: bool, optional (default: False)
            If True, an Enev 2016 conform insulation and window type is chosen
            as given. Nevertheless, the design of the heating system is still
            based on the old values.
        comfortT_lb: float, optional (default: 20.)
            Lower bound of the comfort tepm of the occupants.
        comfortT_ub: float, optional (default: 26.)
            Upper bound of the comfort tepm of the occupants.
        costdata: str, optional (default: default_2016)
            File name where data to the refurbishment measures can be found. 
        elecLoad: np.array, optional (default: None)
            Electricity load profile. Normally it is created by the edemand
            model on parallel with occupancy profile. 
        elecLoadID: str, optional (default: None)
            Identifier of the chosen electricityload. Please set in case
            electricityLoad is not None.
        configurator: dict, optional (default: None)
            Configuration dictionary which includes all parameters required
            for the building optimization
        """

        if configurator is None:
            # Building Physical Parameters

            _kwgs = {}
            _kwgs["buildingYear"] = buildingYear
            _kwgs["buildingType"] = buildingType
            _kwgs["buildingClassification"] = buildingClassification
            _kwgs["eastOrOverall"] = eastOrOverall
            _kwgs["weatherData"] = weatherData
            _kwgs["weatherID"] = weatherID
            if n_persons is not None:
                _kwgs["n_persons"] = n_persons
            if n_apartments is not None:
                _kwgs["n_apartments"] = n_apartments
            _kwgs["latitude"] = latitude
            _kwgs["longitude"] = longitude
            _kwgs["roofTilt"] = roofTilt
            _kwgs["roofOrientation"] = roofOrientation
            _kwgs["a_ref"] = a_ref
            _kwgs["thermalClass"] = thermalClass
            _kwgs["occControl"] = occControl
            _kwgs["nightReduction"] = nightReduction
            _kwgs["capControl"] = capControl
            _kwgs["refurbishment"] = refurbishment
            _kwgs["refurbished"] = isRefurbished
            _kwgs["comfortT_lb"] = comfortT_lb
            _kwgs["comfortT_ub"] = comfortT_ub
            _kwgs["costdata"] = costdata
            if elecLoad is not None:
                _kwgs["elecLoad"] = elecLoad
            if elecLoadID is not None:
                _kwgs["elecLoadID"] = elecLoadID

            self.configurator = config.BuildingConfiguration(_kwgs)

        else:
            if not isinstance(configurator, config.BuildingConfiguration):
                raise ValueError(
                    "configurator needs to be of type "
                    + '"config.BuildingConfiguration"'
                )
            else:
                self.configurator = configurator

        self.cfg = self.configurator.getBdgCfg(includeSupply=False)

        # self.sim5R1C
        self.IDentries = self.configurator.IDentries


        self.cfg = self.get_occupancy_profile(
                self.cfg
            )


        self.ventControl = ventControl

        self.times = self.cfg["weather"].index

        self.thermalmodel = model5R1C.Building5R1C(self.cfg)

        # TODO: kick out
        self.ID = 0  #self.configurator.ID 

        self.times = self.cfg["weather"].index
        return

    @property
    def ID_attr(self):
        """'Returns a dictionary of all attributes which make the building
        unique"""

        return self.IDentries

    def _saveResults(self):
        datapath1 = os.path.join(
            tsib.data.PATH, "results", "buildingprofiles", self.ID + ".csv"
        )
        datapath2 = os.path.join(
            tsib.data.PATH, "results", "buildingstaticresults", self.ID + ".csv"
        )
        self.detailedResults.to_csv(datapath1)
        pd.Series(self.static_results).to_csv(datapath2)
        return

    def _loadResults(self):
        datapath1 = os.path.join(
            tsib.data.PATH, "results", "buildingprofiles", self.ID + ".csv"
        )
        datapath2 = os.path.join(
            tsib.data.PATH, "results", "buildingstaticresults", self.ID + ".csv"
        )
        if os.path.isfile(datapath1):
            self.detailedResults = pd.read_csv(datapath1, index_col=0)
            self.detailedResults.index = pd.to_datetime(
                self.detailedResults.index, utc=True
            )
            self.static_results = pd.read_csv(
                datapath2, index_col=0, squeeze=True, header=None
            ).to_dict()
            return True
        else:
            return False


    def get_occupancy_profile(self, cfg):
        """
        Get all the occupancy related data, like internal heat gain,
        electricity load and occupants activity.
        """

        # get a number of random seeds to generate the profiles
        seeds = np.random.RandomState(cfg['state_seed']).randint(
            0, TOTAL_PROFILE_NUM, size=int(cfg["varyoccupancy"] * cfg["n_apartments"])
        )

        # get the profiles
        hh_profiles = household.get_household_profiles(
            cfg["n_persons"],
            cfg["weather"],
            self.IDentries["weather"],
            seeds=seeds,
            ignore_weather=True,
            mean_load=cfg["mean_load"],
        )

        # get short form apartments
        n_app = int(cfg["n_apartments"])
        # abs number of occupants
        n_occs = int(cfg["n_persons"]) * n_app

        # get from the household profiles buildingprofiles
        bdg_profiles = {}

        # get a profile for each occupancy variation
        for i in range(cfg["varyoccupancy"]):
            if n_app == 1:
                occData = hh_profiles[i]
            else:
                # merge profiles for multiple appartments
                occData = pd.concat(
                    hh_profiles[i * n_app : (i + 1) * n_app],
                    keys=range(i * n_app, (i + 1) * n_app),
                )
                occData = occData.groupby(level=[1]).sum()

            bdg_profiles[i] = {}

            # heat gain values relative (Master thesis cheng feng)
            OccNotActiveHeatGain = 100
            OccActiveHeatGain = 150

            # internal heat gain [kW]
            bdg_profiles[i]["Q_ig"] = (
                occData["AppHeatGain"].values
                + occData["OccActive"].values * OccActiveHeatGain
                + occData["OccNotActive"].values * OccNotActiveHeatGain
            ) / 1000
            #        bdg_profiles[i]['occ_active'] = occData['OccActive'] > 0.0

            # the share of occupants which is not at home
            bdg_profiles[i]["occ_nothome"] = (
                n_occs - occData["OccActive"] - occData["OccNotActive"]
            ) / n_occs

            # the share of occupants which is sleeping
            bdg_profiles[i]["occ_sleeping"] = occData["OccNotActive"].div(n_occs)
            if cfg['tsorb_device_load']:
                bdg_profiles[i]["elecLoad"] = occData["Load"] / 1000
            else:
                bdg_profiles[i]["elecLoad"] = cfg["elecLoad"]

            # get hot water load
            bdg_profiles[i]["hotWaterLoad"] = occData["HotWater"] / 1000

            # overall number of occupants
            bdg_profiles[i]["n_occupants"] = n_occs

            # get fireplace profile
            if cfg["hasFirePlace"]:
                pot_filename = os.path.join(
                    tsib.data.PATH,
                    "results",
                    "fireplaceprofiles",
                    "Profile"
                    + "_apps"
                    + str(int(n_app))
                    + "_occ"
                    + str(int(cfg["n_persons"]))
                    + "_wea"
                    + str(cfg["weatherID"])
                    + "_seed"
                    + str(cfg['state_seed'])
                    + ".csv",
                )
                if os.path.isfile(pot_filename):
                    fireplaceLoad = pd.read_csv(
                        pot_filename,
                        index_col=0,
                        header=None,
                        parse_dates=True,
                        squeeze=True,
                    )
                else:
                    # get oven profile depending on activity and outside temperature
                    fireplaceLoad = tsm.createWoodFireProfile(
                        cfg["weather"]["T"],
                        occData["OccActive"] / n_occs,
                        n_ovens=n_app,
                        T_oven_on=5,
                        t_cool=5.0,
                        fullloadSteps=450,
                        seed=cfg['state_seed'],
                    )
                    fireplaceLoad.to_csv(pot_filename)

                bdg_profiles[i]["fireplaceLoad"] = fireplaceLoad

        # give building config the first profile
        cfg.update(bdg_profiles[i])

        # give the other profiles to the configuration file as well
        if cfg["varyoccupancy"] > 1:
            cfg["vary_profiles"] = bdg_profiles

        return cfg


    def getHeatingSystem(self, T_sup=None):
        """
        Determines the heat transfer coefficient between the heating system
        (e.q. the radiator) and the room based on following assumption:
            The system is designed such that it is able to supply a room
            temperature of 20°C at the outside design temperature with a 
            defined supply temperature. This supply temperature is an 
            assumption for each building.
        """
        # get design heat load
        self.design_Q = self.thermalmodel.calcDesignHeatLoad()

        self.design_T_return = 25
        # derive the heat transfer coefficient of the heating system kW/K

    #        self.design_H_heat = self.design_Q / (self.cfg['T_sup'] - self.design_T_return)

    def sim5R1C(self):
        '''
        Runs the 5R1C simulation.
        '''
        self.thermalmodel.sim5R1C()
        return


    def getLoad(self):
        """
        It manages existing load profiles, and if it there is one, it just
        uses this. Otherwise it simulates a new one.
        Get the relevant heating, cooling and electricity profiles of the
        building.
        
        Returns
        -------
        None, but in Building.results or Building.detailedResults can the
        results be shown as pandas.DataFrame
        """

        resultsExist = self._loadResults()
        if resultsExist:
            return
        else:
            if self.cfg["refurbishment"]:
                befRef = True
                warnings.warn(
                    "For the simulation the refurbishment decisions"
                    + " are deactivated",
                    UserWarning,
                )
                self.cfg["refurbishment"] = False
            else:
                befRef = False

            self.sim5R1C()
            self._saveResults()
            self.cfg["refurbishment"] = befRef
        return

    def toCSV(self, filename="result"):
        """
        Writes the results to 2 csv files:
            filename + "_dynamic"
            filename + "_static"
        
        Parameters
        ----------
        filename: str, optional, default: "result"
            Filepath and name.
        """
        if filename[-4:] is ".csv":
            filename = filename[:-4]
        self.detailedResults.to_csv(filename + "_dynamic.csv", sep=",")
        pd.Series(self.static_results).to_csv(filename + "_static.csv", sep=",")
        return

    def readCSV(self, filename="result"):
        """
        Reads the results from 2 csv files:
            filename + "_dynamic"
            filename + "_static" 
        
        Parameters
        ----------
        filename: str, optional, default: "result"
            Filepath and name.
        """
        if filename[-4:] is ".csv":
            filename = filename[:-4]
        self.detailedResults = pd.read_csv(
            filename + "_dynamic.csv", sep=",", index_col=0, parse_dates=True
        )
        self.results = pd.read_csv(filename + "_static.csv", sep=",").to_dict()
        return

    @property
    def static_results(self):
        return self.thermalmodel.static_results 

    @property
    def detailedResults(self):
        return self.thermalmodel.detailedResults 

    @property
    def detailedRefurbish(self):
        return self.thermalmodel.detailedRefurbish 




if __name__ == "__main__":

    import matplotlib.pyplot as plt

    kwgs = {
        "buildingYear": 1990,
        "latitude": 52.0,
        "longitude": 13.0,
        "comfortT_lb": 21,
        "comfortT_ub": 24,
        "WACC": 0.03,
        "roofTilt": 45.0,
        "surrounding": "Semi",
        "n_apartments": 2,
        "a_ref_app": 100,
        "n_persons": 2,
        "roofOrientation": 135.0,
        "costdata": "default_2016",
        "capControl": True,
    }

    bdg = config.BuildingConfiguration(kwgs)
    #    cfg = bdg.getBdgCfg()
    example = Building(configurator=bdg, refurbishment=False)

    example.sim5R1C()
    example.detailedResults["5R1C Cooling Load T=21-24"] = -example.detailedResults[
        "Cooling Load"
    ]
    example.detailedResults["5R1C Heating Load T=21-24"] = example.detailedResults[
        "Heating Load"
    ]

    example.detailedResults[
        ["5R1C Heating Load T=21-24", "5R1C Cooling Load T=21-24"]
    ].plot()
    plt.show()
    example.detailedResults["T_out"] = example.cfg["weather"]["T"]
    example.detailedResults[["T_air", "T_s", "T_m", "T_out"]].plot()
    plt.show()
    example.thermalmodel.M.exVars.display()

    example.detailedResults[["5R1C Heating Load T=21-24", "T_air"]].plot(subplots=True)
    plt.show()
