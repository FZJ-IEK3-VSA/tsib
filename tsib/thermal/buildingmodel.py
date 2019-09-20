# -*- coding: utf-8 -*-
"""
Created on Sat Dec 10 12:40:17 2016

@author: Leander
"""

import pandas as pd
import os
import numpy as np
import pvlib
import pvlib.irradiance as irr
import time
import pyomo.environ as pyomo
import pyomo.opt as opt
import tsib.thermal.utils as utils
import tsib.buildingconfig as config
import tsam.timeseriesaggregation as tsam
import warnings


DATA_PATH = os.path.join(os.path.dirname(os.path.realpath(__file__)), "..", "data")


class Building(object):
    # list of all attributes which are used to identify the building
    ID_ATTR = [
        "codeBuildingType",
        "weatherID",
        "design_T_min",
        "refurbishment",
        "n_persons",
        "thermalClass",
        "roofTilt",
        "roofOrientation",
        "comfortT_lb",
        "comfortT_ub",
        "ventControl",
        "nightReduction",
        "capControl",
        "costData",
        "elecLoadID",
    ]

    CONST = {
        # Constants for calculation of A_m, dependent of building class
        # (DIN EN ISO 13790, section 12.3.1.2, page 81, table 12)
        "f_Am": [2.5, 2.5, 2.5, 3.0, 3.5],
        # specific heat transfer coefficient between internal air and surface [kW/m^2/K]
        # (DIN EN ISO 13790, section 7.2.2.2, page 35)
        "h_is": 3.45 / 1000,
        # non-dimensional relation between the area of all indoor surfaces
        # and the effective floor area A["f"]
        # (DIN EN ISO 13790, section 7.2.2.2, page 36)
        "lambda_at": 4.5,
        # specific heat transfer coefficient thermal capacity [kW/m^2/K]
        # (DIN EN ISO 13790, section 12.2.2, page 79)
        "h_ms": 9.1 / 1000,
        # ISO 6946 Table 1, Heat transfer resistances for opaque components
        "R_se": 0.04 * 1000,  # external heat transfer coefficient m²K/W
        # ASHRAE 140 : 2011, Table 5.3, page 18 (infrared emittance) (unused --> look at h_r)
        "epsilon": 0.9,
        # external specific radiative heat transfer [kW/m^2/K] (ISO 13790, Schuetz et al. 2017, 2.3.4)
        "h_r": 0.9 * 5.0 / 1000.0,
        # ASHRAE 140 : 2011, Table 5.3, page 18 (absorption opaque comps)
        "alpha": 0.6,
        # average difference external air temperature and sky temperature
        "delta_T_sky": 11.0,  # K
        # density air
        "rho_air": 1.2,  # kg/m^3
        # heat capacity air
        "C_air": 1.006,  # kJ/kg/K
    }

    def __init__(
        self,
        buildingYear=1976,
        buildingType="SFH",
        buildingClassification="Gen",
        weatherData=None,
        weatherID=None,
        eastOrOverall="N",
        annualElecDemand=None,
        n_persons=None,
        n_apartments=None,
        latitude=50.0,
        longitude=8.0,
        roofOrientation=135.0,
        roofTilt=45.0,
        a_ref=None,
        thermalClass=None,
        maxLoad=None,
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
        the physical building of a model. While doing so it runs also an electricity
        load simulation to Inherited from energysystem.Vertex        
        
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
            it is chosen from the IWU-Data.
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
        maxLoad: float, optional (default: DesignHeatLoad)
            Maximal load of the heating system.
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

        self.maxLoad = maxLoad

        if configurator is None:
            # Building Physical Parameters

            kwgs = {}
            kwgs["buildingYear"] = buildingYear
            kwgs["buildingType"] = buildingType
            kwgs["buildingClassification"] = buildingClassification
            kwgs["eastOrOverall"] = eastOrOverall
            kwgs["weatherData"] = weatherData
            kwgs["weatherID"] = weatherID
            if n_persons is not None:
                kwgs["n_persons"] = n_persons
            if n_apartments is not None:
                kwgs["n_apartments"] = n_apartments
            kwgs["latitude"] = latitude
            kwgs["longitude"] = longitude
            kwgs["roofTilt"] = roofTilt
            kwgs["roofOrientation"] = roofOrientation
            kwgs["a_ref"] = a_ref
            kwgs["thermalClass"] = thermalClass
            kwgs["occControl"] = occControl
            kwgs["nightReduction"] = nightReduction
            kwgs["capControl"] = capControl
            kwgs["refurbishment"] = refurbishment
            kwgs["refurbished"] = isRefurbished
            kwgs["comfortT_lb"] = comfortT_lb
            kwgs["comfortT_ub"] = comfortT_ub
            kwgs["costdata"] = costdata
            if elecLoad is not None:
                kwgs["elecLoad"] = elecLoad
            if elecLoadID is not None:
                kwgs["elecLoadID"] = elecLoadID

            self.configurator = config.BuildingConfiguration(kwgs)

        else:
            if not isinstance(configurator, config.BuildingConfiguration):
                raise ValueError(
                    "configurator needs to be of type "
                    + '"config.BuildingConfiguration"'
                )
            else:
                self.configurator = configurator

        self.cfg = self.configurator.getBdgCfg(includeSupply=False)

        self.IDentries = self.configurator.IDentries

        self.costdatapath = os.path.join(
            DATA_PATH, "costdata", self.cfg["costdata"] + ".xlsx"
        )
        self.ventControl = ventControl

        self.times = self.cfg["weather"].index

        # initialize dataframe for irradiance on all surface directions
        self.irrad_surf = pd.DataFrame(index=self.times)

        # initialize result dictionary and dataframes
        self.static_results = {}
        self.detailedResults = pd.DataFrame(index=self.times)
        self.detailedRefurbish = pd.DataFrame()

        # TODO: kick out
        self.ID = (
            self.IDentries["Fabric"].replace(".", "_")
            + "_apps"
            + str(self.IDentries["n_apartments"])
            + "_occ"
            + str(self.IDentries["n_persons"])
            + "_lat"
            + str(self.IDentries["weather"])
            + "_til"
            + str(round(self.IDentries["roofTilt"] * 1e3) / 1e3).replace(".", "-")
            + "_ori"
            + str(round(self.IDentries["roofOrientation"] * 1e3) / 1e3).replace(
                ".", "-"
            )
            + "_area"
            + str(round(self.IDentries["A_ref"] * 1e3) / 1e3).replace(".", "-")
        )
        self.times = self.cfg["weather"].index
        return

    @property
    def ID_attr(self):
        """'Returns a dictionary of all attributes which make the optimization
        unique"""

        return self.IDentries

    def _saveResults(self):
        datapath1 = os.path.join(
            DATA_PATH, "results", "buildingprofiles", self.ID + ".csv"
        )
        datapath2 = os.path.join(
            DATA_PATH, "results", "buildingstaticresults", self.ID + ".csv"
        )
        self.detailedResults.to_csv(datapath1)
        pd.Series(self.static_results).to_csv(datapath2)
        return

    def _loadResults(self):
        datapath1 = os.path.join(
            DATA_PATH, "results", "buildingprofiles", self.ID + ".csv"
        )
        datapath2 = os.path.join(
            DATA_PATH, "results", "buildingstaticresults", self.ID + ".csv"
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

    def _calcUVal(self, listOfU):
        if 0.0 in listOfU:
            return 0.0
        else:
            return 1.0 / sum(1.0 / uVal for uVal in listOfU)

    def _initEnvelop(self, M):
        """
        Initialize the the parameter set which is required for the 
        refurbishment decisions and adds them the pyomo.ConcreteModel instance.
        
        Parameters
        ----------
        M: pyomo.ConcreteModel, required
        """

        # decisions
        M.bComponents = ["Walls", "Roof", "Floor", "Windows", "Ventilation"]

        # for each component, which refurbishment is chosen
        M.bInsul = {}
        # specific heat transfer coefficients for each component for each refurbishment decision
        M.bH = {}
        # specific u values for each component for each refurbishment decision
        M.bU = {}

        # raw insulation/refurbishment data
        refRaw = {}

        # read in potential refurbishment measures
        for comp in M.bComponents:
            refRaw[comp] = pd.read_excel(
                self.costdatapath, sheet_name=comp, skiprows=[1], index_col=0
            )
            refRaw[comp] = refRaw[comp].dropna(how="all")

            # derive u values of each layer
            if comp in ["Walls", "Roof", "Floor"]:
                refRaw[comp]["U_Value"] = (
                    refRaw[comp]["Lambda"] / refRaw[comp]["Thickness"]
                )
            # refurbishment measure or not
            if M.bRefurbishment:
                if self.cfg["force_refurbishment"]:
                    M.bInsul[comp] = refRaw[comp].index.unique()[1:]
                else:
                    M.bInsul[comp] = refRaw[comp].index.unique()
                for dec in M.bInsul[comp]:
                    M.exVarIx.append((comp, dec))
                    M.insulIx.append((comp, dec))
            else:
                ix = refRaw[comp].index.unique()[0]
                M.bInsul[comp] = [ix]
                M.exVarIx.append((comp, ix))
                M.insulIx.append((comp, ix))

            # heat transfer coefficient wall [kW/K] specific for refurbishment decision
            M.bH[comp] = {}

        # TODO cost construction
        # TODO unheated ceiling
        # TODO thermal bridging
        # TODO check for only exterior walls
        # WALL
        # iterate all refurbishment options
        for var in M.bInsul["Walls"]:
            # init heat transfer and capital expenditure
            M.bH["Walls"][var] = 0.0
            M.exVarCAPEX[("Walls", var)] = 0.0
            for wall in ["Wall_1", "Wall_2", "Wall_3"]:
                # get u values
                U_vals_insul = refRaw["Walls"].loc[var, "U_Value"]
                # add existing layer to all refurbishment measures (only for wall)
                if isinstance(U_vals_insul, (np.float64)):
                    U_vals_insul = [U_vals_insul, self.cfg["U_" + wall]]
                else:
                    U_vals_insul = np.append(U_vals_insul.values, self.cfg["U_" + wall])
                # add heat resistance for each wall [kW/K] ()
                M.bH["Walls"][var] += (
                    self._calcUVal(U_vals_insul)
                    * self.cfg["A_" + wall]
                    * self.cfg["b_Transmission_" + wall]
                ) / 1000

                # investment for all walls
                if self.cfg["onlyEnergyInvest"]:
                    M.exVarCAPEX[("Walls", var)] += self.cfg["A_" + wall] * float(
                        refRaw["Walls"].loc[var, "Investment only energy"].sum()
                    )
                else:
                    M.exVarCAPEX[("Walls", var)] += self.cfg["A_" + wall] * float(
                        refRaw["Walls"].loc[var, "Investment"].sum()
                    )

            # no operational expenditure for walls
            M.exVarOPEX[("Walls", var)] = 0.0
            # life time equivalent to default lifetime
            M.exVarLifetime[("Walls", var)] = self.lifetime

        # ROOF
        # iterate all refurbishment options
        for var in M.bInsul["Roof"]:
            # init heat transfer and capital expenditure
            M.bH["Roof"][var] = 0
            M.exVarCAPEX[("Roof", var)] = 0
            for roof in ["Roof_1", "Roof_2"]:
                # get U values
                U_vals_insul = refRaw["Roof"].loc[var, "U_Value"]
                # add existing layer to all refurbishment measures (only for roof)
                if isinstance(U_vals_insul, (np.float64)):
                    U_vals_insul = [U_vals_insul, self.cfg["U_" + roof]]
                else:
                    U_vals_insul = np.append(U_vals_insul.values, self.cfg["U_" + roof])

                # add heat resistance for each roof [kW/K]
                M.bH["Roof"][var] += (
                    self._calcUVal(U_vals_insul)
                    * self.cfg["A_" + roof]
                    * self.cfg["b_Transmission_" + roof]
                    / 1000
                )

                # investment for all roofs
                if self.cfg["onlyEnergyInvest"]:
                    M.exVarCAPEX[("Roof", var)] += self.cfg["A_" + roof] * float(
                        refRaw["Roof"].loc[var, "Investment only energy"].sum()
                    )
                else:
                    M.exVarCAPEX[("Roof", var)] += self.cfg["A_" + roof] * float(
                        refRaw["Roof"].loc[var, "Investment"].sum()
                    )

            M.exVarOPEX[("Roof", var)] = 0
            M.exVarLifetime[("Roof", var)] = self.lifetime

        # FLOOR
        # iterate all refurbishment options
        for var in M.bInsul["Floor"]:
            # init heat transfer and capital expenditure
            M.bH["Floor"][var] = 0
            M.exVarCAPEX[("Floor", var)] = 0
            for floor in ["Floor_1", "Floor_2"]:
                # get U values
                U_vals_insul = refRaw["Floor"].loc[var, "U_Value"]
                # add existing layer to all refurbishment measures (only for floor)
                if isinstance(U_vals_insul, (np.float64)):
                    U_vals_insul = [U_vals_insul, self.cfg["U_" + floor]]
                else:
                    U_vals_insul = np.append(
                        U_vals_insul.values, self.cfg["U_" + floor]
                    )

                # add heat resistance for each floor [kW/K]
                M.bH["Floor"][var] += (
                    self._calcUVal(U_vals_insul)
                    * self.cfg["A_" + floor]
                    * self.cfg["b_Transmission_" + floor]
                    / 1000
                )

                # investment for all Floor
                M.exVarCAPEX[("Floor", var)] += self.cfg["A_" + floor] * float(
                    refRaw["Floor"].loc[var, "Investment"].sum()
                )

            M.exVarOPEX[("Floor", var)] = 0.0
            M.exVarLifetime[("Floor", var)] = self.lifetime

        # WINDOWS
        # add original windows
        refRaw["Windows"].loc["Nothing", "g_gl"] = self.cfg["g_gl_n_Window"]
        refRaw["Windows"].loc["Nothing", "U_Value"] = self.cfg["U_Window"]
        refRaw["Windows"].loc["Nothing", "Investment"] = 0

        M.bg_gl = {}
        M.bU["Windows"] = {}
        for var in M.bInsul["Windows"]:
            M.bH["Windows"][var] = (
                (self.cfg["A_Window"])
                * float(refRaw["Windows"].loc[var, "U_Value"])
                / 1000
            )
            M.bU["Windows"][var] = float(refRaw["Windows"].loc[var, "U_Value"]) / 1000
            M.bg_gl[var] = float(refRaw["Windows"].loc[var, "g_gl"])
            M.exVarCAPEX[("Windows", var)] = self.cfg["A_Window"] * float(
                refRaw["Windows"].loc[var, "Investment"]
            )
            M.exVarOPEX[("Windows", var)] = 0.0
            M.exVarLifetime[("Windows", var)] = self.lifetime

        # VENTILATION
        # heat capacity of the whole air
        C_air = (
            self.cfg["A_ref"]
            * self.cfg["h_room"]
            * self.CONST["rho_air"]
            * self.CONST["C_air"]
        )
        # loop over verntilation options
        for var in M.bInsul["Ventilation"]:
            # ventilation heat flow corrected by the recovery rate for usable air
            M.bH["Ventilation"][var] = (
                C_air
                * (
                    self.cfg["n_air_use"]
                    * (1 - float(refRaw["Ventilation"].loc[var, "Recovery rate"]))
                    + self.cfg["n_air_infiltration"]
                )
                / 3600
            )  # [kW/K]
            M.exVarCAPEX[("Ventilation", var)] = self.cfg["A_ref"] * float(
                refRaw["Ventilation"].loc[var, "Investment"]
            )
            M.exVarOPEX[("Ventilation", var)] = float(
                refRaw["Ventilation"].loc[var, "OPEX-Fix"]
            )
            M.exVarLifetime[("Ventilation", var)] = float(
                refRaw["Ventilation"].loc[var, "Lifetime"]
            )

        return M

    def _initControl(self, M):
        """
        Initialize the the parameter set which is required for the 
        controller decisions and adds them the pyomo.ConcreteModel instance.
        
        Parameters
        ----------
        M: pyomo.ConcreteModel, required
        """
        # read in and add controller
        costRaw = pd.read_excel(
            self.costdatapath, sheet_name="Control", skiprows=[1], index_col=0
        )
        costRaw = costRaw.dropna(how="all")

        # init all cost and vars
        for var in costRaw.index:
            M.exVarCAPEX[("Control", var)] = self.cfg["A_ref"] * float(
                costRaw.loc[var, "Investment spec"]
            ) + float(costRaw.loc[var, "Investment fix"])
            M.exVarLifetime[("Control", var)] = float(costRaw.loc[var, "Lifetime"])
            M.exVarOPEX[("Control", var)] = float(costRaw.loc[var, "OPEX-Fix"])
            M.exVarIx.append(("Control", var))
        # adapt controller depending on previous condition
        if self.cfg["occControl"]:
            M.exVarActive.append(("Control", "Occupancy"))
            M.exVarCAPEX[("Control", "Occupancy")] = 0
        else:
            if not self.cfg["refurbishment"]:
                M.exVarInActive.append(("Control", "Occupancy"))
        if self.cfg["nightReduction"]:
            M.exVarActive.append(("Control", "NightReduction"))
            M.exVarCAPEX[("Control", "NightReduction")] = 0
        else:
            if not self.cfg["refurbishment"]:
                M.exVarInActive.append(("Control", "NightReduction"))
        if self.cfg["capControl"]:
            M.exVarActive.append(("Control", "SmartThermostat"))
            M.exVarCAPEX[("Control", "SmartThermostat")] = 0
        else:
            if not self.cfg["refurbishment"]:
                M.exVarInActive.append(("Control", "SmartThermostat"))
        return M

    def _init5R1C(self, M):
        """
        Initialize all required parameters required for the 5R1C model defined 
        by the DIN EN ISO 13790 and adds them to a pyomo.ConcreteModel
        """
        # Constants
        M.bConst = self.CONST

        # Thermal capacity class definition of the building --> also algebraic
        # variant possible
        bClass_f_lb = {
            "very light": 0.0,
            "light": 95.0,
            "medium": 137.5,
            "heavy": 212.5,
            "very heavy": 313.5,
        }
        bClass_f_ub = {
            "very light": 95.0,
            "light": 137.5,
            "medium": 212.5,
            "heavy": 313.5,
            "very heavy": 313.5 * 2,
        }
        bClass_f_a = {
            "very light": 2.5,
            "light": 2.5,
            "medium": 2.5,
            "heavy": 3.0,
            "very heavy": 3.5,
        }

        # heated floor area
        M.bA_f = self.cfg["A_ref"]

        # calculate effective heat transfer of heat capacity
        M.bA_m = M.bA_f * bClass_f_a[self.cfg["thermalClass"]]
        # effective transfer coefficient [kW/K]
        M.bH_ms = M.bA_m * M.bConst["h_ms"]

        # specific heat [kJ/m^2/K] (Note: The standard value is overwritten!)
        self.cfg["c_m"] = (
            bClass_f_lb[self.cfg["thermalClass"]]
            + bClass_f_ub[self.cfg["thermalClass"]]
        ) / 2.0
        # internal heat capacity [kWh/K]
        M.bC_m = self.cfg["A_ref"] * self.cfg["c_m"] / 3600.0

        # through door [kW/K]
        M.bH_door = (self.cfg["A_Door_1"] * self.cfg["U_Door_1"]) / 1000

        # internal surface area
        M.bA_tot = self.cfg["A_ref"] * M.bConst["lambda_at"]
        # heat transfer between surface and air node [kW/k] (Schuetz et al. 2017 - eq. 11)
        M.bH_is = M.bA_tot * M.bConst["h_is"]

        # comfort temperature
        M.bT_comf_lb = self.cfg["comfortT_lb"]
        M.bT_comf_ub = self.cfg["comfortT_ub"]

        # shorten code
        rt = self.cfg["roofTilt"]
        ro = self.cfg["roofOrientation"]

        # set relevant geometrical data
        surf_az = {
            "North": 0.0 + ro,
            "East": 90.0 + ro,
            "South": 180.0 + ro,
            "West": 270.0 + ro,
            "Horizontal": 180.0,
            "Roof 1": 90.0 + ro,
            "Roof 2": 270.0 + ro,
        }
        surf_tilt = {
            "North": 90.0,
            "East": 90.0,
            "South": 90.0,
            "West": 90.0,
            "Horizontal": 0.0,
            "Roof 1": rt,
            "Roof 2": rt,
        }
        F_r = {"North": 0.5, "East": 0.5, "South": 0.5, "West": 0.5, "Horizontal": 1.0}

        self._calcRadiation(surf_az, surf_tilt)

        # solar radiation dict [kW]
        M.bQ_sol = {}

        # WINDOWS
        # get relevant window area for solar irradiance
        window_area_s = {
            di: self.cfg["A_Window_" + di] * self.cfg["F_sh_vert"]
            for di in ["North", "East", "South", "West"]
        }
        window_area_s["Horizontal"] = (
            self.cfg["A_Window_Horizontal"] * self.cfg["F_sh_hor"]
        )

        # get relevant window area for thermal irradiance
        window_area_t = sum(self.cfg["A_Window_" + key] * F_r[key] for key in F_r)

        # get solar radiation on window area
        irrad_on_windows = (
            self.irrad_surf.mul(pd.Series(window_area_s)).sum(axis=1).values
        )

        # solar gain time series depending on the investment decision
        for var in M.bInsul["Windows"]:
            # thermal radiation [kW] (Schuetz et al. 2017 - eq. 13)
            thermal_rad = (
                window_area_t
                * M.bConst["h_r"]
                * M.bU["Windows"][var]
                * M.bConst["R_se"]
                * M.bConst["delta_T_sky"]
            )
            M.profiles["bQ_sol_Windows" + var] = (
                irrad_on_windows
                * (1.0 - self.cfg["F_f"])
                * self.cfg["F_w"]
                * M.bg_gl[var]
                - thermal_rad
            )

        # WALLS
        mean_ver_irr = (
            self.irrad_surf.loc[:, ["North", "East", "South", "West"]]
            .mean(axis=1)
            .values
        )
        for var in M.bInsul["Walls"]:
            # thermal radiation [kW] (Schuetz et al. 2017 - eq. 13)
            thermal_rad = (
                M.bH["Walls"][var]
                * M.bConst["h_r"]
                * M.bConst["R_se"]
                * M.bConst["delta_T_sky"]
            )
            solar_radiation = (
                mean_ver_irr
                * M.bH["Walls"][var]
                * self.cfg["F_sh_vert"]
                * M.bConst["R_se"]
                * M.bConst["alpha"]
            )
            M.profiles["bQ_sol_Walls" + var] = solar_radiation - thermal_rad

        # ROOF
        mean_roof_irr = self.irrad_surf.loc[:, ["Roof 1", "Roof 2"]].mean(axis=1).values
        for var in M.bInsul["Roof"]:
            # thermal radiation (Schuetz et al. 2017 - eq. 13)
            thermal_rad = (
                M.bH["Roof"][var]
                * M.bConst["h_r"]
                * M.bConst["R_se"]
                * M.bConst["delta_T_sky"]
            )
            solar_radiation = (
                mean_roof_irr
                * M.bH["Roof"][var]
                * self.cfg["F_sh_hor"]
                * M.bConst["R_se"]
                * M.bConst["alpha"]
            )
            M.profiles["bQ_sol_Roof" + var] = solar_radiation - thermal_rad
        return M

    def _calcRadiation(self, surf_az, surf_tilt):
        """
        Calculates the radiation to all considered walls and roofs.
        """
        # init required time series
        SOL_POS = pvlib.solarposition.get_solarposition(
            self.cfg["weather"].index, self.cfg["latitude"], self.cfg["longitude"]
        )
        AM = pvlib.atmosphere.get_relative_airmass(SOL_POS["apparent_zenith"])
        DNI_ET = irr.get_extra_radiation(self.cfg["weather"].index.dayofyear)

        for key in surf_az:
            # calculate total irradiance depending on surface tilt and azimuth
            total = irr.get_total_irradiance(
                surf_tilt[key],
                surf_az[key],
                SOL_POS["apparent_zenith"],
                SOL_POS["azimuth"],
                dni=self.cfg["weather"]["DNI"],
                ghi=self.cfg["weather"]["GHI"],
                dhi=self.cfg["weather"]["DHI"],
                dni_extra=DNI_ET,
                airmass=AM,
                model="perez",
                surface_type="urban",
            )
            # get plane of array (POA) irradiance and replace nan
            self.irrad_surf[key] = total["poa_global"].fillna(0)

        # W to kW
        self.irrad_surf = self.irrad_surf / 1000
        return

    @staticmethod
    def _initOpti(M):
        """
        Initilizes all required list and dicts for indices and params of 
        building model for optmization
        """
        # get profiles dict for all considered time series
        if not hasattr(M, "profiles"):
            M.profiles = utils.PowerDict()
        if not hasattr(M, "profilesEval"):
            M.profilesEval = utils.PowerDict()
        if not hasattr(M, "exVarIx"):
            M.exVarIx = []
        if not hasattr(M, "exVarActive"):
            M.exVarActive = []
        if not hasattr(M, "exVarInActive"):
            M.exVarInActive = []
        if not hasattr(M, "insulIx"):
            M.insulIx = []
        if not hasattr(M, "exVarCost"):
            M.exVarCost = utils.PowerDict()
        if not hasattr(M, "exVarCAPEX"):
            M.exVarCAPEX = utils.PowerDict()
        if not hasattr(M, "exVarOPEX"):
            M.exVarOPEX = utils.PowerDict()
        if not hasattr(M, "exVarLifetime"):
            M.exVarLifetime = utils.PowerDict()
        if not hasattr(M, "bConnectedHeat"):
            M.bConnectedHeat = []
        if not hasattr(M, "bConnectedElec"):
            M.bConnectedElec = []
        if not hasattr(M, "bConnectedCool"):
            M.bConnectedCool = []
        return M

    def _addOpti(self, M):
        """
        Add indicies and parameters of the building to the Optimization model
        M.
        """
        M = self._initOpti(M)

        if hasattr(M, "bMaxLoad"):
            raise NotImplementedError(
                "At the moment only single building"
                + " initialization for enercore network"
                + " possible."
            )

        # limit maximal heat or cooling load
        if self.maxLoad is None:
            designLoad = self.calcDesignHeatLoad()
            M.bMaxLoad = designLoad
            self.maxLoad = designLoad
        else:
            M.bMaxLoad = self.maxLoad

        # define refurbishment
        M.bRefurbishment = self.cfg["refurbishment"]

        # get all envelop values
        M = self._initEnvelop(M)

        # get control decisions
        M = self._initControl(M)

        # get all radiation values
        M = self._init5R1C(M)

        # get occupancy data and internal heat gains and electricity load
        # evaluate all profiles with 0 in case of time series aggregation
        M.profiles["bQ_ig"] = self.cfg["Q_ig"]
        M.profilesEval["bQ_ig"] = tsam.MIN_WEIGHT
        M.profiles["bOccNotHome"] = self.cfg["occ_nothome"].values
        M.profilesEval["bOccNotHome"] = tsam.MIN_WEIGHT
        M.profiles["bOccSleeping"] = self.cfg["occ_sleeping"].values
        M.profilesEval["bOccSleeping"] = tsam.MIN_WEIGHT
        M.profiles["bElecLoad"] = self.cfg["elecLoad"].values
        M.profilesEval["bElecLoad"] = 1.0

        # environment temperature
        M.profiles["T_e"] = self.cfg["weather"]["T"].values
        M.profilesEval["T_e"] = 1.0

        # get WACC and cost
        if self.WACC is None:
            if hasattr(M, "WACC"):
                WACC = M.WACC
            else:
                WACC = 0.08
        else:
            WACC = self.WACC
        for dec in M.exVarIx:
            aFactor = (
                ((1 + WACC) ** M.exVarLifetime[dec])
                * WACC
                / (((1 + WACC) ** M.exVarLifetime[dec]) - 1)
                if WACC > 0
                else 1.0 / M.exVarLifetime[dec]
            )
            M.exVarCost[dec] = (
                M.exVarCAPEX[dec] * (aFactor + M.exVarOPEX[dec]) * M.optInterval / 8760
            )

        # subsets
        M.bX_opaque = []
        for comp in ["Walls", "Roof", "Floor"]:
            for dec in M.bInsul[comp]:
                M.bX_opaque.append((comp, dec))
        M.bX_windows = []
        for comp in ["Windows"]:
            for dec in M.bInsul[comp]:
                M.bX_windows.append((comp, dec))
        M.bX_vent = []
        for comp in ["Ventilation"]:
            for dec in M.bInsul[comp]:
                M.bX_vent.append((comp, dec))
        M.bX_solar = []
        for comp in ["Windows", "Walls", "Roof"]:
            for dec in M.bInsul[comp]:
                M.bX_solar.append((comp, dec))

        # activate or deactivate ventilation control
        if self.ventControl:
            warnings.warn("Ventilation control is not validated", UserWarning)
            M.bVentControl = True
        else:
            M.bVentControl = False

        # calculate big M and small m representing max and min heat flow
        M.bM_q = {}
        M.bm_q = {}
        for comp in M.bInsul:
            M.bM_q[comp] = {}
            M.bm_q[comp] = {}
            for dec in M.bInsul[comp]:
                M.bM_q[comp][dec] = M.bH[comp][dec] * (
                    M.bT_comf_ub - (self.cfg["weather"]["T"].min() - 10)
                )
                M.bm_q[comp][dec] = M.bH[comp][dec] * (
                    M.bT_comf_lb - (self.cfg["weather"]["T"].max() + 10)
                )

        # define design temperature
        M.bT_des = self.cfg["design_T_min"]
        return M

    @staticmethod
    def _addVars(M):
        """
        Static method which adds all nodes as dimensioning variables to the 
        optimization
        """

        if M.bRefurbishment:
            # decision variables which refurbishment measure should be chosen
            M.exVars = pyomo.Var(M.exVarIx, within=pyomo.Binary)
            # auxiliary variables for modelling heat flow on thermal mass surface
            M.bP_X = pyomo.Var(M.bX_windows, M.bX_solar, within=pyomo.Binary)
        else:
            # in case refurbishment is activated, those variables are dropped in the presolve and
            # can be set continuous
            M.exVars = pyomo.Var(M.exVarIx, within=pyomo.NonNegativeReals)
            M.bP_X = pyomo.Var(M.bX_windows, M.bX_solar, within=pyomo.NonNegativeReals)

        # temperature variables
        M.bT_m = pyomo.Var(M.timeIndex)
        M.bT_air = pyomo.Var(M.timeIndex)
        M.bT_s = pyomo.Var(M.timeIndex)

        # heat flows directly into the nodes [kW]
        M.bQ_ia = pyomo.Var(M.timeIndex)  # direct to air node
        M.bQ_m = pyomo.Var(M.timeIndex)  # thermal mass
        M.bQ_st = pyomo.Var(M.timeIndex)  # surface of the thermal mass

        # add ventilation heat flow as variable
        M.bQ_ve = pyomo.Var(M.timeIndex)

        # external heat losses including heat exchange
        M.bQ_comp = pyomo.Var(M.insulIx, M.timeIndex)

        # design heat load
        M.bQ_des = pyomo.Var(within=pyomo.NonNegativeReals)

        return M

    @staticmethod
    def _addCons(M):
        """
        Static method which adds all constraints of the 5R1C-model
        """
        # activate or deactivate existing decisions
        for accVar in M.exVarActive:
            M.exVars[accVar].setlb(1.0)
        for accVar in M.exVarInActive:
            M.exVars[accVar].setub(0.0)

        # set that one component for each insulation type has to be chosen
        # (Schuetz et al. 2017 - eq. 2)
        def chooseComp(M, comp):
            return sum(M.exVars[(comp, dec)] for dec in M.bInsul[comp]) == 1

        M.chooseWall = pyomo.Constraint(["Walls"], rule=chooseComp)
        M.chooseRoof = pyomo.Constraint(["Roof"], rule=chooseComp)
        M.chooseWindows = pyomo.Constraint(["Windows"], rule=chooseComp)
        M.chooseVentilation = pyomo.Constraint(["Ventilation"], rule=chooseComp)
        M.chooseFloor = pyomo.Constraint(["Floor"], rule=chooseComp)

        # set auxiliarty equations for binary * binary to linear
        # (Schuetz et al. 2017 - eq. 16)
        def auxiEq1(M, dec_w, comp_w, dec_x, comp_x):
            return M.bP_X[(dec_w, comp_w, dec_x, comp_x)] <= M.exVars[(dec_x, comp_x)]

        M.auxi1Con = pyomo.Constraint(M.bX_windows, M.bX_solar, rule=auxiEq1)
        # (Schuetz et al. 2017 - eq. 17)
        def auxiEq2(M, dec_w, comp_w, dec_x, comp_x):
            return M.bP_X[(dec_w, comp_w, dec_x, comp_x)] <= M.exVars[(dec_w, comp_w)]

        M.auxi2Con = pyomo.Constraint(M.bX_windows, M.bX_solar, rule=auxiEq2)
        # (Schuetz et al. 2017 - eq. 18)
        def auxiEq3(M, dec_w, comp_w, dec_x, comp_x):
            return (
                M.bP_X[(dec_w, comp_w, dec_x, comp_x)]
                >= M.exVars[(dec_w, comp_w)] + M.exVars[(dec_x, comp_x)] - 1
            )

        M.auxi3Con = pyomo.Constraint(M.bX_windows, M.bX_solar, rule=auxiEq3)

        # calculate nodal gains
        # (Schuetz et al. 2017 - eq. 14)
        def gainAirNode(M, t1, t2):
            return M.bQ_ia[t1, t2] == 0.5 * M.profiles["bQ_ig"][t1, t2]

        M.bGainAirNodeCon = pyomo.Constraint(M.timeIndex, rule=gainAirNode)
        # (Schuetz et al. 2017 - eq. 15)
        def gainMassNode(M, t1, t2):
            return M.bQ_m[t1, t2] == 0.5 * M.bA_m / M.bA_tot * M.profiles["bQ_ig"][
                t1, t2
            ] + M.bA_f / M.bA_tot * sum(
                M.profiles["bQ_sol_" + comp + dec][t1, t2] * M.exVars[(comp, dec)]
                for (comp, dec) in M.bX_solar
            )

        M.bGainMassNodeCon = pyomo.Constraint(M.timeIndex, rule=gainMassNode)
        # (Schuetz et al. 2017 - eq. 16)
        def gainSurNode(M, t1, t2):
            return (
                M.bQ_st[t1, t2]
                == (
                    1
                    - 1
                    / M.bConst["h_ms"]
                    / M.bA_tot
                    * sum(
                        M.bU["Windows"][dec] * M.exVars[("Windows", dec)]
                        for (comp, dec) in M.bX_windows
                    )
                )
                * (0.5 * M.profiles["bQ_ig"][t1, t2])
                - 1
                / M.bConst["h_ms"]
                / M.bA_tot
                * sum(
                    M.bH["Windows"][dec_w]
                    * sum(
                        M.profiles["bQ_sol_" + comp_x + dec_x][t1, t2]
                        * M.bP_X[(comp_w, dec_w, comp_x, dec_x)]
                        for (comp_x, dec_x) in M.bX_solar
                    )
                    for (comp_w, dec_w) in M.bX_windows
                )
                + sum(
                    M.profiles["bQ_sol_" + comp + dec][t1, t2] * M.exVars[(comp, dec)]
                    for (comp, dec) in M.bX_solar
                )
                - M.bQ_m[t1, t2]
            )

        M.bGainSurNodeCon = pyomo.Constraint(M.timeIndex, rule=gainSurNode)

        # if a component is not chosen, force heat flow to be zero
        # (Schuetz et al. 2017 - eq. 23-25)
        def heatFlowZeroUB(M, t1, t2, comp, dec):
            return (
                M.bQ_comp[comp, dec, t1, t2] <= M.exVars[comp, dec] * M.bM_q[comp][dec]
            )

        M.heatFlowZeroUBCon = pyomo.Constraint(
            M.timeIndex, M.insulIx, rule=heatFlowZeroUB
        )
        # (Schuetz et al. 2017 - eq. 23-25)
        def heatFlowZeroLB(M, t1, t2, comp, dec):
            return (
                M.bQ_comp[comp, dec, t1, t2] >= M.exVars[comp, dec] * M.bm_q[comp][dec]
            )

        M.heatFlowZeroLBCon = pyomo.Constraint(
            M.timeIndex, M.insulIx, rule=heatFlowZeroLB
        )

        # if the component is chosen, calculate the heat flow for opaque components
        # (Schuetz et al. 2017 - eq. 8, 23-25)
        def heatFlowActiveUB(M, t1, t2, comp, dec):
            return (
                M.bH[comp][dec] * (M.bT_m[t1, t2] - M.profiles["T_e"][t1, t2])
                - M.bQ_comp[comp, dec, t1, t2]
                <= (1 - M.exVars[comp, dec]) * M.bM_q[comp][dec]
            )

        M.heatFlowActiveUBCon = pyomo.Constraint(
            M.timeIndex, M.insulIx, rule=heatFlowActiveUB
        )
        # (Schuetz et al. 2017 - eq. 8, 23-25)
        def heatFlowActiveLB(M, t1, t2, comp, dec):
            return (
                M.bH[comp][dec] * (M.bT_m[t1, t2] - M.profiles["T_e"][t1, t2])
                - M.bQ_comp[comp, dec, t1, t2]
                >= (1 - M.exVars[comp, dec]) * M.bm_q[comp][dec]
            )

        M.heatFlowActiveLBCon = pyomo.Constraint(
            M.timeIndex, M.insulIx, rule=heatFlowActiveLB
        )

        # TODO: update ventilation flow
        if not M.bVentControl:
            pass
        else:
            warnings.warn("Ventilation control model is not validated", UserWarning)
            # add ventilation flow as control variable
            M.bH_ve_con = pyomo.Var(M.timeIndex)

            # define the working point of the temperature
            M.bT_air_WP = 20.0
            # working point of the ventilation flow is the average ventilation
            # flow, which is pre-defined

            # add the lineraized ventilation flow as constraint
            def ventFlowLin(M, t1, t2):
                return M.bQ_ve[t1, t2] == (
                    M.bH_ve_con[t1, t2] * (M.bT_air_WP - M.profiles["T_e"][t1, t2])
                    + M.bH_ve * (M.bT_air[t1, t2] - M.bT_air_WP)
                )

            M.ventFlowLinCon = pyomo.Constraint(M.timeIndex, rule=ventFlowLin)
            #
            # minimal ventilation flow is the infiltration rate
            def ventMin(M, t1, t2):
                return M.bH_ve_con[t1, t2] >= M.bH_air_inf

            M.ventMinCon = pyomo.Constraint(M.timeIndex, rule=ventMin)

            # maximal is three times the average flow (ASSUMPTION)
            # minimal ventilation flow is the infiltration rate
            def ventmax(M, t1, t2):
                return M.bH_ve_con[t1, t2] <= M.bH_ve * 10

            M.ventMaxCon = pyomo.Constraint(M.timeIndex, rule=ventmax)

            # average ventilation flow per day
            M.stepsPerDay = 24
            M.dayIndex = range(0, int(len(M.timeIndex) / M.stepsPerDay))
            # TODO depending on time resolution
            # TODO for typical days
            def ventMean(M, d):
                return (
                    sum(
                        M.bH_ve_con[1, t]
                        for t in range(d * M.stepsPerDay, (d + 1) * M.stepsPerDay)
                    )
                    >= M.stepsPerDay * M.bH_ve
                )

            M.ventMeanCon = pyomo.Constraint(M.dayIndex, rule=ventMean)
        #        M.bH_air_use = V_air * rho_air * C_air * self.cfg['n_air_use']  / 3600
        #        M.bH_air_inf = V_air * rho_air * C_air * self.cfg['n_air_infiltration']  / 3600

        # define the three energy balances

        # 1 -storage capacity of the building --> include typday behavior
        if M.hasTypPeriods:

            def balMassNode(M, t1, t2):
                return (
                    M.bH_ms * (M.bT_m[t1, t2] - M.bT_s[t1, t2])
                    + sum(M.bQ_comp[comp, dec, t1, t2] for (comp, dec) in M.bX_opaque)
                    + M.bH_door * (M.bT_m[t1, t2] - M.profiles["T_e"][t1, t2])
                    == M.bQ_m[t1, t2]
                    - M.bC_m * (M.bT_m[t1, t2 + 1] - M.bT_m[t1, t2]) / M.stepSize
                )

            M.balMassNodeCon = pyomo.Constraint(
                M.typPeriodIx, M.intraIx[:-1], rule=balMassNode
            )

            def balMassNodeStartEnd(M, period):
                return (
                    M.bH_ms
                    * (M.bT_m[period, M.intraIx[-1]] - M.bT_s[period, M.intraIx[-1]])
                    + sum(
                        M.bQ_comp[comp, dec, period, M.intraIx[-1]]
                        for (comp, dec) in M.bX_opaque
                    )
                    + M.bH_door
                    * (
                        M.bT_m[period, M.intraIx[-1]]
                        - M.profiles["T_e"][period, M.intraIx[-1]]
                    )
                    == M.bQ_m[period, M.intraIx[-1]]
                    - M.bC_m
                    * (M.bT_m[period, M.intraIx[0]] - M.bT_m[period, M.intraIx[-1]])
                    / M.stepSize
                )

            M.balMassNodeStartEndCon = pyomo.Constraint(
                M.typPeriodIx, rule=balMassNodeStartEnd
            )
        else:

            def balMassNode(M, t1, t2):
                return (
                    M.bH_ms * (M.bT_m[t1, t2] - M.bT_s[t1, t2])
                    + sum(M.bQ_comp[comp, dec, t1, t2] for (comp, dec) in M.bX_opaque)
                    + M.bH_door * (M.bT_m[t1, t2] - M.profiles["T_e"][t1, t2])
                    == M.bQ_m[t1, t2]
                    - M.bC_m * (M.bT_m[t1, t2 + 1] - M.bT_m[t1, t2]) / M.stepSize
                )

            M.balMassNodeCon = pyomo.Constraint(M.timeIndex[:-1], rule=balMassNode)

            def balMassNodeStartEnd(M):
                return (
                    M.bH_ms * (M.bT_m[M.timeIndex[-1]] - M.bT_s[M.timeIndex[-1]])
                    + sum(
                        M.bQ_comp[comp, dec, M.timeIndex[-1]]
                        for (comp, dec) in M.bX_opaque
                    )
                    + M.bH_door
                    * (M.bT_m[M.timeIndex[-1]] - M.profiles["T_e"][M.timeIndex[-1]])
                    == M.bQ_m[M.timeIndex[-1]]
                    - M.bC_m
                    * (M.bT_m[M.timeIndex[0]] - M.bT_m[M.timeIndex[-1]])
                    / M.stepSize
                )

            M.balMassNodeStartEndCon = pyomo.Constraint(rule=balMassNodeStartEnd)

        def balSurfaceNode(M, t1, t2):
            return (
                M.bH_ms * (M.bT_s[t1, t2] - M.bT_m[t1, t2])
                + M.bH_is * (M.bT_s[t1, t2] - M.bT_air[t1, t2])
                + sum(M.bQ_comp[comp, dec, t1, t2] for (comp, dec) in M.bX_windows)
                == M.bQ_st[t1, t2]
            )

        M.balSurfaceNodeCon = pyomo.Constraint(M.timeIndex, rule=balSurfaceNode)

        def balAirNode(M, t1, t2):
            return sum(
                M.bQ_comp[comp, dec, t1, t2] for (comp, dec) in M.bX_vent
            ) + M.bH_is * (M.bT_air[t1, t2] - M.bT_s[t1, t2]) == M.bQ_st[t1, t2] - sum(
                M.connectVars[c, t1, t2] for c in M.bConnectedCool
            ) + sum(
                M.connectVars[c, t1, t2] for c in M.bConnectedHeat
            )

        M.balAirNodeCon = pyomo.Constraint(M.timeIndex, rule=balAirNode)

        # add calculation of design heat load including refurbishment decisions
        # adjustment factor for floor
        M.bDesignAdjust = {
            "Floor": 1.45,
            "Walls": 1.0,
            "Roof": 1.0,
            "Windows": 1.0,
            "Ventilation": 1.0,
        }

        def designLoad(M):
            return M.bQ_des == (
                (
                    M.bH_door
                    + sum(
                        M.bH[comp][dec] * M.bDesignAdjust[comp] * M.exVars[comp, dec]
                        for comp, dec in M.insulIx
                    )
                )
                * (22.917 - M.bT_des)
            )

        M.designLoadCon = pyomo.Constraint(rule=designLoad)

        if M.bRefurbishment:
            # define that a controller has to be installed for the installation of another one
            def controlOrder1(M):
                return (
                    M.exVars[("Control", "Occupancy")]
                    <= M.exVars[("Control", "SmartThermostat")]
                )

            M.controlOrder1 = pyomo.Constraint(rule=controlOrder1)

            def controlOrder2(M):
                return (
                    M.exVars[("Control", "Occupancy")]
                    <= M.exVars[("Control", "NightReduction")]
                )

            M.controlOrder2 = pyomo.Constraint(rule=controlOrder2)

        # define upper bound of the temperature
        def temperatureUB(M, t1, t2):
            return M.bT_air[t1, t2] <= (
                M.bT_comf_lb
                + (M.bT_comf_ub - M.bT_comf_lb)
                * M.exVars[("Control", "SmartThermostat")]
                - (M.bT_comf_ub - 30.0)
                * M.profiles["bOccNotHome"][t1, t2]
                * M.exVars[("Control", "Occupancy")]
            )

        M.temperatureUBCon = pyomo.Constraint(M.timeIndex, rule=temperatureUB)

        # define lower bound of the temperature
        def temperatureLB(M, t1, t2):
            return M.bT_air[t1, t2] >= (
                M.bT_comf_lb
                - (M.bT_comf_lb - 18.0)
                * M.profiles["bOccSleeping"][t1, t2]
                * M.exVars[("Control", "NightReduction")]
                - (M.bT_comf_lb - 14.0)
                * M.profiles["bOccNotHome"][t1, t2]
                * M.exVars[("Control", "Occupancy")]
            )

        M.temperatureLBCon = pyomo.Constraint(M.timeIndex, rule=temperatureLB)

        # add electricity constraint
        def electricityBalanceRule(M, t1, t2):
            return (
                sum(M.connectVars[c, t1, t2] for c in M.bConnectedElec)
                == M.profiles["bElecLoad"][t1, t2]
            )

        M.electricityBalanceCon = pyomo.Constraint(
            M.timeIndex, rule=electricityBalanceRule
        )

        return M

    def _readResults(self, M):
        """
        Reads the results from the optimization model M and stores them in a
        result dictionary.
        
        """
        self.detailedResults["Heating Load"] = np.array(
            [
                sum(M.connectVars[c, t].value for c in M.bConnectedHeat)
                for t in M.fullTimeIndex
            ]
        )
        self.detailedResults["Cooling Load"] = np.array(
            [
                sum(M.connectVars[c, t].value for c in M.bConnectedCool)
                for t in M.fullTimeIndex
            ]
        )
        self.detailedResults["T_air"] = np.array(
            [M.bT_air[t].value for t in M.fullTimeIndex]
        )
        self.detailedResults["T_s"] = np.array(
            [M.bT_s[t].value for t in M.fullTimeIndex]
        )
        self.detailedResults["T_m"] = np.array(
            [M.bT_m[t].value for t in M.fullTimeIndex]
        )
        self.detailedResults["T_e"] = np.array(
            [M.profiles["T_e"][t] for t in M.fullTimeIndex]
        )
        self.detailedResults["Electricity Load"] = self.cfg["elecLoad"].values

        for dec in M.exVarIx:
            if M.exVars[dec].stale:
                if M.exVars[dec].lb == M.exVars[dec].ub:
                    M.exVars[dec].value = M.exVars[dec].lb
                else:
                    warnings.warn(
                        "Stale value in result of "
                        + str(dec)
                        + " detected. Result is set to the lb "
                        "of the variable",
                        UserWarning,
                    )
                    M.exVars[dec].value = M.exVars[dec].lb
            self.detailedRefurbish[dec] = {}
            self.detailedRefurbish[dec]["Capacity"] = M.exVars[dec].value
            self.detailedRefurbish[dec]["FixCost"] = (
                M.exVarCost[dec] * M.exVars[dec].value
            )
            self.detailedRefurbish[dec]["CAPEX"] = (
                M.exVarCAPEX[dec] * M.exVars[dec].value
            )
        self.static_results["Capacity"] = M.bQ_des.value
        self.static_results["FixCost"] = 0
        self.static_results["CAPEX"] = 0
        self.static_results["OPEX fix"] = 0.0
        self.static_results["VarCost"] = 0.0
        self.static_results["OPEX var"] = 0.0
        self.static_results["OPEX"] = 0.0
        self.result_load = self.detailedResults["Heating Load"]

        return

    def calcDesignHeatLoad(self):
        """
        Calculates the design heat load which is needed to satisfy the
        nominal outdoor temperature.
        
        Returns
        -------
        designHeatLoad [kW]
        """

        b = self.cfg
        designHeatLoad = (
            b["A_Roof_1"] * b["U_Roof_1"] * b["b_Transmission_Roof_1"]
            + b["A_Roof_2"] * b["U_Roof_2"] * b["b_Transmission_Roof_2"]
            + b["A_Wall_1"] * b["U_Wall_1"] * b["b_Transmission_Wall_1"]
            + b["A_Wall_2"] * b["U_Wall_2"] * b["b_Transmission_Wall_2"]
            + b["A_Wall_3"] * b["U_Wall_3"] * b["b_Transmission_Wall_3"]
            + b["A_Window"] * b["U_Window"]
            + b["A_Door_1"] * b["U_Door_1"]
            + (
                b["A_ref"]
                * b["h_room"]
                * 1.2
                * 1006
                * (b["n_air_infiltration"] + b["n_air_use"])
                / 3600
            )
        ) * (22.917 - self.cfg["design_T_min"]) + (
            (
                b["A_Floor_1"] * b["U_Floor_1"] * b["b_Transmission_Floor_1"]
                + b["A_Floor_2"] * b["U_Floor_2"] * b["b_Transmission_Floor_2"]
            )
            * (22.917 - self.cfg["design_T_min"])
            * 1.45
        )
        return designHeatLoad / 1000

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
        self.design_Q = self.calcDesignHeatLoad()

        self.design_T_return = 25
        # derive the heat transfer coefficient of the heating system kW/K

    #        self.design_H_heat = self.design_Q / (self.cfg['T_sup'] - self.design_T_return)

    def sim5R1C(self, solver=None, tee=True):
        """
        Simulate the building model independently from any other optimization
        based on the 5R1C model of the DIN 137900/Thomas Schuetz 2017 with
        pyomo.
        
        Returns
        -------
        None, but in Building.results or Building.detailedResults can the
        results be shown as pandas.DataFrame
        """

        # take default solver from environment variables
        if solver is None:
            try:
                solver = os.environ["SOLVER"]
            except KeyError:
                solver = "gurobi"

        self.WACC = None
        self.lifetime = 40

        self.heatCost = 0.08
        self.coolCost = 0.02
        self.elecCost = 0.25

        # init concrete model
        M = pyomo.ConcreteModel()

        # set time index
        M.timeIndex = [(1, t) for t in range(len(self.times))]
        M.fullTimeIndex = M.timeIndex
        timediff = self.times[1] - self.times[0]
        M.stepSize = timediff.total_seconds() / 3600  # step size in hours
        M.hasTypPeriods = False
        M.WACC = 0.08
        M.optInterval = 8760

        # create heat and electricity input (PLEASE DO NOT CHANGE THE
        # SYNTAX OF THIS BLOCK - equivalent to enercore)
        M.connectIndex = [
            ("Heat", "Building"),
            ("Electricity", "Building"),
            ("Cool", "Building"),
        ]
        M.connectVars = pyomo.Var(
            M.connectIndex, M.timeIndex, within=pyomo.NonNegativeReals
        )
        M.bConnectedHeat = [("Heat", "Building")]
        M.bConnectedElec = [("Electricity", "Building")]
        M.bConnectedCool = [("Cool", "Building")]

        # init all indices and parameters
        M = self._initOpti(M)
        M = self._addOpti(M)

        # all profiles have been initialized here --> divide to constraint part
        if not M.hasTypPeriods:
            for profile_name in M.profiles:
                M.profiles[profile_name] = {
                    timeIndex: M.profiles[profile_name][i]
                    for i, timeIndex in enumerate(M.timeIndex)
                }

        # add other variables
        M = self._addVars(M)

        # add constraints
        M = self._addCons(M)

        # add cost
        M.bHeatCost = self.heatCost  # eur/kWh
        M.bCoolCost = self.coolCost  # eur/kWh
        M.bEectricityCost = self.elecCost  # eur/kWh

        # add the maximal heat load as soft constraint in order to avoid infeasibility
        M.bMaxLoadViolation = pyomo.Var(within=pyomo.NonNegativeReals)

        # create a maximal heat load as soft constraint in order to avoid infeasibility
        def maxHeatingLoad(M, t1, t2):
            return (
                sum(M.connectVars[c, t1, t2] for c in M.bConnectedHeat)
                - M.bMaxLoadViolation
                <= M.bMaxLoad
            )

        M.maxHeatingLoadCon = pyomo.Constraint(M.timeIndex, rule=maxHeatingLoad)

        def maxCoolingLoad(M, t1, t2):
            return (
                sum(M.connectVars[c, t1, t2] for c in M.bConnectedCool)
                - M.bMaxLoadViolation
                <= M.bMaxLoad
            )

        M.maxCoolingLoadCon = pyomo.Constraint(M.timeIndex, rule=maxCoolingLoad)

        # objective function
        def minLoad(M):
            return (
                sum(
                    sum(M.connectVars[c, t] for c in M.bConnectedCool)
                    for t in M.timeIndex
                )
                * M.bCoolCost
                + sum(
                    sum(M.connectVars[c, t] for c in M.bConnectedHeat)
                    for t in M.timeIndex
                )
                * M.bHeatCost
                + sum(
                    sum(M.connectVars[c, t] for c in M.bConnectedElec)
                    for t in M.timeIndex
                )
                * M.bEectricityCost
                + sum(M.exVarCost[dec] * M.exVars[dec] for dec in M.exVarIx)
                + M.bMaxLoadViolation * 1e5  # penalty term
            )

        M.obj = pyomo.Objective(rule=minLoad)

        # optimize
        self.M = M
        optprob = opt.SolverFactory(solver)
        optprob.options = utils.manageSolverOpts(
            solver, {"Threads": 1, "LogFile": ""}
        )  # , 'LogToConsole':0

        # solve
        self.opt_results = optprob.solve(M, tee=tee)

        if not M.bMaxLoadViolation.value is None:
            if M.bMaxLoadViolation.value > 0.0:
                warnings.warn(
                    "Maximal heat load exceeded by "
                    + str(round(M.bMaxLoadViolation.value, 3))
                    + " kW",
                    UserWarning,
                )

        # read results
        self._readResults(M)

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

    #    example= Building( a_ref = 300.,
    #        refurbishment = False, isRefurbished = False,
    #        nightReduction = True,
    #        comfortT_lb = 20, comfortT_ub = 26,
    #    )
    designLoad = example.calcDesignHeatLoad()
    #    buildingYear = 1990, buildingType = 'SFH',
    #                 buildingClassification = 'Gen', weatherData = None,
    #                 eastOrOverall = 'N', annualHeatDemand = None,occControl = True,
    #                 n_persons = 4, refurbishment = False, ventControl = False,

    # example.sim5R1C(solver="gurobi", tee=True)
    example.sim5R1C(solver="gurobi", tee=True)
    example.detailedResults["5R1C Cooling Load T=20-26"] = -example.detailedResults[
        "Cooling Load"
    ]
    example.detailedResults["5R1C Heating Load T=20-26"] = example.detailedResults[
        "Heating Load"
    ]

    example.detailedResults[
        ["5R1C Heating Load T=20-26", "5R1C Cooling Load T=20-26"]
    ].plot()
    plt.show()
    example.detailedResults["T_out"] = example.cfg["weather"]["T"]
    example.detailedResults[["T_air", "T_s", "T_m", "T_out"]].plot()
    plt.show()
    example.M.exVars.display()

    example.detailedResults[["5R1C Heating Load T=20-26", "T_air"]].plot(subplots=True)
    plt.show()
