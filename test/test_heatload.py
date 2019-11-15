# -*- coding: utf-8 -*-
"""
Created on Fri Apr 08 11:33:01 2016

@author: Leander Kotzur
"""

import time
import os

import pandas as pd

import tsib
import tsib.data

def test_heatload():

    starttime = time.time()

    # get raw building data set
    buildingSet = pd.read_csv(
        os.path.join(tsib.data.PATH, "episcope", "tabula_DE_wPersons.csv"), header=0, index_col=0
    )

    # get a random building ID
    ix = 24
    ID = buildingSet.index[ix]

    # get time series data
    try_data, loc = tsib.readTRY(try_num=4)

    # parameterize a building
    bdgcfg = tsib.BuildingConfiguration(
        {
            "ID": ID,
            "weatherData": try_data,
            "weatherID": "TRY_4",
            "refurbishment": False,
            "nightReduction": False,
            "occControl": False,
            "capControl": True,
            "n_persons": 2,
            "comfortT_lb": 20.0,
            "comfortT_ub": 26.0,
            "roofOrientation": 0.0,
            "n_apartments": 1,
            "longitude": loc["longitude"],
            "latitude": loc["latitude"],
        }
    )
    # setup a building with the configuration
    bdgObj = tsib.Building(configurator=bdgcfg)

    # get the occupancy profiles to manipulate them
    bdgObj._get_occupancy_profile( bdgObj.cfg)

    # manipulate internal gains with tabula mean value
    bdgObj.cfg["Q_ig"] = (
        bdgObj.cfg["Q_ig"] * 15.552 / (bdgObj.cfg["Q_ig"].sum() / bdgObj.cfg["A_ref"])
    )

    # run simulation
    bdgObj.getHeatLoad()  # take solver from environment variable

    # get specific heat demand
    q_sim = bdgObj.timeseries["Heating Load"].sum() / bdgObj.cfg["A_ref"]

    # get calculated heat demand by IWU
    q_iwu = buildingSet.loc[ID, "q_h_nd"]

    print("Profile generation took " + str(time.time() - starttime))

    print("Spec. heat demand IWU [kWh/m²/a]: " + str(round(q_iwu)))
    print("Spec. heat demand 5R1C [kWh/m²/a]: " + str(round(q_sim)))

    if abs(q_sim - q_iwu) > 20:
        raise ValueError(
            "The difference between simulation and the values listed by the IWU is too high."
        )

    if ix == 24:
        if not round(q_sim) == 197.0:
            raise ValueError("Different result for mean heat load than expected.")
    return




if __name__ == "__main__":
    test_heatload()
