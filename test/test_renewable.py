# -*- coding: utf-8 -*-
"""
Created on Fri Apr 08 11:33:01 2016

@author: Leander Kotzur
"""

import time
import os

import pandas as pd

import tsib
import tsib
import tsib.data

def test_renewable():

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
            "roofOrientation": 0.0,
            "longitude": loc["longitude"],
            "latitude": loc["latitude"],
        }
    )

    # setup a building with the configuration
    bdgObj = tsib.Building(configurator=bdgcfg)

    # get the renewable profiles to manipulate them
    bdgObj.getRenewables()

    #TODO
