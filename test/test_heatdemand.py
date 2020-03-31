# -*- coding: utf-8 -*-
"""
Created on Fri Apr 08 11:33:01 2016

@author: Leander Kotzur
"""

import time
import os

import pandas as pd
import numpy as np

import tsib
import tsib.data


def test_scaleHeatDemand():
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
        }
    )
    bdgObj = tsib.Building(configurator=bdgcfg)

    # get the design heat load
    origDesignLoad = bdgObj.thermalmodel.calcDesignHeatLoad()

    # scale to a reduced value
    bdgObj.thermalmodel.scaleHeatLoad(scale=0.5)

    # get updated load
    reducedDesignLoad = bdgObj.thermalmodel.calcDesignHeatLoad()

    np.testing.assert_almost_equal(reducedDesignLoad/0.5, origDesignLoad, decimal=2)

if __name__ == "__main__":
    test_scaleHeatDemand()
