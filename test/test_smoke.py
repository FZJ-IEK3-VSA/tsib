# -*- coding: utf-8 -*-
"""
Created on Tue Oct 29 16:17:06 2019

@author: Kevin Knosala 
"""

import os
import pytest

import pandas as pd

import tsib



@pytest.mark.skip(reason="Seperation of weather and building configuration yet not implemented.")
def test_smoke():

    # Read buildings from episcope database
    bdgs_df = pd.read_csv(os.path.join("tsib/data/episcope/episcope.csv"), index_col=1)
    # Select random building by building code
    bdg_dict = bdgs_df.loc['DE.N.SFH.10.Gen.ReEx.001.001'].to_dict()

    # Create building configuration object
    # NOTE: The right kwargs have to be created
    bdg_cfg = tsib.BuildingConfiguration(bdg_dict)

    # Get weather data
    try_data, loc = tsib.readTRY()

    # Initialize a building object with this configuration for given weather
    # NOTE: Weather data seperated from building configuration
    bdg_obj = tsib.Building(configurator=bdg_cfg, weather=try_data)

    # Calculate loads
    bdg_obj.getLoad()
