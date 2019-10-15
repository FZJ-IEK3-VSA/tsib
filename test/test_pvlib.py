# -*- coding: utf-8 -*-
"""
Created on  Tue Sep 24 09:33:07 2019

@author: Kevin
"""
import time

import tsib
import pandas as pd


def test_pvlib():

    starttime = time.time()

    # Run PV simulation with Pv_lib from Sandia (default values)
    try_data, loc = tsib.readTRY()
    tmy_data = tsib.TRY2TMY(try_data)

    specific_load, space_coverage = tsib.simPhotovoltaic(tmy_data)

    print("Summary of specific load profile [kW/kWp]:")
    print(specific_load.describe())
    print("Space coverage [m^2/kWp]: " + str(space_coverage))
    print("Profile generation took " + str(time.time() - starttime) + "s")


if __name__ == "__main__":
    test_pvlib()