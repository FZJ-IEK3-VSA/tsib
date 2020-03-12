# -*- coding: utf-8 -*-
"""
Created on Fri Apr 08 11:33:01 2016

@author: Leander Kotzur
"""

import time
import os

import tsib

def test_heatpump_cop():

    # get weather data
    try_data, loc = tsib.readTRY()

    # get cop profile
    CoP = tsib.simHeatpump(try_data["T"], COP_limit = 6.)

    assert CoP.max() <= 6.0
    assert CoP.min() >= 0.0

