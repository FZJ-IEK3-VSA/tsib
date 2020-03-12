# -*- coding: utf-8 -*-
"""
Created on Fri Apr 08 11:33:01 2016

@author: Leander Kotzur
"""

import time
import os

import pandas as pd

import tsib

def test_fireplace():

    # get weather data
    try_data, loc = tsib.readTRY()
    
    n_occs = 3
    occupancy_data = tsib.simSingleHousehold(n_occs, 2010, get_hot_water=True, resample_mean=True)

    # test for different seeds
    fullloadhours = 450

    # test for different seeds
    for seed in range(5):
        fireplaceprofile = tsib.simFireplace(
            try_data["T"],
            occupancy_data["OccActive"] / n_occs,
            n_ovens=1,
            T_oven_on=5,
            t_cool=5.0,
            fullloadSteps=fullloadhours,
            seed=seed,
        )
        assert fireplaceprofile.sum() < fullloadhours + 100
        assert fireplaceprofile.sum() > fullloadhours - 100

