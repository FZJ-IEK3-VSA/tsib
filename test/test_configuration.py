# -*- coding: utf-8 -*-
"""
Created on Fri Apr 08 11:33:01 2016

@author: Leander Kotzur
"""


import tsib.buildingconfig as cfg

def test_configuration_1():
    # parameterize a building
    bdgcfg = cfg.BuildingConfiguration(
        {
            "refurbishment": False,
            "nightReduction": False,
            "occControl": False,
            "capControl": True,
            "n_persons": 2,
            "roofOrientation": 0.0,
            "n_apartments": 1,
            "longitude": 58.,
            "latitude": 8.,
        }
    )
    test = bdgcfg.getBdgCfg()
    return


def test_configuration_2():
    # parameterize a building
    bdgcfg = cfg.BuildingConfiguration(
        {
            "buildingYear": 1980,
            "n_floors": 3,
            "n_persons": 2,
            "roofOrientation": 0.0,
            "n_apartments": 1,
            "longitude": 58.,
            "latitude": 8.,
        }
    )
    test = bdgcfg.getBdgCfg()
    return
