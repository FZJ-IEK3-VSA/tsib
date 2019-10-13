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
            "latitude": 49.,
            "longitude": 12.,
        }
    )
    test = bdgcfg.getBdgCfg()
    return


def test_configuration_2():
    # parameterize a building
    bdgcfg = cfg.BuildingConfiguration(
        {
            "buildingYear": 1980,
            "n_persons": 2,
            "roofOrientation": 0.0,
            "n_apartments": 2,
            "a_ref": 300.,
            "surrounding": "Detached",
            "latitude": 52.,
            "longitude": 13.,
        }
    )
    test = bdgcfg.getBdgCfg()
    return


def test_configuration_3():
    # parameterize a building
    kwgs = {
        "buildingYear": 1990,
        "latitude": 52.0,
        "longitude": 13.0,
        "comfortT_lb": 21.,
        "comfortT_ub": 24.,
        "WACC": 0.03,
        "roofTilt": 45.0,
        "surrounding": "Semi",
        "n_apartments": 2,
        "a_ref_app": 100.,
        "n_persons": 2,
        "roofOrientation": 135.0,
        "costdata": "default_2016",
        "capControl": True,
    }

    bdgcfg = cfg.BuildingConfiguration(kwgs)

    test = bdgcfg.getBdgCfg()

    return