# -*- coding: utf-8 -*-
"""
Created on Fri Apr 08 11:33:01 2016

@author: Leander Kotzur
"""


import tsib

def test_get_ID():
    # parameterize a building
    bdgcfg = tsib.BuildingConfiguration(
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

    bdgObj = tsib.Building(configurator=bdgcfg)
    print('ID is : ' + str(bdgObj.ID))

    return


def test_set_ID():
    # parameterize a building
    bdgcfg = tsib.BuildingConfiguration(
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
    bdgObj = tsib.Building(configurator=bdgcfg)
    
    bdgObj.ID = 'custom'

    if not bdgObj.ID == 'custom':
        raise ValueError()

    return
