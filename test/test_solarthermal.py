"""
Created on 13.10.2019

@author: Leander Kotzur
"""


import numpy as np

import tsib
import tsib.renewables.solar as sol 


def test_solarthermal():
    # read weather data
    try_data, loc = tsib.readTRY(try_num=4)
    
    # format to tmy format
    tmy_data = tsib.TRY2TMY(try_data)
    
    # get solar thermal potential
    spec_load_st = sol.simSolarThermal(tmy_data,     
                                latitude=loc['latitude'],
                                longitude=loc['longitude'],
                                surface_tilt=30,
                                surface_azimuth=180,
                                )
    
    # expected yield in kWh/sqm/a
    np.testing.assert_array_almost_equal(spec_load_st.sum(), 584.956523695128, decimal=0
    )



