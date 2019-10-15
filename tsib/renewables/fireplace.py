"""
Created on Sat Dec 10 12:40:17 2017

@author: Leander Kotzur
"""

import warnings

import pandas as pd
import numpy as np


def simFireplace(
    temperature,
    occ_act,
    n_ovens=1,
    T_oven_on=5,
    t_cool=5.0,
    fullloadSteps=450,
    seed=None,
):
    """
    Creates the profile of the heating of wood ovens based on the outside
    temperature and the activity of the occupants. The profile is generated
    based on stochastics.
    
    Parameters
    ----------
    temperature: pandas.Series(), required
        Outside temperature profile of the location.
    occ_act: pandas.Series(), required
        Series of values between 0 and 1 desribing the share of overall active
        occpants at every time step.
    n_ovens: int, optional (default:1)
        Number of ovens in the building.
    T_oven_on: int or float, optional (default:5.)
        Threeshold outside temperature [°C] when the oven is turned on.
    t_cool: int, optional (default:5)
        Number of timesteps till when the oven is cooled down again.
    fullloadSteps: int or float, optional (default:450)
        Resulting number of full load timesteps. Attention: This value is not
        exact since it is a stochastic profile generation.
    seed: int, optional (default:None)
        Seed required for reproduceability. 
        If none, it is completely random.
        
    Returns
    ----------
    load_profile: pandas.Series()
        Relative load profile of the ovens in kW/kWp.
    """

    # Oven is only turned under a temperature threeshold
    tempBool = temperature < T_oven_on

    # Increase probability that the oven is turned on in case it is colder outside
    relCold = (T_oven_on - temperature) / (T_oven_on - temperature.min())

    # Caclulate fire activation probability
    prob = occ_act.values * tempBool.values * relCold.values

    # avoid rounding errors
    prob[prob < 0] = 0

    load_profile = pd.Series(0, index=temperature.index)
    for n in range(int(n_ovens)):

        # Modifier to reduce probability in order to fit the full load hours
        p_mod = fullloadSteps / (prob.sum() * t_cool / 2)

        overallProb = prob * p_mod
        overallProb[overallProb > 1.0] = 1.0

        # Binary decision if an oven can be activated
        initLogArr = pd.Series(np.random.RandomState(seed).binomial(1, overallProb))

        # create the profile
        heatLoad = []
        loadBefore = 0
        for initLog in initLogArr:
            if initLog:
                load = 1.0
            else:
                if loadBefore > 0.001:
                    load = loadBefore - 1.0 / t_cool
                else:
                    load = 0
            heatLoad.append(load)
            loadBefore = load
        load_profile += pd.Series(heatLoad, index=temperature.index)
    profile = load_profile / n_ovens
    if abs(profile.sum() - fullloadSteps) > 100:
        warnings.warn(
            "Fullload hour deviation is higher than 100. "
            + "Input parameters make it difficult or impossible "
            + "to generate the expected profile"
        )
    return load_profile / n_ovens
