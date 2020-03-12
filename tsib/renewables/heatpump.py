"""
Created on Sat Dec 10 12:40:17 2017

@author: Leander Kotzur
"""

import copy

import numpy as np

def simHeatpump(T_cold, T_hot = 50., efficiency = 0.45, T_limit = -20., COP_limit = 7.):
    '''
    Creates a timedepent Coefficient of Performance (COP) based on the potential carnot
    efficiency and a quality grade/efficiency of the system.

    Parameters
    -----------
    T_cold: float, np.array or list, required
        Outside temperature or in case of a ground source heat pump the ground temperature in degree C.
    T_hot: float,np.array or list, optional (default: 50.0)
        Temperature level in degree C at which the heat pump has to provide the heat.
    efficiency: float, optional (default: 0.45) (Lauinger,2016)
        Factor to multiplicate the carnot efficiency with.
    T_limit: float, optional (default: -5)
        Temperature at which the heat pump is shut off.
    COP_limit: float, optional (default: 6.)
        Maximum Coefficient of Performance that can be reached due to the 
        technical limitations of the system.

    Returns
    ----------
    CoP: pandas.Series()
        Time series of the Coefficient of Performance.
    '''
    
    
    # calculate timedependet COP      
    CoP = efficiency * (T_hot + 273.15) / (T_hot - T_cold) 
    
    # cutt of temperature
    CoP[T_cold < T_limit ] = 0.0

    # limit too high CoPS
    if len(CoP) >  1:
        CoP[CoP > COP_limit] = COP_limit
    else:
        CoP = min(CoP, COP_limit)

    return CoP