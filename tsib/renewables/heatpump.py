"""
Created on Sat Dec 10 12:40:17 2017

@author: Leander Kotzur
"""

import copy

import numpy as np

def simHeatpump(T_cold, T_hot = 50., efficiency = 0.45, T_limit = -5.):
    '''
    Creates a timedepent Coefficient of Performance (Cop)

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

    Returns
    ----------
    CoP: pandas.Series()
        Time series of the Coefficient of Performance.
    '''
    
    # limit outside temperature to a maximum value
    T_cold_limit = copy.deepcopy(T_cold)
    if isinstance(T_cold_limit, (np.ndarray, np.generic) ):
        T_cold_limit[ T_cold_limit > 15.0 ] = 15.0
    
    # calculate timedependet COP      
    CoP = efficiency * (T_hot + 273.15) / (T_hot - T_cold_limit) 

    # cutt of temperature
    CoP[T_cold<T_limit ] = 0.0

    return CoP