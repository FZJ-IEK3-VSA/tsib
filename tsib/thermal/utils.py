# -*- coding: utf-8 -*-
"""
Created on Sat Dec 10 12:40:17 2016

@author: Leander Kotzur
"""

class PowerDict(dict):  
    '''
    Dictionary with additional functions
    '''
    def __init__(self, parent=None, key=None):
        self.parent = parent
        self.key = key

    def __missing__(self, key): 
        '''
        Creation of subdictionaries on fly
        '''
        self[key] = PowerDict(self, key)
        return self[key]

    def append(self, item):
        '''
        additional append function for lists in dict
        '''
        self.parent[self.key] = [item]

    def __setitem__(self, key, val):
        dict.__setitem__(self, key, val)
        try:
            val.parent = self
            val.key = key
        except AttributeError:
            pass


def manageSolverOpts(solver, solverOpts):
    '''
    Adds solver specific options.

    Parameters
    ----------
    solver: str, required
        Solver name used by pyomo
    solverOpts: dict, required
        Other solveroptions
    
    Returns
    -------
    solverOpts: dict
    '''

    defaultOpts = {}

    defaultOpts_gurobi = {
        "Threads": 3,
        "OptimalityTol": 1e-8,
        "Method": 2,  # interior point/barrier
        "Cuts": 0,  # no precut of solution spce
        "NodeMethod": 2,  # interior points
        "IntFeasTol": 1e-9,  # small values in order to avoid errors with BigM
    }

    # note: scip needs a param file, which is located in the dir of he model - so, no direct options can be given
    defaultOpts_scip = {}

    defaultOpts_cbc = {"primalT": 1e-3}

    defaultOpts_glpk = {}

    defaultOpts_cplex = {"threads": 3}

    # append default options
    if solver == "gurobi":
        defaultOpts.update(defaultOpts_gurobi)
    elif solver == "scip":
        # note: scip needs a param file, which is located in the dir of he model - so, no direct options can be given
        defaultOpts.update(defaultOpts_scip)
    elif solver == "cplex":
        defaultOpts.update(defaultOpts_cplex)
    elif solver == "cbc":
        defaultOpts.update(defaultOpts_cbc)
    elif solver == "glpk":
        defaultOpts.update(defaultOpts_glpk)
        if "Threads" in solverOpts:
            solverOpts.pop("Threads")
        if "LogFile" in solverOpts:
            solverOpts.pop("LogFile")
    else:
        raise ValueError(
            'Solver name unknown. Please use one of "gurobi", "scip", "cbc", "glpk" or "cplex".'
        )

    # just add default options if not defined in solverOpts
    for option in defaultOpts:
        if not option in solverOpts:
            solverOpts[option] = defaultOpts[option]

    return solverOpts
