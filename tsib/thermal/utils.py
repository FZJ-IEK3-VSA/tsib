#!/usr/bin/env python
""" Test Script only
.. module:: utils
    :platform Unix, Windows
    :version $Id: utils.py 68 2014-03-02 20:31:16Z leko $
.. codeauthor:: Leander Kotzur
.. codeauthor:: Roland Doll

"""
import os
import sys
import numpy as np
import pandas as pd
import datetime


class PowerDict(dict):  # Power Dictionary with additional functions
    def __init__(self, parent=None, key=None):
        self.parent = parent
        self.key = key

    def __missing__(self, key):  # creation of subdictionaries on fly
        self[key] = PowerDict(self, key)
        return self[key]

    def append(self, item):  # additional append function for lists in dict
        self.parent[self.key] = [item]

    def __setitem__(self, key, val):
        dict.__setitem__(self, key, val)
        try:
            val.parent = self
            val.key = key
        except AttributeError:
            pass


def calcFixCost(WACC, optInterval, Lifetime, CAPEX, OPEXfix):
    """
    Annualizes the capital expenditure CAPEX [eur or eur/kW] and the fix operation expenditure OPEXfix [-/a] as relative value of the capital expenditure. 
    This is done for a horizon or lifetime [years] with an interest rate or WACC [-/a].
    The annualized cost are then linearly amortized to the optInterval [hours]
    """
    raise Warning('Function "calcFixCost" is deprecated')
    aFactor = (
        ((1 + WACC) ** Lifetime) * WACC / (((1 + WACC) ** Lifetime) - 1)
        if WACC > 0
        else 1.0 / Lifetime
    )
    if CAPEX > 0:
        annCosts = CAPEX * (aFactor + OPEXfix)
    else:
        annCosts = OPEXfix
    return annCosts * optInterval / 8760


def readProfile(filePath):
    """
    Reads a profile from a .txt file, where each profile entry is in a new line and returns it as a list
    """
    data = []
    f = open(os.path.dirname(os.path.realpath(__file__)) + filePath)
    for line in f.readlines():
        data.append(float(line.replace(",", ".")))
    f.close()
    return data


def writeProfile(data, filePath):
    """
    writes a profile t a .txt file, where each profile entry is in a new line and returns it as a list
    """
    f = open(os.path.dirname(os.path.realpath(__file__)) + filePath, "w")
    for datapoint in data:
        f.write(str(datapoint) + "\n")
    f.close()
    return


def profilesToHourlyProfiles(sourceFolderPath, aimFolderPath):
    files = os.listdir(os.path.dirname(os.path.realpath(__file__)) + sourceFolderPath)
    for file in files:
        if file[-4:] == ".txt":
            data = readProfile(sourceFolderPath + "\\" + file)
            if (
                len(data) > 8760 and len(data) % 8760 == 0
            ):  # more timesteps then a year + ganzzahliges vielfaches
                hourlydata = data[0 :: len(data) / 8760]
            # normalize to average value of 1
            sum_val = sum(hourlydata)
            final_data = [cor_val / sum_val * 8760 for cor_val in hourlydata]
            writeProfile(final_data, aimFolderPath + "\\" + file[:-4] + "_hourly.txt")
    return




def manageSolverOpts(solver, solverOpts):

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

    for option in defaultOpts:
        if not option in solverOpts:
            solverOpts[option] = defaultOpts[option]

    return solverOpts
