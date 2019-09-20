#!/usr/bin/env python
''' Test Script only
.. module:: utils
    :platform Unix, Windows
    :version $Id: utils.py 68 2014-03-02 20:31:16Z leko $
.. codeauthor:: Leander Kotzur
.. codeauthor:: Roland Doll

'''
import os
import sys
import numpy as np
import pandas as pd
import datetime

        
class PowerDict(dict): #Power Dictionary with additional functions
    def __init__(self, parent = None, key = None):
        self.parent = parent
        self.key = key
    def __missing__(self, key): #creation of subdictionaries on fly
        self[key] = PowerDict(self, key)
        return self[key]
    def append(self, item): #additional append function for lists in dict
        self.parent[self.key] = [item]
    def __setitem__(self, key, val):
        dict.__setitem__(self, key, val)
        try:
            val.parent = self
            val.key = key
        except AttributeError:
            pass



         
def calcFixCost(WACC,optInterval,Lifetime,CAPEX,OPEXfix):
    '''
    Annualizes the capital expenditure CAPEX [eur or eur/kW] and the fix operation expenditure OPEXfix [-/a] as relative value of the capital expenditure. 
    This is done for a horizon or lifetime [years] with an interest rate or WACC [-/a].
    The annualized cost are then linearly amortized to the optInterval [hours]
    '''
    raise Warning('Function "calcFixCost" is deprecated')
    aFactor = ((1+WACC)**Lifetime)*WACC/(((1+WACC)**Lifetime)-1) if WACC > 0 else 1./Lifetime;
    if CAPEX>0:
        annCosts=CAPEX*(aFactor+OPEXfix);
    else:
        annCosts = OPEXfix
    return annCosts*optInterval/8760;      

def readProfile(filePath):
    '''
    Reads a profile from a .txt file, where each profile entry is in a new line and returns it as a list
    '''
    data=[]
    f = open(os.path.dirname(os.path.realpath(__file__)) + filePath)
    for line in f.readlines():
        data.append(float(line.replace(',','.')))
    f.close()    
    return data    


def writeProfile(data,filePath):
    '''
    writes a profile t a .txt file, where each profile entry is in a new line and returns it as a list
    '''
    f = open(os.path.dirname(os.path.realpath(__file__)) + filePath, 'w')
    for datapoint in data:
        f.write(str(datapoint)+'\n')
    f.close()
    return


def profilesToHourlyProfiles(sourceFolderPath,aimFolderPath):
    files=os.listdir(os.path.dirname(os.path.realpath(__file__)) + sourceFolderPath)
    for file in files:
        if file[-4:]=='.txt':
            data=readProfile(sourceFolderPath + '\\' + file)
            if len(data)>8760 and len(data)%8760==0: # more timesteps then a year + ganzzahliges vielfaches
                hourlydata=data[0::len(data)/8760]
            # normalize to average value of 1
            sum_val=sum(hourlydata)
            final_data=[cor_val/sum_val*8760 for cor_val in hourlydata ]
            writeProfile(final_data,aimFolderPath+'\\'+file[:-4]+'_hourly.txt')
    return



def getKwarg(kwargs,name,default,obligatory):
    'OPEN'
    return

def isString(var):
    '''
    Python 2 and 3 compatible string testing.
    '''
    PY3 = sys.version_info[0] == 3
    if PY3:
        string_type = str,
    else:
        string_type = basestring,
    return isinstance(var,string_type)
        
        
def LCOE(WACC,CAPEX,OPEX_fix,OPEX_var,Lifetime,fullloadhours):
    aFactor = ((1+WACC)**Lifetime)*WACC/(((1+WACC)**Lifetime)-1) if WACC > 0 else 1./Lifetime;
    return (CAPEX*aFactor+OPEX_fix)/fullloadhours
    
def write2gif(orgi_path_list, aim_path, duration = 5.0):
    import imageio
    frames = [imageio.imread(path) for path in orgi_path_list]
#    print('read')
    kargs = { 'duration': duration }
    imageio.mimsave(aim_path, frames, 'GIF', **kargs)
    return
#    with imageio.get_writer(PACKAGE_PATH +'/__enersysopt_cache__/movie.gif', mode='I') as writer:
#        for filename in orgi_path_list:
#            image = imageio.imread(filename)
#            writer.append_data(image)
    
def write2gif_alt(orgi_path_list, aim_path, duration = 5.0):
    from .images2gif import writeGif
    from PIL import Image
    frames = [Image.open(path) for path in orgi_path_list]
#    print('read')
    kargs = { 'duration': duration }
    writeGif(aim_path, frames,  **kargs)#'GIF',
    return
    
    

def seasonValue(dates):
    '''
    Recieves a pd.DateIndex array and calculates the "seasonal" value of each
    entry from 0 (Winter 21.12) to 1 (Summer 21.06).
    '''
    years = np.unique(dates.year)
    years = np.insert(years,0,years[0]-1)
    years = np.append(years,years[-1]+1)
    suwiDiff = datetime.datetime(years[0],6,21) - datetime.datetime(years[0],12,21)
    suwiDiff = suwiDiff.total_seconds()
    diffs = pd.DataFrame(index = dates)
    for year in years:     
        wiDiff = dates - datetime.datetime(year,12,21)
        diffs[str(year) + 'w'] = abs(wiDiff.total_seconds() / suwiDiff)
    seasonVals = diffs.min(axis=1)
    return seasonVals.values
    
    
def getElectroChemicalEff(voltage, current_density, elecMode = True):
    '''
    Calculates the efficiency of an electrolysis or fuel cell derived from
    the voltage (V) and current_density (A/cm^2)
    '''
    LHV_H2 = 119972 # J/g
    molmass_H2 = 2.01588  #g/mol
    faraday_eff = 0.98
    hydrogen_rate = current_density * faraday_eff / 96485.332 / 2.0
    LHV_rate = hydrogen_rate * molmass_H2 * LHV_H2
    if elecMode:
        efficiency = (LHV_rate * faraday_eff)/ (voltage * current_density)
    else:
        efficiency = (voltage * current_density)/ (LHV_rate / faraday_eff)
    return efficiency

def manageSolverOpts(solver, solverOpts):

    defaultOpts = {}

    defaultOpts_gurobi = {
        'Threads': 3,
        'OptimalityTol': 1e-8,
        'Method': 2,  # interior point/barrier
        'Cuts': 0,  # no precut of solution spce
        'NodeMethod': 2,  # interior points
        'IntFeasTol': 1e-9,  # small values in order to avoid errors with BigM
        }

    # note: scip needs a param file, which is located in the dir of he model - so, no direct options can be given
    defaultOpts_scip = {}

    defaultOpts_cbc = {'primalT': 1e-3}

    defaultOpts_glpk = {}

    defaultOpts_cplex = {'threads': 3,
                        }

    # append default options
    if solver == 'gurobi':
        defaultOpts.update(defaultOpts_gurobi)
    elif solver == 'scip':
        # note: scip needs a param file, which is located in the dir of he model - so, no direct options can be given
        defaultOpts.update(defaultOpts_scip)
    elif solver == 'cplex':
        defaultOpts.update(defaultOpts_cplex)
    elif solver == 'cbc':
        defaultOpts.update(defaultOpts_cbc)
    elif solver == 'glpk':
        defaultOpts.update(defaultOpts_glpk)
        if 'Threads' in solverOpts:
            solverOpts.pop('Threads')
        if 'LogFile' in solverOpts:
            solverOpts.pop('LogFile')
    else:
        raise ValueError('Solver name unknown. Please use one of "gurobi", "scip", "cbc", "glpk" or "cplex".')

    for option in defaultOpts:
        if not option in solverOpts:
            solverOpts[option] = defaultOpts[option]

    return solverOpts


def sumSubComponents(resultManager, components, onlyNodes=False):
    '''
    Sums the results of technologies which are initialized several times
    (Photovoltaic 1 and 2 etc.) to a single result and cleans the heat pump capacity
    '''


    for node in components:
        # avoid capacity error for multiple heat pump intialization - set capacity of virtual heat pumps to zero
        if node == 'Heat pump':
            for column in resultManager.nodes.columns:
                if 'Heat pump' in column and not 'Heat pump' == column:
                    if isinstance(resultManager.nodes.index, pd.MultiIndex):
                        resultManager.nodes.loc[pd.IndexSlice[:, 'Capacity'], column] = 0
                    else:
                        resultManager.nodes.loc['Capacity', column] = 0


        if not onlyNodes:
            for column in resultManager.flows.columns:
                FromTo = column.split(' to ')
                if node in FromTo[0] and not node == FromTo[0]:
                    if not node + ' to ' + FromTo[1] in resultManager.flows.columns:
                        resultManager.flows[node + ' to ' + FromTo[1]] = 0
                    resultManager.flows.loc[:, node + ' to ' + FromTo[1]] += resultManager.flows.loc[:, column]
                    resultManager.flows.drop([column], axis=1, inplace=True)
                if node in FromTo[1] and not node == FromTo[1]:
                    if not FromTo[0] + ' to ' + node in resultManager.flows.columns:
                        resultManager.flows[FromTo[0] + ' to ' + node] = 0
                    resultManager.flows.loc[:, FromTo[0] + ' to ' + node] += resultManager.flows.loc[:, column]
                    resultManager.flows.drop([column], axis=1, inplace=True)

            for column in resultManager.loads.columns:
                if node in column and not node == column:
                    if not node in resultManager.loads.columns:
                        resultManager.loads[node] = 0
                    resultManager.loads.loc[:, node] += resultManager.loads.loc[:, column]
                    resultManager.loads.drop([column], axis=1, inplace=True)

        for column in resultManager.nodes.columns:
            if node in column and not node == column:
                if not node in resultManager.nodes.columns:
                    resultManager.nodes[node] = 0
                resultManager.nodes.loc[:, node] += resultManager.nodes.loc[:, column]
                resultManager.nodes.drop([column], axis=1, inplace=True)

    return resultManager


def dropComponents(resultManager, components, onlyNodes=False):
    '''
    Sums the results of technologies which are initialized several times
    (Photovoltaic 1 and 2 etc.) to a single result
    '''

    for node in components:
        if not onlyNodes:
            for column in resultManager.flows.columns:
                FromTo = column.split(' to ')
                if node == FromTo[0] or node == FromTo[1]:
                    resultManager.flows.drop([column], axis=1, inplace=True)
            for column in resultManager.loads.columns:
                if node == column:
                    resultManager.loads.drop([column], axis=1, inplace=True)
        for column in resultManager.nodes.columns:
            if node == column:
                resultManager.nodes.drop([column], axis=1, inplace=True)

    return resultManager

def renameComponent(resultManager, oldName, newName, onlyNodes = False):
    '''
    Renames a component in the result manager
    '''

    resultManager.nodes.rename(columns = {oldName:newName}, inplace=True)
    if not onlyNodes:
        resultManager.loads.rename(columns = {oldName:newName}, inplace=True)
        for column in resultManager.flows.columns:
            FromTo = column.split(' to ')
            if oldName == FromTo[0]:
                resultManager.loads.rename(columns = {str(column): str(column).replace(oldName,newName)},
                                           inplace = True)
    return resultManager
    
def _aggregateColumns(df, column_names):
    '''
    Sums all columns which columns names are containt by the invidividual entries
    of column_names.

    Parameters
    ----------
    df: pandas.DataFrame()
    column_names: list

    Returns
    -------
    pandas.DataFrame()
    '''
    for node in column_names:
        for column in df.columns:
            if node in column and not node == column:
                if not node in df.columns:
                    df[node] = 0
                df.loc[:, node] += df.loc[:, column]
                df.drop([column], axis=1, inplace=True)
    return df



def getGridExchange(resultManager):
    '''
    Returns the grid exchange from the flow data of the result manager 
    '''
#    
#    SUB_COMPS = ['Heat pump', 'Heat storage', 'HNode', 'Photovoltaic', 'Control',
#                                 'Solar thermal','Walls','Windows','Roof']
#    
#    # reduce the input data set
#    resultManager = sumSubComponents(resultManager, SUB_COMPS)
    
    gridExchange = pd.DataFrame()
    gridExchange = -resultManager.flows[['Photovoltaic 1 to FiTPV',
                                         'Photovoltaic 2 to FiTPV',
                                   'CHP to FiTCHP',
                                   'Fuel cell to FiTCHP']].sum(axis=1)
    gridExchange = (gridExchange +resultManager.flows['Electricity supply to AC Node']
                    + resultManager.loads['HP Tarif'])
    gridExchange = gridExchange.unstack(level=0)
    return gridExchange

import matplotlib.colors as mcolors
def make_colormap(seq):
    """
    Return a LinearSegmentedColormap
    seq: a sequence of floats and RGB-tuples. The floats should be increasing
    and in the interval (0,1).
    """
    seq = [(None,) * 3, 0.0] + list(seq) + [1.0, (None,) * 3]
    cdict = {'red': [], 'green': [], 'blue': []}
    for i, item in enumerate(seq):
        if isinstance(item, float):
            r1, g1, b1 = seq[i - 1]
            r2, g2, b2 = seq[i + 1]
            cdict['red'].append([item, r1, r2])
            cdict['green'].append([item, g1, g2])
            cdict['blue'].append([item, b1, b2])
    return mcolors.LinearSegmentedColormap('CustomMap', cdict)


OWN_CMAP = make_colormap([ np.array([2., 61., 100.])/255*0.8,
                     np.array([2., 61., 100.])/255,0.1,
                     np.array([2., 61., 100.])/255,
                     np.array([2., 61., 100.])/255*1.3,0.25,
                     np.array([2., 61., 100.])/255*1.3,
                     np.array([ 30., 150., 81.])/255*1.3,0.7,
                     np.array([ 30., 150., 81.])/255*1.3,
                     np.array([250., 235., 90.])/255, ])
    
OWN_CMAP_2 = make_colormap([np.array([2., 61., 100.])/255,
                         np.array([2., 61., 100.])/255*1.3,0.25,
                         np.array([2., 61., 100.])/255*1.3,
                         np.array([ 30., 150., 81.])/255*1.3,0.7,
                         np.array([ 30., 150., 81.])/255*1.3,
                         np.array([250., 235., 90.])/255, ])
        
CMAP_DISS = make_colormap([np.array([2., 61., 100.])/255,
                     np.array([20., 129., 129.])/255,0.25,
                     np.array([20., 129., 129.])/255,
                     np.array([108., 139., 70.])/255,0.5,
                     np.array([108., 139., 70.])/255,
                     np.array([250.,150.,90.])/255,0.75,
                     np.array([250.,150.,90.])/255,
                     np.array([255.,192.,0.])/255, ])


if __name__ == '__main__':
    print(getElectroChemicalEff(1.1,0.5, elecMode = True))