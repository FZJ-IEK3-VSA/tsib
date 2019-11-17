# -*- coding: utf-8 -*-
"""
Created on Fri Apr 08 11:33:01 2016

@author: Leander Kotzur
"""

import time
import os
import traceback
import logging
import warnings

import multiprocessing as mp

import pandas as pd
import numpy as np

from tsorb.utils.InputData import DataExchangeCsv
from tsorb.ElectricalLoadProfile import ElectricalLoadProfile
import tsib.data


def simSingleHousehold(residents, year, **elp_kwargs):
    """
    Function to build an Electrical Profile for the given number of residents
    in the given year and run it for exactly one household.
    ---------------------------------------------------------------------------
    Parameters:
        residents: int, required
            set the number of residents per household
        year: int, required
            set the calendar year the calculation will take place
        elp_kwargs: keyword arguments from parametrizing the 
            ElectricalLoadProfile class.
    ---------------------------------------------------------------------------
    Returns:
        time_series as numpy array with the total electric energy consumption
        for one year with a minutewise resolution
    """
    data_ex_main = DataExchangeCsv()
    elp = ElectricalLoadProfile(data_ex_main, residents, **elp_kwargs)

    result = elp.get_rescheduled_profiles(year)
    return result


def _run_household_year(queue, residents, year, work_seed, **elp_kwargs):
    """
    Parallel Worker Function to build an Electrical Profile for the given
    number of residents in the given year and run it. It is given a queue to
    save the result of different processes while multiprocessing.
    ---------------------------------------------------------------------------
    Parameters:
        q: Queue-Object from Multiprocessing lib, required
            Object to save data generated from processes run on parallel
            processors.
        residents: int, required
            set the number of residents per household
        year: int, required
            set the calendar year the calculation will take place        
        elp_kwargs: keyword arguments from parametrizing the 
            ElectricalLoadProfile class.
    ---------------------------------------------------------------------------
    Returns:
        Saved the time_series of one year inside the given queue
    """
    try:
        np.random.seed(work_seed)
        data_ex_main = DataExchangeCsv()
        elp = ElectricalLoadProfile(data_ex_main, residents, **elp_kwargs)
        result = elp.get_rescheduled_profiles(year)
        queue.put((None, result))
        return
    except Exception as e:
        queue.put((e, traceback.format_exc()))
        #        traceback.print_exc(file=sys.stdout)
        return


def simHouseholdsParallel(
    residents,
    year,
    no_of_households,
    singleProfiles=False,
    cores=mp.cpu_count() - 1,
    **elp_kwargs
):
    """
    Use to get a calculation of a chosen number of households with a full year
    electrical LoadProfile.
    Build n-1 (n = number of cores) parallel Processes, start them on different
    cores with the ElectricalLoadProfile function "run_for_year". If number of
    Households calculated is higher than the number of cores, the cores will be
    started repetetive. Wait for program to finish for correct closing of tasks
    on the different cores.
    ---------------------------------------------------------------------------
    Parameters:
        residents: int, required
            set the number of residents per household
        year: int, required
            set the calendar year the calculation will take place
        no_of_households: int, required
            set the number of households that will be calculated
        singleProfiles: bool, optional (default: False)
            If set, a list of single profiles is returned instead of the 
            aggregated ones.
        cores: int, optional (default: CPU-Count -1)
            Number of cores/threads used for parallel profile generation.
        elp_kwargs: keyword arguments from parametrizing the 
            ElectricalLoadProfile class.
    ---------------------------------------------------------------------------
    Returns:
        return the sum of the ElectricLoadProfiles results
    """
    # setup queue for outputs
    output = mp.Queue()
    # Number of runs (rounded) with usage of n-1 cores
    no_of_full_loops = int(no_of_households / cores)
    # last run with remainding calculations
    last_loop_calcs = int(no_of_households - (cores * no_of_full_loops))

    # initialize empty aggregated results
    if singleProfiles:
        agg_results = []
    else:
        agg_results = None

    masterseed = 7
    #    ws = [2]*no_of_households  # Test if workerseeds are all set the same all Households will be the same
    ws = [masterseed + i for i in range(no_of_households)]

    if last_loop_calcs > 0:
        no_loops = no_of_full_loops + 1
    else:
        no_loops = no_of_full_loops

    # Loop with parallel calculation for the number of households
    for runs in range(no_loops):
        # number of parallel processes
        if runs == no_of_full_loops:
            no_parallel = last_loop_calcs
        else:
            no_parallel = cores
        # Setup a list of processes that we want to run
        processes = [
            mp.Process(
                target=_run_household_year,
                args=(output, residents, year, ws[x + (cores * runs)]),
                kwargs=(elp_kwargs),
            )
            for x in range(no_parallel)
        ]
        # Run processes
        print(
            "Start "
            + str(no_parallel)
            + " processes, round "
            + str(runs + 1)
            + " of "
            + str(no_loops)
        )
        for i, p in enumerate(processes):
            print("  Process " + str(i + 1))
            p.start()
        # Get process results from the output queue
        result_list = [output.get() for p in processes]
        # check if the process failed
        for result in result_list:
            if isinstance(result[0], Exception):
                print(result[1])
                raise result[0]
            else:
                if singleProfiles:
                    agg_results.append(result[1])
                else:
                    if agg_results is None:
                        agg_results = result[1]
                    else:
                        agg_results += result[1]
        
        # Exit/Close the completed processes
        for p in processes:
            p.join()
        print("- closed the round")
        mp.active_children()

    # Transponate the Array to Columns = Number of Households
    return agg_results



def getHouseholdProfiles(
    n_persons,
    weather_data,
    weatherID,
    seeds=[0],
    ignore_weather=True,
    mean_load=True,
    cores=mp.cpu_count() - 1,
):
    """
    Gets or creates the relevant occupancy profiles for a building
    simulation or optimization.
    
    
    Parameters
    ----------
    n_persons: integer, required
        Number of persons living in a single appartment.
    weather_data: pd.DataFrame(), required
        A time indexed pandas dataframe containing weather data with 
        the GHI as a column.
    weatherID: str, required
        Giving an ID to the weather data to identify the resulting profile.
    seeds: list, optional (default: [0])
        List of integer seeds to create a number of profiles which have
        similar input parameters, but a varying output. Default, no seed is
        chosen. 
    ignore_weather: bool, optional (default: False)
        Since atm only the GHI is required for the electricity load profile,
        the weather plays a minor role and can be ignored by the identificaiton
        of profiles.
    mean_load: bool, optional (default: True)
        Decides if the created load profiles on 1-minute basis shall be 
        downsampled by taking the mean of 60 minutes or the first value in
        every 60 minutes.
    cores: int, optiona(default: mp.cpu_count() - 1)
        Number of cores used for profile generation.
    """

    # get the potential profile names
    filenames = {}
    for seed in seeds:
        profile_ID = "Profile" + "_occ" + str(int(n_persons)) + "_seed" + str(seed)
        if not ignore_weather:
            profile_ID = profile_ID + "_wea" + str(weatherID)

        if mean_load:
            profile_ID = profile_ID + "_mean"

        filenames[seed] = os.path.join(
            tsib.data.PATH, "results", "occupantprofiles", profile_ID + ".csv"
        )

    # check how many profiles do not exist#
    not_existing_profiles = {}
    for seed in seeds:
        if not os.path.isfile(filenames[seed]):
            not_existing_profiles[seed] = filenames[seed]

    # info about runtime
    if cores < 1:
        warnings.warn('Recognized cores are less than one. The code will behave as the number is one.')
        cores = 1

    _runtime = np.floor(float(len(not_existing_profiles))/cores)
    _log_str = str(len(not_existing_profiles)) + " household profiles need to get calculated. \n"
    _log_str += "With " + str(cores) + " threads, the estimated runtime is " + str(_runtime) + " minutes."
    logging.info(_log_str)

    # run in parallel all profiles
    if len(not_existing_profiles) > 1:
        new_profiles = simHouseholdsParallel(
            int(n_persons),
            2010,
            len(not_existing_profiles),
            singleProfiles=True,
            weather_data=weather_data,
            get_hot_water=True,
            resample_mean=mean_load,
            cores=cores,
        )
    # if single profile just create one profile and avoid multiprocessing
    elif len(not_existing_profiles) > 0:
        one_profile = simSingleHousehold(
            int(n_persons),
            2010,
            weather_data=weather_data,
            get_hot_water=True,
            resample_mean=mean_load,
        )
        new_profiles = [one_profile]

    # write results to csv files
    for i, seed in enumerate(not_existing_profiles):
        new_profiles[i].to_csv(not_existing_profiles[seed])

    # load all profiles
    profiles = []
    for seed in seeds:
        profile = pd.read_csv(filenames[seed], index_col=0)
        # TODO get a proper indexing in tsorb based on the weather data
        profile.index = weather_data.index

        profiles.append(profile)

    return profiles



if __name__ == "__main__":
    
    start_time = time.time()
    load_5_unique = simHouseholdsParallel(
        3, 2010, 5, singleProfiles=True
    )  # , resolved_load = True


    print("--- %s seconds ---" % (time.time() - start_time))
