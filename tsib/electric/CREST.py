import multiprocessing as mp
from tsib.electric.utils.InputData import DataExchangeCsv
from tsib.electric.ElectricalLoadProfile import ElectricalLoadProfile
import time
import pandas as pd
import numpy as np
import traceback


def one_household_year(residents, year, **elp_kwargs):
    '''
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
    '''
    data_ex_main = DataExchangeCsv()
    elp = ElectricalLoadProfile(data_ex_main, residents,
                                **elp_kwargs)
        
    result = elp.get_rescheduled_profiles(year)
    return result
    
def _run_household_year(queue, residents, year, work_seed, **elp_kwargs):
    '''
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
    '''
    try: 
        np.random.seed(work_seed)
        data_ex_main = DataExchangeCsv()
        elp = ElectricalLoadProfile(data_ex_main, residents, **elp_kwargs)  
        result = elp.get_rescheduled_profiles(year)
        queue.put((None,result))
        return
    except Exception as e:
        queue.put((e,traceback.format_exc()))
#        traceback.print_exc(file=sys.stdout)
        return


def run_district_year(residents, year, no_of_households, 
                      singleProfiles = False, cores = mp.cpu_count() - 1,
                      **elp_kwargs):
    '''
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
        elp_kwargs: keyword arguments from parametrizing the 
            ElectricalLoadProfile class.
        singleProfiles: bool, optional (default: False)
            If set, a list of single profiles is returned instead of the 
            aggregated ones.
        cores: int, optional (default: CPU-Count -1)
            Number of cores/threads used for parallel profile generation.
    ---------------------------------------------------------------------------
    Returns:
        return the sum of the ElectricLoadProfiles results
    '''
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

    if last_loop_calcs >0:
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
        processes = [mp.Process(target=_run_household_year, 
                                  args=(output, residents, year, 
                                        ws[x+(cores*runs)], ),
                                  kwargs = (elp_kwargs),) for x in range(no_parallel)]
        # Run processes
        print('Start '+str(no_parallel)+ ' processes, round ' + str(runs+1) + ' of ' 
                + str(no_loops))
        for i,p in enumerate(processes):
            print('  Process ' + str(i+1))
            p.start()
        # Get process results from the output queue
        result_list = [output.get() for p in processes]
        # check if the process failed
        for result in result_list:
            if isinstance(result[0],Exception):
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
#        # otherwise create list of all results
#        results = [res[1] for res in result_list]
#        results_all_temp1 = pd.DataFrame.from_records(results)
#        df = pd.concat([df, results_all_temp1])
        # Exit/Close the completed processes
        for p in processes:
            p.join()
        print('- closed the round')
        mp.active_children()

    # Transponate the Array to Columns = Number of Households
    return agg_results
#    results_sum = results_all.sum(axis=1)
#
#
#    leap_year = (year % 4)
#    if leap_year == 0:
#        index = pd.date_range(start=pd.datetime(year, 1, 1), periods=1440 * 366, freq='1min')
#    elif leap_year != 0:
#        index = pd.date_range(start=pd.datetime(year, 1, 1), periods=1440 * 365, freq='1min')
#
#    results_sum.index = index
#    results_all.index = index
#    results_sum_hour = results_sum.resample('60min', label='left').mean()
#    # check whether all output was extracted from queue
#    if output.empty():
#        print('All Outputs retrieved from Queue')
#    print('start Export to csv')
##    export_data('time_series_'+ str(no_of_households) + '_households' + str(year) + '.csv', pd.DataFrame(data=results_all))
#    export_data('time_series_' + str(no_of_households) + '_sum_houses_' +
#                str(year) + '_occ_' + str(residents) + '.csv', pd.DataFrame(data=results_sum))
#    export_data('time_series_' + str(no_of_households) + '_sum_houses_hourly_' +
#                str(year) + '_occ_' + str(residents) + '.csv', pd.DataFrame(data=results_sum_hour))
#    print('Finished Export to .csv')


if __name__ == "__main__":
    start_time = time.time()
    load_5_unique = run_district_year(3, 2010, 5, singleProfiles = True) # , resolved_load = True
    #one_household_year(2, 2015, irradiance_dir='TRY2010_06_Jahr.dat', TRY_data=True)
   # _run_household_year(True,4, 2015, 1, irradiance_dir='TRY2010_06_Jahr.dat', TRY_data=True)
    print("--- %s seconds ---" % (time.time() - start_time))