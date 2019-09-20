import pandas as pd
import numpy as np
import math
from tsib.electric.utils.InputData import DataExchangeCsv

# from tsib.electric.utils.ExportData import export_data
from tsib.electric.utils.DPD import DPD


class OccupancyModel(object):
    def __init__(self, data_ex, residents, four_states=True):
        """
        Initialize OccupancyModel Object
        -----------------------------------------------------------------------
        Parameters:
            data_ex: DataExchange Object, required
                DataExchange Object to exchange data between the Models and 
                save runspecific Data
            residents: int, required
                Number of residents living in the household
            four_states: bool, optional
                Decides if the simulation should be run for four state 
                data
        -----------------------------------------------------------------------
        """
        self._data_ex = data_ex
        self._residents = residents
        self._four_states = four_states
        self._has_run = False
        # Preload the matrices with the activity propabilities of the residents
        if self._four_states:
            self._start_state = self._data_ex.get_four_state_start_state
            self._occ_act_trans_data = self._data_ex.get_four_state_trans_data
        else:
            self._start_state = self._data_ex.get_occ_start_states_data
            self._occ_act_trans_data = self._data_ex.get_occ_act_trans_data
        # Initialize an empty numpy array for the resulting activity
        self.occ_activity = np.zeros(144, int)

    @property
    def get_residents(self):
        return self._residents

    @property
    def get_occ_no_activity(self):
        """
        Returns the number of occupants when they are at home, but not active.
        """
        if not self._has_run:
            self.run()
        if not self._four_states:
            raise ValueError(
                "Just a four state occupancy model can return"
                + "these states. Please set four_states parameter to True."
            )
        return self.occ_no_activity

    @property
    def get_occ_activity(self):
        """
        Return the number of active occupants at home a 10-minutes interval
        """
        if not self._has_run:
            self.run()
        return self.occ_activity

    def run(self, day_of_week):
        """
        Runs the Object of OccupancyModel
        -----------------------------------------------------------------------
        Parameters:
            day_of_week: str, required
                String either it is a weekday "wd" or weekend "we"
        -----------------------------------------------------------------------
        Returns:
            Calculates the residents activity for one day with a resolution of
            10 minutes. The result is saved in a numpy array self.occ_activity
            Then set the variable _has_run to True
        -----------------------------------------------------------------------
        """
        # reset occ_activity for the purpose of full year calculations
        states = np.zeros(144, int)
        # Get the Data for the given number of residents and type of day
        start_state = self._start_state[day_of_week][:, self._residents - 1]
        act_transition = self._occ_act_trans_data[day_of_week + str(self._residents)]
        dpd = DPD(start_state)

        # get start states and index jumper per interval
        if self._four_states:
            num_row_jump = int((self._residents + 1) ** 2)
            start_ix = dpd.get_random_interval()
            # state index for 6 person household to index of no. residents
            state = int(
                (
                    (self._residents + 1) * math.floor(start_ix / 7)
                    + int(math.fmod(start_ix, 7))
                )
            )
        else:
            num_row_jump = 7
            state = dpd.get_random_interval()

        # in case of multiple calls, get the state from the last run as start state
        if self._has_run:
            state = self._states[-1]
        states[0] = state

        # loop over whole day
        for interval in range(0, 144):
            row = act_transition[num_row_jump * interval + state, 2:]
            dpd = DPD(row)
            state = dpd.get_random_interval()
            states[interval] = state

        # save results
        if self._four_states:
            home = np.floor(states / (self._residents + 1)).astype(int)
            active = np.fmod(states, self._residents + 1)
            self.occ_activity = np.minimum(home, active)
            self.occ_no_activity = np.maximum(home - active, 0)
        else:
            self.occ_activity = states

        self._states = states
        #        self.export_occ_activity()
        self._has_run = True

    # def export_occ_activity(self):
    #     '''
    #     Export the occ_activity as a .csv file in the Data Folder
    #     '''
    #     data = self.occ_activity
    #     index = pd.date_range(start='00:00', end='23:50', freq='10min').time
    #     df = pd.DataFrame(data=data, index=index)
    #     export_data('occ_activity.csv', df)

    @staticmethod
    def set_seed(seed):
        """
        Sets the seed for the random number generation
        :param seed:
        """
        np.random.seed(seed)


if __name__ == "__main__":
    data_ex_main = DataExchangeCsv()
    occ_model = OccupancyModel(data_ex_main, 1, True)
    occ_model.set_seed(1)
    occ_model.run("wd")
    # occ_model.export_occ_activity()
    # occ_model.run_n_times_and_plot(10)
#    occ_model.average_activity()2
# for i in range(1,6):
#     occ_model = OccupancyModel(data_ex_main, i, 'wd')
#     occ_model.run_n_times_and_plot(10000)
# for i in range(1,6):
#     occ_model = OccupancyModel(data_ex_main, i, 'wd')
#     occ_model.run_n_times_and_plot(10000)
