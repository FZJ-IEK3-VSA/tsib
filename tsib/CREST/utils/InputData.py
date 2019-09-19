import os
import pandas as pd
import numpy as np

DATA_PATH = os.path.join(os.path.dirname(__file__),'..','data')

class DataExchange(object):
    def __init__(self):
        '''
        Init an Object to save data used in the different models, so the read
        functions from .csv file won't be called in the repetetively used code.
        '''
        # Data for the OccupancyModel states of user activity and change probs
        self._occ_start_states_data = {}
        self._occ_act_trans_data = {}
        # Data for the LightingModel
        self._bulbs = np.array([])
        self._bulbs_duration_config = np.array([])
        self._irradiance = np.array([])
        self._occ_sharing_factor = np.array([])
        self._bulbs_number = np.array([])
        self._bulbs_type = np.array([])
        self._bulbs_power = None
        # Data for the ApplianceModel
        self._app_data = np.array([])
        self._app_activity = {}
        self._rel_temperature_modifier = np.array([])
        self._activities = np.array([])
        self._four_state_start_state = {}
        self._four_state_trans_data = {}
        self._four_state_24_hour_occ = np.array([])
        self._app_pre_setup = np.array([])

    def create():
        return DataExchangeCsv()


# noinspection PyTypeChecker
class DataExchangeCsv(DataExchange):
    def __init__(self):
        self._path = DATA_PATH
        super(DataExchangeCsv, self).__init__()

    @property
    def get_occ_act_trans_data(self):
        if not bool(self._occ_act_trans_data):
            for residents in range(1, 6):
                self._occ_act_trans_data['wd' + str(residents)] \
                    = self.import_data('occ_act_wd_' + str(residents) + '.csv').values
                self._occ_act_trans_data['we' + str(residents)] \
                    = self.import_data('occ_act_we_' + str(residents) + '.csv').values
        return self._occ_act_trans_data

    @property
    def get_occ_start_states_data(self):
        if not bool(self._occ_start_states_data):
            self._occ_start_states_data['we'] = self.import_data('occ_start_states_weekend.csv').values
            self._occ_start_states_data['wd'] = self.import_data('occ_start_states_weekday.csv').values
        return self._occ_start_states_data

    @property
    def get_bulbs(self):
        if not self._bulbs.any():
            self._bulbs = self.import_data('bulbs.csv', header=None).values
        return self._bulbs
#
#    @property
#    def get_bulbs_number(self):
#        if not self._bulbs_number.any():
#            self._bulbs_number =  self.import_data('bulb_number.csv',header = None)
#        return  self._bulbs_number

    @property
    def get_bulbs_type(self):
        if not self._bulbs_type.any():
            self._bulbs_type = self.import_data('bulb_type.csv',header=None).values
        return  self._bulbs_type
    
    @property
    def get_bulbs_power(self):
        if self._bulbs_power is None:
            self._bulbs_power = self.import_data('bulbs_power.csv',header=0)
        return self._bulbs_power
    
    @property
    def get_bulbs_duration_config(self):
        if not self._bulbs_duration_config.any():
            self._bulbs_duration_config = self.import_data('bulbs_duration_config.csv').values
        return self._bulbs_duration_config

    @property
    def get_irradiance(self, directory=None):
        if not self._irradiance.any():
            self._irradiance = self.import_data('irradiance.csv').values
        return self._irradiance

    @property
    def get_occ_sharing_factor(self):
        if not self._occ_sharing_factor.any():
            self._occ_sharing_factor = self.import_data('occ_sharing_factor.csv', decimal = ',').values
        return self._occ_sharing_factor

    @property
    def get_app_data(self):
        if not self._app_data.any():
            self._app_data = self.import_data('appliances.csv',decimal = ',').values
        return self._app_data

    @property
    def get_app_pre_setup(self):
        self.app_pre_setup = self.import_data('app_pre_setup.csv').values
        return self.app_pre_setup
        
#    @property
#    def get_app_activity(self) -> dict:
#        if not bool(self._app_activity):
#            for residents in range(1, 6):
#                self._app_activity['wd' + str(residents)] \
#                    = self.import_data('app_act_wd_' + str(residents) + '.csv').values
#                self._app_activity['we' + str(residents)] \
#                    = self.import_data('app_act_we_' + str(residents) + '.csv').values
#        return self._app_activity

    @property
    def get_activities(self):
        if not self._activities.any():
            self._activities = self.import_data('activities.csv').values
        return self._activities

#    @property  # Used for Electric Heating purposes now not used anymore
#    def get_rel_temperature_modifier(self) -> np.array:
#        if not self._rel_temperature_modifier.any():
#            self._rel_temperature_modifier = self.import_data('rel_temperature_modifier.csv').values
#        return self._rel_temperature_modifier

    @property
    def get_four_state_trans_data(self):
        if not bool(self._four_state_trans_data):
            for residents in range(1, 7):
                self._four_state_trans_data['wd' + str(residents)] \
                    = self.import_data('four_state_wd_' + str(residents) + '.csv', skip_rows=9) \
                          .dropna(axis=1).values
                self._four_state_trans_data['we' + str(residents)] \
                    = self.import_data('four_state_we_' + str(residents) + '.csv', skip_rows=9) \
                          .dropna(axis=1).values
        return self._four_state_trans_data

    @property
    def get_four_state_start_state(self):
        if not bool(self._four_state_start_state):
            self._four_state_start_state['wd'] = \
                self.import_data('four_state_starting_state_weekday.csv', skip_rows=2).dropna(axis=1).values[:,1:]
            self._four_state_start_state['we'] = \
                self.import_data('four_state_starting_state_weekend.csv', skip_rows=2).dropna(axis=1).values[:,1:]
        return self._four_state_start_state

    @property
    def get_four_state_24_hour_occ(self):
        if not self._four_state_24_hour_occ.any():
            self._four_state_24_hour_occ = self.import_data('four_state_24_hour_occupancy.csv', skip_rows=2) \
                                               .dropna(axis=1).values
        return self._four_state_24_hour_occ

    def import_data(self, filename, header=0, sep=';', skip_rows=0,
                    decimal = '.'):
        try:
            if filename == '':
                raise ValueError('Import from CSV-file: Filename is empty')
            combined_path = os.path.normpath(os.path.join(self._path, filename))
            return pd.read_csv(filepath_or_buffer=combined_path, header=header, 
                               sep=sep, skiprows=skip_rows, decimal = decimal)
        except ValueError as e:
            print(e)

if __name__ == "__main__":
    data_ex_main = DataExchangeCsv()
    #a = data_ex_main.get_four_state_trans_data
#    a=data_ex_main.get_bulbs_number
#    b=data_ex_main.get_bulbs_type
#    c=data_ex_main.get_bulbs_power
#    print(a)
#    print(b)
#    print(c)