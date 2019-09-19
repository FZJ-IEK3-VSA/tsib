from __future__ import division
from tsib.CREST.Load import Load
from tsib.CREST.utils.Calculations import Calculations
import numpy as np
import math


class Appliance(Load):
    def __init__(self, key, long_name, type_name, act_use_profile,
                 ownership_probability,
                 standby_power, mean_cycles_power, cycles_per_year, mean_cycle_length,
                 restart_delay, calibration, total_energy, active_occ_dependent,
                 average_activity_profile, mean_power_factor, heat_gain):
        '''
        Object of class Load used to store the data connected to the
        applications in the household. 
        '''
        """
        :param key: Upper case name
        :param long_name: complete name
        :param type_name: appliance type name
        :param act_use_profile: describing the use profile
        :param ownership_probability: ownership probability
        :param standby_power: standby power demand in Watts
        :param mean_cycles_power: mean power demand of a cycle in Watts
        :param cycles_per_year: average cycle of cycles per year
        :param mean_cycle_length: length of an average cycle in minutes
        :param restart_delay: delay between cycles in minutes
        :param calibration: calibration constant
        <link>http://www.sciencedirect.com/science/article/pii/S0378778802002414</link>
        :param total_energy:
        :param active_occ_dependent:
        :param average_activity_profile
        :param mean_power_factor
        :param heat_gain # ratio how much of the light is transformed to heat
        :return:
        """
        super(Appliance, self).__init__(key.upper())
        # input set from appliances.csv
        self.long_name = long_name
        self.type_name = type_name
        # use_profiles = LEVEL, ACTIVE_OCC, ACT_TV, ACT_COOKING, ACT_LAUNDRY,
        # ACT_WASHDRESS, ACT_IRON, ACT_HOUSECLEAN, CUSTOM (Storage Heater)
        self.act_use_profile = act_use_profile
        self._ownership_probability = ownership_probability
        self._standby_power = standby_power
        self._mean_cycles_power = mean_cycles_power
        self._cycles_per_year = cycles_per_year
        self._mean_cycle_length = mean_cycle_length
        self._restart_delay = restart_delay
        self._calibration = calibration
        # not needed input from appliances.csv
        self._total_energy = total_energy
        self._active_occ_dependent = active_occ_dependent
        self._average_activity_profile = average_activity_profile
        self._mean_power_factor = mean_power_factor
        self.heat_gain = heat_gain
        # set inside the class
        self.owned = False
        self._power = float
        self._cycle_time_left = 0
        self._restart_delay_time_left = self._set_restart_delay()
        self._rated_power = self._set_rated_power()
        self.switched_on_time_series = np.zeros(1440)
        self.consumption = np.zeros(1440)
        self.heat_prod = np.zeros(1440)

    def __str__(self):
        return 'Key: {0}, activity use profile: {1}, ownership: {2}, standby power: {3}, mean cycles power: {4},' \
               'base cycles/y: {5}, mean cycles length (min): {6}, delay restart (min): {7}, calibration scalar: {8}' \
            .format(self._key, self.act_use_profile, self._ownership_probability, self._standby_power,
                    self._mean_cycles_power, self._cycles_per_year, self._mean_cycle_length, self._restart_delay,
                    self._calibration)

    def _set_restart_delay(self):
        """
        Assumes that the true restart delay is evenly distributed about the restart delay
        """
        return int(2 * np.random.random() * self._restart_delay)

    def _set_rated_power(self):
        """
        Assumes that the true rated_power can be calculated from a normal distribution
        """
        return Calculations.calc_from_normal_distr(self._mean_cycles_power, self._mean_cycles_power / 10)

    @staticmethod
    def _is_between(x, lower, upper):
        """
        Convenience function to test if a value is within specified limits.
        :param x: value to test
        :param lower: lower limit
        :param upper: upper limit
        :return: True if x is between lower and upper
        """
        return lower <= x <= upper

    def start(self, minute):
        # Determine how long this appliance is going to be on for
        self._cycle_time_left = self._calc_cycle_length()
        # Determine if this appliance has a delay after the cycle before it can restart
        self._restart_delay_time_left = self._restart_delay
        # Set the power
        self._power = self._get_power_usage(self._cycle_time_left)
        self.switched_on_time_series[minute] = 2
        # Decrement the cycle time left
        self._cycle_time_left -= 1

    def keep_running(self, minute):
        # Set the power
        self._power = self._get_power_usage(self._cycle_time_left)
        self.switched_on_time_series[minute] = 2
        # Decrement the cycle time left
        self._cycle_time_left -= 1

    def _calc_cycle_length(self):
        length = self._mean_cycle_length
        # Use the TV watching length data approximation, derived from the TUS data
        if self._key == 'TV1' or self._key == 'TV2' or self._key == 'TV3':
            # The cycle length is approximated by the following function
            # The average viewing time is approximately 73 minutes
            length = 70 * math.pow((-np.log(1 - np.random.random())), 1.1)
        return int(length)

    def _get_power_usage(self, cycle_time_left):
        """
        Gets the power consumption of this appliance in a time interval
        :param cycle_time_left: minutes left in the running cycle of the appliance
        :return: power consumption of this appliance in this time interval
        """
        # Set the working power to the rated power
        power = self._rated_power
        # Some appliances have a custom (variable) power profile depending
        # on the time left
        if self._key == 'WASHING_MACHINE' or self._key == 'WASHER_DRYER':
            total_cycle_time = 0

            if self._key == 'WASHING_MACHINE':
                total_cycle_time = 138
            if self._key == 'WASHER_DRYER':
                total_cycle_time = 198

            # This is an example power profile for an example washing
            # machine. This simplistic model is based upon data from personal
            # communication with a major washing machine manufacturer
            # todo check if the one is really needed
            time = total_cycle_time - cycle_time_left + 1
            if self._is_between(time, 1, 8):
                power = 73  # start-up and fill
            elif self._is_between(time, 9, 29):
                power = 2056  # heating
            elif self._is_between(time, 30, 81):
                power = 73  # wash and drain
            elif self._is_between(time, 82, 92):
                power = 73  # spin
            elif self._is_between(time, 93, 94):
                power = 250  # rinse
            elif self._is_between(time, 95, 105):
                power = 73  # spin
            elif self._is_between(time, 106, 107):
                power = 205  # rinse
            elif self._is_between(time, 108, 118):
                power = 73  # spin
            elif self._is_between(time, 119, 120):
                power = 250  # rinse
            elif self._is_between(time, 121, 131):
                power = 73  # spin
            elif self._is_between(time, 132, 133):
                power = 250  # rinse
            elif self._is_between(time, 134, 138):
                power = 568  # fast spin
            elif self._is_between(time, 139, 198):
                power = 2500  # drying cycle
            else:
                power = self._standby_power
        return power

    def set_ownership(self, random_app_per_run=True):
#        np.random.seed(5)       # Test fuer fixe App Ausstattung
        if random_app_per_run is True:
#            np.random.seed()         # Test fuer fixe App Ausstattung
            self.owned = np.random.random() < self._ownership_probability
        else:
            # set fixed seed for app-ownership (debug purposes)
            np.random.seed(random_app_per_run)
            self.owned = np.random.random() < self._ownership_probability
            # set seed back to another random seed
            np.random.seed()  # for test debug commented out

    def is_owned(self):
        return self.owned

    def is_off(self):
        return self._cycle_time_left <= 0

    def waiting_for_restart(self):
        return self._restart_delay_time_left > 0



if __name__ == "__main__":
    apps = [Appliance] * 10
    for i in range(10):
        if np.random.random() < 0.5:
            act_occ_dependent = False
        else:
            act_occ_dependent = True
        b = Appliance('app_' + str(i), 'app_' + str(i), 'type_name_' + str(i), 'profile_' + str(1), np.random.random(),
                      10 * np.random.random(), 50 * np.random.random(), 1000 * np.random.random(), 100 * np.random.random(),
                      50 * np.random.random(), 0.03 * np.random.random(), 100000000 * np.random.random(), act_occ_dependent,
                      np.random.random(), 0.8 * np.random.random())
        apps[i] = b
    
    for app in apps:
        print(app)
        app.set_ownership()
        if app.owned:
            print('owned')
        else:
            print('not owned')
