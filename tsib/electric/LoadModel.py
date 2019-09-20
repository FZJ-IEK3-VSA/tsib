from __future__ import division
from tsib.electric.Load import Load
from tsib.electric.Appliance import Appliance
import numpy as np
import pandas as pd
import datetime
from collections import Counter
from tsib.electric.Bulb import Bulb
from tsib.electric.utils.Calculations import Calculations
from tsib.electric.utils.ApplianceProbabilityModifier import (
    ApplianceProbabilityModifier,
)


app_act_profiles = [
    "LEVEL",
    "ACTIVE_OCC",
    "CUSTOM",
    "ACT_TV",
    "ACT_COOKING",
    "ACT_LAUNDRY",
    "ACT_WASHDRESS",
    "ACT_IRON",
    "ACT_HOUSECLEAN",
]


class AppliancesModel(object):
    def __init__(
        self,
        data_ex,
        occ_model,
        random_app_seed_per_run=True,
        pre_setup=False,
        get_hot_water=False,
    ):
        """
        Initializes an Appliance Model Object. When initialized a randomized
        application setup is build for the household. If pre_setup is set then
        a pre chosen setup from .csv data will be loaded.
        -----------------------------------------------------------------------
        Parameters:
            data_ex: DataExchange Object, required
                DataExchange Object to exchange data between the Models and
                save runspecific Data
            occ_model: OccupancyModel Object, required
                Object containing the behaviour of the residents in the house-
                hold. Used to calculate the usage of the applications based on
                the activity of the residents.
            random_app_per_run: int, optional (default=True)
                If set the setting of the Appliances owned in the household
                will be seeded. (Same "random" distribution for each
                Initilization with this seed)
            pre_setup: string, optional (default=False)
                If set to a directory Path as a string, a pre chosen setup for
                the Appliance distribution in the household will be loaded
            get_hot_water: bool, optional (default=False)
                If set, an extra hot water profile is created independent
                from the electricity load of the other devices.
        -----------------------------------------------------------------------
        """
        self._data_ex = data_ex
        self._occ_model = occ_model
        self.loads = [Load]

        # Variables for saving Data:
        self.pd_app_type_loads = pd.DataFrame(
            data=np.zeros((1440, 9)),
            columns=app_act_profiles,
            index=pd.date_range(start="00:00", periods=1440, freq="1min").time,
        )
        self.random_app_seed_per_run = random_app_seed_per_run
        self.pre_setup = pre_setup
        # if hot water shall
        self._get_hot_water = get_hot_water
        # all possible appliances
        self.loads = self._get_appliances()
        # assigns appliances to the household (random or pre_setup)
        self._configure_appliances(self.loads)
        # Load List of Activities
        self.activities = self._data_ex.get_activities

    def run(self, year, day_of_week, day_in_year):
        """
        Runs the Object of Type LoadModel
        -----------------------------------------------------------------------
        Parameters:
            year: int, required
                Sets the year in which the Calculations will take place
            day_of_week: str, required
                String either it is a weekday "wd" or weekend "we"----
        Returns:
            Saves the total_consumption of the Appliances and Lighting of the
            household in the specific variable of the Model
        -----------------------------------------------------------------------
        """
        self.total_consumption = np.zeros(1440)
        self.total_heat_gain = np.zeros(1440)
        if self._get_hot_water:
            self.total_hot_water = np.zeros(1440)
        self.run_model(year, day_of_week, day_in_year)
        for load in self.loads:
            if self._get_hot_water and load.type_name == "Water heating":
                self.total_hot_water += load.consumption
            else:
                self.total_consumption += load.consumption

            self.total_heat_gain += np.multiply(load.consumption, load.heat_gain)

    def run_model(self, year, day_of_week, day_in_year):
        """
        Run a ApplianceModel for the given year and type of day. The total
        electric consumption of all applications in the household is calculated
        -----------------------------------------------------------------------
        Parameters:
            year: int, required
                Sets the year in which the Calculations will take place
            day_of_week: str, required
                String either it is a weekday "wd" or weekend "we"
        -----------------------------------------------------------------------
        Return
            saves the appliance electricity consumption in the different
            appliance objects. Also saves the total consumption of all
            appliances in the Appliance Model
        """
        # active occupants over the day
        active_occ_data = self._occ_model.get_occ_activity
        # probability of a certain kind of activity over the day
        activities = self._load_activity_statistics()
        # Simulation of each appliance
        for app in self.loads:
            if app.is_owned():
                minute = 0
                #                print('True')   # Test fuer App-Ausstattung
                while minute < 1440:
                    # ten minute period
                    ten_minute_count = int(minute / 10)
                    # set the standby power at this minute
                    app._power = app._standby_power
                    # number of active occupants
                    # adjusted for a 10 minute interval
                    active_occ = active_occ_data.item(ten_minute_count)
                    # if the appliance is off,
                    # has complete a cycle and a restart delay
                    if app.is_off() and app.waiting_for_restart():
                        app._restart_delay_time_left -= 1
                    # if the appliance is off but able to restart
                    elif app.is_off():
                        if app.act_use_profile == "LEVEL":
                            self._determine_starting_event(minute, app, 1)
                        elif app.act_use_profile == "ACTIVE_OCC":
                            if active_occ > 0:
                                self._determine_starting_event(minute, app, 1)
                        elif app.act_use_profile in [
                            "ACTIVE_OCC",
                            "ACT_TV",
                            "ACT_COOKING",
                            "ACT_LAUNDRY",
                            "ACT_WASHDRESS",
                            "ACT_IRON",
                            "ACT_HOUSECLEAN",
                        ]:
                            if active_occ > 0:
                                # Get the probability of this use profile
                                pm = self.get_probability_modifier(
                                    activities,
                                    active_occ,
                                    app.act_use_profile,
                                    day_of_week,
                                )
                                # Get the activity statistics
                                # for this profile at this time step
                                self._determine_starting_event(
                                    minute, app, pm.get_modifiers[ten_minute_count]
                                )
                    else:
                        # Appliance is on, if the occupants become inactive,
                        # switch off the appliance
                        if (
                            active_occ == 0
                            and app.act_use_profile != "LEVEL"
                            and app.act_use_profile != "ACT_LAUNDRY"
                            and app.act_use_profile != "CUSTOM"
                        ):
                            # Do nothing. The activity will be completed upon the return of the active occupancy.
                            # Note that LEVEL means that the appliance use is not related to active occupancy.
                            # Note also that laundry appliances do not switch off upon a transition to inactive
                            # occupancy.
                            pass
                        else:
                            app.keep_running(minute)
                    app.consumption[minute] = app._power

                    if app._power != 0 and app.switched_on_time_series[minute] != 2:
                        app.switched_on_time_series[minute] = 1
                    minute += 1

        #            else:               # Test Ausstattung Apps
        #                print('False')  # Test Ausstattung Apps
        self._group_loads_by_type()

    def _determine_starting_event(self, minute, app, activity_probability):
        """
        Calculates whether the Appliance is started in the given minute
        """
        if np.random.random() < app._calibration * activity_probability:
            # starts event for the appliance
            app.start(minute)

    def get_probability_modifier(
        self, activities, active_occ, app_use_profile, day_of_week
    ):
        return [
            x
            for x in activities
            if x.day_of_week == day_of_week
            and x.count_active_occupants == active_occ
            and x.key.upper() == app_use_profile.upper()
        ][0]

    def set_seed(self, seed):
        """
        Sets the seed for the random number generation
        :param seed:
        """
        np.random.seed(seed)

    def _load_activity_statistics(self):
        """
        Gets the activity statistics from self.activities where the activity
        propabilities were saved when initializing the Object. To calculate the
        probability the function ApplianceProbabilityModifier is used.
        """
        #        activities = self._data_ex.get_activities
        rows = self.activities.shape[0]
        act_stat_list = [None] * rows
        for row in range(rows):
            if self.activities.item(row, 0) == 1:
                week_of_day = "we"
            else:
                week_of_day = "wd"
            act_stat_list[row] = ApplianceProbabilityModifier(
                week_of_day,
                self.activities.item(row, 1),
                self.activities.item(row, 2).upper(),
                self.activities[row, 3:],
            )
        return act_stat_list

    def _configure_appliances(self, apps):
        """
        Set the "owned" variable inside an appliance object for all appliances
        in the given List. If an appliance is owned by a household is
        calculated by a random generator. This random generator can be seeded
        by the random_app_per_run variable when initializing the Object.
        The other Option is to load a pre chosen setup from a csv file, which
        this option is also defined when init is called.
        """
        # without pre setup possessing of app is randomly chosen
        if self.pre_setup is False:
            for app in apps:
                # set the existance of water heating technologies to 1 if it
                # is an extra profile
                if self._get_hot_water and app.type_name == "Water heating":
                    if not app.get_key == "E_INST".upper():
                        app.owned = True
                    else:
                        app.owned = False
                else:
                    app.set_ownership(random_app_per_run=self.random_app_seed_per_run)
        # other option is to get the app pre_setup via a pre_setup getter
        else:
            i = 0
            for app in apps:
                pre_setup_bool = self._data_ex.get_app_pre_setup
                app.owned = pre_setup_bool[i, 2]
                i += 1

    def _get_appliances(self):
        """
        Returns all the possible applications of a Household
        """
        # get all appliances via app getter
        appliances_data = self._data_ex.get_app_data
        num_of_appliances = appliances_data.shape[0]
        appliances = [Appliance] * num_of_appliances
        for i in range(num_of_appliances):
            # save the data from the loaded csv file inside Appliance Objects
            appliances[i] = Appliance(
                appliances_data[i, 0],
                appliances_data[i, 1],
                appliances_data[i, 2],
                appliances_data[i, 3],
                appliances_data[i, 4],
                appliances_data[i, 5],
                appliances_data[i, 6],
                appliances_data[i, 7],
                appliances_data[i, 8],
                appliances_data[i, 9],
                appliances_data[i, 10],
                appliances_data[i, 11],
                appliances_data[i, 12],
                appliances_data[i, 13],
                appliances_data[i, 14],
                appliances_data[i, 15],
            )
        return appliances

    def _group_loads_by_type(self):
        self.pd_app_type_loads = pd.DataFrame(
            data=np.zeros((1440, 9)),
            columns=app_act_profiles,
            index=pd.date_range(start="00:00", periods=1440, freq="1min").time,
        )
        """
        """
        for a in (a for a in self.loads if a.owned):
            if a.act_use_profile == "LEVEL":
                self.pd_app_type_loads["LEVEL"] += a.consumption
            elif a.act_use_profile == "ACTIVE_OCC":
                self.pd_app_type_loads["ACTIVE_OCC"] += a.consumption
            elif a.act_use_profile == "CUSTOM":
                self.pd_app_type_loads["CUSTOM"] += a.consumption
            elif a.act_use_profile == "ACT_TV":
                self.pd_app_type_loads["ACT_TV"] += a.consumption
            elif a.act_use_profile == "ACT_COOKING":
                self.pd_app_type_loads["ACT_COOKING"] += a.consumption
            elif a.act_use_profile == "ACT_LAUNDRY":
                self.pd_app_type_loads["ACT_LAUNDRY"] += a.consumption
            elif a.act_use_profile == "ACT_WASHDRESS":
                self.pd_app_type_loads["ACT_WASHDRESS"] += a.consumption
            elif a.act_use_profile == "ACT_IRON":
                self.pd_app_type_loads["ACT_IRON"] += a.consumption
            elif a.act_use_profile == "ACT_HOUSECLEAN":
                self.pd_app_type_loads["ACT_HOUSECLEAN"] += a.consumption


class LightingModel(object):
    def __init__(self, data_ex, occ_model, weather_data=None):
        """
        Initializes an Appliance Model Object. When initialized a randomized
        application setup is build for the household. If pre_setup is set then
        a pre chosen setup from .csv data will be loaded.
        -----------------------------------------------------------------------
        Parameters:
            data_ex: DataExchange Object, required
                DataExchange Object to exchange data between the Models and
                save runspecific Data
            occ_model: OccupancyModel Object, required
                Object containing the behaviour of the residents in the house-
                hold. Used to calculate the usage of the applications based on
                the activity of the residents.
            weather_data: pd.DataFrame, optional (default=None)
                A pandas dataframe containing weather data with the GHI as a 
                column.
        -----------------------------------------------------------------------        
        """
        self._data_ex = data_ex
        self._occ_model = occ_model
        self.loads = [Load]

        # House external global irradiance threshold
        self._mean_irradiance = 60
        # Standard deviation
        self._sd_irradiance = 10
        self._has_run = False
        # self.irradiance_dir = irradiance_dir # Keyword argument for directory
        # variable with the lightbulbs of the household
        self.loads = self._get_bulbs()
        self.weather_data = weather_data
        if weather_data is None:
            # irradiance data for one year with typical day for each month
            self._irradiance = self._data_ex.get_irradiance
        else:
            # irradiance TRY-Data
            _GHI = weather_data["GHI"].tz_localize(None)
            freq = _GHI.index[1] - _GHI.index[0]
            # add values in the beginning and end to extrapolate
            startVal = pd.Series(_GHI[0], index=[_GHI.index[0] - freq])
            endVal = pd.Series(_GHI[-1], index=[_GHI.index[-1] + freq])
            _GHI_extended = pd.concat([startVal, _GHI, endVal])
            _GHI_minute = (
                _GHI_extended.resample("1min", label="left", closed="left")
                .mean()
                .interpolate()
            )
            self._irradiance = _GHI_minute
        # factor for collective usage of lighting with more than one resident
        self._occ_sharing_factor_data = self._data_ex.get_occ_sharing_factor
        self._bulb_duration = self._data_ex.get_bulbs_duration_config
        # print(_GHI_year_minute) #test

    @property
    def get_has_run(self):
        return self._has_run

    def run(self, year, day_of_week, day_in_year):
        """
        Runs the Object of Type LoadModel
        -----------------------------------------------------------------------
        Parameters:
            year: int, required
                Sets the year in which the Calculations will take place
            day_of_week: str, required
                String either it is a weekday "wd" or weekend "we"----
        Returns:
            Saves the total_consumption of the Appliances and Lighting of the
            household in the specific variable of the Model
        -----------------------------------------------------------------------
        """
        self.total_consumption = np.zeros(1440)
        self.total_heat_gain = np.zeros(1440)
        self.run_model(year, day_of_week, day_in_year)
        for load in self.loads:
            self.total_consumption += load.consumption
            self.total_heat_gain += np.multiply(load.consumption, load.heat_gain)

    def run_model(self, year, day_of_week, day_in_year):
        """
        Run a LightingModel for the given year and type of day. The total
        electric consumption of all lighting in the household is calculated
        -----------------------------------------------------------------------
        Parameters:
            year: int, required
                Sets the year in which the Calculations will take place
            day_of_week: str, required
                String either it is a weekday "wd" or weekend "we"
            day_in_year: int, required
                numerical day in year counted from 01.01.xxxx counted as 1. 
                Is generated automatically when run in full year calculations
        -----------------------------------------------------------------------
        Return
            saves the lighting electricity consumption in the different
            Bulb objects. Also saves the total consumption of all
            Lighting in the Lighting Model
        """
        # calculate the irradiance threshold for the house
        irr_threshold = Calculations.calc_from_normal_distr(
            self._mean_irradiance, self._sd_irradiance
        )
        if self.weather_data is None:
            # get the irradiance data for the given month
            date = datetime.datetime(year, 1, 1) + datetime.timedelta(day_in_year - 1)
            irradiance_data = self._irradiance[:, date.month - 1]
        else:
            # Get TRY-irradiance Data for the given day
            irradiance_data = self._irradiance[
                (self._irradiance.index.dayofyear == day_in_year)
                & (self._irradiance.index.year == year)
            ].values
        if not irradiance_data.size:
            raise ValueError(
                "No irradiance data available for year "
                + str(year)
                + " and day in year "
                + str(day_in_year)
            )
        # get the activity data from the OccupancyModel
        active_occ_data = self._occ_model.get_occ_activity
        # simulation loop for each bulb
        for bulb in self.loads:
            minute = 0
            # Consumption set to zero before each run per day
            bulb.consumption = np.zeros(1440)
            # start Loop for one day (1440=24h*60min)
            while minute < 1440:
                # irradiance for a day in this month at this minute
                irradiance = irradiance_data.item(minute)
                irradiance = round(irradiance)
                # number of active occupants, adjusted for a 10 minute interval
                active_occ = active_occ_data.item(int(minute / 10))
                # inefficient irradiance plus \
                # a small chance to turn on the bulb anyway
                irr_below_threshold = (
                    irradiance < irr_threshold or np.random.random() < 0.05
                )
                # Sharing Factor for collective usage
                occ_sharing_factor = self._occ_sharing_factor_data[active_occ, 1]
                # Combine sharing and bulb_weight factor
                sharing_and_weight_factor = occ_sharing_factor * float(bulb.weight)
                if irr_below_threshold and (
                    np.random.random() < sharing_and_weight_factor
                ):
                    duration = self._calc_bulb_duration(self._bulb_duration)
                    for i in range(0, duration):
                        if minute >= 1440:
                            break
                        if active_occ_data.item(int(minute / 10)) != 0:
                            bulb.switch_on(minute)
                        minute += 1
                else:
                    minute += 1
        self._has_run = True

    def sum_total_consumption(self):
        """
        Calculates the aggregated consumption over a day
        """
        return self.total_consumption.sum(0) / 60 / 1000

    def _get_bulbs(self):
        """
       
        Returns the installed bulbs in the household
        based on 100 sample bulb configurations
        """

        x = 0
        bulbs_number_rand = np.random.normal(25, 4.426778, 1)
        bulbs_number_rand = bulbs_number_rand.round(0)
        bulbs_number_rand_i = int(bulbs_number_rand)
        # create List with i numbers of bulbs to put on configuration
        bulb_config = list(range(bulbs_number_rand_i))

        types = list()
        types = self._data_ex.get_bulbs_type[:, 0]
        random = np.random.choice(types, len(bulb_config), replace=False)
        # counts the the number of a certain type of bulb
        counts_type = Counter(random)
        # Create DataFrame out of .Counter
        counts_type_pd = pd.DataFrame(counts_type, index=[x])
        # counts_type_a=counts_type.append(counts_type)
        # x_flats.append(counts_type_pd)
        counts_type_pd_i = counts_type_pd.fillna(0).astype(int)

        i = 0
        # Lists with the length of the number of each type of bulb
        if "Glueh" in counts_type_pd_i.columns:
            Glueh = counts_type_pd_i.loc[i, "Glueh"]
        else:
            Glueh = 0
        if "Halogen" in counts_type_pd_i.columns:
            Halogen = counts_type_pd_i.loc[i, "Halogen"]
        else:
            Halogen = 0
        if "Led" in counts_type_pd_i.columns:
            Led = counts_type_pd_i.loc[i, "Led"]
        else:
            Led = 0
        if "Spar" in counts_type_pd_i.columns:
            Spar = counts_type_pd_i.loc[i, "Spar"]
        else:
            Spar = 0
        if "Leucht" in counts_type_pd_i.columns:
            Leucht = counts_type_pd_i.loc[i, "Leucht"]
        else:
            Leucht = 0
        # preparing DataFrames to select a power
        power_glueh = self._data_ex.get_bulbs_power["Glueh"].dropna()
        power_spar = self._data_ex.get_bulbs_power["Spar"].dropna()
        power_led = self._data_ex.get_bulbs_power["Led"].dropna()
        power_halogen = self._data_ex.get_bulbs_power["Halogen"].dropna()
        power_leucht = self._data_ex.get_bulbs_power["Leucht"].dropna()

        p_glueh_random = np.random.choice(power_glueh, Glueh, replace=True)
        p_spar_random = np.random.choice(power_spar, Spar, replace=True)
        p_led_random = np.random.choice(power_led, Led, replace=True)
        p_halogen_random = np.random.choice(power_halogen, Halogen, replace=True)
        p_leucht_random = np.random.choice(power_leucht, Leucht, replace=True)

        bulb_rating = np.concatenate(
            (
                p_glueh_random,
                p_spar_random,
                p_halogen_random,
                p_led_random,
                p_leucht_random,
            ),
            axis=0,
        )
        bulbs = [Bulb] * len(bulb_rating)
        for i in range(len(bulb_rating)):
            bulbs[i] = Bulb("BULB_" + str(i), bulb_rating[i])
        return bulbs

    def _get_bulbs_old(self):
        """
       
        Returns the installed bulbs in the household
        based on 100 
        """

        # get all given bulb configs from csv file
        bulbs_data = self._data_ex.get_bulbs
        num_of_bulb_configs = bulbs_data.shape[0]
        # house number from 0..num_of_bulb_configs-1
        rand_house = np.random.randint(0, num_of_bulb_configs)
        num_of_bulbs = int(bulbs_data[rand_house, 0])
        bulbs = [Bulb] * num_of_bulbs
        for i in range(num_of_bulbs):
            bulbs[i] = Bulb("BULB_" + str(i), bulbs_data[rand_house, i + 1])
        return bulbs

    def _calc_bulb_duration(self, bulb_duration):
        """
        Calculates randomly the duration the lightbulb is turned on
        Returns integer with duration as int
        """
        rand = np.random.randint(9)
        lower = bulb_duration[rand, 1]
        upper = bulb_duration[rand, 2]
        return int(np.random.random() * (upper - lower) + lower)

    def set_seed(self, seed):
        """
        Sets the seed for the random number generation
        :param seed:
        """
        np.random.seed(seed)


#    def export_total_consumption_time_series(self):
#        data = self.total_consumption
#        index = pd.date_range(start='00:00:00', periods=1440, freq='1min').time
#        df = pd.DataFrame(data=data, index=index)
#        self._data_ex.export_data('lig_total_consumption.csv', df)
#
# def run_lig_model_n_times(num_runs: 10):
#    data_ex = DataExchangeFactory(dbtype=DbType.csv).create()
#    num_occ = 5
#    day = DayOfWeek.weekday
#    month = 1
#    average_bulb_consumption = 0
#    for i in range(num_runs):
#        lig_model = LightingModel(data_ex, OccupancyModel(data_ex, num_occ, day), month)
#        lig_model.run()
#        average_bulb_consumption += lig_model.total_consumption / 60
#        print(i)
#    average_bulb_consumption /= num_runs
#    print(average_bulb_consumption.sum(0))
#    data = average_bulb_consumption
#    index = pd.date_range(start='00:00:00', periods=1440, freq='1min').time
#    df = pd.DataFrame(data=data, index=index)
#    file_str = 'bulb_loads_average_' + str(num_runs) + 'runs_' + str(num_occ) + 'pers_' + str(day) + '_' + str(
#        month) + '_month' + '.csv'
#    data_ex_main.export_data(file_str, df)
