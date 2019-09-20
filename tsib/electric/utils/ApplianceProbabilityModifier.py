import numpy as np


class ApplianceProbabilityModifier(object):
    '''
    Describes a factor used to modify the probability of an appliance
    running. Each modifier contains a vector of doubles which describes the
    proportion of households where at least one occupant is engaged in a
    particular activity during a particular ten minute period
    '''

    def __init__(self, day_of_week, count_active_occupants, key, modifiers):
        self.day_of_week = day_of_week
        self.count_active_occupants = count_active_occupants
        self.key = key
        self._modifiers = modifiers

    @property
    def get_modifiers(self):
        return self._modifiers