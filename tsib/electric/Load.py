import numpy as np


class Load(object):
    def __init__(self, key):
        self._key = key
        self.consumption = np.zeros(
            1440
        )  # pd.DataFrame(np.zeros((1440, 1), dtype=int), columns=['Load'])

    @property
    def get_key(self):
        return self._key

    def is_on(self, minute):
        if self.consumption.item(minute) != 0:
            return True
        else:
            return False

    def __str__(self):
        return "Key: {0}".format(self._key)
