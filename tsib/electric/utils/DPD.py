import numpy as np


class DPD(object):
    """
    Discrete Probability Distribution
    """

    def __init__(self, values):
        self._values = values.astype(np.float)
        self._normalize_values()

    @property
    def get_values(self):
        return self._values

    def _normalize_values(self):
        """
        Normalizes values to one
        :return:
        """
        sum = self._sum()
        for i in range(self._values.size):
            self._values[i] /= sum

    def _sum(self):
        """
        Sum of the Value-array
        """
        sum = 0
        for v in self._values:
            sum = sum + v
        return float(sum)

    def _get_cumu_dist_of_values(self):
        """
        Returns the cumulative distribution of this probability distribution
        :return:
        """
        cumu_dist = np.zeros(self._values.size)
        previous = 0
        for i in range(cumu_dist.size):
            cumu_dist[i] = previous + self._values[i]
            previous = cumu_dist[i]
        return cumu_dist

    def get_random_interval(self):
        """
        Pulls a random number from the cumulative distribution
        :return:
        """
        if self._sum() != 1.0:
            self._normalize_values()
        # draw from Uniform Distribution
        rand = np.random.random()
        interval = 0
        cumu_dist = self._get_cumu_dist_of_values()
        for interval in range(self._values.size):
            if rand <= cumu_dist[interval]:
                break
        return interval

    def set_seed(self, seed):
        """
        Sets the seed for the random number generation
        :param seed:
        """
        np.random.seed(seed)
