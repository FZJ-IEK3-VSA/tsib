from __future__ import division
import numpy as np


class Calculations(object):
    @staticmethod
    def calc_from_normal_distr(mean, sigma):
        # guess number in range of [mean-4*sigma;mean+4*sigma]
        guess = 0.0
        if sigma != 0.0:
            while True:
                guess = (2 * np.random.random() - 1) * 4 * sigma + mean
                # probability of the guess tested against the normal distribution
                # the first factor doesn't make sense IMO
                # 1/(self._sigma * (2*np.pi)**(1/2)) * np.exp(-((irr_threshold-self._mean)**2)/(2*self._sigma**2))
                p_guess = np.exp(-((guess - mean) ** 2) / (2 * sigma ** 2))
                if np.random.random() <= p_guess:
                    break
            return guess
        else:
            return mean
