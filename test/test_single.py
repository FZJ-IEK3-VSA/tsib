import os
import time

import pandas as pd

import tsib


def test_single():

    starttime = time.time()

    # init empty
    profile = tsib.simSingleHousehold(3, 2010, get_hot_water=True, resample_mean=True)

    print("Profile generation took " + str(time.time() - starttime))


if __name__ == "__main__":
    test_single()
