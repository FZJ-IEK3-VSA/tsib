import os
import time

import pandas as pd

import tsib


def test_parallel():

    starttime = time.time()

    res = tsib.simHouseholdsParallel(3, 2010, 2, cores=2, get_hot_water=True)

    print("Profile generation took " + str(time.time() - starttime))


if __name__ == "__main__":
    test_parallel()
