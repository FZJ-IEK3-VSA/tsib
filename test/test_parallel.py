import os
import time

import pandas as pd

import tsib.household.profiles as CREST


def test_parallel():

    starttime = time.time()

    res = CREST.run_households_parallel(3, 2010, 2, cores=2, get_hot_water=True)

    print("Profile generation took " + str(time.time() - starttime))


if __name__ == "__main__":
    test_parallel()
