
import cdsapi

def test_get_era5_example():
    '''
    Gets an era5 test set from the cds via api.
    '''
    c = cdsapi.Client()
    c.retrieve("reanalysis-era5-pressure-levels",
    {
    "variable": "temperature",
    "pressure_level": "1000",
    "product_type": "reanalysis",
    "year": "2008",
    "month": "01",
    "day": "01",
    "time": "12:00",
    "format": "grib"
    },
    "download.grib")

