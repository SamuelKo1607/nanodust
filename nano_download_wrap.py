import sys
import datetime as dt
import os
import sunpy_soar
from sunpy.net import Fido
import sunpy.net.attrs as a

from conversions import YYYYMMDD2date
from conversions import date2YYYYMMDD
from paths import cdf_stat_location
from paths import cdf_tswf_e_location
from paths import cdf_mamp_location



# Disable print
def blockPrint():
    sys.stdout = open(os.devnull, 'w')


# Restore print
def enablePrint():
    sys.stdout = sys.__stdout__


def soar_download(dt_from,
                  dt_to,
                  data_instrument,
                  data_level,
                  data_product,
                  folder):
    """
    A fetching function for data files from SOAR.

    Parameters
    ----------
    dt_from : datetime.datetime()
        Starting date of the selection.
    dt_to : datetime.datetime()
        The last downloaded day.
    data_instrument : str
        Such as "EUI" or "SPICE".
    data_level : int
        0 : raw, 1 : engineering, 2 : scientific, 3 : higher
    data_product : str
        Such as "EUI-FSI304-IMAGE" 
            or "rpw-tds-surv-tswf-e" 
            or "rpw-tds-surv-stat".
    folder : str
        Where to place the downloaded files.

    Returns
    -------
    files : parfive.Results
        The download report.

    """
    
    time_from = str(dt_from.date())
    time_to = str(dt_to.date())
    
    instrument = a.Instrument(data_instrument)
    time = a.Time(time_from, time_to)
    level = a.Level(data_level)
    product = a.soar.Product(data_product)
    
    result = Fido.search( instrument & time & level & product )
    files = Fido.fetch(result,path=folder)
    return files


def download_update_cdf(start_date,
                        end_date,
                        data_products = ["rpw-tds-surv-stat",
                                         "rpw-tds-surv-tswf-e",
                                         "rpw-tds-surv-mamp"]):
    """
    The function to download all the SOAR files 
    that are nneded fro the nanodust analysis.

    Parameters
    ----------
    strat_date : datetime.datetime()
        Starting date of the selection.
    end_date : datetime.datetime()
        The last downloaded day (inclusive).

    data_products: list of string
        The list of data products to download. Defaults to: 
            ["rpw-tds-surv-stat",
             "rpw-tds-surv-tswf-e",
             "rpw-tds-surv-mamp"].

    Returns
    -------
    None.

    """

    folders = [cdf_stat_location,
               cdf_tswf_e_location,
               cdf_mamp_location]

    for product in data_products:
        if product == "rpw-tds-surv-stat":
            folder = cdf_stat_location
        elif product == "rpw-tds-surv-tswf-e":
            folder = cdf_tswf_e_location
        elif product == "rpw-tds-surv-mamp":
            folder = cdf_mamp_location
        else:
            raise Exception("Unexpected product")
        print("Downloading "+product+" to "+folder)
        
        soar_download(start_date,
                      end_date,
                      "RPW",
                      2,
                      product,
                      folder)


def main(YYYYMMDD_from,
         YYYYMMDD_to,
         data_products = ["rpw-tds-surv-stat",
                          "rpw-tds-surv-tswf-e",
                          "rpw-tds-surv-mamp"]):
    """
    The function to go through all the days in the range and chechk whether 
    the data for the day si downlaoded, downloads if it was not downloaded 
    beforehand. 

    Parameters
    ----------
    YYYYMMDD_from : str / int
        Date from (inclusive).
    YYYYMMDD_to : str / int
        Date to (inclusive).

    data_products: list of string
        The list of data products to download. Defaults to: 
            ["rpw-tds-surv-stat",
             "rpw-tds-surv-tswf-e",
             "rpw-tds-surv-mamp"].

    Returns
    -------
    None.

    """

    start_date = YYYYMMDD2date(YYYYMMDD_from)
    end_date = YYYYMMDD2date(YYYYMMDD_to)
    download_update_cdf(start_date,end_date,data_products)

    
    
if __name__ == "__main__":
    YYYYMMDD_from = sys.argv[1]
    YYYYMMDD_to = sys.argv[2]
    if len(sys.argv)>3:
        prod = sys.argv[3]
        main(YYYYMMDD_from, YYYYMMDD_to, [prod])
    else:
        main(YYYYMMDD_from, YYYYMMDD_to)

