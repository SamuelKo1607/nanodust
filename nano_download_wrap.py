import sys
sys.path.insert(0, 'C:\\Users\\skoci\\Documents\\dust\\000_commons')
sys.path.insert(0, 'C:\\Users\\skoci\\Documents\\dust\\003_solar_orbiter')
import datetime as dt
import os
import sunpy_soar
from sunpy.net import Fido
import sunpy.net.attrs as a

from conversions import YYYYMMDD2date
from conversions import date2YYYYMMDD
from solo_download import fetch
from keys import cdf_stat_location
from keys import cdf_tswf_e_location
from keys import cdf_mamp_location



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
                        end_date):
    """
    The function to download all the SOAR files 
    that are nneded fro the nanodust analysis.

    Parameters
    ----------
    strat_date : datetime.datetime()
        Starting date of the selection.
    end_date : datetime.datetime()
        The last downloaded day (inclusive).

    Returns
    -------
    None.

    """

    folders = [cdf_stat_location,
               cdf_tswf_e_location,
               cdf_mamp_location]

    data_products = ["rpw-tds-surv-stat",
                     "rpw-tds-surv-tswf-e",
                     "rpw-tds-surv-mamp"]

    for i in range(len(folders)):
        soar_download(start_date,
                      end_date,
                      "RPW",
                      2,
                      data_products[i],
                      folders[i])


def main(YYYYMMDD_from,
         YYYYMMDD_to):
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

    Returns
    -------
    None.

    """
    idate = YYYYMMDD2date(YYYYMMDD_from)
    date_to = YYYYMMDD2date(YYYYMMDD_to)
    while idate <= date_to:
        try:
            blockPrint()
            cdf_file_e = fetch(date2YYYYMMDD(idate),
                         'C:\\Users\\skoci\\Disk Google\\000 Škola\\UIT\\getting data\\solo\\rpw\\tds_wf_e',
                         "tds_wf_e",
                         "_rpw-tds-surv-tswf-e_",
                         ["V06.cdf","V05.cdf","V04.cdf","V03.cdf","V02.cdf","V01.cdf"],    
                         True)
        except:
            enablePrint()
            print("no data online for "+date2YYYYMMDD(idate))
        else:
            enablePrint()
            print(date2YYYYMMDD(idate)+" OK")
        finally:
            idate += dt.timedelta(days=1)
    
    
if __name__ == "__main__":
    YYYYMMDD_from = sys.argv[1]
    YYYYMMDD_to = sys.argv[2]
    main(YYYYMMDD_from, YYYYMMDD_to)

