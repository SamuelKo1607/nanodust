import sys
sys.path.insert(0, 'C:\\Users\\skoci\\Documents\\dust\\000_commons')
sys.path.insert(0, 'C:\\Users\\skoci\\Documents\\dust\\003_solar_orbiter')
import datetime as dt
import os
from conversions import YYYYMMDD2date
from conversions import date2YYYYMMDD
from solo_download import fetch


# Disable print
def blockPrint():
    sys.stdout = open(os.devnull, 'w')


# Restore print
def enablePrint():
    sys.stdout = sys.__stdout__


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
    
    
main(sys.argv[1],sys.argv[2])

