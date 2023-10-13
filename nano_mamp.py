import numpy as np
import cdflib
import matplotlib.pyplot as plt
import datetime as dt

from keys import cdf_mamp_location

from nano_load_days import get_cdfs_to_analyze
from nano_load_days import load_all_days
from nano_load_days import Impact
from nano_load_days import Day
from conversions import tt2000_to_date


def dt_isclose(datetime1,
               datetime2_list,
               thershold = dt.timedelta(milliseconds=65)):
    """
    Will decide if two datetimes are close to each other.

    Parameters
    ----------
    datetime1 : dt.datetime
        The first time of interest.
    datetime2_list : dt.datetime or a list of dt.datetime
        The second time of interest, possibly a list.
    thershold : dt.timedelta, optional
        The maximum difference, inclusive. If datetime2_list is a list,
        then true is returned if there is a match between datetime_1 and
        any of the datetime2_list's elements.
        The default is dt.timedelta(milliseconds=65).

    Returns
    -------
    isclose : bool
        True if they are close.

    """
    if max(datetime1,datetime2_list) - min(datetime1,datetime2_list) > thershold:
        isclose = False
    else:
        isclose = True

    #TBD make it compare a single time1 with a list2

    return isclose


def is_confirmed(days,datetime):
    """
    The function to return the information whether the
    waveform was classified as dust or no-dust.
 
    Parameters
    ----------
    days : list of Day object
        These days will be checked when confirming whether 
        the time is a known impact or non-impact.
    datetime : 
        The time of interest.

    Returns
    -------
    confirmed : int
        1, -1 or 0, where 0 means unknown, 
        1 means dust and -1 means non-dust.

    """
    dust_times = [day.impact_times for day in days]
    flat_dust_times = [item for sublist in dust_times for item in sublist]

    non_dust_times = [day.non_impact_times for day in days]
    flat_non_dust_times = [item for sublist in non_dust_times for item in sublist]

    #TBD take the datetime and compare it to both dusts and non-dusts,
    # return 1 if dust, return -1 if non-dust and return 0 if non classfied
    # or in the unlikely case that both dust and non-dust are a match

    classification = 0

    return classification





cdfs = get_cdfs_to_analyze(cdf_mamp_location,"*mamp*.cdf")

file = 'C:\\Users\\skoci\\Disk Google\\000 Škola\\UIT\\getting data\\solo\\rpw\\mamp\\solo_L2_rpw-tds-surv-mamp_20220219_V02.cdf'


cdf_file = cdflib.CDF(file)
print(cdf_file.file)
wfs = cdf_file.varget("WAVEFORM_DATA")[:,0:3]

#floor at -0.1V
wfs = np.maximum(wfs,np.zeros_like(wfs)-0.1)

#remove median but remember it for later
medians = np.median(wfs,axis=0)
wfs -= medians

#get dust suspects
threshold = 0.001
suspects = np.arange(len(wfs[:,0]))[(wfs[:,0]>threshold)+
                                    (wfs[:,1]>threshold)+
                                    (wfs[:,2]>threshold)]

#flag as sure dust or sure no dust or unknown based on the cnn results
# TBD save no-dusts from CNN


epochs = cdf_file.varget("Epoch")
dates = tt2000_to_date(epochs)
channel_ref = cdf_file.varget("CHANNEL_REF")[0]
print(cdf_file.varget("sampling_rate")[0])
print(dates[2]-dates[1])





plot_from, plot_to = 0, len(wfs[:,0]) #249230, 249270
print(dates[plot_from])
for i in range(3):
    plt.plot(dates[plot_from:plot_to],
             wfs[plot_from:plot_to,i],
             label = channel_ref[i])
plt.xticks(rotation=90)
plt.legend()
plt.show()