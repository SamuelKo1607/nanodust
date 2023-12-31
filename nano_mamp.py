import numpy as np
import pandas as pd
import cdflib
import matplotlib.pyplot as plt
import datetime as dt
import sys
import glob
import os

from paths import cdf_mamp_location

from nano_load_days import Impact
from nano_load_days import Day
from nano_load_days import get_cdfs_to_analyze
from nano_load_days import load_all_days
from nano_load_days import load_list
from nano_load_days import save_list
from nano_load_days import moving_average
from nano_load_days import moving_stdev
from conversions import tt2000_to_date



class ImpactSuspect:
    """
    Each instance holds the attributes of 
    one suspected dust impact as per MAMP. 
    """
    def __init__(self,
                 datetime,
                 sampling_rate,
                 coverage_hours,
                 period,
                 index,
                 classification,
                 channel_ref,
                 amplitudes,
                 threshold):
        self.datetime = datetime
        self.YYYYMMDD = datetime.strftime('%Y%m%d')
        self.sampling_rate = sampling_rate
        self.coverage_hours = coverage_hours
        self.period = period
        self.index = index
        self.classification = classification
        if classification == 1:
            self.classification_readable = "Confirmed dust"
        elif classification == -1:
            self.classification_readable = "Confirmed non-dust"
        else:
            self.classification_readable = "Not confirmed"
        self.channel_ref = channel_ref
        self.amplitudes = amplitudes
        self.produced = dt.datetime.now()

    def info(self):
        print(self.datetime.strftime("%d.%m.%Y %H:%M:%S, ")+
              self.classification_readable + ","
              " \n sampling rate: "+ str(self.sampling_rate)+
              " \n amplitude: " + str(self.amplitudes))

    def show(self):
        print(vars(self))







def load_all_suspects_particular(location = os.path.join("998_generated","mamp_processed",""),
                                 date_from = dt.datetime(2010,1,1),
                                 date_to = dt.datetime(2050,1,1)):
    """
    The function to load all the suspected impacts as per MAMP, loads all the
    pickles from the given directory.

    Parameters
    ----------
    location : str, optional
        The data directory. Default is "998_generated\\mamp_processed\\".
    date_from : dt.datetime
        The first relevant moment (filtering the instances by >=).
    date_to : dt.datetime
        The last relevant moment (filtering the instances by <=).

    Returns
    -------
    suspects : list of ImpactSuspect object
        Impact suspects, class ImpactSuspect.
    """

    location = os.path.join(os.path.normpath( location ), '')

    files = glob.glob(location+"*.pkl")
    suspects = []
    for file in files:
        suspect = load_list(file,"")
        suspects.append(suspect)

    flat_list_of_impacts = [item for sublist in suspects for item in sublist]

    filtered = [i for i in flat_list_of_impacts if date_from <= i.datetime <= date_to ]

    return filtered


def load_all_suspects(location = os.path.join("998_generated","mamp_processed",""),
                      aggregate_file = os.path.join("data_synced","all_suspects.pkl"),
                      date_from = dt.datetime(2010,1,1),
                      date_to = dt.datetime(2050,1,1)):
    """
    Loads all suspects into a list of ImpactSuspect.

    Parameters
    ----------
    location : str, optional
        The data directory. Default is "998_generated\\mamp_processed\\".
    aggregate_file : TYPE, optional
        The aggregate (synced) data directory. 
        The default is os.path.join("data_synced","all_suspects.pkl").
    date_from : dt.datetime
        The first relevant moment (filtering the instances by >=).
    date_to : dt.datetime
        The last relevant moment (filtering the instances by <=).

    Returns
    -------
    all_suspects : list of ImpactSuspect
        All the impact suspects.

    """

    try:
        with open(aggregate_file, "rb") as f:
            print("loaded aggregate file:")
            print(f)
            all_suspects = load_list(aggregate_file,location=None)

    except:
        print("loading non-aggregate files:")
        print(location)
        all_suspects = load_all_suspects_particular(location)

    finally:
        all_suspects = [i for i in all_suspects if date_from <= i.datetime <= date_to ]

    return all_suspects


def is_confirmed(dust_datetimes,
                 non_dust_datetimes,
                 suspect_datetimes,
                 batch = 200,
                 threshold = dt.timedelta(milliseconds=65)):
    """
    The function to return the information whether the
    waveform was classified as dust or no-dust.
 
    Parameters
    ----------
    dust_datetimes : np.array of dt.datetimes
        These are the confirmed dust times to check against.
    non_dust_datetimes : np.array of dt.datetimes
        These are the confirmed non-dust times to check against.
    suspect_datetimes : np.array of dt.datetime
        The times of interest.
    batch : int
        how many suspect_datetimes to process at once. Done due to ram, with
        200 batches, one instance takes up about 2GiB
    threshold : dt.timedelta
        The maximum difference, inclusive. If datetime2_list is a list,
        then true is returned if there is a match between datetime_1 and
        any of the datetime2_list's elements.
        The default is dt.timedelta(milliseconds=65).
    

    Returns
    -------
    classification : np.array of int, length = len(suspect_datetimes)
        1, -1 or 0, where 0 means unknown, 
        1 means dust and -1 means non-dust.

    """

    confirmed_dust = np.zeros(0,dtype=int)
    confirmed_non_dust = np.zeros(0,dtype=int)

    for i in range(len(suspect_datetimes)//batch+1):
        #get a slice of all the suspects
        suspect_datetimes_chunk = suspect_datetimes[batch*i:batch*(i+1)]

        #the check if the suspects in this chunk are confirmed dusts
        dust_diffs = np.abs(np.subtract(suspect_datetimes_chunk.reshape((len(suspect_datetimes_chunk),1)),
                                 dust_datetimes.reshape((1,len(dust_datetimes)))))
        min_dust_diffs = np.min(dust_diffs,axis=1)
        confirmed_dust_chunk = min_dust_diffs < threshold

        #the check if the suspects in this chunk are confirmed non-dusts
        non_dust_diffs = np.abs(np.subtract(suspect_datetimes_chunk.reshape((len(suspect_datetimes_chunk),1)),
                                     non_dust_datetimes.reshape((1,len(non_dust_datetimes)))))
        min_non_dust_diffs = np.min(non_dust_diffs,axis=1)
        confirmed_non_dust_chunk = min_non_dust_diffs < threshold

        #append
        confirmed_dust = np.append(confirmed_dust,
                                   confirmed_dust_chunk)
        confirmed_non_dust = np.append(confirmed_non_dust,
                                   confirmed_non_dust_chunk)

    classification = confirmed_dust - confirmed_non_dust

    return classification


def get_coverage(epoch,threshold = 0.95):
    """
    A function to get find out the data coverage for a given day, in 
    case of intermittent coverage.

    Parameters
    ----------
    epoch : np.array of float
        The times of datapoints.
    threshold : float, optional
        What fraction of the diffs in time we want to cover. 
        Too low value gives only approximate result, 
        but too high may include the long waiting times between the bursts,
        which would be completely wrong. This is a workaroundsince I did not
        find the vale of data coverage or something similar in the cdf data.
        The default is 0.95.

    Returns
    -------
    time_covered : float
        Coverage in hours, can't be much hugher than 24.

    """
    values, counts = np.unique(np.diff(epoch), return_counts=True)
    inds = (-counts).argsort()
    sorted_values = values[inds]
    sorted_counts = counts[inds]

    sorted_counts_cumsum = np.cumsum(sorted_counts)
    elements_needed = np.min(
        np.arange(1,1+len(sorted_counts_cumsum))[
        sorted_counts_cumsum > threshold*len(epoch)
        ])
    weighted_mean_period = np.sum(
        sorted_counts[0:elements_needed]*sorted_values[0:elements_needed]
        / np.sum(sorted_counts[0:elements_needed]) )

    time_covered = weighted_mean_period * len(epoch) / 3.6e12
    return time_covered


def get_mamp_suspects(wfs,
                      threshold=4,
                      window=10**5-1,
                      consecutive_allowed=True):
    """
    The function to return the indices of suspected dust impacts
    among the electrical waveforms of MAMP.

    Parameters
    ----------
    wfs : np.array of float, shape (n, 3)
        The three waveform channels as recorded with MAMP in XLD1.

    threshold : float
        The voltage that triggers a suspect, in sigmas (stdevs).

    Returns
    -------
    suspects : np.array of int, 1D
        The indices of suspected dust impacts.

    thresholds : np.array of float, 1D
        The thresholds applicable to every 
    """

    #floor at 0V
    wfs_floor = np.maximum(wfs[:,2],np.zeros_like(wfs[:,2]))
    length = len(wfs_floor)

    #check that the window is not too long
    if window>length//100:
        window = max(length//100-1+(length//100)%2, 3)

    #compute mean
    wfs_mean = moving_average(wfs_floor,n=window)

    #compute moving stdev
    wfs_stdev = moving_stdev(wfs_floor,n=window)

    #get effective threshold
    verge = wfs_mean + threshold*wfs_stdev

    #get dust suspects
    suspects = np.arange(len(wfs_floor))[(wfs[:,2]>verge)]

    if not consecutive_allowed:
        #abbandon the suspects that are consecutive, both of them
        neighbouring = np.concatenate( (np.diff(suspects)==1,
                                       (np.zeros(1,dtype=bool)) ))
        neighbouring += np.concatenate((np.zeros(1,dtype=bool),neighbouring[:-1]))
        suspects = suspects[neighbouring==False]

    #get thresholds
    thresholds = verge[suspects]

    return suspects, thresholds


def evaluate_thresholds(cdf_file_path,
                        sigmas=np.arange(1,40,0.4),
                        figures_location=os.path.join("998_generated","figures","threshlod_eval","")):
    """
    To decide what should be the right sigmas threshold when looking for dust.

    Call as >>> [evaluate_thresholds(cdf) for cdf in 
                get_cdfs_to_analyze(cdf_mamp_location,"*mamp*.cdf")]

    Parameters
    ----------
    cdf_file_path : str
        The path to a MAMP cdf file.
    sigmas : np.array of float, optional
        The sigmas to test. The default is np.arange(1,10,0.2).
    figures_location : str, optional
        The path there to save the figure. 
        The default is os.path.join("998_generated","figures","").

    Raises
    ------
    Exception
        DESCRIPTION.

    Returns
    -------
    None.

    """

    cdf_file = cdflib.CDF(cdf_file_path)
    wfs = cdf_file.varget("WAVEFORM_DATA")[:,0:3]
    channel_refs = cdf_file.varget("CHANNEL_REF")
    if np.allclose(channel_refs[:,2],2): #XLD1
        pass
    else:
        raise Exception("they are not all XLD1")

    suspects = [np.shape(get_mamp_suspects(wfs,threshold=sigma))[1] for sigma in sigmas]
    fig,ax = plt.subplots()

    ax.semilogy(sigmas, suspects)
    ax.set_xlabel("sigmas")
    ax.set_ylabel("threshold crossings")

    fig.tight_layout()
    fig.savefig(figures_location+"thresholds_"+cdf_file_path[-21:-4]+'.png',
                format='png', dpi=600)
    fig.show()


def suspects_stat(date_from = dt.datetime(2010,1,1),
                  date_to = dt.datetime(2050,1,1)):
    """
    To get a quick overview of the suspects extracted from MAMP 
    and how they compare to the data from Days.

    Parameters
    ----------
    date_from : dt.datetime
        The first relevant moment (filtering the instances by >=).
    date_to : dt.datetime
        The last relevant moment (filtering the instances by <=).

    Returns
    -------
    None.

    """

    all_suspects = load_all_suspects(date_from = date_from,
                                     date_to = date_to)
    all_days = load_all_days()
    YYYYMMDDs = list(set([suspect.YYYYMMDD for suspect in all_suspects]))

    for YYYYMMDD in YYYYMMDDs:
        print(YYYYMMDD)
        day_match = [day for day in all_days if day.YYYYMMDD == YYYYMMDD]
        if len(day_match)>0:
            day = day_match[0]
            suspects = [suspect for suspect in all_suspects if suspect.YYYYMMDD == YYYYMMDD]
            confirmed_positive = [suspect for suspect in suspects if suspect.classification == 1]
            confirmed_negative = [suspect for suspect in suspects if suspect.classification == -1]
            print("impacts on Day: "+str(len(day.impact_times)))
            print(f"conf. pos/neg {len(confirmed_positive)}/{len(confirmed_negative)} of {len(suspects)} MAMP suspects")
        else:
            print("Day data not available")


def main(target_input_cdf,
         target_output_pkl,
         threshold=12):
    """
    The main routine, takes the given MAMP cdf, finds all the dust suspects 
    based on the amplitudes and the threshold. 

    Parameters
    ----------
    target_input_cdf : TYPE
        DESCRIPTION.
    target_output_pkl : TYPE
        DESCRIPTION.
    threshold : float, optional
        The sigmas that trigger a suspect, as in get_mamp_suspects().
        The default is 4.

    Returns
    -------
    None.

    """

    cdf_file = cdflib.CDF(target_input_cdf)
    wfs = cdf_file.varget("WAVEFORM_DATA")[:,0:3]
    channel_refs = cdf_file.varget("CHANNEL_REF")
    sampling_rates = cdf_file.varget("sampling_rate")

    suspect_events = []
    if np.allclose(channel_refs[:,2],2): #XLD1
        suspects,thresholds = get_mamp_suspects(wfs,threshold=threshold)
    else:
        suspects,thresholds = [],[]
    epochs = cdf_file.varget("Epoch")
    suspect_epochs = epochs[suspects]
    coverage = get_coverage(epochs)

    if len(suspects):
        suspect_datetimes = tt2000_to_date(suspect_epochs)
    
        days = load_all_days()
        dust_times = [day.impact_times for day in days]
        flat_dust_times = np.array([item for sublist in dust_times for item in sublist])
        non_dust_times = [day.non_impact_times for day in days]
        flat_non_dust_times = np.array([item for sublist in non_dust_times for item in sublist])
    
        confirmed_classification = is_confirmed(flat_dust_times,
                                                flat_non_dust_times,
                                                suspect_datetimes)

        for i, datetime in enumerate(suspect_datetimes):
            channel_ref = channel_refs[i]
            sampling_rate = sampling_rates[i]
            period = min([epochs[i]-epochs[i-1],epochs[i+1]-epochs[i]])
            classification = confirmed_classification[i]
            amplitudes = wfs[i,:]
    
            suspect = ImpactSuspect(datetime,
                                    sampling_rate,
                                    coverage,
                                    period,
                                    i,
                                    classification,
                                    channel_ref,
                                    amplitudes,
                                    thresholds[i])
    
            suspect_events.append(suspect)
    else:
        pass

    save_list(suspect_events, target_output_pkl)
    


"""

cdfs = get_cdfs_to_analyze(cdf_mamp_location,"*mamp*.cdf")
#file = 'C:\\Users\\skoci\\Disk Google\\000 Škola\\UIT\\getting data\\solo\\rpw\\mamp\\solo_L2_rpw-tds-surv-mamp_20220219_V02.cdf'

for cdf in cdfs:
    main(cdf,
         ("C:\\Users\\skoci\\Documents\\nanodust\\998_generated\\mamp_processed\\"+
         cdf[cdf.find("solo_L2_rpw-tds-surv-mamp"):-4]+
         ".pkl"))

"""

#%%
if __name__ == "__main__":
    if len(sys.argv)>1:
        input_file = sys.argv[1]
        output_file = sys.argv[2]
        main(input_file, output_file)
    else:
        suspects_stat()


