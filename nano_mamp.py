import numpy as np
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
from conversions import tt2000_to_date



class ImpactSuspect:
    """
    Each instance holds the attributes of 
    one suspected dust impact as per MAMP. 
    """
    def __init__(self,
                 datetime,
                 sampling_rate,
                 period,
                 index,
                 classification,
                 channel_ref,
                 amplitudes):
        self.datetime = datetime
        self.YYYYMMDD = datetime.strftime('%Y%m%d')
        self.sampling_rate = sampling_rate
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








def load_all_suspects(location = os.path.join("998_generated","mamp_processed",""),
                      date_from = dt.datetime(2010,1,1),
                      date_to = dt.datetime(2050,1,1)):
    """
    The function to load all the suspected impacts as per MAMP.

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



def get_mamp_suspects(wfs,threshold = 0.002):
    """
    The function to return the indices of suspected dust impacts
    among the electrical waveforms of MAMP.

    Parameters
    ----------
    wfs : np.array of float, shape (n, 3)
        The three waveform channels as recorded with MAMP in XLD1.

    threshold : float
        The voltage that triggers a suspect, in Volts.

    Returns
    -------
    suspects : np.array of in, 1D
        The indices of suspected dust impacts.
    """

    #floor at -0.1V
    wfs = np.maximum(wfs,np.zeros_like(wfs)-0.1)

    #remove median but remember it for later
    medians = np.median(wfs,axis=0)
    wfs -= medians

    #get dust suspects

    suspects = np.arange(len(wfs[:,0]))[(wfs[:,0]>threshold)+
                                        (wfs[:,1]>threshold)+
                                        (wfs[:,2]>threshold)]

    return suspects


def suspects_stat(location = os.path.join("998_generated","mamp_processed",""),
                  target_location = os.path.join("data_synced",""),
                  date_from = dt.datetime(2010,1,1),
                  date_to = dt.datetime(2050,1,1)):
    """
    To get a quick overview of the suspects extracted from MAMP 
    and how they compare to the data from Days.

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
    None.

    """

    all_suspects = load_all_suspects(location,
                                     date_from,
                                     date_to)

    save_list(all_suspects, "all_suspects.pkl" ,target_location)

    all_days = load_all_days()

    YYYYMMDDs = list(set([suspect.YYYYMMDD for suspect in all_suspects]))

    for YYYYMMDD in YYYYMMDDs:
        print(YYYYMMDD)
        day = [day for day in all_days if day.YYYYMMDD == YYYYMMDD][0]
        suspects = [suspect for suspect in all_suspects if suspect.YYYYMMDD == YYYYMMDD]
        confirmed_positive = [suspect for suspect in suspects if suspect.classification == 1]
        confirmed_negative = [suspect for suspect in suspects if suspect.classification == -1]
        print("impacts on Day: "+str(len(day.impact_times)))
        print(f"conf. pos/neg {len(confirmed_positive)}/{len(confirmed_negative)} positive of {len(suspects)} MAMP suspects")



def main(target_input_cdf,
         target_output_pkl,
         threshold = 0.05):
    """
    The main routine, takes the given MAMP cdf, finds all the dust suspects 
    based on the amplitudes and the threshold. 

    Parameters
    ----------
    target_input_cdf : TYPE
        DESCRIPTION.
    target_output_pkl : TYPE
        DESCRIPTION.
    threshold : float
        The voltage that triggers a suspect, as in get_mamp_suspects().

    Returns
    -------
    None.

    """

    cdf_file = cdflib.CDF(target_input_cdf)
    wfs = cdf_file.varget("WAVEFORM_DATA")[:,0:3]
    channel_refs = cdf_file.varget("CHANNEL_REF")
    sampling_rates = cdf_file.varget("sampling_rate")

    suspects = get_mamp_suspects(wfs,threshold=threshold)
    epochs = cdf_file.varget("Epoch")
    suspect_epochs = epochs[suspects]
    suspect_datetimes = tt2000_to_date(suspect_epochs)

    days = load_all_days()
    dust_times = [day.impact_times for day in days]
    flat_dust_times = np.array([item for sublist in dust_times for item in sublist])
    non_dust_times = [day.non_impact_times for day in days]
    flat_non_dust_times = np.array([item for sublist in non_dust_times for item in sublist])

    confirmed_classification = is_confirmed(flat_dust_times,
                                            flat_non_dust_times,
                                            suspect_datetimes)

    suspect_events = []

    for i, datetime in enumerate(suspect_datetimes):
        channel_ref = channel_refs[i]
        sampling_rate = sampling_rates[i]
        period = min([epochs[i]-epochs[i-1],epochs[i+1]-epochs[i]])
        classification = confirmed_classification[i]
        amplitudes = wfs[i,:]

        suspect = ImpactSuspect(datetime,
                                sampling_rate,
                                period,
                                i,
                                classification,
                                channel_ref,
                                amplitudes)

        suspect_events.append(suspect)

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


