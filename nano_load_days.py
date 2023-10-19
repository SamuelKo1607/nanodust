import numpy as np
import pandas as pd
import datetime as dt
import glob
import cdflib
import pickle
from scipy.signal import butter
from scipy.signal import sosfilt
from scipy import stats

from conversions import tt2000_to_date
from conversions import YYYYMMDD2jd
from conversions import YYYYMMDD2date
from nano_ephemeris import fetch_heliocentric_solo

from paths import cdf_stat_location
from paths import cdf_tswf_e_location
from paths import lo_f_cat_location
from paths import lo_f_cat_new_location
from paths import hi_f_cat_location
from paths import solo_ephemeris_file



class Impact:
    """
    Pretty self explanatory, each instance holds the attributes of one assumed 
    dust impact. 
    """
    def __init__(self, datetime, sampling_rate, amplitude, index):
        self.datetime = datetime
        self.YYYYMMDD = datetime.strftime('%Y%m%d')
        self.sampling_rate = sampling_rate
        self.amplitude = amplitude
        self.index = index
        self.produced = dt.datetime.now()


    def info(self):
        print(self.datetime.strftime("%d.%m.%Y %H:%M:%S")+
              " \n sampling rate: "+ str(self.sampling_rate)+
              " \n amplitude: " + str(self.amplitude))

    def show(self):
        print(vars(self))


class Day:
    """
    Aggreaget results of the measurement day, plus SolO and RPW status.
    """
    def __init__(self,
                 date,
                 impact_times,
                 non_impact_times,
                 duty_hours,
                 sampling_rate,
                 heliocentric_distance,
                 spacecraft_speed,
                 heliocentric_radial_speed):
        self.date = date
        self.YYYYMMDD = date.strftime('%Y%m%d')
        self.impact_count = len(impact_times)
        self.impact_times = impact_times
        self.non_impact_times = non_impact_times
        self.duty_hours = duty_hours
        self.sampling_rate = sampling_rate
        self.heliocentric_distance = heliocentric_distance
        self.spacecraft_speed = spacecraft_speed
        self.heliocentric_radial_speed = heliocentric_radial_speed
        self.produced = dt.datetime.now()


    def info(self):
        print(self.YYYYMMDD+
              " \n impact count: "+ str(self.impact_count)+
              " \n duty hours: " + str(self.duty_hours))






def save_list(data,name,location=""):
    """
    A simple function to save a given list to a specific location using pickle.

    Parameters
    ----------
    data : list
        The data to be saved. In our context: mostly a list of Impact objects.
    
    name : str
        The name of the file to be written. 
        
    location : str, optional
        The relative path to the data folder. May not be used if the name
        contains the folder as well. Default is therefore empty.

    Returns
    -------
    none
    """
    
    with open(location+name, "wb") as f:  
        pickle.dump(data, f)
        
        
def load_list(name,location):
    """
    A simple function to load a saved list from a specific location using pickle.

    Parameters
    ----------    
    name : str
        The name of the file to load. 
        
    location : str, optional
        The relative path to the data folder.

    Returns
    -------
    data : list
        The data to be loaded. In our context: mostly a list of Impact objects.
    """
    
    with open(location+name, "rb") as f:
        data = pickle.load(f)
    return data


def load_all_days(days_location = "998_generated\\days\\"):
    """
    The function to load all the measurement days data.

    Parameters
    ----------
    days_location : str, optional
        The data directory. Default is "998_generated\\days\\".

    Returns
    -------
    days : list of Day object
        Measurement days, class Day from nano_load_days.
    """
    files = glob.glob(days_location+"*.pkl")
    days = []
    for file in files:
        day = load_list(file,"")
        days.append(day)

    return days


def load_all_impacts(impacts_location = "998_generated\\impacts\\",
                     date_from = dt.datetime(2010,1,1),
                     date_to = dt.datetime(2050,1,1)):
    """
    The function to load all the impacts data.

    Parameters
    ----------
    impacts_location : str, optional
        The data directory. Default is "998_generated\\impacts\\".
    date_from : dt.datetime
        The first relevant moment (filtering the instances by >=).
    date_to : dt.datetime
        The last relevant moment (filtering the instances by <=).

    Returns
    -------
    impacts : list of Impact object
        Impact data points, class Impact from nano_load_days.
    """
    files = glob.glob(impacts_location+"*.pkl")
    impacts = []
    for file in files:
        impact = load_list(file,"")
        impacts.append(impact)

    flat_list_of_impacts = [item for sublist in impacts for item in sublist]

    filtered = [i for i in flat_list_of_impacts if date_from <= i.datetime <= date_to ]

    return filtered




def extract_variables_from_days(days):
    """
    The function to extract the relevant attributes from the provided days.

    Parameters
    ----------
    days : list of Day object
        Measurement days, class Day from nano_load_days.

    Returns
    -------
    dates : np.array of float
        The measuremetn days.
    counts : np.array of int
        The dust counts, enocutered on the dates.
    duty_hours : np.array of float
        The duty cycle in hours, per day. Often around 1.5.
    sampling_rates : np.array of float
        The sampling rate on that day. No dates with a vraiable 
        sampling rate were processed in the first place.
    """
    dates = np.zeros(0,dtype=dt.datetime)
    counts = np.zeros(0,dtype=int)
    duty_hours = np.zeros(0,dtype=float)
    sampling_rates = np.zeros(0,dtype=float)
    for day in days:
        dates = np.append(dates,day.date)
        counts = np.append(counts,day.impact_count)
        duty_hours = np.append(duty_hours, day.duty_hours)
        sampling_rates = np.append(sampling_rates, day.sampling_rate)

    return dates, counts, duty_hours, sampling_rates


def get_errors(days, prob_coverage = 0.9):
    """
    The function to calculate the errorbars for flux 
    assuming Poisson distribution and taking into account
    the number of detection.

    Parameters
    ----------
    days : list of Day object
        Measurement days, class Day from nano_load_days.
    prob_coverage : float, optional
        The coverage of the errobar interval. The default is 0.9.

    Returns
    -------
    err_plusminus_flux : np.array of float
        The errorbars, lower and upper bound, shape (2, n).

    """
    dates, counts, duty_hours, sampling_rates = extract_variables_from_days(days)

    counts_err_minus = -stats.poisson.ppf(0.5-prob_coverage/2, mu=counts)+counts
    counts_err_plus  = +stats.poisson.ppf(0.5+prob_coverage/2, mu=counts)-counts
    err_plusminus_flux = np.array([counts_err_minus,counts_err_plus]) / (duty_hours/(24))

    return err_plusminus_flux


def get_cdfs_to_analyze(cdf_files_directory,
                        search_mask = "*rpw-tds-surv-tswf-e*.cdf"):
    """
    Returns the list of triggered snapshot E-field 
    cdfs from the given directory.

    Parameters
    ----------
    cdf_files_directory : str
        The directory.

    Returns
    -------
    cdfs : list of str
        List of filepaths to the files.

    """
    cdfs = glob.glob(cdf_files_directory+search_mask)
    return cdfs


def hipass_correction(signal,epochstep,cutoff=370):
    """
    Correct the time-series for a high-pass filter of X Hz 
    through Laplace domain as Amalia did.

    Parameters
    ----------
    signal : np.array
        The wavefrom to be corrected.
    epochstep : float
        A quantum of time in s.
    cutoff : float, optional
        The cutoff frequency in Hz. The default is 370.

    Returns
    -------
    corrected : np.array
        The corrected wavefrom.
    """
    cumsignal = np.cumsum(signal)
    corrected = signal + (2*np.pi*cutoff)*(cumsignal*epochstep)
    corrected -= np.mean(corrected)
    return corrected


def event_filter(wf,time_step):
    """
    Psacks all the filtering that is done to 
    the unprocessed "WAVEFORM_DATA_VOLTAGE" signal before evaluating
    properties of the dust impact.

    Parameters
    ----------
    wf : np.array
        The wavefrom to be corrected.
    time_step : float
        A quantum of time in ns.

    Returns
    -------
    corrected : np.array
        The corrected wavefrom.
    """
    bg = np.mean(wf)
    sos = butter(32,7e4,btype="lowpass",fs=(1/(time_step/1e9)),output="sos")
    smooth = sosfilt(sos, wf)-bg
    #smooth = np.append(smooth[10:],10*[smooth[-1]])
    corrected = hipass_correction(smooth,time_step/1e9)
    return corrected


def get_duty_hours(YYYYMMDD,location):
    """
    Gets the duty hours of RPW triggered TDS on a given day

    Parameters
    ----------
    YYYYMMDD : str
        The date of interest.
    location : str
        folder of the _stats_ files.

    Raises
    ------
    Exception
        If no file is found.

    Returns
    -------
    duty_hours : float
        Number of duty hours, usually around 1.5.

    """
    stats_file_name_pattern = (location+
                             "solo_L2_rpw-tds-surv-stat_"+
                             YYYYMMDD+"*.cdf")
    stats_file_name = glob.glob(stats_file_name_pattern)
    if len(stats_file_name) == 0:
        raise Exception("no stats file for "+YYYYMMDD)
    else:
        stats_cdf_file = cdflib.CDF(stats_file_name[0])
        snapshot_len = np.average(stats_cdf_file.varget("SNAPSHOT_LEN"))
        nr_events = stats_cdf_file.varget("SN_NR_EVENTS")
        sampling = stats_cdf_file.varget("SAMPLING_RATE")
        duty_hours = snapshot_len*sum(nr_events/sampling/3600)
        return duty_hours


def find_amplitude(smooth_1,
                   smooth_2,
                   smooth_3):
    """
    Find the primary peak (= body peak) of the impact 
    based on the three waveforms.

    Parameters
    ----------
    smooth1 : np.array
        Filtered waveform V1.
    smooth2 : np.array
        Filtered waveform V2.
    smooth3 : np.array
        Filtered waveform V3.

    Returns
    -------
    amplitude : float
        The primary peak amplitude.
    """

    max1 = np.max(smooth_1)
    max2 = np.max(smooth_2)
    max3 = np.max(smooth_3)

    amplitude = np.min([max1,max2,max3])

    return amplitude


def process_impact(cdf_file, i):
    """
    Takes a CDF and an index of a TDS event that is believed to be dust 
    and constructs the Impact class object accordingly.

    Parameters
    ----------
    cdf_file : cdflib.cdfread.CDF
        A loaded CDF file containing TDS event.
    i : int
        An index of a dust impact.

    Returns
    -------
    impact : Impact
        In instance of Impact class.
    """

    e = cdf_file.varget("WAVEFORM_DATA_VOLTAGE")[i,:,:]
    sampling_rate = cdf_file.varget("SAMPLING_RATE")[i]
    date_time = tt2000_to_date(cdf_file.varget("EPOCH")[i])
    epoch_offset = cdf_file.varget("EPOCH_OFFSET")

    time_step = epoch_offset[i][1]-epoch_offset[i][0]

    monopole_1 = e[2,:]-e[1,:]
    monopole_2 = e[2,:]
    monopole_3 = e[2,:]-e[1,:]-e[0,:]

    smooth_1 = event_filter(monopole_1, time_step)
    smooth_2 = event_filter(monopole_2, time_step)
    smooth_3 = event_filter(monopole_3, time_step)

    amplitude = find_amplitude(smooth_1,smooth_2,smooth_3)

    impact = Impact(date_time, sampling_rate, amplitude, i)

    return impact


def get_solo_state(YYYYMMDD, ephem_file):
    """
    The function to get the heliocentric distance and the velocity of SolO.

    Parameters
    ----------
    YYYYMMDD : str
        The date of interest.
    ephem_file : str
        The location of the ephemereis file.

    Returns
    -------
    r : float
        heliocentric distance in AU.
    v : float
        total speed in km/s.
    v_rad : float
        radial speed in km/s.

    """
    f_hel_r, f_hel_phi, f_rad_v, f_tan_v = fetch_heliocentric_solo(
                                                solo_ephemeris_file)
    jd = YYYYMMDD2jd(YYYYMMDD)
    r = float(f_hel_r(jd))
    v_rad = float(f_rad_v(jd))
    v_tan = float(f_tan_v(jd))
    v = float(np.sqrt(v_rad**2+v_tan**2))

    return r, v, v_rad


def process_cdf(cdf_file):
    """
    The function to process one CDF file. Depending on whther it is 
    the higher or the lower sampling rate, the correct decision dust/nodust 
    is done. Different approach is taken for CNN files 
    before 2022 11 24 (incl) and after. List of Impact instances 
    and the list of missing files are returned.

    Parameters
    ----------
    cdf_file : cdf_file : cdflib.cdfread.CDF
        A loaded CDF file to be processed.

    Raises
    ------
    Exception
        - If the CDF file firmat is not acceptable.
        - If the CNN classification file is empty

    Returns
    -------
    impacts : list of Impact 
        The impacts from that day.
    missing_files : list of string
        The files that were expected but were not found.
    """
    e = cdf_file.varget("WAVEFORM_DATA_VOLTAGE")
    epochs = cdf_file.varget("EPOCH")
    quality_fact = cdf_file.varget("QUALITY_FACT")
    channel_ref = cdf_file.varget("CHANNEL_REF")
    sampling_rate = cdf_file.varget("SAMPLING_RATE")
    sw = cdf_file.attget("Software_version",entry=0).Data
    YYYYMMDD = str(cdf_file.file)[-16:-8]

    events_count = np.shape(e)[0]
    impacts = []
    missing_files = []

    all_close = np.all(np.isclose(sampling_rate, sampling_rate[0]))
    is_xld = np.all(np.isclose(channel_ref[:,:3],[13, 21, 20]))
    is_current = sw in ["2.1.1"]

    if all_close and is_xld and is_current and np.isclose(sampling_rate[0],
                                                          262137.5):
        #lower sampling rate, process and append impacts with Impact object

        if YYYYMMDD2date(YYYYMMDD) <= YYYYMMDD2date("20221124"):
            #legacy
            cnn_file_name_pattern = (lo_f_cat_location+
                                     "solo_L2_rpw-tds-surv-tswf-e_"+
                                     str(cdf_file.file)[-16:-6]+"*.txt")
            cnn_file_name = glob.glob(cnn_file_name_pattern)
            if len(cnn_file_name) == 0:
                missing_files += cnn_file_name_pattern
            else:
                cnn_classification = pd.read_csv(cnn_file_name[0])
                plain_index = np.arange(events_count)
                flagged_dust = quality_fact==65535
                reordered_index = np.append(plain_index[flagged_dust==1],plain_index[flagged_dust==0])
                dust = np.array(cnn_classification["Label"])
                dust = np.array(dust, dtype=bool)
                indices_analyzed = np.array(cnn_classification["Index"])[dust]-1
                try:
                    indices = reordered_index[indices_analyzed]
                except:
                    raise Exception("non-classified data @ "+
                                    str(cdf_file.file)[-16:-4])
                else:
                    for i in indices:
                        impact = process_impact(cdf_file, i)
                        impacts.append(impact)

        else:
            #new version as of oct/2023
            cnn_files_pattern = (lo_f_cat_new_location+
                                     "solo_L2_rpw-tds-surv-tswf-e_"+
                                     YYYYMMDD+"*.txt")
            cnn_files = glob.glob(cnn_files_pattern)
            if len(cnn_files) == 0:
                missing_files += cnn_files_pattern
            indices_analyzed = np.zeros(0,dtype=int)
            for cnn_file in cnn_files:
                indices_analyzed = np.append(indices_analyzed, int(cnn_file[-7:-4]))
            indices_analyzed -= 1 #because matlab

            plain_index = np.arange(events_count)
            flagged_dust = quality_fact==65535
            reordered_index = np.append(plain_index[flagged_dust==1],plain_index[flagged_dust==0])
            try:
                indices = reordered_index[indices_analyzed]
            except:
                raise Exception("non-classified data @ "+
                                str(cdf_file.file)[-16:-4])
            else:
                for i in indices:
                    impact = process_impact(cdf_file, i)
                    impacts.append(impact)


    elif all_close and is_xld and is_current and np.isclose(sampling_rate[0],
                                                            524275):
        #higher sampling rate, process and append impacts with Impact object

        for i in range(events_count):
            cnn_cat_filename = YYYYMMDD+"_"+str(i).zfill(4)+".txt"
            try:
                cnn_cat_data = pd.read_csv(hi_f_cat_location+cnn_cat_filename)
            except:
                missing_files.append(cnn_cat_filename)
            else:
                cnn_cat_dust = cnn_cat_data["Label"]>0.5
                if cnn_cat_dust[0]:
                    impact = process_impact(cdf_file, i)
                    impacts.append(impact)
                else:
                    pass #not an impact


    else:
        raise Exception("Undefined data @ "+
                        str(cdf_file.file)[-16:-4])

    #construct day
    try:
        date = dt.datetime.strptime(YYYYMMDD, "%Y%m%d").date()
        all_suspect_times = [tt2000_to_date(epoch) for epoch in epochs]
        impact_times=[i.datetime for i in impacts]
        non_impact_times = list(set(all_suspect_times) - set(impact_times))
        duty_hours = get_duty_hours(YYYYMMDD, cdf_stat_location)
        r, v, v_rad = get_solo_state(YYYYMMDD, solo_ephemeris_file)


    except Exception as err:
        print("day construction failed @ "+YYYYMMDD)
        raise err
    else:
        day = Day(date,
                  impact_times,
                  non_impact_times,
                  duty_hours,
                  sampling_rate[0],
                  r,
                  v,
                  v_rad)

    return day, impacts, missing_files


def main():
    all_missing_files = []
    for file in get_cdfs_to_analyze(cdf_tswf_e_location):
        cdf_file = cdflib.CDF(file)
        short_name = file[-44:-4]
        try:
            day, impacts, missing_files = process_cdf(cdf_file)
        except Exception as err:
            print(short_name)
            print(err)
        else:
            save_list(impacts,short_name+"_impacts.pkl","998_generated\\impacts\\")
            save_list(day,short_name+"_day.pkl","998_generated\\days\\")
            all_missing_files.append(missing_files)
    return all_missing_files



#missing_files = main()
