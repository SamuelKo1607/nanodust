import numpy as np
import pandas as pd
import glob
import cdflib
import pickle

from conversions import tt2000_to_date

from keys import cdf_tswf_e_location
from keys import lo_f_cat_location
from keys import hi_f_cat_location




class Impact:
    def __init__(self, datetime, sampling_rate, amplitude):
        self.datetime = datetime
        self.sampling_rate = sampling_rate
        self.amplitude = amplitude


    def info(self):
        print(self.datetime.strftime("%d.%m.%Y %H:%M:%S")+
              " @ sampling rate "+ self.sampling_rate+
              " - amplitude "*self.amplitude)







def save_list(data,name,location):
    """
    A simple function to save a given list to a specific location using pickle.

    Parameters
    ----------
    data : list
        The data to be saved. In our context: mostly a list of Impact objects.
    
    name : str
        The name of the file to be written. 
        
    location : str
        The relative path to the data folder.

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


def download_available_cdf(strat_date, end_date):
    #TODO, use Soar this time, please
    pass


def get_cdfs_to_analyze(cdf_files_directory):
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
    cdfs = glob.glob(cdf_files_directory+"*rpw-tds-surv-tswf-e*.cdf")
    return cdfs


def process_impact(cdf_file, i):

    e = cdf_file.varget("WAVEFORM_DATA_VOLTAGE")[i,:,:]
    sampling_rate = cdf_file.varget("SAMPLING_RATE")[0]
    date_time = tt2000_to_date(cdf_file.varget("EPOCH")[0])

    monopole_1 = e[2,:]-e[1,:]
    monopole_2 = e[2,:]
    monopole_3 = e[2,:]-e[1,:]-e[0,:]

    amplitude = 0 #TODO

    impact = Impact(date_time, sampling_rate, amplitude)

    return impact



def process_cdf(cdf_file):
    e = cdf_file.varget("WAVEFORM_DATA_VOLTAGE")
    quality_fact = cdf_file.varget("QUALITY_FACT")
    channel_ref = cdf_file.varget("CHANNEL_REF")
    sampling_rate = cdf_file.varget("SAMPLING_RATE")

    events_count = np.shape(e)[0]
    impacts = []

    all_close = np.all(np.isclose(sampling_rate, sampling_rate[0]))
    is_xld = np.all(np.isclose(channel_ref,[13, 21, 20]))

    if all_close and is_xld and np.isclose(262137.5, sampling_rate[0]):
        #lower sampling rate, process and append impacts with Impact object


        pass

    elif all_close and is_xld and np.isclose(524275, sampling_rate[0]):
        #higher sampling rate, process and append impacts with Impact object
        for i in range(events_count):
            TDS_cat_dust = quality_fact[i] == 65535
            if TDS_cat_dust:
                YYYYMMDD = str(cdf_file.file)[-16:-8]
                cnn_cat_filename = YYYYMMDD+"_"+str(i).zfill(4)+".txt"
                cnn_cat_data = pd.read_csv(hi_f_cat_location+cnn_cat_filename)
                cnn_cat_dust = cnn_cat_data["Label"]>0.5
                if cnn_cat_dust:
                    impact = process_impact(cdf_file, i)
                else:
                    pass
            else:
                pass #TODO when we have data we can do those as well

    else:
        raise Exception("Undefined data @ "+
                        str(cdf_file.file)[-16:-4])

    return impacts


def main():
    for file in get_cdfs_to_analyze(cdf_tswf_e_location):
        cdf_file = cdflib.CDF(file)
        short_name = file[-44:-4]
        try:
            impacts = process_cdf(cdf_file)
        except Exception as err:
            print(err)
        else:
            save_list(impacts,short_name+".pkl","998_generated\\impacts\\")





#main()


"""
The rough plan:
    - go through the downloaded cdfs
    - for each decide if it is low or high sampling rate
        - if low sampling rate, use legacy classification files
        - if high sampling rate, use new classification files
    - for each impact, fill the class impact with data 
"""