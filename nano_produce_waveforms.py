import numpy as np
import sys
sys.path.insert(0, 'C:\\Users\\skoci\\Documents\\dust\\000_commons')
sys.path.insert(0, 'C:\\Users\\skoci\\Documents\\dust\\003_solar_orbiter')
import cdflib
import csv
import glob
import matplotlib.pyplot as plt
from scipy import signal as signal



def stats(targetFileName,
          message):
    """
    Save a single string into the txt file of the given name. 

    Parameters
    ----------
    targetFileName : str
        filepath and the name of the target file
    message : str
        the debug message to be written

    Returns
    -------
    None.

    """
    with open(targetFileName, 'w') as f:
        f.write(message)
        f.close()


def pad(wf,where_to_start):
    """
    A function to pad the given array and return an array of double the lenght. 

    Parameters
    ----------
    wf : array of number
        An array to be padded.
    where_to_start : int
        There in the array of double length will the original data start.

    Returns
    -------
    wf_padded : array of number
        The padded array.
    """
    
    #noise background
    wf_padded = np.random.normal(np.mean(wf),np.std(wf),2*len(wf))
    
    #softened signal
    softening_mask = np.hstack((np.arange(0,1,0.1),
                                np.ones(len(wf)-20),
                                np.arange(0.9,-0.1,-0.1)))
    wf_soft = wf*softening_mask
    
    #inverse softening mask to cut out a part of the noise
    mask = range(where_to_start,where_to_start+len(wf))
    wf_padded[mask] *= 1-softening_mask
    
    #sum
    wf_padded[mask] += wf_soft
    
    return wf_padded
    
    
def subsample(wf):
    return signal.resample(wf,16384)
    

def main(target_input_cdf,
         target_output_stats):
    
    #path = "C:\\Users\\skoci\\Disk Google\\000 škola\\UIT\\getting data\\solo\\rpw\\tds_wf_e"
    #cdfs = glob.glob(path+"\\*L2*tswf-e*.cdf")
    
    YYYYMMDD = target_input_cdf[-16:-8]
    
    try:            
        cdf_file_e = cdflib.CDF(target_input_cdf)
        
        e = cdf_file_e.varget("WAVEFORM_DATA_VOLTAGE") #[event,channel,time], channel = 2 is monopole in XLD1
        epoch = cdf_file_e.varget("EPOCH")              #start of each impact
        epoch_offset = cdf_file_e.varget("EPOCH_OFFSET")     #time variable for each sample of each impact
        sw = cdf_file_e.attget("Software_version",entry=0)["Data"]    #should be 2.1.1 or higher
        sw = int(sw.replace(".",""))
        channel_ref = cdf_file_e.varget("CHANNEL_REF")  #should be xld1-like
        sampling_rate = cdf_file_e.varget("SAMPLING_RATE")    
        quality_fact = cdf_file_e.varget("QUALITY_FACT")
        
    except:
        stats(target_output_stats,"error at "+target_input_cdf)
        
    else:
        file_report = ""
        #for each waveform
        for i in range(len(e)):
            #only if data is XLD1
            if sw>=211 and min(channel_ref[i]==[13, 21, 20])==1:
                
                #assigned dust / no dust:
                if quality_fact[i] == 65535:
                    folder = "C:\\Users\\skoci\\Documents\\nanodust\\998_generated\\waveforms\\TDS_dust\\"
                else:
                    folder = "C:\\Users\\skoci\\Documents\\nanodust\\998_generated\\waveforms\\TDS_other\\"
                
                #decide on the treatment
                if len(e[i][0])==16384 and np.isclose(262137.5,
                                                      sampling_rate[i],
                                                      rtol=1e-6,atol=5e-1):
                    treatment = "none"
                elif len(e[i][0])==16384 and np.isclose(524275.0,
                                                      sampling_rate[i],
                                                      rtol=1e-6,atol=5e-1):
                    treatment = "pad, subsample"
                elif len(e[i][0])==32768 and np.isclose(524275.0,
                                                      sampling_rate[i],
                                                      rtol=1e-6,atol=5e-1):
                    treatment = "subsample"
                else:
                    treatment = "report: " + str(len(e[i][0])) + str(sampling_rate[i])
            
                ch0 = e[i][0]
                ch1 = e[i][1]
                ch2 = e[i][2]
                
                #padding with noise
                if "pad" in treatment:
                    # the three channels are supposed to be padded symmetrically
                    start_orig_data = np.random.choice(np.arange(0,len(e[i][0])))
                    
                    ch0 = pad(ch0,start_orig_data)
                    ch1 = pad(ch1,start_orig_data)
                    ch2 = pad(ch2,start_orig_data)
                
                #resampling to 16384 samples
                if "subsample" in treatment:
                    
                    ch0 = subsample(ch0)
                    ch1 = subsample(ch1)
                    ch2 = subsample(ch2)
    
                #save the processed data as .csv
                if "subsample" in treatment:
                    with open(folder+YYYYMMDD+"_"+str(i).zfill(4)+".csv", 'w', encoding='UTF8') as f:
                        writer = csv.writer(f)
                
                        # write the data row by row
                        for j in range(len(ch0)):
                            writer.writerow((ch0[j],ch1[j],ch2[j]))
            
            else:
                treatment = "report: " + str(sw) + str(channel_ref[i])
                
            #add line of the stats file
            file_report += str(i).zfill(4)+" "+treatment+"\n"
        
        stats(target_output_stats,file_report)
        
        # would bbe nice to add a total report, but for that we 
        # need to aggregate the report files, I need that to search globally


main(sys.argv[1],sys.argv[2])    