import sys
sys.path.insert(0, 'C:\\Users\\skoci\\Documents\\dust\\000_commons')
sys.path.insert(0, 'C:\\Users\\skoci\\Documents\\dust\\003_solar_orbiter')
import cdflib


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


def main(target_input_cdf,
         target_output_stats):
    
    
    try:            
        cdf_file_e = cdflib.CDF(target_input_cdf)
        
        e = cdf_file_e.varget("WAVEFORM_DATA_VOLTAGE") #[event,channel,time], channel = 2 is monopole in XLD1
        epoch = cdf_file_e.varget("EPOCH")              #start of each impact
        epoch_offset = cdf_file_e.varget("EPOCH_OFFSET")     #time variable for each sample of each impact
        sw = cdf_file_e.attget("Software_version",entry=0)["Data"]    #should be 2.1.1 or higher
        
        #tbd check the SW version
        #tbd check the XLD1 'CHANNEL_REF'
        #tbd check the length of E ~= 16384 or double that
        #tbd check 'SAMPLING_RATE' /1e5 ~= 2*2.62
        #tbd store which version do we have and what has to be done
        
        
    except:
        stats(target_output_stats,"error at "+target_input_cdf)
    else:
        
        #TBD padding of the data, according to what has to be done
        #TBD save the data as .txt, resampled to 16384
        
        
        
        stats(target_output_stats,"TBD")
   









   
#main(sys.argv[1],sys.argv[2])    