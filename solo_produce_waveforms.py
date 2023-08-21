import numpy as np
import sys
sys.path.insert(0, 'C:\\Users\\skoci\\Documents\\dust\\000_commons')
sys.path.insert(0, 'C:\\Users\\skoci\\Documents\\dust\\003_solar_orbiter')
import cdflib
import csv


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
        for i in len(e):
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
            else:
                treatment = "report: " + str(sw) + str(channel_ref[i])
            
            #add line of the stats file
            file_report += str(i).zfill(4)+" "+treatment+"\n"    
            
            #padding with noise
            if "pad" in treatment:
                pass #define function
            
            #resampling to 16384 samples
            if "subsample" in treatment:
                pass #define function
            

            #TBD save the processed data as .csv
            with open(folder+YYYYMMDD+"_"+str(i).zfill(4).csv, 'w', encoding='UTF8') as f:
                writer = csv.writer(f)
        
                # write the data row by row
                # for something in something:
                writer.writerow(row)
        
        
        stats(target_output_stats,file_report)
   









   
#main(sys.argv[1],sys.argv[2])    