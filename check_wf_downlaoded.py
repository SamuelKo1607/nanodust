import numpy as np
import cdflib
import datetime as dt
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.signal import savgol_filter

from nano_load_days import get_cdfs_to_analyze
from conversions import YYYYMMDD2date
import figure_standards as figstd

from paths import cdf_stat_location
from paths import cdf_tswf_e_location

axes_size = figstd.set_rcparams_dynamo(mpl.rcParams, num_cols=1, ls='thin')
mpl.rcParams['figure.dpi']= 600




if __name__ == "__main__":
    cdfs_stat = get_cdfs_to_analyze(cdf_stat_location,
                                    search_mask = "*rpw*stat*.cdf")
    cdfs_tswfe = get_cdfs_to_analyze(cdf_tswf_e_location)
    
    matching = 0
    all_dust_downloaded = 0
    hope_to_find_total = np.zeros(0)
    found_total = np.zeros(0)
    dtimes = np.zeros(0,dtype=dt.datetime)
    for cdf_stat in cdfs_stat:
        YYYYMMDD = cdf_stat[cdf_stat.find("surv-stat_")+10:-8]
        matching_tswfe = [s for s in cdfs_tswfe if YYYYMMDD in s]
        if len(matching_tswfe)>0:
            try:
                matching += 1
                cdf_tswfe = matching_tswfe[0]
        
                loaded_stat = cdflib.CDF(cdf_stat)
                stat_cdf_epochs = loaded_stat.varget("epoch")
                dust_events = loaded_stat.varget("DU_NR_IMPACT")
                wave_events = loaded_stat.varget("WA_NR_EVENTS")
                dust_count_stat = sum(dust_events)
                wa_count_stat = sum(wave_events)
                dust_epochs_16s = stat_cdf_epochs[dust_events>0]

        
                loaded_tswfe = cdflib.CDF(cdf_tswfe)
                tswfe_cdf_epochs = loaded_tswfe.varget("epoch")
                dust_count_tswfe = sum(loaded_tswfe.varget("quality_fact")==65535)
                wa_count_tswfe = sum(loaded_tswfe.varget("quality_fact")!=65535)

                #check if for every dust approx epoch there is a downloaded waveform
                hope_to_find = len(dust_epochs_16s)
                found = 0
                for i, epoch in enumerate(dust_epochs_16s):
                    epoch_delays = tswfe_cdf_epochs - epoch
                    pos_epoch_delays = epoch_delays[epoch_delays>0]
                    nearest = min(pos_epoch_delays)
                    if nearest/1000000000 < 16.5:
                        #found
                        found += 1
                    else:
                        #not found
                        pass
                hope_to_find_total = np.append(hope_to_find_total,hope_to_find)
                found_total = np.append(found_total,found)
                dtimes = np.append(dtimes,YYYYMMDD2date(YYYYMMDD))

                print("#################")
                print(YYYYMMDD)
                print("wfs downloaded total:"+str(len(tswfe_cdf_epochs)))
                print("stat dust, waves: "+str(dust_count_stat)+", "+str(wa_count_stat))
                print("  wf dust, waves: "+str(dust_count_tswfe)+", "+str(wa_count_tswfe))
                print(f"{found} dusts out of {hope_to_find} were matched")
                if found == hope_to_find:
                    all_dust_downloaded+=1

            except:
                print("#################")
                print(YYYYMMDD)
                print("exception")
        else:
            print("#################")
            print(YYYYMMDD)
            print("no tswfe match for stat")
    
    print("#################")
    print("numer of stat files: "+str(len(cdfs_stat)))
    print("numer of tswf files: "+str(len(cdfs_tswfe)))
    print("all dust downloaded "+str(all_dust_downloaded)+
          " out of "+str(matching))
    print("dusts downloaded "+str(found_total)+
          " out of "+str(hope_to_find_total))

    fig,ax = plt.subplots()
    ax.plot(dtimes,hope_to_find_total,label="STAT")
    ax.plot(dtimes,found_total,label="TSWF_E")
    ratio_smoothed = savgol_filter(100*found_total/hope_to_find_total, 11, 2)
    ax.plot(dtimes,ratio_smoothed,label="ratio %")
    ax.legend(fontsize="small")
    ax.tick_params(axis='x',labelrotation=60)
    ax.set_ylim(0,100)
    fig.show()


    
    

    
    
    
    
    
