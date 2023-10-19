import numpy
import cdflib

from nano_load_days import get_cdfs_to_analyze

from paths import cdf_stat_location
from paths import cdf_tswf_e_location





cdfs_stat = get_cdfs_to_analyze(cdf_stat_location,
                                search_mask = "*rpw*stat*.cdf")
cdfs_tswfe = get_cdfs_to_analyze(cdf_tswf_e_location)

matching = 0
all_dust_downloaded = 0
too_many_dust_downloaded = 0
too_few_dust_downloaded = 0
for cdf_stat in cdfs_stat:
    YYYYMMDD = cdf_stat[cdf_stat.find("surv-stat_")+10:-8]
    matching_tswfe = [s for s in cdfs_tswfe if YYYYMMDD in s]
    if len(matching_tswfe)>0:
        try:
            matching += 1
            cdf_tswfe = matching_tswfe[0]
    
            loaded_stat = cdflib.CDF(cdf_stat)
            stat_cdf_epochs = loaded_stat.varget("epoch")
            dust_count_stat = sum(loaded_stat.varget("DU_NR_IMPACT"))
            wa_count_stat = sum(loaded_stat.varget("WA_NR_EVENTS"))

    
            loaded_tswfe = cdflib.CDF(cdf_tswfe)
            tswfe_cdf_epochs = loaded_tswfe.varget("epoch")
            dust_count_tswfe = sum(loaded_tswfe.varget("quality_fact")==65535)
            wa_count_tswfe = sum(loaded_tswfe.varget("quality_fact")!=65535)


    
            print("#################")
            print(YYYYMMDD)
            print("wfs downloaded total:"+str(len(tswfe_cdf_epochs)))
            print("stat dust, waves: "+str(dust_count_stat)+", "+str(wa_count_stat))
            print("  wf dust, waves: "+str(dust_count_tswfe)+", "+str(wa_count_tswfe))
            if dust_count_stat == dust_count_tswfe:
                all_dust_downloaded += 1
            elif dust_count_stat < dust_count_tswfe:
                too_many_dust_downloaded += 1
            else:
                too_few_dust_downloaded += 1
        except:
            print("#################")
            print(YYYYMMDD)
            print("exception")
    else:
        print("#################")
        print(YYYYMMDD)
        print("no tswfe par for stat")

print("#################")
print("numer of stat files: "+str(len(cdfs_stat)))
print("numer of tswf files: "+str(len(cdfs_tswfe)))
print("all dust downloaded:")
print(str(all_dust_downloaded)+" out of "+str(matching))
print("too few dust downloaded:")
print(str(too_few_dust_downloaded)+" out of "+str(matching))
print("too many dust downloaded:")
print(str(too_many_dust_downloaded)+" out of "+str(matching))









