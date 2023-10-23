import os

#cdf_stat_location = os.path.join(os.path.normpath(      "..\\data\\solo\\rpw\\tds_stat\\"), '')
cdf_stat_location = os.path.join("..","data","solo","rpw","tds_stat","")
#cdf_tswf_e_location = os.path.join(os.path.normpath(    "..\\data\\solo\\rpw\\tds_wf_e\\"), '')
cdf_tswf_e_location = os.path.join("..","data","solo","rpw","tds_wf_e","")
#cdf_mamp_location = os.path.join(os.path.normpath(      "..\\data\\solo\\rpw\\mamp\\"), '')
cdf_mamp_location = os.path.join("..","data","solo","rpw","mamp","")
#lo_f_cat_location = os.path.join(os.path.normpath(      "..\\data\\solo\\rpw\\Amplitude_data\\Index-Label-MaxAmplitude\\"), '')
lo_f_cat_location = os.path.join("..","data","solo","rpw","Amplitude_data","Index-Label-MaxAmplitude","")
#lo_f_cat_new_location = os.path.join(os.path.normpath(  "..\\data\\solo\\rpw\\Amplitude_data\\low_frequency_new\\Label-MaxAmplitude\\"), '')
lo_f_cat_new_location = os.path.join("..","data","solo","rpw","Amplitude_data","low_frequency_new","Label-MaxAmplitude","")
#hi_f_cat_location = os.path.join(os.path.normpath(      "..\\data\\solo\\rpw\\Amplitude_data\\nanodust_classification\\all_events\\"), '')
hi_f_cat_location = os.path.join("..","data","solo","rpw","Amplitude_data","nanodust_classification","all_events","")
#solo_ephemeris_file = os.path.normpath(    "..\\data\\ephemerides\\solo_fine_ephemeris_noheader.txt")
solo_ephemeris_file = os.path.join("..","data","ephemerides","solo_fine_ephemeris_noheader.txt")
#venus_ephemeris_file = os.path.normpath(   "..\\data\\ephemerides\\venus_fine_ephemeris_noheader.txt")
venus_ephemeris_file = os.path.join("..","data","ephemerides","venus_fine_ephemeris_noheader.txt")


if __name__ == "__main__":
    print(cdf_stat_location)
    print(cdf_mamp_location)
    print(cdf_tswf_e_location)