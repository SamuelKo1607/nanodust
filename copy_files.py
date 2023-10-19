import shutil
import cdflib
import glob
import numpy as np
from paths import cdf_tswf_e_location
from nano_load_days import get_cdfs_to_analyze


def main(destination, sampling_rate = 524275):
    files = get_cdfs_to_analyze(cdf_tswf_e_location)
    for file in files:
        cdf_file = cdflib.CDF(file)
        try:
            file_sampling_rate = cdf_file.varget("sampling_rate")[0]
            if np.isclose(file_sampling_rate, sampling_rate):
                if len(glob.glob(destination+"*"+file[-16:]))>0:
                    pass
                else:
                    shutil.copy(file, destination)
        except:
            pass


#main("D:\\dust_cdfs_high_sampling_rate\\")