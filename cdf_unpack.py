import os
import numpy as np
import glob
import cdflib
import matplotlib.pyplot as plt
import matplotlib as mpl
from paths import cdf_tswf_e_location
from paths import unpacked_waveforms
from conversions import YYYYMMDD2date
from conversions import tt2000_to_date
from nano_load_days import event_filter

import figure_standards as figstd
axes_size = figstd.set_rcparams_dynamo(mpl.rcParams, num_cols=1, ls='thin')
mpl.rcParams['figure.dpi']= 200


def plot(epoch,
         e1,e2,e3,
         dust=False,
         name="filename",
         save=True):

    #fig = plt.figure(constrained_layout=True)
    #gs = fig.add_gridspec(3,1,hspace=0.1)
    #ax = gs.subplots()

    fig = plt.figure(figsize=(6,4))
    ax = fig.subplots(3, 1, sharex=True, gridspec_kw=dict(hspace=0.1))

    time_step = epoch[1]-epoch[0]
    time = epoch/1e6
    e1 -= np.mean(e1)
    e2 -= np.mean(e2)
    e3 -= np.mean(e3)
    dynamic_range = np.max(np.abs(np.vstack((e1,e2,e3))))
    smooth1 = event_filter(e1, time_step)
    smooth2 = event_filter(e2, time_step)
    smooth3 = event_filter(e3, time_step)

    if dust:
        color = "black"
    else:
        color = "red"

    ax[0].plot(time,e1,color="red",lw=0.5,alpha=0.3)
    ax[1].plot(time,e2,color="red",lw=0.5,alpha=0.3)
    ax[2].plot(time,e3,color="red",lw=0.5,alpha=0.3)
    ax[0].plot(time,smooth1,color=color,lw=0.5,alpha=1)
    ax[1].plot(time,smooth2,color=color,lw=0.5,alpha=1)
    ax[2].plot(time,smooth3,color=color,lw=0.5,alpha=1)

    for a in ax:
        a.set_ylim(-dynamic_range,dynamic_range)
    ax[2].set_xlabel(r"Time [$ms$]")
    ax[1].set_ylabel(r"Voltage [$V$]")

    fig.tight_layout()
    if dust:
        filename = os.path.join(unpacked_waveforms,
                                "dust_suspect","")+name+".png"
    else:
        filename = os.path.join(unpacked_waveforms,
                                "non_suspect","")+name+".png"
    if save:
        plt.savefig(filename,dpi=200)
    plt.close()


def unpack(cdf_file):

    YYYYMMDD = str(cdf_file.file)[-16:-8]
    e = cdf_file.varget("WAVEFORM_DATA_VOLTAGE")
    #[event,channel,time], channel = 2 is monopole in XLD1

    epoch = cdf_file.varget("EPOCH")
    #start of each impact

    epoch_offset = cdf_file.varget("EPOCH_OFFSET")
    #time variable for each sample of each impact

    sw = cdf_file.attget("Software_version",entry=0).Data
    #should be 2.1.1 or higher
    sw = int(sw.replace(".",""))

    channel_ref = cdf_file.varget("CHANNEL_REF")
    #should be xld1-like

    sampling_rate = cdf_file.varget("SAMPLING_RATE")
    quality_fact = cdf_file.varget("QUALITY_FACT")

    print("__________")
    print(YYYYMMDD2date(YYYYMMDD).strftime("%B %d, %Y"))
    print(f"{len(quality_fact)} events")
    print(f"{sum(quality_fact==65535)} dust suspects")

    for i in range(len(quality_fact)):
        if min(channel_ref[i]==np.array([13,21,20])):
            mono_1 = e[i,2,:] - e[i,1,:]
            mono_2 = e[i,2,:]
            mono_3 = e[i,2,:] - e[i,1,:] - e[i,0,:]
        elif min(channel_ref[i]==np.array([10,20,30])):
            mono_1 = e[i,0,:]
            mono_2 = e[i,1,:]
            mono_3 = e[i,2,:]
        plot(epoch_offset[i,:],
             mono_1,
             mono_2,
             mono_3,
             quality_fact[i]==65535,
             f"{YYYYMMDD}_{i}")


def main(files=None):

    if files is None:
        files = glob.glob(cdf_tswf_e_location+"\\*L2*tswf-e*.cdf")

    for file in files:
        cdf_file_e = cdflib.CDF(file)
        unpack(cdf_file_e)




#%%
if __name__ == "__main__":
    files_of_interest = [
        '..\\data\\solo\\rpw\\tds_wf_e\\solo_L2_rpw-tds-surv-tswf-e-cdag_20231004_V02.cdf',
        '..\\data\\solo\\rpw\\tds_wf_e\\solo_L2_rpw-tds-surv-tswf-e-cdag_20231005_V02.cdf',
        '..\\data\\solo\\rpw\\tds_wf_e\\solo_L2_rpw-tds-surv-tswf-e-cdag_20231006_V01.cdf']
    main(files = files_of_interest)