import matplotlib.pyplot as plt
import matplotlib as mpl
import glob
import numpy as np
import datetime as dt
from scipy import interpolate
from scipy.signal import argrelextrema
from scipy.signal import savgol_filter

from nano_load_days import load_all_impacts
from nano_load_days import load_all_days
from nano_load_days import load_list
from nano_load_days import save_list
from nano_load_days import Impact
from nano_load_days import Day
from nano_ephemeris import load_hae
from conversions import hae2vse
from conversions import date2jd
from conversions import jd2date
from conversions import YYYYMMDD2jd
import figure_standards as figstd

axes_size = figstd.set_rcparams_dynamo(mpl.rcParams, num_cols=1, ls='thin')
mpl.rcParams['figure.dpi'] = 600

from keys import venus_ephemeris_file
from keys import solo_ephemeris_file
from nano_ephemeris import au



def solo_vse(jd,
             venus_ephemeris_file=venus_ephemeris_file,
             solo_ephemeris_file=solo_ephemeris_file,
             location = "998_generated\\assets\\"):
    """
    A function to get the Solar Orbiter's position in VSE coordinates as
    a function of time. If the assets are not ready, 
    build them and then use them.

    This requiers that Venus ephemeris file spans grater range than 
    Solo ephemeris file. Ideally, the ephemeris data points agree on time,
    alternatively the sampling of Venus's location should be higher that that 
    of Solo's location.

    Parameters
    ----------
    jd : float or np.array of float
        Julian date of interest.

    Returns
    -------
    vse : np.array of float
        Is either np.shape(solo_vse(jd))==(3,) if input was float or 
        np.shape(solo_vse(jd))==(n,3) if the input was array of len(jd)==n.
    """
    #try looking for a pre-computed function
    #if found, use it
    #if not found, build it and then use it
    try:
        f_vse_x = load_list("vse_x.pkl",location)[0]
        f_vse_y = load_list("vse_y.pkl",location)[0]
        f_vse_z = load_list("vse_z.pkl",location)[0]
    except:
        venus_jd, venus_hae = load_hae(venus_ephemeris_file)
        solo_jd, solo_hae = load_hae(solo_ephemeris_file)
        vses = np.zeros((0,3))
        for i in range(len(solo_jd)):
            venus_i = (np.abs(venus_jd - solo_jd[i])). argmin()
            venus_hae_i = venus_hae[venus_i]
            vse = hae2vse(solo_hae[i],venus_hae_i)
            vses = np.vstack((vses,vse))
        f_vse_x = interpolate.interp1d(solo_jd, vses[:,0], kind = "cubic")
        f_vse_y = interpolate.interp1d(solo_jd, vses[:,1], kind = "cubic")
        f_vse_z = interpolate.interp1d(solo_jd, vses[:,2], kind = "cubic")
        save_list([f_vse_x],"vse_x.pkl",location)
        save_list([f_vse_y],"vse_y.pkl",location)
        save_list([f_vse_z],"vse_z.pkl",location)
    else:
        pass
    finally:
        vse_x = f_vse_x(jd)
        vse_y = f_vse_y(jd)
        vse_z = f_vse_z(jd)
        vse = np.vstack((vse_x,vse_y,vse_z)).transpose()
        if type(jd) in [np.float64,float]:
            return vse[0]
        else:
            return vse


def get_venus_approaches(distance=0.1):
    """
    Gets an array of all the Venus approaches at which SolO got closer than 
    the specified distance.

    Parameters
    ----------
    distance : float, optional
        The approach defining distance in AU. The default is 0.1.

    Returns
    -------
    approaches : np.array of floats
        Julian dates of the approaches.

    """
    solo_jd, solo_hae = load_hae(solo_ephemeris_file)
    jds = np.arange(min(solo_jd),max(solo_jd),1/1440)
    vses = solo_vse(jds)
    venus_distances = (vses[:,0]**2 + vses[:,1]**2 + vses[:,2]**2)**0.5
    local_minima = argrelextrema(venus_distances, np.less)[0]
    approaches = jds[local_minima[venus_distances[local_minima]<distance]]
    return approaches


def oversampling_flux(time,flux,window=5,order=2):
    """
    The function to return smoothed spline of the flux using Savitzky Golay.

    Parameters
    ----------
    time : np.array of float
        The time (x-axis of the plot).
    flux : np.array of float
        The flux (y-axis of the plot).
    window : int, optional
        The smoothing window, has to be odd. The default is 5.
    order : int, optional
        The polynomial order. The default is 2.

    Returns
    -------
    oversampled_time : np.array of float
        The time (x-axis of the plot), oversampled.
    oversampled_flux : np.array of float
        The flux (y-axis of the plot), oversampled.

    """
    savgoled = savgol_filter(flux,window,order)
    spline = interpolate.interp1d(time, savgoled, kind = "cubic")
    oversampled_time = np.arange(min(time),max(time),0.1)
    oversampled_flux = spline(oversampled_time)
    return oversampled_time, oversampled_flux


def plot_venus_approach_profiles(deltadays=14):
    """
    A procedure to show daily impac counts and how they evolve close to Venus,
    one panel for each Venus encounter.

    Parameters
    ----------
    deltadays : int, optional
        How many days before and after encounter to show. The default is 14.

    Returns
    -------
    None.

    """
    days = load_all_days()
    jds = []
    for day in days:
        jds.append(YYYYMMDD2jd(day.YYYYMMDD))
    jds = np.array(jds)
    approaches = get_venus_approaches(distance=0.1)
    approaches = approaches[(approaches>min(jds))*(max(jds)>approaches)]
    fig = plt.figure(figsize=(3,2))
    gs = fig.add_gridspec(3,1,hspace=0.05)
    ax = gs.subplots(sharex=True)
    for i in range(len(approaches)):
        filtered_days = list(filter(lambda x:
                                    abs(x.date-jd2date(approaches[i]).date())
                                        <
                                        dt.timedelta(days=deltadays), days))
        delta = []
        flux = []
        for day in filtered_days:
            delta.append((day.date-jd2date(approaches[i]).date()).days)
            flux.append(day.impact_count/day.duty_hours*24)
        oversampled_time, oversampled_flux = oversampling_flux(delta,flux)
        ax[i].plot(oversampled_time, oversampled_flux,lw=0.5,color="red")
        ax[i].scatter(delta,flux)
        ax[i].text(.05, .75, str(jd2date(approaches[i]))[:16],
                   fontsize="small", ha='left',
                   transform=ax[i].transAxes)
    for a in ax:
        a.set_ylim(0,650)
    ax[-1].set_xlabel("Days to Venus encounter [day]")
    ax[len(ax)//2].set_ylabel("Impact rate [$day^{-1}$]")

    fig.show()


def plot_venus_impacts(zoom=0.03):
    """
    The procedure to show the zoom on Venus and all the dust impacts 
    encountered near Venus.

    Parameters
    ----------
    zoom : float, optional
        Captured are, in AU. The default is 0.03.

    Returns
    -------
    None.

    """
    solo_jd, solo_hae = load_hae(solo_ephemeris_file)
    trajectory_x = np.zeros(0)
    trajectory_y = np.zeros(0)
    for jd in solo_jd:
        if jd<2460104.5:
            vse = solo_vse(jd)
            trajectory_x = np.append(trajectory_x,vse[0])
            trajectory_y = np.append(trajectory_y,vse[1])

    impacts = load_all_impacts()
    impact_times = []
    vse_xs = np.zeros(0)
    vse_ys = np.zeros(0)
    for i in range(len(impacts)):
        impact = impacts[i]
        impact_times.append(impact.datetime)
        jd = date2jd(impact.datetime)
        vse = solo_vse(jd)
        vse_xs = np.append(vse_xs,vse[0])
        vse_ys = np.append(vse_ys,vse[1])
    
    fig,ax  = plt.subplots()
    ax.plot(trajectory_x,trajectory_y,color="red",lw=0.2)
    ax.scatter(vse_xs,vse_ys,s=5,alpha=0.1,lw=0,color="blue")
    ax.set_xlim(-zoom,zoom)
    ax.set_ylim(-zoom,zoom)
    ax.set_aspect(1)
    ax.set_xlabel("VSE X [AU]")
    ax.set_ylabel("VSE Y [AU]")
    fig.show()




plot_venus_approach_profiles()

plot_venus_impacts()


















