import matplotlib.pyplot as plt
import matplotlib as mpl
import glob
import numpy as np
import os
import datetime as dt
from scipy import interpolate
from scipy.signal import argrelextrema
from scipy.signal import savgol_filter

from nano_load_days import load_all_impacts
from nano_load_days import load_all_days
from nano_load_days import load_list
from nano_load_days import save_list
from nano_load_days import get_errors
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

from paths import venus_ephemeris_file
from paths import solo_ephemeris_file
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

    location = os.path.join(os.path.normpath( location ), '')

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


def get_venus_distances(jds):
    """
    Gives the distance of SolO from Venus in times jds. 

    Parameters
    ----------
    jds : float or np.array of float
        Julian dates of interest.

    Returns
    -------
    venus_distances : float or np.array of float
        Distances.

    """
    vses = solo_vse(jds)
    if vses.ndim == 1:
        venus_distances = (vses[0]**2 + vses[1]**2 + vses[2]**2)**0.5
    elif vses.ndim == 2:
        venus_distances = (vses[:,0]**2 + vses[:,1]**2 + vses[:,2]**2)**0.5
    else:
        raise ValueError("unsupported jds dim "+str(vses.ndim-1))
    return venus_distances


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
    venus_distances = get_venus_distances(jds)
    local_minima = argrelextrema(venus_distances, np.less)[0]
    approaches = jds[local_minima[venus_distances[local_minima]<distance]]
    return approaches


def oversampling_flux(time,flux,window=7,order=2):
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


def plot_approach_profiles(approaches,
                           deltadays=14,
                           target = "Venus encounter",
                           distance_measure = None,
                           errors_on_flux = None,
                           force_ylim = None,
                           spline = True,
                           cscheme = "blue",
                           figures_location = "998_generated\\figures\\",
                           name = 'venus_approach_profiles.png'):
    """
    A procedure to show daily impact counts and how they evolve close to 
    a point, like an approach. One panel for each encounter is shown.

    Parameters
    ----------
    approaches : list of float
        Julin date of the approaches in question.

    deltadays : int, optional
        How many days before and after encounter to show. The default is 14.

    target : str, optional
        Name of the body of interest, such as Venus. The default is "Venus".

    distance_measure : function : jd -> float, optional
        A function to measure the distance on the day of interes.
        Julian date of the approach is the input, while distance in AU
        is the output. The default in None.

    errors_on_flux : function Day -> np.array of float, optional
        A function that can return the errorbars on the day. 
        The output must be tThe errorbars, 
        lower and upper bound, shape (2, n). 
        The default in None.

    force_ylim : float, optional
        If used, it will force the force_ylim value as the value for the 
        common upper y limit. Bottom is 0. The default is None.
    
    spline : boolean, optional
        If True, the spline is drawn. The default is True.

    cscheme : str, optional
        Color scheme. The default is "blue".

    figures_location : str
        The directory where to put the plot. 
        The default is "998_generated\\figures\\"

    name : str
        Name of the plot to be saved. 
        The default is 'venus_approach_profiles.png'.
        

    Returns
    -------
    None.

    """

    figures_location = os.path.join(os.path.normpath( figures_location ), '')

    days = load_all_days()
    jds = []
    for day in days:
        jds.append(YYYYMMDD2jd(day.YYYYMMDD))
    jds = np.array(jds)

    approaches = approaches[(approaches>min(jds))*(max(jds)>approaches)]
    if distance_measure != None:
        approach_distances = distance_measure(approaches)
    fig = plt.figure(figsize=(3,2*len(approaches)/3))
    gs = fig.add_gridspec(len(approaches),1,hspace=0.05)
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
        if spline:
            oversampled_time, oversampled_flux = oversampling_flux(delta,flux)
            ax[i].plot(oversampled_time, oversampled_flux,
                       lw=0.5,color="red")
        if errors_on_flux != None:
            errors = errors_on_flux(filtered_days)
            ax[i].errorbar(delta,flux, errors,
                           color=cscheme,alpha=0.2,elinewidth=1,lw=0)
        ax[i].scatter(delta,flux,color=cscheme)
        ax[i].text(.05, .75, str(jd2date(approaches[i]))[:16],
                   fontsize="x-small", ha='left',
                   transform=ax[i].transAxes)
        if distance_measure != None:
            ax[i].text(.97, .75, str(int(approach_distances[i]*au))+" km",
                       fontsize="x-small", ha='right',
                       transform=ax[i].transAxes)

    for a in ax:
        a.set_ylim(bottom=0)
        if force_ylim != None:
            a.set_ylim(top=force_ylim)
        top = a.get_ylim()[1]
        a.set_xlim(-deltadays,deltadays)
        a.vlines(0,0,top,color="gray",lw=0.5,alpha=0.3,zorder=1,ls="dotted")
        a.set_ylim(bottom=0,top=top)
    ax[-1].set_xlabel("Days since "+target+" [day]")
    ax[len(ax)//2].set_ylabel("Impact rate [$day^{-1}$]")
    fig.savefig(figures_location+name,
                format='png', dpi=600)
    fig.show()


def plot_venus_impacts(zoom=0.0005,
                       figures_location = "998_generated\\figures\\",
                       date_from = dt.datetime(2010,1,1),
                       date_to = dt.datetime(2050,1,1)):
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

    figures_location = os.path.join(os.path.normpath( figures_location ), '')

    solo_jd, solo_hae = load_hae(solo_ephemeris_file)

    solo_jd_oversampled = np.arange(min(solo_jd[solo_jd<2460104.5]),
                                    max(solo_jd[solo_jd<2460104.5]),
                                    0.015625)

    vse = solo_vse(solo_jd_oversampled)
    trajectory_x = vse[:,0]
    trajectory_y = vse[:,1]

    impacts = load_all_impacts(date_from = date_from, date_to = date_to)
    indices = np.array([impact.index for impact in impacts])
    impact_times = []
    jds = []
    for i in range(len(impacts)):
        impact = impacts[i]
        impact_times.append(impact.datetime)
        jds.append(date2jd(impact.datetime))
    jds = np.array(jds)
    vse = solo_vse(jds)
    vse_xs = vse[:,0]
    vse_ys = vse[:,1]
    
    fig,ax  = plt.subplots()
    venus = plt.Circle(( 0 , 0), 6000, facecolor="brown", edgecolor="none")
    ax.plot(trajectory_x*au,trajectory_y*au,color="red",lw=0.2)
    ax.scatter(vse_xs*au,vse_ys*au,s=5,alpha=0.1,lw=0,color="blue")
    for i,index in enumerate(indices.astype(str)):
        ax.annotate(index,(vse_xs[i]*au,vse_ys[i]*au),
                    (vse_xs[i]*au+np.random.uniform(-1.5,0)*20000,
                     vse_ys[i]*au+np.random.uniform(0,1.5)*20000),
                    fontsize="xx-small")
    ax.add_artist(venus)
    ax.text(.05, .92, date_from.strftime("%Y-%m-%d %H:%M"),
               fontsize="small", ha='left',
               transform=ax.transAxes)
    ax.text(.05, .85, date_to.strftime("%Y-%m-%d %H:%M"),
               fontsize="small", ha='left',
               transform=ax.transAxes)
    ax.set_xlim(-zoom*au,zoom*au)
    ax.set_ylim(-zoom*au,zoom*au)
    ax.set_aspect(1)
    ax.set_xlabel("VSE X [km]")
    ax.set_ylabel("VSE Y [km]")
    fig.savefig(figures_location+'venus_approach_impacts.png',
                format='png', dpi=600)
    fig.show()



if __name__ == "__main__":
    plot_approach_profiles(get_venus_approaches(distance=0.1),
                           distance_measure=get_venus_distances,
                           errors_on_flux = get_errors)
    
    plot_venus_impacts(date_from = dt.datetime(2022,9,4,0,0),
                       date_to = dt.datetime(2022,9,4,6,0))


















