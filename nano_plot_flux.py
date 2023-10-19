import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib as mpl
import numpy as np
import datetime as dt
from scipy import interpolate
from scipy.signal import argrelextrema
from scipy.signal import savgol_filter
import os

from paths import solo_ephemeris_file

from nano_load_days import load_all_days
from nano_load_days import load_list
from nano_load_days import save_list
from nano_load_days import Impact
from nano_load_days import Day
from nano_load_days import get_errors
from nano_load_days import extract_variables_from_days
from nano_ephemeris import load_hae
from conversions import jd2date
import figure_standards as figstd
from venus_view_impacts import plot_approach_profiles


axes_size = figstd.set_rcparams_dynamo(mpl.rcParams, num_cols=1, ls='thin')
mpl.rcParams['figure.dpi']= 600


apo_jd = 2459000.5+np.array([134,368,566,765,954,1134])
peri_jd = 2459000.5+np.array([255,470,666,865,1045,1224])








def plot_flux(days,
              figures_location = "998_generated\\figures\\"):
    """
    A plot of daily flux is made with the data from the provided days files.

    Parameters
    ----------
    days : list of Day object
        Measurement days, class Day from nano_load_days.

    figures_location : str, optional
        Where to put the wrawn figure. Default is "figures_location".

    Returns
    -------
    None.
    """

    figures_location = os.path.join(os.path.normpath( figures_location ), '')

    dates, counts, duty_hours, sampling_rates = extract_variables_from_days(days)

    err_plusminus_flux = get_errors(days)

    colorcodes = np.zeros(0,dtype=str)
    for sampling_rate in sampling_rates:
        if sampling_rate < 263000:
            colorcodes = np.append(colorcodes,"firebrick")
        else:
            colorcodes = np.append(colorcodes,"teal")

    fig, ax = plt.subplots(figsize=(3.8, 3))
    ax.set_ylabel("Impact rate (duty-cycle corrected) [$day^{-1}$]"
                  , fontsize="medium")
    ax.set_title('CNN Dust Impacts: '+str(np.round(sum(counts)))
                 , fontsize="medium", fontweight="bold")
    ax.tick_params(axis='x',labeltop=False,labelbottom=True)
    ax.tick_params(axis='y',labelleft=True,labelright=False)
    ax.tick_params(axis='y',which="minor",left=True,right=False)
    ax.scatter(dates,counts/duty_hours*24,
               c=colorcodes, s=1,zorder=100)
    for color in np.unique(colorcodes):
        mask = colorcodes==color
        ax.errorbar(dates[mask], counts[mask]/duty_hours[mask]*24, err_plusminus_flux[:,mask],
                    c=color, lw=0, elinewidth=0.4,alpha=0.35,zorder=100)
    ax.xaxis.set_major_locator(mdates.MonthLocator(bymonth=(1, 3, 5, 7, 9, 11)))
    ax.xaxis.set_minor_locator(mdates.MonthLocator())
    ax.set_xlim(left = min(dates), right = max(dates))
    ax.tick_params(axis='x',which="minor",bottom=True,top=True)
    ax.tick_params(axis='x',labelrotation=60)
    ax.tick_params(labelsize="medium")
    ax.set_ylim(0,1400)
    ax.tick_params(labelleft=True,
                   labelright=True,
                   labelbottom = True,
                   labeltop = False)
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))
    for jd in apo_jd:
        ax.axvline(x=jd2date(jd),color="darkgray",lw=0.7)
    for jd in peri_jd:
        ax.axvline(x=jd2date(jd),color="darkgray",ls="--",lw=0.7)
    ax.text(.82, .96, 'Aphelion', ha='left', va='top', rotation=90, color="gray", fontsize="small", transform=ax.transAxes)
    ax.text(.74, .96, 'Perihelion', ha='left', va='top', rotation=90, color="gray", fontsize="small", transform=ax.transAxes)
    ax.text(.07, .92, r'$f_s = 262 \, ksps$', ha='left', va='top', color="firebrick", backgroundcolor="white", transform=ax.transAxes)
    ax.text(.07, .83, r'$f_s = 524 \, ksps$', ha='left', va='top', color="teal", backgroundcolor="white", transform=ax.transAxes)
    fig.tight_layout()
    fig.savefig(figures_location+'cnn_flux.png', format='png', dpi=600)
    fig.show()


def solo_hae(jd,
             solo_ephemeris_file=solo_ephemeris_file,
             location = "998_generated\\assets\\"):
    """
    A function to get the Solar Orbiter's position in HAE coordinates as
    a function of time. If the assets are not ready, 
    build them and then use them.

    Parameters
    ----------
    jd : float or np.array of float
        Julian date of interest.

    Returns
    -------
    hae : np.array of float
        Is either np.shape(solo_hae(jd))==(3,) if input was float or 
        np.shape(solo_hae(jd))==(n,3) if the input was array of len(jd)==n.
    """

    location = os.path.join(os.path.normpath( location ), '')

    try:
        f_hae_x = load_list("hae_x.pkl",location)[0]
        f_hae_y = load_list("hae_y.pkl",location)[0]
        f_hae_z = load_list("hae_z.pkl",location)[0]
    except:
        solo_jd, solo_hae = load_hae(solo_ephemeris_file)
        f_hae_x = interpolate.interp1d(solo_jd, solo_hae[:,0], kind = "cubic")
        f_hae_y = interpolate.interp1d(solo_jd, solo_hae[:,1], kind = "cubic")
        f_hae_z = interpolate.interp1d(solo_jd, solo_hae[:,2], kind = "cubic")
        save_list([f_hae_x],"hae_x.pkl",location)
        save_list([f_hae_y],"hae_y.pkl",location)
        save_list([f_hae_z],"hae_z.pkl",location)
    else:
        pass
    finally:
        hae_x = f_hae_x(jd)
        hae_y = f_hae_y(jd)
        hae_z = f_hae_z(jd)
        hae = np.vstack((hae_x,hae_y,hae_z)).transpose()
        if type(jd) in [np.float64,float]:
            return hae[0]
        else:
            return hae


def get_heliocentric_distances(jds):
    """
    Gives the distance of SolO from the Sun in times jds. 

    Parameters
    ----------
    jds : float or np.array of float
        Julian dates of interest.

    Returns
    -------
    heliocentric_distances : float or np.array of float
        Distances in AU.

    """
    hae = solo_hae(jds)
    if hae.ndim == 1:
        heliocentric_distances = (hae[0]**2 + hae[1]**2 + hae[2]**2)**0.5
    elif hae.ndim == 2:
        heliocentric_distances = (hae[:,0]**2 + hae[:,1]**2 + hae[:,2]**2)**0.5
    else:
        raise ValueError("unsupported jds dim "+str(hae.ndim-1))
    return heliocentric_distances


def get_sun_approaches(distance=0.6):
    """
    Gets an array of all the Sun approaches at which SolO got closer than 
    the specified distance.

    Parameters
    ----------
    distance : float, optional
        The approach defining distance in AU. The default is 0.6.

    Returns
    -------
    approaches : np.array of floats
        Julian dates of the approaches.

    """
    solo_jd, solo_hae = load_hae(solo_ephemeris_file)
    jds = np.arange(min(solo_jd),max(solo_jd),1/1440)
    heliocentric_distances = get_heliocentric_distances(jds)
    local_minima = argrelextrema(heliocentric_distances, np.less)[0]
    approaches = jds[local_minima[heliocentric_distances[local_minima]<distance]]
    return approaches







if __name__ == "__main__":
    plot_flux(load_all_days())

    plot_approach_profiles(get_sun_approaches(distance=0.6),
                           deltadays = 14,
                           target = "perihelion",
                           force_ylim = 1400,
                           spline = False,
                           cscheme = "forestgreen",
                           distance_measure = get_heliocentric_distances,
                           errors_on_flux = get_errors,
                           figures_location = "998_generated\\figures\\",
                           name = 'perihelia_approach_profiles.png')





