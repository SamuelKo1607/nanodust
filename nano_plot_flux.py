import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib as mpl
import numpy as np
import datetime as dt
from scipy import interpolate
from scipy.signal import argrelextrema
from scipy.signal import savgol_filter
import os
import glob

from paths import solo_ephemeris_file

from nano_load_days import load_all_days
from nano_load_days import load_list
from nano_load_days import save_list
from nano_load_days import Impact
from nano_load_days import Day
from nano_mamp import ImpactSuspect
from nano_load_days import moving_average
from nano_load_days import get_errors
from nano_load_days import extract_variables_from_days
from nano_ephemeris import load_hae
from nano_mamp import load_all_suspects
from conversions import jd2date
from conversions import date2jd
from conversions import YYYYMMDD2date
import figure_standards as figstd
from venus_view_impacts import plot_approach_profiles


axes_size = figstd.set_rcparams_dynamo(mpl.rcParams, num_cols=1, ls='thin')
mpl.rcParams['figure.dpi']= 600


apo_jd = 2459000.5+np.array([134,368,566,765,954,1134])
peri_jd = 2459000.5+np.array([255,470,666,865,1045,1224])





def fill_nan(A):
    """
    Interpolate to fill NaNs. Found on stack overflow: 
    https://stackoverflow.com/questions/6518811/interpolate-nan-values-in-a-numpy-array

    Parameters
    ----------
    A : np.array of float
        The array with some nans.

    Returns
    -------
    B : np.array of float
        The array with no nans, they are interpolated over.

    """
    inds = np.arange(A.shape[0])
    good = np.where(np.isfinite(A))
    f = interpolate.interp1d(inds[good], A[good],bounds_error=False)
    B = np.where(np.isfinite(A),A,f(inds))
    return B


def plot_flux(days,
              figures_location = os.path.join("998_generated","figures",""),
              overplot = None,
              styles = None,
              aspect = 1.2666):
    """
    A plot of daily flux is made with the data from the provided days files.

    Parameters
    ----------
    days : list of Day object
        Measurement days, class Day from nano_load_days.

    figures_location : str, optional
        Where to put the wrawn figure. Default is "figures_location".

    overplot : list of functions: np.array of dt.datetime -> np.array of float, optional
        A list of functons that will be used to overplot over the data. 
        Default is None, in which case nothing is overplotted.

    styles : list of str, optional
        list of fmt string, such as "g:2" or something, this assigns 
        a style to the overplot lines. Default is None, in which case "b-" 
        is used.

    Returns
    -------
    None.
    """

    if styles == None and overplot != None:
        styles = ["b-"]*len(overplot)

    figures_location = os.path.join(os.path.normpath( figures_location ), '')

    dates, counts, duty_hours, sampling_rates = extract_variables_from_days(days)

    err_plusminus_flux = get_errors(days)

    colorcodes = np.zeros(0,dtype=str)
    for sampling_rate in sampling_rates:
        if sampling_rate < 263000:
            colorcodes = np.append(colorcodes,"firebrick")
        else:
            colorcodes = np.append(colorcodes,"teal")

    fig, ax = plt.subplots(figsize=(3*aspect, 3))
    ax.set_ylabel("Impact rate (duty-cycle corrected) [$day^{-1}$]"
                  , fontsize="medium")
    ax.set_title('CNN Dust Impacts: '+str(np.round(sum(counts)))
                 , fontsize="medium", fontweight="bold")
    ax.tick_params(axis='x',labeltop=False,labelbottom=True)
    ax.tick_params(axis='y',labelleft=True,labelright=False)
    ax.tick_params(axis='y',which="minor",left=True,right=False)
    if overplot!= None:
        for i, line in enumerate(overplot):
            ax.plot(dates,line(dates),styles[i],lw=1,ms=0)
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


def plot_mamp_suspects(all_suspects,
                       compare_days = None,
                       overplot = None,
                       styles = None,
                       ylim = 2000,
                       figures_location = os.path.join("998_generated",
                                                       "figures","")):
    """
    A plot of MAMP suspect flux, compared to CNN impacts.
    
    Parameters
    ----------
    all_suspects : list of ImpactSuspect
        Impact suspects, class ImpactSuspect from nano_mamp.
    
    compare_days : list of Day, optional
        Measurement days, class Day from nano_load_days. Default is None, in 
        which case nothing is shown.

    overplot : list of functions: np.array of dt.datetime -> 
               np.array of float, optional
        A list of functons that will be used to overplot over the data. 
        Default is None, in which case nothing is overplotted.

    styles : list of str, optional
        list of fmt string, such as "g:2" or something, this assigns 
        a style to the overplot lines. Default is None, in which case "b-" 
        is used.

    ylim : float, optional
        The upper bound on y-axis. The default is 1500.        

    figures_location : str, optional
        Where to put the wrawn figure. Default is "998_generated\\figures\\".
    
    Returns
    -------
    None.
    """

    if styles == None and overplot != None:
         styles = ["b-"]*len(overplot)

    YYYYMMDDs = [suspect.YYYYMMDD for suspect in all_suspects]
    if compare_days != None:
        YYYYMMDDs += [day.YYYYMMDD for day in compare_days]
    unique_YYYYMMDDs = np.unique(np.array(YYYYMMDDs))
    if compare_days != None:
        counts_days = np.zeros(len(unique_YYYYMMDDs),dtype=float)
    counts_suspects = np.zeros(len(unique_YYYYMMDDs),dtype=float)
    coverage_hours = np.zeros(len(unique_YYYYMMDDs),dtype=float)
    uniqe_dates = np.zeros(len(unique_YYYYMMDDs),dtype=dt.datetime)
    for i, YYYYMMDD in enumerate(unique_YYYYMMDDs):
        counts_suspects[i] = len([suspect for suspect in all_suspects
                                  if suspect.YYYYMMDD == YYYYMMDD])
        coverage_hours[i] = np.mean([suspect.coverage_hours
                                     for suspect in all_suspects
                                     if suspect.YYYYMMDD == YYYYMMDD])
        uniqe_dates[i] = YYYYMMDD2date(YYYYMMDD)
        if counts_suspects[i] == 0:
            counts_suspects[i] = np.nan
        if compare_days != None:
            counts_days[i] += sum([day.impact_count * 24 / day.duty_hours
                                   for day in compare_days
                                   if day.YYYYMMDD == YYYYMMDD])
            if counts_days[i] == 0:
                counts_days[i] = np.nan

    counts_suspects_normalized = counts_suspects*24/coverage_hours

    fig, ax = plt.subplots()

    #plot MAMP suspects
    ax.plot(uniqe_dates,counts_suspects_normalized,
            alpha=.6,color="tab:blue",
            label=r"MAMP (>12$\sigma$)")
    #add exrtremes
    crossings = [[d, c] for d, c in zip(uniqe_dates,counts_suspects_normalized)
                 if c>ylim]
    for crossing in crossings:
        ax.text(crossing[0]-dt.timedelta(weeks=2), ylim*0.5,
                "{:.0f}".format(crossing[1]),
                rotation = 90, horizontalalignment='right',
                fontsize="x-small", color="tab:blue",)

    #compare to CNN data
    if compare_days != None:
        ax.plot(uniqe_dates,counts_days,
                alpha=.6, color="tab:orange",
                label="CNN (/24h equiv.)")
        #add exrtremes
        crossings = [[d, c] for d, c in zip(uniqe_dates,counts_days)
                     if c>ylim]
        for crossing in crossings:
            ax.text(crossing[0]-dt.timedelta(weeks=2), ylim*0.8,
                    "{:.0f}".format(crossing[1]),
                    rotation = 90, horizontalalignment='right',
                    fontsize="x-small", color="tab:orange")

    #overplot models
    if overplot!= None:
        for i, line in enumerate(overplot):
            ax.plot(uniqe_dates,line(uniqe_dates),
                    styles[i],ms=0,label = "Bayesian fit")


    ax.legend(fontsize="x-small")
    ax.set_ylim(0,ylim)
    ax.set_ylabel("Impact rate [$day^{-1}$]", fontsize="medium")
    ax.tick_params(axis='x',labelrotation=60)
    fig.tight_layout()
    fig.savefig(figures_location+'suspects_vs_cnn.png', format='png', dpi=600)
    fig.show()


def solo_hae(jd,
             solo_ephemeris_file=solo_ephemeris_file,
             location = os.path.join("998_generated","assets","")):
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


def build_bayesian_fit(fit_data):
    """
    A wrapper function to give the fit result of the Bayesian fitting. 
    Requires loading the data in the form of X, Y1, Y2, Y3; handily provided 
    by a different script of the "dust" project. 

    Parameters
    ----------
    array_of_datetimes : np.array of dt.datetime
        Usually something corresponding to the X-axis that is to be plotted.
    fit_data : list [0:3] of array of float, identical array lengths.
        Underlying data for the curves to be plotted. Usually comes from 
        a pickle provided by a different script, shaping the output of a
        Bayesian fitting routine. The first array is the underlying 
        julian dates, while the other 3 are the curves.

    Returns
    -------
    v_f_mean_dt : function : np.array of dt.datetime -> np.array of float
        The flux in /day, bottom 5% quantile as provided by the Baysian fitting.
    v_f_q5_dt : function : np.array of dt.datetime -> np.array of float
        The flux in /day, mean as provided by the Baysian fitting.
    v_f_q95_dt : function : np.array of dt.datetime -> np.array of float
        The flux in /day, top 5% quantile as provided by the Baysian fitting.

    """
    jd_span = fit_data[0]
    mean_s = fit_data[1]
    q5_s = fit_data[2]
    q95_s = fit_data[3]

    f_mean = interpolate.interp1d(jd_span, mean_s, fill_value="extrapolate", kind=3)
    f_q5 = interpolate.interp1d(jd_span, q5_s, fill_value="extrapolate", kind=3)
    f_q95 = interpolate.interp1d(jd_span, q95_s, fill_value="extrapolate", kind=3)

    v_f_mean = np.vectorize(f_mean)
    v_f_q5 = np.vectorize(f_q5)
    v_f_q95 = np.vectorize(f_q95)

    v_f_mean_dt = lambda x : v_f_mean(date2jd(x)) #x is dt.datetime
    v_f_q5_dt = lambda x : v_f_q5(date2jd(x)) #x is dt.datetime
    v_f_q95_dt = lambda x : v_f_q95(date2jd(x)) #x is dt.datetime

    return v_f_mean_dt,v_f_q5_dt,v_f_q95_dt


def plot_asymmetric_flux(days,
                         asymmetry=0.2,
                         figures_location=os.path.join("998_generated","figures",""),
                         impacts_location=os.path.join("998_generated","impacts",""),
                         perihelia=[]):
    """
    To study the evolution of the flux of the very asymmetric impacts.

    Parameters
    ----------
    days : list of Day object
        Measurement days, class Day from nano_load_days.

    asymmetry : float, optional
        The threshold of symmetry to consider something asymmetric. 
        The default is 0.2.

    figures_location : str, optional
        Where to put the plots. 
        The default is os.path.join("998_generated","figures","").

    impacts_location : str, optional
        Where the Impact pickles are. 
        The default is os.path.join("998_generated","impacts","").

    perihelia : list of float
        List of perihelia dates in JD. If provided, vlines are shown.    

    Returns
    -------
    None.

    """

    dates = []
    total = np.zeros(0, dtype=int)
    asymmetric = np.zeros(0, dtype=int)
    duty_hours = np.zeros(0, dtype=float)
    heliocentric_distances = np.zeros(0, dtype=float)

    for i, day in enumerate(days):
        if not i%100:
            print(f"{i+1}/{len(days)}")
        file = glob.glob(impacts_location+"*"+day.YYYYMMDD+"*.pkl")
        if len(file) == 1:
            impacts = load_list(file[0],"")

            dates.append(day.date)
            total = np.append(total,len(impacts))
            asymmetric = np.append(asymmetric,
                                   len([impact for impact in impacts if impact.symmetry<0.2]))
            duty_hours = np.append(duty_hours,day.duty_hours)
            heliocentric_distances = np.append(heliocentric_distances, day.heliocentric_distance)
        else:
            pass

    fig = plt.figure()
    gs = fig.add_gridspec(2,1,wspace=0.1)
    ax = gs.subplots(sharex=True)

    ax[0].plot(dates,total*24/duty_hours,label="Total")
    ax[0].plot(dates,asymmetric*24/duty_hours,label="Asymmetric")
    ax[0].legend(fontsize="small")
    ax[0].set_ylim(0,1500)
    ax[0].set_ylabel("Impact\n rate [$day^{-1}$]", fontsize="medium")
    ax[0].tick_params(axis='x',labelrotation=60)

    ax[1].plot(dates,100*asymmetric/total,color="dodgerblue",alpha=0.3)
    smooth = moving_average(fill_nan(asymmetric/total),29)
    ax[1].plot(dates[30:-30],100*smooth[30:-30],color="navy")
    ax[1].set_ylim(0,100)
    ax[1].hlines(0,min(dates),max(dates),color="gray")
    ax[1].set_ylabel("Asymmetric\n fraction [\%]", fontsize="medium")
    ax[1].tick_params(axis='x',labelrotation=60)

    ax_alt = ax[1].twinx()
    ax_alt.set_ylim(0,1.2)
    ax_alt.plot(dates,heliocentric_distances,color="darkred")
    ax_alt.set_ylabel("Heliocentric\n distnace [AU]", fontsize="medium")

    if len(perihelia):
        for a in ax:
            for i,aproach in enumerate(perihelia):
                left,right = a.get_xlim()
                bottom,top = a.get_ylim()
                a.vlines(jd2date(aproach),
                         bottom,top,
                         color="gray",zorder=1,ls="dotted")
                a.set_ylim(bottom=bottom,top=top)
                a.set_xlim(left=left,right=right)


    fig.tight_layout()
    fig.savefig(figures_location+'asymmetric_flux.png', format='png', dpi=600)
    fig.show()







#%%
if __name__ == "__main__":

    perihelia = get_sun_approaches(distance=0.6)

    #loading the bayesian fit data
    mean, bottom5, top5 = build_bayesian_fit(load_list("bayesian_fit_orig.pkl",
                                                       os.path.join("data_synced","")
                                                       ))

    plot_flux(load_all_days(),
              aspect = 2)

    plot_flux(load_all_days(),
              overplot=[bottom5, mean, top5],
              styles=["k:","k-","k:"])

    plot_approach_profiles(perihelia,
                           deltadays = 30,
                           target = "perihelion",
                           force_ylim = 1400,
                           spline = False,
                           cscheme = "forestgreen",
                           distance_measure = get_heliocentric_distances,
                           errors_on_flux = get_errors,
                           overplot=[bottom5, mean, top5],
                           styles=["k:","k-","k:"],
                           figures_location = "998_generated\\figures\\",
                           name = 'perihelia_approach_profiles.png')

    plot_mamp_suspects(load_all_suspects(),
                       compare_days = load_all_days(),
                       overplot=[mean],
                       styles=["k-"])

    plot_asymmetric_flux(load_all_days(),
                         perihelia=perihelia)


