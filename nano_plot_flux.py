import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib as mpl
import glob
import numpy as np
import datetime as dt
from scipy import stats

from nano_load_days import load_list
from nano_load_days import Impact
from nano_load_days import Day
from conversions import jd2date
import figure_standards as figstd

axes_size = figstd.set_rcparams_dynamo(mpl.rcParams, num_cols=1, ls='thin')
mpl.rcParams['figure.dpi']= 600


apo_jd = 2459000.5+np.array([134,368,566,765,954,1134])
peri_jd = 2459000.5+np.array([255,470,666,865,1045,1224])


def load_all_days(days_location = "998_generated\\days\\"):
    """
    The function to load all the measurement days data.

    Parameters
    ----------
    days_location : str, optional
        The data directory. Default is "998_generated\\days\\".

    Returns
    -------
    days : list of Day object
        Measurement days, class Day from nano_load_days.
    """
    files = glob.glob(days_location+"*.pkl")
    days = []
    for file in files:
        day = load_list(file,"")
        days.append(day)

    return days


def extract_variables_from_days(days):
    """
    The function to extract the relevant attributes from the provided days.

    Parameters
    ----------
    days : list of Day object
        Measurement days, class Day from nano_load_days.

    Returns
    -------
    dates : np.array of float
        The measuremetn days.
    counts : np.array of int
        The dust counts, enocutered on the dates.
    duty_hours : np.array of float
        The duty cycle in hours, per day. Often around 1.5.
    sampling_rates : np.array of float
        The sampling rate on that day. No dates with a vraiable 
        sampling rate were processed in the first place.
    """
    dates = np.zeros(0,dtype=dt.datetime)
    counts = np.zeros(0,dtype=int)
    duty_hours = np.zeros(0,dtype=float)
    sampling_rates = np.zeros(0,dtype=float)
    for day in days:
        dates = np.append(dates,day.date)
        counts = np.append(counts,day.impact_count)
        duty_hours = np.append(duty_hours, day.duty_hours)
        sampling_rates = np.append(sampling_rates, day.sampling_rate)

    return dates, counts, duty_hours, sampling_rates


def get_errors(days, prob_coverage = 0.9):
    """
    The function to calculate the errorbars for flux 
    assuming Poisson distribution and taking into account
    the number of detection.

    Parameters
    ----------
    days : list of Day object
        Measurement days, class Day from nano_load_days.
    prob_coverage : float, optional
        The coverage of the errobar interval. The default is 0.9.

    Returns
    -------
    err_plusminus_flux : np.array of float
        The errorbars, lower and upper bound, shape (2, n).

    """
    dates, counts, duty_hours, sampling_rates = extract_variables_from_days(days)

    counts_err_minus = -stats.poisson.ppf(0.5-prob_coverage/2, mu=counts)+counts
    counts_err_plus  = +stats.poisson.ppf(0.5+prob_coverage/2, mu=counts)-counts
    err_plusminus_flux = np.array([counts_err_minus,counts_err_plus]) / (duty_hours/(24))

    return err_plusminus_flux


def plot_flux(days, figures_location = "998_generated\\figures\\"):
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
    ax.text(.85, .96, 'Aphelion', ha='left', va='top', rotation=90, color="gray", fontsize="small", transform=ax.transAxes)
    ax.text(.74, .96, 'Perihelion', ha='left', va='top', rotation=90, color="gray", fontsize="small", transform=ax.transAxes)
    fig.tight_layout()
    fig.savefig(figures_location+'cnn_flux.png', format='png', dpi=600)
    fig.show()








plot_flux(load_all_days())





