import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib as mpl
import numpy as np
import datetime as dt


from nano_load_days import load_all_days
from nano_load_days import Impact
from nano_load_days import Day
from nano_load_days import get_errors
from nano_load_days import extract_variables_from_days
from conversions import jd2date
import figure_standards as figstd

axes_size = figstd.set_rcparams_dynamo(mpl.rcParams, num_cols=1, ls='thin')
mpl.rcParams['figure.dpi']= 600


apo_jd = 2459000.5+np.array([134,368,566,765,954,1134])
peri_jd = 2459000.5+np.array([255,470,666,865,1045,1224])








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
    ax.text(.82, .96, 'Aphelion', ha='left', va='top', rotation=90, color="gray", fontsize="small", transform=ax.transAxes)
    ax.text(.74, .96, 'Perihelion', ha='left', va='top', rotation=90, color="gray", fontsize="small", transform=ax.transAxes)
    ax.text(.07, .92, r'$f_s = 262 \, ksps$', ha='left', va='top', color="firebrick", backgroundcolor="white", transform=ax.transAxes)
    ax.text(.07, .83, r'$f_s = 524 \, ksps$', ha='left', va='top', color="teal", backgroundcolor="white", transform=ax.transAxes)
    fig.tight_layout()
    fig.savefig(figures_location+'cnn_flux.png', format='png', dpi=600)
    fig.show()








plot_flux(load_all_days())





