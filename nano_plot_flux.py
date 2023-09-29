import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib as mpl
import glob
import numpy as np
import datetime as dt

from solo_load_days import load_list
from solo_load_days import Impact
from solo_load_days import Day
from conversions import jd2date
import figure_standards as figstd

axes_size = figstd.set_rcparams_dynamo(mpl.rcParams, num_cols=1, ls='thin')
mpl.rcParams['figure.dpi']= 600


apo_jd = 2459000.5+np.array([134,368,566,765,954,1134])
peri_jd = 2459000.5+np.array([255,470,666,865,1045,1224])


def load_all_days(days_location = "998_generated\\days\\"):
    files = glob.glob(days_location+"*.pkl")
    days = []
    for file in files:
        day = load_list(file,"")
        days.append(day)

    return days


def plot_flux(days):
    dates = np.zeros(0,dtype=dt.datetime)
    counts = np.zeros(0)
    duty_hours = np.zeros(0)
    for day in days:
        dates = np.append(dates,day.date)
        counts = np.append(counts,day.impact_count)
        duty_hours = np.append(duty_hours, day.duty_hours)

    fig, ax = plt.subplots(figsize=(3.8, 3))
    ax.set_ylabel("Impact rate (duty-cycle corrected) [$day^{-1}$]"
                  , fontsize="medium")
    ax.set_title('CNN Dust Impacts: '+str(np.round(sum(counts)))[:-2]
                 , fontsize="medium", fontweight="bold")
    ax.tick_params(axis='x',labeltop=False,labelbottom=True)
    ax.tick_params(axis='y',labelleft=True,labelright=False)
    ax.tick_params(axis='y',which="minor",left=True,right=False)
    ax.scatter(dates,counts/duty_hours*24,
               color="maroon", s=1,zorder=100)
    #ax.errorbar(flux_date, flux_cnn, err_plusminus_flux_cnn,
    #            color="firebrick", lw=0, elinewidth=0.4,alpha=0.35,zorder=100)
    ax.xaxis.set_major_locator(mdates.MonthLocator(bymonth=(1, 3, 5, 7, 9, 11)))
    ax.xaxis.set_minor_locator(mdates.MonthLocator())
    ax.set_ylim(bottom = 0, top = 1000)
    ax.set_xlim(left = min(dates), right = max(dates))
    ax.tick_params(axis='x',which="minor",bottom=True,top=True)
    ax.tick_params(axis='x',labelrotation=60)
    ax.tick_params(labelsize="medium")
    ax.set_ylim(0,1000)
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))
    for jd in apo_jd:
        ax.axvline(x=jd2date(jd),color="darkgray",lw=0.7)
    for jd in peri_jd:
        ax.axvline(x=jd2date(jd),color="darkgray",ls="--",lw=0.7)
    ax.text(.85, .96, 'Aphelion', ha='left', va='top', rotation=90, color="gray", fontsize="small", transform=ax.transAxes)
    ax.text(.74, .96, 'Perihelion', ha='left', va='top', rotation=90, color="gray", fontsize="small", transform=ax.transAxes)
    fig.tight_layout()
    fig.savefig('cnn_flux.png', format='png', dpi=600)
    fig.show()



plot_flux(load_all_days())



"""

fig, ax = plt.subplots(figsize=(3.8, 3))
ax.set_ylabel("Impact rate (duty-cycle corrected) [$day^{-1}$]", fontsize=mylabelsize)
ax.set_title('CNN Dust Impacts: '+str(np.round(sum(flux_cnn_imported)))[:-2], fontsize=mytitlesize, fontweight="bold")
ax.tick_params(axis='x',labeltop=False,labelbottom=True)
ax.tick_params(axis='y',labelleft=True,labelright=False)
ax.tick_params(axis='y',which="minor",left=True,right=False)
ax.scatter(flux_date, flux_cnn, color="maroon", s=1,zorder=100)
ax.errorbar(flux_date, flux_cnn, err_plusminus_flux_cnn, color="firebrick", lw=0, elinewidth=0.4,alpha=0.35,zorder=100)
ax.xaxis.set_major_locator(mdates.MonthLocator(bymonth=(1, 3, 5, 7, 9, 11)))
ax.xaxis.set_minor_locator(mdates.MonthLocator())
ax.set_ylim(bottom = 0, top = 1000)
ax.set_xlim(left = min(flux_date), right = max(flux_date))
ax.tick_params(axis='x',which="minor",bottom=True,top=True)
ax.tick_params(axis='x',labelrotation=60)
ax.tick_params(labelsize="medium")
ax.set_ylim(0,1000)
ax2 = ax.twinx()
ax2.set_ylim(0,1000*0.0625)
ax2.set_ylabel('Count [1]')
ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))
for jd in apo_jd:
    ax.axvline(x=jd2date(jd),color="darkgray",lw=0.7)
for jd in peri_jd:
    ax.axvline(x=jd2date(jd),color="darkgray",ls="--",lw=0.7)
ax.text(.85, .96, 'Aphelion', ha='left', va='top', rotation=90, color="gray", fontsize="small", transform=ax.transAxes)
ax.text(.74, .96, 'Perihelion', ha='left', va='top', rotation=90, color="gray", fontsize="small", transform=ax.transAxes)
fig.tight_layout()
fig.savefig('C:\\Users\\skoci\\Disk Google\\000 škola\\UIT\\presentations\\Venus 2023\\cnn_flux.png', format='png', dpi=600)
fig.show()


"""