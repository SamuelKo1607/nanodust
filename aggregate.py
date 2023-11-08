import os
import csv

from nano_load_days import Day
from nano_mamp import ImpactSuspect

from nano_mamp import load_all_suspects
from nano_load_days import save_list
from conversions import date2jd
from nano_load_days import load_all_days




def make_flux_to_fit_inla(days,
                          name = "flux_readable.csv",
                          location = os.path.join("data_synced","")):
    """
    Creates a CSV in a reasonably readable format, to use in an external
    application, such as RStudio. 

    Parameters
    ----------
    days : list of Day
        Impact data to include.
    name : str, optional
        the name of the file. The default is "flux_readable.csv".
    location : str, optional
        The path where to put the result. 
        The default is os.path.join("data_synced","").

    Returns
    -------
    None.

    """

    with open(location+name, 'w', newline='') as f:
        writer = csv.writer(f, delimiter=',')
        writer.writerow(["Julian date",
                         "Fluxes [/day]",
                         "Radial velocity [km/s]",
                         "Tangential velocity [km/s]",
                         "Radial distance [au]",
                         "Detection time [hours]",
                         "Velocity phase angle [deg]",
                         "Velocity inclination [deg]",
                         "V_X (HAE) [km/s]",
                         "V_Y (HAE) [km/s]",
                         "V_Z (HAE) [km/s]"])
        for day in days:
            writer.writerow([str(date2jd(day.date)),
                             str(day.impact_count),
                             str(day.heliocentric_radial_speed),
                             str(day.heliocentric_tangential_speed),
                             str(day.heliocentric_distance),
                             str(day.duty_hours),
                             str(day.velocity_phase),
                             str(day.velocity_inclination),
                             str(day.velocity_HAE_x),
                             str(day.velocity_HAE_y),
                             str(day.velocity_HAE_z)])



def aggregate_flux_readable(location = os.path.join("998_generated","days",""),
                            target_location = os.path.join("data_synced",""),
                            force = False):
    """
    A wrapper function to aggregate all Day objects creted with 
    nano_make_file into one human redable file, this checks if it exists 
    and if not, then saves it.

    Parameters
    ----------
    location : str, optional
        The location of the source data. 
        The default is os.path.join("998_generated","days","").
    target_location : str, optional
        The location where the aggregated data will be saved. 
        The default is os.path.join("data_synced","").
    force : bool, optional
        Whether to overwrite any existing file. The default is False.

    Returns
    -------
    None.

    """
    if force:
        print("saving readable file forcedly")
        make_flux_to_fit_inla(load_all_days(days_location = location),
                              location = target_location)
    else:
        try:
            with open(target_location+"flux_readable.csv") as f:
                print("readable file OK:")
                print(f)
        except:
            print("saving readable file")
        make_flux_to_fit_inla(load_all_days(days_location = location),
                              location = target_location)


def aggregate_mamp_suspects(location = os.path.join("998_generated","mamp_processed",""),
                            target_location = os.path.join("data_synced",""),
                            force = False):
    """
    A wrapper function to aggregate all ImpactSuspect objects creted with 
    nano_mamp into one pickle, this checks if it exists 
    and if not, then saves it.

    Parameters
    ----------
    location : str, optional
        The location of the source data. 
        The default is os.path.join("998_generated","mamp_processed","").
    target_location : str, optional
        The location where the aggregated data will be saved. 
        The default is os.path.join("data_synced","").
    force : bool, optional
        Whether to overwrite any existing file. The default is False.

    Returns
    -------
    None.

    """
    if force:
        print("saving all mamp suspects forcedly")
        all_suspects = load_all_suspects(location)
        save_list(all_suspects, "all_suspects.pkl", target_location)
    else:
        try:
            with open(target_location+"all_suspects.pkl") as f:
                print("mamp suspects OK:")
                print(f)
        except:
            print("saving all mamp suspects")
            all_suspects = load_all_suspects(location)
            save_list(all_suspects, "all_suspects.pkl", target_location)


#%%

if __name__ == "__main__":
    aggregate_flux_readable()
    aggregate_mamp_suspects()