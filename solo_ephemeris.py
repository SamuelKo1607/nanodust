import numpy as np
import csv
from scipy import interpolate

au = 149597870.7 #astronomical unit, km


def load_ephemeris(ephemeris_file):
    """
    Parameters
    ----------
    ephemeris_file : str
        The ephemeris file to access.

    Raises
    ------
    LookupError
        in case the file is not loaded correctly.

    Returns
    -------
    time : numpy.ndarray(1,:) of float
        Julian date of an ephemeris table record.
    hae_r : numpy.ndarray(3,:) of float
        3D vector of HAE coordinates for every record in the epehemeris file.
    hae_v : numpy.ndarray(3,:) of float
        3D vector of HAE velocities for every record in the epehemeris file.
    hae_phi : numpy.ndarray(1,:) of float
        Phase angle for every ephemeris table record.
    radial_v : numpy.ndarray(1,:) of float
        Radial velocity in km/s for every ephemeris table record.
    tangential_v : numpy.ndarray(1,:) of float
        Tangential velocity in km/s for every ephemeris table record.
    """

    time = np.zeros(0)
    hae_r = np.zeros(0)
    hae_v = np.zeros(0)

    try:
        with open(ephemeris_file) as file:
            reader = csv.reader(file, delimiter=',')
            for row in reader:
                time = np.append(time,float(row[0]))
                hae_r = np.append(hae_r,[float(row[2]),float(row[3]),float(row[4])])
                hae_v = np.append(hae_v,[float(row[5]),float(row[6]),float(row[7])])
            hae_r = np.reshape(hae_r,((len(hae_r)//3,3)))
            hae_v = np.reshape(hae_v,((len(hae_v)//3,3)))
            r = (hae_r[:,0]**2+hae_r[:,1]**2+hae_r[:,2]**2)**0.5
            hae_phi = np.degrees(np.arctan2(hae_r[:,1],hae_r[:,0]))
            hae_theta = np.degrees(np.arccos((hae_r[:,2]/r[:])))

    except:
        raise LookupError("Unable to load file "+ephemeris_file)

    else:
        #compute radial and tangential velocities
        radial_v = np.zeros(len(hae_r[:,0]))
        tangential_v = np.zeros(len(hae_r[:,0]))
        for i in range(len(hae_r[:,0])):
            unit_radial = hae_r[i,:]/np.linalg.norm(hae_r[i,:])
            radial_v[i] = np.inner(unit_radial,hae_v[i,:])
            tangential_v[i] = np.linalg.norm(hae_v[i,:]-radial_v[i]*unit_radial)     
        return time, hae_r, hae_v, hae_phi, radial_v, tangential_v, hae_theta


def fetch_heliocentric(file):
    """
    This function returns helicoentric distance & phase, in addition to 
    radial and tangential velocities in a form of 1D functions. Be careful, 
    there is extrapolation.

    Parameters
    ----------
    file : str
       The ephemeris file to access.

    Returns
    -------
    f_hel_r : 1D function: float -> float
        heliocentric distance in AU.
    f_hel_phi : 1D function: float -> float
        heliocentric phase angle, measured from the first point of Aries.
    f_rad_v : 1D function: float -> float
        heliocentric radial velocity in km/s.
    f_tan_v : 1D function: float -> float
        heliocentric tangential veloctiy in km/s.
    """

    jd_ephem, hae_r, hae_v, hae_phi, radial_v, tangential_v, hae_theta = load_ephemeris(file)
    heliocentric_distance = np.sqrt(hae_r[:,0]**2+hae_r[:,1]**2+hae_r[:,2]**2)/au #in au
    f_hel_r = interpolate.interp1d(jd_ephem,heliocentric_distance,fill_value="extrapolate",kind=3)
    f_hel_phi = interpolate.interp1d(jd_ephem,hae_phi,fill_value="extrapolate",kind=3)
    f_rad_v = interpolate.interp1d(jd_ephem,radial_v,fill_value="extrapolate",kind=3)
    f_tan_v = interpolate.interp1d(jd_ephem,tangential_v,fill_value="extrapolate",kind=3)
    
    return f_hel_r, f_hel_phi, f_rad_v, f_tan_v














