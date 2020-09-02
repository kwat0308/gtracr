'''
Utility class for gtracr.
Contains conversions between geodesic coordinates and debuggers for the package.
may also contain other things (not known as of now)
'''

import os
import sys
import numpy as np
import pickle

CURRENT_DIR = os.path.dirname(os.path.realpath(__file__))
DATA_DIR = os.path.join(CURRENT_DIR, "data")


def dec_to_dms(lat_dec, lng_dec):
    '''
    Converts geodesic coordinates expressed in decimal notation 
    to DMS notation (degree, minutes, seconds).

    Parameters
    ----------
    - lat_dec (float): the latitude in decimal notation ([-90, 90])
    - lng_dec (float): the longitude in decimal notation ([-180,180])

    Returns
    --------
    - lat_dms (str) : latitude in DMS notation
    - lng_dms (str) : longitude in DMS notation

    The evaluation is performed with reference to the DMS to decimal calculator:
    https://www.rapidtables.com/convert/number/degrees-to-degrees-minutes-seconds.html

    The convention follows the WGS 84 convention as posed here 
    : https://gisgeography.com/latitude-longitude-coordinates/

    The convention shows that:
    - A positive latitude (latitude > 0) is towards the North Pole
    - A negative latitude (latitude < 0) is towards the South Pole
    - A positive longitude (longitude > 0) is towards the East relative to one's location
    - A negative longitude (longitude < 0) is towards the West relative to one's location

    So with this, we see that:
    - latitude goes counter-clockwise from the equator
    - longitude goes counter-clockwise from the Prime Meridian
    '''

    lat_deg = int(np.floor(lat_dec))
    lat_min = int(np.floor((lat_dec - lat_deg) * 60.))
    lat_sec = int(np.floor((lat_dec - lat_deg - (lat_min / 60.)) * 60.))

    # add north or south depending on sign of lat_dec
    lat_symb = "N" if lat_dec >= 0 else "S"

    lat_dms = "{:d}°{:d}\'{:d}\"{:s}".format(lat_deg, lat_min, lat_sec,
                                             lat_symb)

    lng_deg = int(np.floor(lng_dec))
    lng_min = int(np.floor((lng_dec - lng_deg) * 60.))
    lng_sec = int(np.floor((lng_dec - lng_deg - (lng_min / 60.)) * 60.))

    # add east or west depending on sign of lng_dec
    lng_symb = "E" if lat_dec >= 0 else "W"

    lng_dms = "{:d}°{:d}\'{:d}\"{:s}".format(lng_deg, lng_min, lng_sec,
                                             lng_symb)

    return lat_dms, lng_dms

def import_dict(fname):
    '''
    Import the dictionary with filepath fname.
    '''
    with open(fname, "rb") as f:
        the_dict = pickle.load(f)
    return the_dict


def get_particledict():
    '''
    Get the particle dictionary from the .pkl file
    '''
    fpath = os.path.join(DATA_DIR, "particle_dict.pkl")
    particle_dict = import_dict(fpath)
    return particle_dict


def get_locationdict():
    '''
    Get the location dictionary from the .pkl file
    '''
    fpath = os.path.join(DATA_DIR, "location_dict.pkl")
    location_dict = import_dict(fpath)
    return location_dict