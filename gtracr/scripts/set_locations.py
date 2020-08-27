from gtracr.lib.location import Location
'''
Creates the dictionary that contains the latitude, longitude, and altitude 
of different locations, and exports them as a .pkl file for other files
to use.
'''
import os
import sys
import numpy
import pickle

CURRENT_DIR = os.path.dirname(os.path.realpath(__file__))
PARENT_DIR = os.path.dirname(CURRENT_DIR)
# sys.path.append(CURRENT_DIR)

DATA_DIR = os.path.join(PARENT_DIR, "data")


def export_dict(location_dict, fname):
    with open(fname, "wb") as f:
        pickle.dump(location_dict, f, protocol=-1)


def import_dict(fname):
    with open(fname, "rb") as f:
        location_dict = pickle.load(f)
    return location_dict


if __name__ == "__main__":
    fname = "location_dict.pkl"

    fpath = os.path.join(DATA_DIR, fname)

    # some important locations of interest
    locations = [
        Location("Kamioka", 36.434800, 137.276599, -1.),
        Location("IceCube", -89.99, 0., -1.),
        Location("SNOLAB", 46.471983, -81.186701, -2.),
        Location("UofA", 53.523230, -113.526319, 0.)
    ]

    # check if file exists first
    if os.path.exists(fpath):
        location_dict = import_dict(fpath)
    else:
        location_dict = {}

    # add location to location_dict if it does not exist
    for loc in locations:
        if loc.name not in list(location_dict.keys()):
            location_dict[loc.name] = loc
        else:
            continue

    # export back to data directory
    export_dict(location_dict, fpath)
