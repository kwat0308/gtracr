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
import argparse

CURRENT_DIR = os.path.dirname(os.path.realpath(__file__))
PARENT_DIR = os.path.dirname(CURRENT_DIR)
# sys.path.append(CURRENT_DIR)

DATA_DIR = os.path.join(PARENT_DIR, "data")

FILE_NAME = "location_dict.pkl"

FILE_PATH = os.path.join(DATA_DIR, FILE_NAME)


def export_dict(location_dict, fname):
    # set protocol=3 for < py3.8 compatibility
    protocol = 3

    with open(fname, "wb") as f:
        pickle.dump(location_dict, f, protocol=protocol)


def import_dict(fname):
    with open(fname, "rb") as f:
        location_dict = pickle.load(f)
    return location_dict


def set_locations(locations, args):
    # check if file exists first and
    # reset dictionary if we want to clean it
    if os.path.exists(FILE_PATH) and not args.reset_dict:
        location_dict = import_dict(FILE_PATH)
    else:
        location_dict = {}

    # if clean_dict:
    #     location_dict = {}

    # add location to location_dict if it does not exist
    for loc in locations:
        if loc.name not in list(location_dict.keys()):
            location_dict[loc.name] = loc
        else:
            continue

    # export back to data directory
    export_dict(location_dict, FILE_PATH)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=
        'Set the data file that contains geodesic coordinate data of different locations.'
    )
    parser.add_argument('-r',
                        '--reset',
                        dest="reset_dict",
                        action="store_true",
                        help='Clean the dataset and create from scratch.')

    args = parser.parse_args()

    # some important locations of interest
    # mostly neutrino telescopes, detectors, and telescope arrays

    # below is alpha version with depth support
    # locations = [
    #     Location("Kamioka", 36.434800, 137.276599, -1.),
    #     Location("IceCube", -89.99, 0., -1.),
    #     Location("SNOLAB", 46.471983, -81.186701, -2.),
    #     Location("UofA", 53.523230, -113.526319, 0.),
    #     Location("CTA-North", 28.76216, -17.89201, 2.2),
    #     Location("CTA-South", -24.68342778, -24.68342778, 0.),
    #     Location("ORCA", 42.80000000, 6.0333333, -2.475),
    #     Location("ANTARES", 42.8, 6.1666666, -0.35),
    #     Location("Baikal-GVD", 51.77139, 104.3978, -1.366),
    #     Location("TA", 39.208099, -113.121285, 1.4)
    # # ]

    # below is working version without depth support
    locations = [
        Location("Kamioka", 36.434800, 137.276599, 0.),
        Location("IceCube", -89.99, 0., 0.),
        Location("SNOLAB", 46.471983, -81.186701, 0.),
        Location("UofA", 53.523230, -113.526319, 0.),
        Location("CTA-North", 28.76216, -17.89201, 0.),
        Location("CTA-South", -24.68342778, -24.68342778, 0.),
        Location("ORCA", 42.80000000, 6.0333333, 0.),
        Location("ANTARES", 42.8, 6.1666666, 0.),
        Location("Baikal-GVD", 51.77139, 104.3978, 0.),
        Location("TA", 39.208099, -113.121285, 0.)
    ]

    set_locations(locations, args)
