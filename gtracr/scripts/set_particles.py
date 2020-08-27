from gtracr.lib.particle import Particle
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


def export_dict(particle_dict, fname):
    with open(fname, "wb") as f:
        pickle.dump(particle_dict, f, protocol=-1)


def import_dict(fname):
    with open(fname, "rb") as f:
        particle_dict = pickle.load(f)
    return particle_dict


if __name__ == "__main__":
    fname = "particle_dict.pkl"

    fpath = os.path.join(DATA_DIR, fname)

    # some important locations of interest
    particles = [
        Particle("positron", -11, 0.5109 * (1e-3), 1, "e+"),
        Particle("electron", 11, 0.5109 * (1e-3), -1, "e-"),
        Particle("proton", 2212, 0.937272, 1, "p+"),
        Particle("anti-proton", -2212, 0.937272, -1, "p-")
    ]

    # check if file exists first
    if os.path.exists(fpath):
        particle_dict = import_dict(fpath)
    else:
        particle_dict = {}

    # add location to particle_dict if it does not exist
    for particle in particles:
        if particle.name not in list(particle_dict.keys()):
            particle_dict[particle.label] = particle
        else:
            continue

    # export back to data directory
    export_dict(particle_dict, fpath)
