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
import argparse

CURRENT_DIR = os.path.dirname(os.path.realpath(__file__))
PARENT_DIR = os.path.dirname(CURRENT_DIR)
# sys.path.append(CURRENT_DIR)

DATA_DIR = os.path.join(PARENT_DIR, "data")
FILE_NAME = "particle_dict.pkl"
FILE_PATH = os.path.join(DATA_DIR, FILE_NAME)


def export_dict(particle_dict):
    # set protocol=3 for < py3.8 compatibility
    protocol = 3

    with open(FILE_PATH, "wb") as f:
        pickle.dump(particle_dict, f, protocol=protocol)


def import_dict():
    with open(FILE_PATH, "rb") as f:
        particle_dict = pickle.load(f)
    return particle_dict

def set_particles(particles, args):

    # check if file exists first
    if os.path.exists(FILE_PATH) and not args.reset_dict:
        particle_dict = import_dict()
    else:
        particle_dict = {}

    # add location to particle_dict if it does not exist
    for particle in particles:
        if particle.name not in list(particle_dict.keys()):
            particle_dict[particle.label] = particle
        else:
            continue

    # export back to data directory
    export_dict(particle_dict)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
    description=
    'Set the data file that contains different particles for the cosmic rays.'
    )
    parser.add_argument('-r',
                        '--reset',
                        dest="reset_dict",
                        action="store_true",
                        help='Clean the dataset and create from scratch.')
    
    
    # different cosmic ray particles
    particles = [
        Particle("positron", -11, 0.5109 * (1e-3), 1, "e+"),
        Particle("electron", 11, 0.5109 * (1e-3), -1, "e-"),
        Particle("proton", 2212, 0.937272, 1, "p+"),
        Particle("anti-proton", -2212, 0.937272, -1, "p-")
    ]

    args = parser.parse_args()

    set_particles(particles, args)
    
