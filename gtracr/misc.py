'''
File for miscellaneous tasks
'''

import os
import sys
import numpy
import pickle

CURRENT_DIR = os.path.dirname(os.path.realpath(__file__))
PARENT_DIR = os.path.dirname(CURRENT_DIR)
sys.path.append(CURRENT_DIR)

DATA_DIR = os.path.join(CURRENT_DIR, "data")


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
