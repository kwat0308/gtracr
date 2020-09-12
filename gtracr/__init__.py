import os
import sys

CURRENT_DIR = os.path.dirname(os.path.realpath(__file__))
LIB_DIR = os.path.join(CURRENT_DIR, "lib")
# PARENT_DIR = os.path.dirname(CURRENT_DIR)
sys.path.append(CURRENT_DIR)
sys.path.append(LIB_DIR)

# define global dictionaries that will contain the
# location and particle information that will be used
# throughout the code

# from gtracr.utils import set_locationdict, set_particledict

# global location_dict = set_locationdict()
# global particle_dict = set_particledict()
