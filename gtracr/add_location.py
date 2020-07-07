'''
Creates dictionary so that we can add locations to global dictionary
Essentially the same thing with add_particle.py

As with add_particle.py, there should be a more systematic way to do this.
'''
import os, sys
import numpy as np

sys.path.append(os.getcwd())
sys.path.append(os.path.join(os.getcwd(), "gtracr"))

# from gtracr.lib.location import Location
from Location import Location

location_dict = {}

def add_to_dict(Loc):
    location_dict[loc.name] = Loc

# some basic locations of interest
kamioka = Location("Kamioka", 36.434800, 137.276599, 0.)
icecube = Location("IceCube",-89.99, 0., 0.)
uofa = Location("UofA", 53.523230, -113.526319, 0.)

# print(kamioka.name())

for loc in [kamioka, icecube, uofa]:
    # print(loc)
    # print(loc.name())
    add_to_dict(loc)
    # print(loc)

# print(location_dict)