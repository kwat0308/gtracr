'''
Adds particles and creates a particle dictionary 
'''

import os, sys
import numpy as np

sys.path.append(os.getcwd())
sys.path.append(os.path.join(os.getcwd(), "gtracr"))

# from gtracr.lib.particle import Particle
from Particle import Particle

particle_dict = {}

def add_to_dict(part):
    particle_dict[part.label] = part


# create particles
ep = Particle("positron", -11, 0.5109*(1e-3), 1, "e+")  
em = Particle("electron", 11, 0.5109*(1e-3), -1, "e-")
pp = Particle("Proton", 2122, 0.937272, 1, "p+")
pm = Particle("anti-proton", -2122, 0.937272, -1, "p-")

for part in [ep, em, pp, pm]:
    add_to_dict(part)

# print(Particle())
# print(particle_dict)