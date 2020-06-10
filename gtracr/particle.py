import os, sys
import numpy as np

sys.path.append(os.getcwd())
sys.path.append(os.path.join(os.getcwd(), "gtracr"))

from gtracr.constants import *

class Particle:

    """
    Utility class for cosmic ray particles
    Members:
        - name: the name of the particle (string)
        - pid: the particle id as in pdg (int)
        - mass: the particle rest mass (float) [units of GeV]
        - charge: particle's charge Z (int) [units of elementary charge]
        - label: the shorthand name for the particle (string)

    Notes:
        - PDGID obtained from here: http://pdg.lbl.gov/2007/reviews/montecarlorpp.pdf
        - The mass of the particles are also obtained from PDG

    Example:
    proton: proton = Particle("Proton", 2212, 0.938272, "p+")
    """

    def __init__(self, name, pid, mass, charge, label):
        self.name = name
        self.pid = pid
        self.mass = mass
        self.charge = charge
        self.label = label
    
    # rest energy 
    def energy(self):
        return self.mass*(SPEED_OF_LIGHT**2.)

    # momentum
    def momentum(self, energy):
        return np.sqrt((self.mass*(SPEED_OF_LIGHT**2.))**2. + energy**2.) / SPEED_OF_LIGHT

    # rigidity (R = pc / Ze)
    def rigidity(self, energy):
        return np.sqrt((self.mass*(SPEED_OF_LIGHT**2.))**2. + energy**2.) / (np.abs(self.charge)*ELEMENTARY_CHARGE)
    
    # string represetation for print output
    def __str__(self):
        return "{:s}: PID = {:d}, m = {:f}GeV, Z = {:d}e".format(self.name, self.pid, self.mass, self.charge)

# example using proton
if __name__ == '__main__':
    proton = Particle("Proton", 2122, 0.937272, 1, "p+")
    print(proton)
    print(proton.energy())
    