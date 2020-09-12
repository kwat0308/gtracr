from gtracr.lib.constants import SPEED_OF_LIGHT
import os
import sys
import numpy as np


class Particle:
    """
    Utility class for cosmic ray particles

    Parameters
    ------------

    - name : str
        the name of the particle
    - pid : int
        the particle id as in the PDG database
    - mass : float
        the particle's rest mass in GeV
    - charge : int
        the particle's charge in e
    - label : str
        the shorthand name for the particle

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
        self.momentum = 0.
        self.rigidity = 0.

    def set_from_energy(self, energy):
        self.momentum = np.sqrt(energy**2. - self.mass**2.)
        self.rigidity = self.momentum / np.abs(self.charge)

    def set_from_rigidity(self, rigidity):
        self.momentum = rigidity * np.abs(self.charge)
        self.rigidity = rigidity

    def set_from_momentum(self, momentum):
        self.momentum = momentum
        self.rigidity = self.momentum / np.abs(self.charge)

    def get_energy_rigidity(self):
        return np.sqrt((self.rigidity * np.abs(self.charge))**2. +
                       self.mass**2.) + self.mass

    # string represetation for print output
    def __str__(self):
        return "{:s}: PID = {:d}, m = {:.6f}GeV, Z = {:d}e \n Momentum = {:.6e}, Rigidity = {:.6e}".format(
            self.name, self.pid, self.mass, self.charge, self.momentum,
            self.rigidity)


# example using proton
if __name__ == '__main__':
    proton = Particle("Proton", 2122, 0.937272, 1, "p+")
    print(proton)
    print(proton.momentum)
