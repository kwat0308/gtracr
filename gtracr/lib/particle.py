import os, sys
import numpy as np

sys.path.append(os.getcwd())
sys.path.append(os.path.join(os.getcwd(), "gtracr"))

from gtracr.utils import SPEED_OF_LIGHT


class Particle:
    """
    Utility class for cosmic ray particles
    Members:
    - name: the name of the particle (string)
    - pid: the particle id as in pdg (int)
    - mass: the particle rest mass (float) [units of GeV / c^2]
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
        # set properties as members
        self.momentum = 0.
        self.velocity = 0.
        self.rigidity = self.momentum / np.abs(self.charge)

    # momentum [units GeV/c]
    def set_momentum(self, energy):
        self.momentum = np.sqrt(self.mass**2. + energy**2.)
        return self.momentum

    # velocity [units in m/s]
    def set_velocity(self):
        self.velocity = ((self.momentum * SPEED_OF_LIGHT) /
                         np.sqrt(self.momentum**2. +
                                 (self.mass * SPEED_OF_LIGHT)**2.))
        return self.velocity

    # rigidity (R = pc / Ze) from energy [units GV]
    def set_rigidity_from_energy(self, energy):
        self.rigidity = self.set_momentum(energy) / (np.abs(self.charge))

    # rigidity (R = pc / Ze) from momentum [units GV]
    def set_rigidity_from_momentum(self):
        self.rigidity = self.momentum / (np.abs(self.charge))

    # string represetation for print output
    def __str__(self):
        return "{:s}: PID = {:d}, m = {:f}GeV, Z = {:d}e".format(
            self.name, self.pid, self.mass, self.charge)


# example using proton
if __name__ == '__main__':
    proton = Particle("Proton", 2122, 0.937272, 1, "p+")
    print(proton)
    print(proton.momentum)
