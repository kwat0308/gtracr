'''
Library that controls the equation for the Earth's magnetic field.
'''

import os, sys
import numpy as np

sys.path.append(os.getcwd())
sys.path.append(os.path.join(os.getcwd(), "gtracr"))

from gtracr.constants import EARTH_RADIUS, g10

# magnetic field (ideal dipole)


def B_r(r, theta, phi):
    return 2. * (EARTH_RADIUS / r)**3. * g10 * np.cos(theta)
    # return 2.*(1./r)**3.*g10*np.cos(theta)


#
def B_theta(r, theta, phi):
    return (EARTH_RADIUS / r)**3. * g10 * np.sin(theta)
    # return (1./r)**3.*g10*np.sin(theta)


def B_phi(r, theta, phi):
    return 0.