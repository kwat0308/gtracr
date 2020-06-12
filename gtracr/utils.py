'''
Utility script for various small things (constants, conversions)
'''

import scipy.constants as sc
import numpy as np

SPEED_OF_LIGHT = sc.c

# ELEMENTARY_CHARGE = sc.e

EARTH_RADIUS = 6.371e3  # earth radius (kilometers)
g10 = -29404.8 * (1e-9)  # B-field parameter from IGRF 2020 (in Teslas)

DEG_TO_RAD = np.pi / 180.
RAD_TO_DEG = 180. / np.pi


# converts from spherical coordinates to Cartesian ones
def spherical_to_cartesian(r, theta, phi):
    x = r * np.cos(phi) * np.sin(theta)
    y = r * np.sin(phi) * np.sin(theta)
    z = r * np.cos(theta)
    return np.array([x, y, z])


# get momentum spherical components from magnitude of momentum in 3-D
def get_sphcomp_momentum(p, r, theta):
    pr = p
    ptheta = p / r
    pphi = p / (r * np.sin(theta))
    return (pr, ptheta, pphi)

def get_momentum_magnitude(pr, ptheta, pphi, r, theta):
    return np.sqrt(pr**2. + (r*ptheta)**2. + (r*np.sin(theta)*pphi)**2.)