'''
Utility script for various small things (conversions, info function for debugging)
'''

# import scipy.constants as sc
import os, sys
import numpy as np

sys.path.append(os.getcwd())
sys.path.append(os.path.join(os.getcwd(), "gtracr"))

from gtracr.constants import SPEED_OF_LIGHT

# SPEED_OF_LIGHT = sc.c

# # ELEMENTARY_CHARGE = sc.e

# EARTH_RADIUS = 6.371e3  # earth radius (kilometers)
# g10 = -29404.8 * (1e-9)  # B-field parameter from IGRF 2020 (in Teslas)

# DEG_TO_RAD = np.pi / 180.
# RAD_TO_DEG = 180. / np.pi


def gamma(v):
    return np.reciprocal(np.sqrt(1 - (v / SPEED_OF_LIGHT)**2.))


def SphCoord_to_CarCoord(r, theta, phi):
    x = r * np.cos(phi) * np.sin(theta)
    y = r * np.sin(phi) * np.sin(theta)
    z = r * np.cos(theta)
    return np.array([x, y, z])


def CarCoord_to_SphCoord(x, y, z):
    r = np.sqrt(x**2. + y**2. + z**2.)
    theta = np.arccos(z / r)
    phi = np.arctan2(y, x)
    return np.array([r, theta, phi])


def SphVel_to_CarVel(vr, vtheta, vphi, r, theta, phi):
    vx = vr * np.sin(theta) * np.cos(phi) + r * vtheta * np.cos(
        theta) * np.cos(phi) - r * vphi * np.sin(theta) * np.sin(phi)
    vy = vr * np.sin(theta) * np.sin(phi) + r * vtheta * np.cos(
        theta) * np.sin(phi) + r * vphi * np.sin(theta) * np.cos(phi)
    vz = vr * np.cos(theta) - r * vtheta * np.sin(theta)

    return np.array([vx, vy, vz])


def CarVel_to_SphVel(vx, vy, vz, r, theta, phi):
    vr = vx * np.sin(theta) * np.cos(phi) + vy * np.sin(theta) * np.sin(
        phi) + vz * np.cos(theta)
    vtheta = (vx * np.cos(theta) * np.cos(phi) +
              vy * np.cos(theta) * np.sin(phi) - vz * np.sin(theta)) / r
    vphi = (-vx * np.sin(phi) + vy * np.cos(phi)) / (r * np.sin(theta))

    return np.array([vr, vtheta, vphi])


# # converts from spherical coordinates to Cartesian ones
# def spherical_to_cartesian(r, theta, phi):
#     x = r * np.cos(phi) * np.sin(theta)
#     y = r * np.sin(phi) * np.sin(theta)
#     z = r * np.cos(theta)
#     return np.array([x, y, z])

# # get momentum spherical components from magnitude of momentum in 3-D
# def vp_components_spherical(p, r, theta):
#     pr = p
#     ptheta = p / r
#     pphi = p / (r * np.sin(theta))
#     return (pr, ptheta, pphi)


# magnitude of velocity from spherical components
def vmag_spherical(vr, vtheta, vphi, r, theta):
    return np.sqrt(vr**2. + (r * vtheta)**2. + (r * np.sin(theta) * vphi)**2.)


# # info function
