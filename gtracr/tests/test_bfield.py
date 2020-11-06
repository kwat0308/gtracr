'''
Test the magnetic field models (dipole and igrf) by comparing the 
magnetic field magnitude values
'''
import os
import sys
import numpy as np
import pytest

# from gtracr.lib._libgtracr import IGRF
# import gtracr.lib.igrf_utils as iuf
# from gtracr.lib.magnetic_field import MagneticField, IGRF13
from gtracr.lib.constants import EARTH_RADIUS

CURRENT_DIR = os.path.dirname(os.path.realpath(__file__))
PARENT_DIR = os.path.dirname(CURRENT_DIR)
DATA_DIR = os.path.join(PARENT_DIR, "data")

# import gauss coefficients from shc file
CURRENT_YEAR = 2015

LEAP_YEAR = False
# check if we have a leap year or not
if CURRENT_YEAR % 4 == 0:
    LEAP_YEAR = True

# r, theta, phi values of interest
coord_list = [(EARTH_RADIUS, np.pi / 2., np.pi), (EARTH_RADIUS, 0.5, np.pi),
              (2. * EARTH_RADIUS, np.pi / 2., np.pi),
              (2. * EARTH_RADIUS, np.pi / 2., 2. * np.pi),
              (10. * EARTH_RADIUS, 0., 2. * np.pi),
              (10. * EARTH_RADIUS, 0., np.pi / 4.),
              (10. * EARTH_RADIUS, 0., np.pi / 6.),
              (2. * EARTH_RADIUS, np.pi / 6., np.pi),
              (2.35 * EARTH_RADIUS, (4. * np.pi) / 6., np.pi),
              (5. * EARTH_RADIUS, (4. * np.pi) / 6., (3. * np.pi) / 2.)]


def test_pydipole():
    '''
    Test the dipole model in the Python version.
    '''

    from gtracr.lib.magnetic_field import MagneticField

    expected_bmag = [
        2.94048000e-05, 5.35010091e-05, 3.67560000e-06, 3.67560000e-06,
        5.88096000e-08, 5.88096000e-08, 5.88096000e-08, 6.62628213e-06,
        2.99732384e-06, 3.11191153e-07
    ]

    pydip = MagneticField()

    for iexp, coord in enumerate(coord_list):
        bf_values = pydip.values(*coord)
        bmag = np.linalg.norm(np.array(bf_values))

        assert np.allclose(bmag, expected_bmag[iexp])


def test_pyigrf():
    '''
    Test the IGRF model in the Python version.
    '''
    from gtracr.lib.magnetic_field import IGRF13

    expected_bmag = [
        3.42851920e-05, 5.42888711e-05, 4.07163707e-06, 3.45071901e-06,
        5.97223497e-08, 5.97223497e-08, 5.97223497e-08, 6.83667250e-06,
        3.42108503e-06, 2.67410913e-07
    ]

    pyigrf = IGRF13(CURRENT_YEAR)

    for iexp, coord in enumerate(coord_list):
        bf_values = pyigrf.values(*coord)
        bmag = np.linalg.norm(np.array(bf_values))

        assert np.allclose(bmag, expected_bmag[iexp])


def test_dipole():
    '''
    Test the Python model in the C++ version.
    '''
    from gtracr.lib._libgtracr import MagneticField

    expected_bmag = [
        2.94048000e-05, 5.35010091e-05, 3.67560000e-06, 3.67560000e-06,
        5.88096000e-08, 5.88096000e-08, 5.88096000e-08, 6.62628213e-06,
        2.99732384e-06, 3.11191153e-07
    ]

    dip = MagneticField()

    for iexp, coord in enumerate(coord_list):
        bf_values = dip.values(*coord)
        bmag = np.linalg.norm(np.array(bf_values))

        assert np.allclose(bmag, expected_bmag[iexp])


# def test_igrf():
#     '''
#     Test the IGRF model in the C++ version.
#     '''
#     from gtracr.lib._libgtracr import IGRF

#     coord_list = [
#         # (EARTH_RADIUS, np.pi / 2., np.pi),
#         (EARTH_RADIUS, 0.5, np.pi),
#         (2. * EARTH_RADIUS, np.pi / 2., np.pi),
#         (2. * EARTH_RADIUS, np.pi / 2., 2.*np.pi),
#         # (10. * EARTH_RADIUS, 0., 2.*np.pi),
#         (10. * EARTH_RADIUS, 0., np.pi / 4.),
#         (10. * EARTH_RADIUS, 0., np.pi / 6.),
#         # (2. * EARTH_RADIUS, np.pi / 6., np.pi),
#         (2.35 * EARTH_RADIUS, (4. *np.pi) / 6., np.pi),
#         (5. * EARTH_RADIUS, (4. *np.pi) / 6., (3. * np.pi) / 2.)
#     ]

#     # in nT
#     expected_bmag = [
#         # 18541.77277708073,
#         # 2104.4220463334086,
#         # 1681.845888088114,
#         # 919.3370935769924,
#         # # 5.41646351617013,
#         # 5.538620048328442,
#         # 5.4857208582223365,
#         # # 494.4790949331748,
#         996.3453533330078,
#         64.65244708389457
#     ]

#     expected_bmag = np.array(expected_bmag)

#     DATA_PATH = os.path.join(DATA_DIR, "igrf13.json")

#     igrf = IGRF(DATA_PATH, CURRENT_YEAR)

#     for iexp, coord in enumerate(coord_list):
#         bf_values = igrf.values(*coord)
#         bmag = np.linalg.norm(np.array(bf_values))

#         assert np.allclose(bmag, expected_bmag[iexp])
