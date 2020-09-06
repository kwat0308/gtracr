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
coord_list = [
    (EARTH_RADIUS, np.pi / 2., np.pi),
    (EARTH_RADIUS, 0.5, np.pi),
    (2. * EARTH_RADIUS, np.pi / 2., np.pi),
    (2. * EARTH_RADIUS, np.pi / 2., 2.*np.pi),
    (10. * EARTH_RADIUS, 0., 2.*np.pi),
    (10. * EARTH_RADIUS, 0., np.pi / 4.),
    (10. * EARTH_RADIUS, 0., np.pi / 6.),
    (2. * EARTH_RADIUS, np.pi / 6., np.pi),
    (2.35 * EARTH_RADIUS, (4. *np.pi) / 6., np.pi),
    (5. * EARTH_RADIUS, (4. *np.pi) / 6., (3. * np.pi) / 2.)
]

def test_pydipole():
    '''
    Test the dipole model in the Python version.
    '''

    from gtracr.lib.magnetic_field import MagneticField

    expected_bmag = [
        2.94048e-05,
        5.3501009058777e-05,
        3.6756e-06,
        3.6756e-06,
        5.880960000000001e-08,
        5.880960000000001e-08,
        5.880960000000001e-08,
        6.6262821340477195e-06,
        2.997323835820273e-06,
        3.111911526063683e-07
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
        5.913755366904459e-05,
        5.9118115425034316e-05,
        7.3917929316092776e-06,
        7.39041733659261e-06,
        5.90976751976999e-08,
        5.90976751976999e-08,
        5.90976751976999e-08,
        7.389422687896002e-06,
        4.556819761739493e-06,
        4.7303278534088914e-07
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
        2.94048e-05,
        5.3501009058777e-05,
        3.6756e-06,
        3.6756e-06,
        5.880960000000001e-08,
        5.880960000000001e-08,
        5.880960000000001e-08,
        6.6262821340477195e-06,
        2.997323835820273e-06,
        3.111911526063683e-07
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


