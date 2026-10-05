"""
Numerical regression tests with tight tolerances.
Golden values were captured from the frozen-field RK4 implementation (2026-03-01).
These tests verify that optimizations do not change physical results.

IGRF_DATE pins the evaluation date so golden values don't drift day-to-day.
Trajectory goldens were regenerated after correcting Gregorian date conversion;
the calendar date and numerical tolerances are unchanged.
"""

from pathlib import Path

import numpy as np
import pytest

from gtracr._libgtracr import IGRF
from gtracr.trajectory import Trajectory

DATA_PATH = str(
    Path(__file__).parent.parent / "src" / "gtracr" / "data" / "igrf13.json"
)

IGRF_DATE = "2026-03-01"  # pinned so golden values don't drift with calendar date

# ---------------------------------------------------------------------------
# B-field component tests
# ---------------------------------------------------------------------------

# (r [m], theta [rad], phi [rad])
BFIELD_COORDS = [
    (6.471e6, 1.570796326794897, 0.000000000000000),
    (6.471e6, 0.785398163397448, 1.570796326794897),
    (6.471e6, 2.356194490192345, 3.141592653589793),
    (6.471e6, 0.523598775598299, -1.047197551196598),
    (6.471e6, 2.617993877991494, 0.523598775598299),
    (7.000e6, 1.570796326794897, 0.785398163397448),
    (8.000e6, 1.047197551196598, -0.523598775598299),
    (6.471e6, 0.100000000000000, 0.500000000000000),
    (6.471e6, 3.041592653589793, 1.500000000000000),
    (6.600e6, 1.570796326794897, 3.141592653589793),
]

# (Br [T], Btheta [T], Bphi [T])
EXPECTED_BFIELD = [
    (1.4855695312886571e-05, -2.6275042651584387e-05, -2.1783012782847714e-06),
    (-5.022698287666843e-05, -2.2136022738173563e-05, 7.670865558988388e-07),
    (4.998055401773263e-05, -1.7598037567259555e-05, 8.735265501917444e-06),
    (-5.159088313804225e-05, -1.0660607633772066e-05, -4.726147369878846e-06),
    (2.9451320469060627e-05, -1.23928974362497e-05, -1.116463946819446e-05),
    (7.098511625973886e-06, -2.428599076702391e-05, -7.740370143861933e-07),
    (-1.3781812451778208e-05, -1.3714018905197889e-05, -2.416042461267865e-06),
    (-5.3446224075424234e-05, -3.7948639669908606e-06, 1.5139925094950205e-06),
    (5.042536466649683e-05, 4.981942197366722e-06, -1.4294509205006488e-05),
    (3.086411652001845e-06, -3.0081399087441946e-05, 5.111542735065936e-06),
]


@pytest.mark.parametrize("idx", range(len(BFIELD_COORDS)))
def test_igrf_bfield_components(idx):
    r, theta, phi = BFIELD_COORDS[idx]
    br_exp, btheta_exp, bphi_exp = EXPECTED_BFIELD[idx]
    igrf = IGRF(DATA_PATH, 2020.0)
    br, btheta, bphi = igrf.values(r, theta, phi)
    assert np.isclose(br, br_exp, rtol=1e-10), f"Br mismatch at coord {idx}"
    assert np.isclose(btheta, btheta_exp, rtol=1e-10), f"Btheta mismatch at coord {idx}"
    assert np.isclose(bphi, bphi_exp, rtol=1e-10), f"Bphi mismatch at coord {idx}"


# ---------------------------------------------------------------------------
# IGRF trajectory sixvector tests
# ---------------------------------------------------------------------------

# (plabel, zenith, azimuth, palt, lat, lng, dalt, rig, energy)
INITIAL_VARIABLES = [
    ("p+", 90.0, 90.0, 100.0, 0.0, 0.0, 0.0, 30.0, None),
    ("p+", 120.0, 90.0, 100.0, 0.0, 0.0, -1.0, 30.0, None),
    ("p+", 0.0, 155.0, 100.0, 50.0, 100.0, 0.0, 50.0, None),
    ("p+", 90.0, 175.0, 100.0, 89.0, 20.0, 0.0, 20.0, None),
    ("p+", 90.0, 175.0, 100.0, -90.0, 20.0, 0.0, 20.0, None),
    ("e-", 90.0, 175.0, 100.0, 40.0, 200.0, 0.0, 20.0, None),
    ("p+", 45.0, 275.0, 0.0, 40.0, 200.0, 0.0, 20.0, None),
    ("p+", 45.0, 0.0, 10.0, 40.0, 200.0, 0.0, 20.0, None),
    ("p+", 45.0, 180.0, 0.0, 89.0, 0.0, 0.0, 20.0, None),
    ("p+", 45.0, 180.0, 0.0, 0.0, 180.0, 100.0, 20.0, None),
    ("p+", 45.0, 180.0, 0.0, 0.0, 180.0, 100.0, 5.0, None),
    ("p+", 45.0, 180.0, 0.0, 0.0, 180.0, 100.0, None, 10.0),
    ("p+", 9.0, 100.0, 0.0, 50.0, 260.0, 100.0, None, 50.0),
]

EXPECTED_IGRF_SIXVEC = [
                           [6471199.517110254, 1.5707962229024632, 0.0004630453910298668, -5.167561417742061e-21, -7.194403601679309e-21, 1.6032855531086903e-17],
                           [6394645.886638521, 1.3997609474866473, 1.9748336958845865, -8.553805939763706e-18, -3.059500442273911e-18, 1.3210768181011588e-17],
                           [63712058.12519359, 0.8081886540372633, 2.178093999442155, 2.6652315327407214e-17, 8.686393620383588e-19, 1.7130010611797745e-18],
                           [63712200.20608276, 0.20252942311566047, -2.007063219296211, 1.067769304047924e-17, -4.1268604633239844e-19, -2.492642104543562e-19],
                           [63713785.92970307, 2.6757983740299207, -0.03139613660959638, 1.0669670861043922e-17, -3.940716760807771e-19, -4.983997564699268e-19],
                           [63714464.22404205, 1.459081267639741, -2.8585222373986885, 1.0593744235089239e-17, 2.7714528666840403e-19, -1.3933214481024215e-18],
                           [63713952.99241918, 1.0095830827950918, -2.296389987185877, 1.0573122748200227e-17, 1.0194153354631925e-18, 1.189722157254854e-18],
                           [63712988.451991454, 1.7952143195601666, -1.790178661644296, 1.0484644724397159e-17, 1.668911793934461e-18, 1.237953760601734e-18],
                           [63714144.86715386, 0.35150931147524134, -1.5738660347638471, 1.0685480131342008e-17, 9.633064375278565e-20, -2.3833979088655394e-19],
                           [63712447.55720404, 1.3923767558552427, -0.45185530842511, 1.0242398569969793e-17, -1.5673106271703685e-19, 3.051930207919755e-18],
                           [6470348.094086879, 1.7058454152230478, -3.0124865833331866, -2.346691842293289e-18, 7.380818897647043e-19, 1.043369883693765e-18],
                           [6468905.210524317, 1.7882786005498825, -2.7841826806617416, -4.645995014056802e-18, -1.0411698697217315e-18, 2.3750046091697656e-18],
                           [63713064.911237076, 0.8986707511061053, -1.3190057400505286, 2.66599609056324e-17, 1.0612470315444263e-18, 1.3797162849404555e-18],
                       ]

EXPECTED_IGRF_ESCAPED = [False, False, True, True, True, True, True, True, True, True, False, False, True]


@pytest.mark.parametrize("idx", range(len(INITIAL_VARIABLES)))
def test_igrf_sixvector(idx):
    if EXPECTED_IGRF_ESCAPED[idx]:
        pytest.skip(
            "Exit point on escape sphere is architecture-sensitive; use test_igrf_escaped_flag instead"
        )
    plabel, zenith, azimuth, palt, lat, lng, dalt, rig, en = INITIAL_VARIABLES[idx]
    traj = Trajectory(
        plabel=plabel,
        zenith_angle=zenith,
        azimuth_angle=azimuth,
        particle_altitude=palt,
        latitude=lat,
        longitude=lng,
        detector_altitude=dalt,
        rigidity=rig,
        energy=en,
        bfield_type="igrf",
        date=IGRF_DATE,
    )
    traj.get_trajectory(dt=1e-5, max_time=1.0)
    assert np.allclose(traj.final_sixvector, EXPECTED_IGRF_SIXVEC[idx], rtol=1e-4), (
        f"sixvector mismatch at case {idx}"
    )


@pytest.mark.parametrize("idx", range(len(INITIAL_VARIABLES)))
def test_igrf_escaped_flag(idx):
    plabel, zenith, azimuth, palt, lat, lng, dalt, rig, en = INITIAL_VARIABLES[idx]
    traj = Trajectory(
        plabel=plabel,
        zenith_angle=zenith,
        azimuth_angle=azimuth,
        particle_altitude=palt,
        latitude=lat,
        longitude=lng,
        detector_altitude=dalt,
        rigidity=rig,
        energy=en,
        bfield_type="igrf",
        date=IGRF_DATE,
    )
    traj.get_trajectory(dt=1e-5, max_time=1.0)
    assert traj.particle_escaped == EXPECTED_IGRF_ESCAPED[idx], (
        f"escaped flag mismatch at case {idx}: got {traj.particle_escaped}, expected {EXPECTED_IGRF_ESCAPED[idx]}"
    )
