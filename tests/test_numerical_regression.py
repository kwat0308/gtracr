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
    ("p+", 0.0, 25.0, 100.0, 50.0, 100.0, 0.0, 50.0, None),
    ("p+", 90.0, 5.0, 100.0, 89.0, 20.0, 0.0, 20.0, None),
    ("p+", 90.0, 5.0, 100.0, -90.0, 20.0, 0.0, 20.0, None),
    ("e-", 90.0, 5.0, 100.0, 40.0, 200.0, 0.0, 20.0, None),
    ("p+", 45.0, 265.0, 0.0, 40.0, 200.0, 0.0, 20.0, None),
    ("p+", 45.0, 180.0, 10.0, 40.0, 200.0, 0.0, 20.0, None),
    ("p+", 45.0, 0.0, 0.0, 89.0, 0.0, 0.0, 20.0, None),
    ("p+", 45.0, 0.0, 0.0, 0.0, 180.0, 100.0, 20.0, None),
    ("p+", 45.0, 0.0, 0.0, 0.0, 180.0, 100.0, 5.0, None),
    ("p+", 45.0, 0.0, 0.0, 0.0, 180.0, 100.0, None, 10.0),
    ("p+", 9.0, 80.0, 0.0, 50.0, 260.0, 100.0, None, 50.0),
]

EXPECTED_IGRF_SIXVEC = [
                           [6471199.678073804, 1.570796257531502, 0.00046304539704920204, -5.167500353599133e-21, -7.19463695247204e-21, 1.6032856346671467e-17],
                           [6394420.968217134, 1.3997150358997243, 1.97520655907173, -8.558157967292086e-18, -3.0571297538557243e-18, 1.3210428154331233e-17],
                           [63712046.461093366, 0.8081663834552908, 2.178122292796411, 2.6652474353989612e-17, 8.685688014767103e-19, 1.7130168335439596e-18],
                           [63714830.77024167, 0.20296324677064523, -2.0061260905786527, 1.068183157523584e-17, -4.125245282112378e-19, -2.4900884461461534e-19],
                           [63712428.500871226, 4.610658306371075, 62847097628.365524, 3.577595358218272e-06, 3.651999849196411e-07, 1.7355431537059573e-11],
                           [63713797.91604345, 1.4592302227136231, -2.8584920798550146, 1.0596800038021129e-17, 2.771308604460814e-19, -1.3935541162250434e-18],
                           [63714221.61281417, 1.0094570207744513, -2.296455561787531, 1.0575069924813004e-17, 1.0193268898190514e-18, 1.1898605096713222e-18],
                           [63714114.85374857, 1.7949552464577723, -1.7900522033747195, 1.048864491080826e-17, 1.6693199787203039e-18, 1.2378719513167588e-18],
                           [63713931.92615502, 0.35165768266511555, -1.5734197269913082, 1.068630207924647e-17, 9.638443922627212e-20, -2.3818386395450007e-19],
                           [63713577.573224165, 1.392378137488974, 5.830159593088955, 1.0244434104127667e-17, -1.5820866886328317e-19, 3.0518433807088616e-18],
                           [6469415.14434462, 1.7059264959612939, 3.27097543340976, -2.3542226955696143e-18, 7.336257217006058e-19, 1.040420820197456e-18],
                           [6468776.190936486, 1.7882367744243577, 3.499317186168287, -4.649978166385247e-18, -1.0468766737533179e-18, 2.374906992556875e-18],
                           [63713058.96803148, 0.898651727617785, -1.3189597288015469, 2.6660031354156934e-17, 1.0611519656487583e-18, 1.3797361745242354e-18],
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
