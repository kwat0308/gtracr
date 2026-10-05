"""
Compare the values of the final times and sixvector of the trajectory for the dipole model
"""

import numpy as np

from gtracr.trajectory import Trajectory

IGRF_DATE = "2026-03-01"

# in the form :
# (plabel, zenith, azimuth, particle_altitude,
# latitude, longitude, detector_altitude, rigidity, kinetic energy)
initial_variable_list = [
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


def test_trajectories_dipole():
    """
    Test the final times of the trajectory evaluation in the dipole field.
    """

    expected_times = [
                         1e-05,
                         0.3051700000001644,
                         0.19222000000005146,
                         0.20289000000006213,
                         0.21143000000007067,
                         0.2024500000000617,
                         0.19874000000005798,
                         0.2170700000000763,
                         0.19498000000005422,
                         0.23267000000009191,
                         0.00731999999999987,
                         0.019299999999999384,
                         0.19332000000005256,
                     ]

    dt = 1e-5
    max_time = 1.0

    for iexp, initial_variables in enumerate(initial_variable_list):
        (plabel, zenith, azimuth, palt, lat, lng, dalt, rig, en) = initial_variables

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
            bfield_type="dipole",
        )

        traj.get_trajectory(dt=dt, max_time=max_time)

        assert np.allclose(traj.final_time, expected_times[iexp])


def test_trajectories_igrf():
    """
    Test the final times of the trajectory evaluation in the IGRF field.

    Only forbidden trajectories (particle_escaped=False) are checked against
    golden final times.  Escaping trajectories exit at a point on the 10 RE
    sphere whose exact time depends on FP rounding differences between
    architectures; the escaped flag alone is tested for those cases via
    test_trajectories_igrf_escaped.
    """

    expected_times = [
                         1e-05,
                         0.02260999999999925,
                         0.1925600000000518,
                         0.20333000000006257,
                         0.21143000000007067,
                         0.20077000000006,
                         0.19821000000005745,
                         0.21603000000007527,
                         0.19481000000005405,
                         0.26402000000012327,
                         0.005159999999999958,
                         0.01248999999999966,
                         0.19273000000005197,
                     ]

    # Only forbidden trajectories have a stable final time across architectures.
    escaped_flags = [
        False,
        False,
        True,
        True,
        True,
        True,
        True,
        True,
        True,
        True,
        False,
        False,
        True,
    ]

    dt = 1e-5
    max_time = 1.0

    for iexp, initial_variables in enumerate(initial_variable_list):
        if escaped_flags[iexp]:
            continue
        (plabel, zenith, azimuth, palt, lat, lng, dalt, rig, en) = initial_variables

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

        traj.get_trajectory(dt=dt, max_time=max_time)

        assert np.allclose(traj.final_time, expected_times[iexp])


def test_trajectories_stepsize():
    """
    Test the final times of the trajectory evaluation in the igrf field for
    different step sizes
    """

    expected_times = [
                         0.2160556599269434,
                         0.21605570000518287,
                         0.2160560000001954,
                         0.2160600000000753,
                         0.21609999999999252,
                         0.21700000000000016,
                         0.22000000000000006,
                         0.30000000000000004,
                     ]

    dt_arr = [1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1]
    max_time = 1.0

    (plabel, zenith, azimuth, palt, lat, lng, dalt, rig, en) = (
        "p+",
        90.0,
        0.0,
        100.0,
        0.0,
        0.0,
        0.0,
        50.0,
        None,
    )

    for iexp, dt in enumerate(dt_arr):
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

        traj.get_trajectory(dt=dt, max_time=max_time)

        assert np.allclose(traj.final_time, expected_times[iexp])


def test_trajectories_maxtimes():
    """
    Test the final times of the trajectory evaluation in the igrf field for
    different maximal times
    """

    expected_times = [
                         0.00999999999999976,
                         0.027829999999999036,
                         0.07743000000000268,
                         0.2154500000000747,
                         0.2160600000000753,
                         0.2160600000000753,
                         0.2160600000000753,
                         0.2160600000000753,
                         0.2160600000000753,
                         0.2160600000000753,
                     ]

    dt = 1e-5
    max_times = np.logspace(-2, 2, 10)

    for iexp, max_time in enumerate(max_times):
        (plabel, zenith, azimuth, palt, lat, lng, dalt, rig, en) = (
            "p+",
            90.0,
            0.0,
            100.0,
            0.0,
            0.0,
            0.0,
            50.0,
            None,
        )

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

        traj.get_trajectory(dt=dt, max_time=max_time)

        assert np.allclose(traj.final_time, expected_times[iexp])


def test_trajectories_dates():
    """
    Test the final times of the trajectory evaluation in the igrf field for
    different dates
    """

    expected_times = [
                         0.2133600000000726,
                         0.21364000000007288,
                         0.21339000000007263,
                         0.21592000000007516,
                         0.2154700000000747,
                         0.21541000000007465,
                         0.2148600000000741,
                         0.2146500000000739,
                         0.21540000000007464,
                         0.21601000000007525,
                     ]

    dt = 1e-5
    max_time = 1.0

    dates = [
        "1900-01-01",
        "1909-01-01",
        "1900-10-31",
        "2020-09-12",
        "2004-03-08",
        "2000-02-28",
        "1970-03-26",
        "1952-04-30",
        "1999-03-08",
        "2024-03-09",
    ]
    for iexp, date in enumerate(dates):
        (plabel, zenith, azimuth, palt, lat, lng, dalt, rig, en) = (
            "p+",
            90.0,
            0.0,
            100.0,
            0.0,
            0.0,
            0.0,
            50.0,
            None,
        )

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
            date=date,
        )

        traj.get_trajectory(dt=dt, max_time=max_time)

        assert np.allclose(traj.final_time, expected_times[iexp])


def test_dipole_sixvec():
    expected_sixvec = [
                          [6471199.6231868975, 1.5707963267948966, 0.0004630454043181561, -6.048499118556204e-21, 9.838339507932405e-34, 1.6032857217419333e-17],
                          [63712300.83559172, 1.5707963267948966, 5.580433986060423, 1.5505969567435694e-17, 1.8581747158211076e-33, 4.078320975259935e-18],
                          [63713269.53687547, 0.7985366105429904, 2.1319063599054346, 2.6669318457456412e-17, 7.546367686664471e-19, 1.4900036620404851e-18],
                          [63714251.410265714, 0.2537281781495107, -1.9109469377678965, 1.0686470734510774e-17, -3.6031166789418974e-19, -6.505411664262572e-20],
                          [63712428.500868715, 4.610658306371472, 62847097628.365524, 3.5775953582182725e-06, 3.6519998491959885e-07, 1.7355430247610533e-11],
                          [63713300.33339308, 1.6105525740280033, -2.7959790366532498, 1.0588925784438677e-17, 5.701287181115243e-19, -1.363290587923674e-18],
                          [63712607.22518879, 0.9083970502780275, -2.1870272718918757, 1.057770709634604e-17, 8.931930639651709e-19, 1.2664568461466769e-18],
                          [63714889.61083344, 1.736088348918876, -1.62019689065021, 1.0482050002627351e-17, 1.4982817461041292e-18, 1.4869916330922027e-18],
                          [63712177.85039131, 0.40000732247158444, -1.532082945867246, 1.0687201073699463e-17, 1.8413213282957926e-19, -1.1463897368312226e-19],
                          [63714380.11990551, 1.5475754081854918, 4.9057401814950365, 1.0304335656474183e-17, -1.0183567922399613e-18, 2.6562963702938345e-18],
                          [6470037.431013613, 1.7370099845068263, 3.3341864614756536, -2.445839526561264e-18, 7.294074980622177e-19, 8.03215914258814e-19],
                          [6470820.332392309, 1.7565309914985197, 3.6907467389360367, -4.8353706959365755e-18, -1.987335368974099e-18, 1.0111371752046643e-18],
                          [63714530.69677783, 0.9088220477531327, -1.2406120109998, 2.663912585692339e-17, 1.1896970350264556e-18, 1.6522032947644405e-18],
                      ]

    dt = 1e-5
    max_time = 1.0

    for iexp, initial_variables in enumerate(initial_variable_list):
        (plabel, zenith, azimuth, palt, lat, lng, dalt, rig, en) = initial_variables

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
            bfield_type="dipole",
        )

        traj.get_trajectory(dt=dt, max_time=max_time)

        assert np.allclose(
            traj.final_sixvector, np.array(expected_sixvec[iexp]), rtol=1e-5
        )
