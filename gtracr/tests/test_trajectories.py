'''
Compare the values of the final times and sixvector of the trajectory for the dipole model
'''

import os
import sys
import numpy as np
import pytest

from gtracr.trajectory import Trajectory
from gtracr.lib.trajectorypoint import TrajectoryPoint

# in the form :
# (plabel, zenith, azimuth, particle_altitude, 
# latitude, longitude, detector_altitude, rigidity, kinetic energy)
initial_variable_list = [
    ("p+", 90., 90., 100., 0., 0., 0., 30., None),  
    ("p+", 120., 90., 100., 0., 0., -1., 30., None),
    ("p+", 0., 25., 100., 50., 100., 0., 50., None),
    ("p+", 90., 5., 100., 89., 20., 0., 20., None),
    ("p+", 90., 5., 100., -90., 20., 0., 20., None),
    ("e-", 90., 5., 100., 40., 200., 0., 20., None),
    ("p+", 45., 265., 0., 40., 200., 0., 20., None),
    ("p+", 45., 180., 10., 40., 200., 0., 20., None),
    ("p+", 45., 0., 0., 89., 0., 0., 20., None),
    ("p+", 45., 0., 0., 0., 180., 100., 20., None),
    ("p+", 45., 0., 0., 0., 180., 100., 5., None),
    ("p+", 45., 0., 0., 0., 180., 100., None, 10.),
    ("p+", 9., 80., 0., 50., 260., 100., None, 50.),
]

def test_trajectories_dipole():
    '''
    Test the final times of the trajectory evaluation in the dipole field.
    '''

    expected_times= [
        0.004199999999999997,
        0.29988000000015913,
        0.19221000000005145,
        0.20288000000006212,
        0.21143000000007067,
        0.2024600000000617,
        0.19868000000005792,
        0.2169700000000762,
        0.19499000000005423,
        0.23232000000009156,
        0.0076999999999998545,
        0.019769999999999364,
        0.19331000000005255
    ]

    dt = 1e-5
    max_time = 1.

    for iexp, initial_variables in enumerate(initial_variable_list):

        (plabel, zenith, azimuth, palt, lat, lng, dalt, rig, en) = initial_variables
        
        traj = Trajectory(
            plabel=plabel,
            zenith_angle = zenith,
            azimuth_angle = azimuth,
            particle_altitude = palt,
            latitude = lat,
            longitude = lng,
            detector_altitude = dalt,
            rigidity = rig,
            energy = en,
            bfield_type = "dipole"
        )
        
        traj.get_trajectory(dt = dt, max_time=max_time)

        assert np.allclose(traj.final_time, expected_times[iexp])

def test_trajectories_igrf():
    '''
    Test the final times of the trajectory evaluation in the IGRF field.
    '''
    
    expected_times= [
        0.004199999999999997,
        0.29988000000015913,
        0.19221000000005145,
        0.20288000000006212,
        0.21143000000007067,
        0.2024600000000617,
        0.19868000000005792,
        0.2169700000000762,
        0.19499000000005423,
        0.23232000000009156,
        0.0076999999999998545,
        0.019769999999999364,
        0.19331000000005255
    ]

    dt = 1e-5
    max_time = 1.

    for iexp, initial_variables in enumerate(initial_variable_list):

        (plabel, zenith, azimuth, palt, lat, lng, dalt, rig, en) = initial_variables
        
        traj = Trajectory(
            plabel=plabel,
            zenith_angle = zenith,
            azimuth_angle = azimuth,
            particle_altitude = palt,
            latitude = lat,
            longitude = lng,
            detector_altitude = dalt,
            rigidity = rig,
            energy = en,
            bfield_type="igrf"
        )
        
        traj.get_trajectory(dt = dt, max_time=max_time)

        assert np.allclose(traj.final_time, expected_times[iexp])


def test_trajectories_stepsize():
    '''
    Test the final times of the trajectory evaluation in the igrf field for 
    different step sizes
    '''
    
    expected_times= [
        0.22073792288899635,
        0.22073792269455356,
        0.22073795420430928,
        0.22073818352306868,
        0.22073922285798694,
        0.22074043737707236,
        0.22082955572826515,
        0.22185737144760803,
        0.23247893970267883,
        0.30000000000000004
    ]

    dt_arr = np.logspace(-9, -1, 10)
    max_time = 1.

    for iexp, dt in enumerate(dt_arr):

        (plabel, zenith, azimuth, palt, lat, lng, dalt, rig, en) = ("p+", 90., 0., 100., 0., 0., 0., 50., None)
        
        traj = Trajectory(
            plabel=plabel,
            zenith_angle = zenith,
            azimuth_angle = azimuth,
            particle_altitude = palt,
            latitude = lat,
            longitude = lng,
            detector_altitude = dalt,
            rigidity = rig,
            energy = en,
            bfield_type="igrf"
        )
        
        traj.get_trajectory(dt = dt, max_time=max_time)

        assert np.allclose(traj.final_time, expected_times[iexp])

def test_trajectories_maxtimes():
    '''
    Test the final times of the trajectory evaluation in the igrf field for 
    different maximal times
    '''
    
    expected_times= [
        0.00999999999999976,
        0.027829999999999036,
        0.07743000000000268,
        0.2154500000000747,
        0.22074000000007998,
        0.22074000000007998,
        0.22074000000007998,
        0.22074000000007998,
        0.22074000000007998,
        0.22074000000007998
    ]

    dt = 1e-5
    max_times = np.logspace(-2, 2, 10)

    for iexp, max_time in enumerate(max_times):

        (plabel, zenith, azimuth, palt, lat, lng, dalt, rig, en) = ("p+", 90., 0., 100., 0., 0., 0., 50., None)
        
        traj = Trajectory(
            plabel=plabel,
            zenith_angle = zenith,
            azimuth_angle = azimuth,
            particle_altitude = palt,
            latitude = lat,
            longitude = lng,
            detector_altitude = dalt,
            rigidity = rig,
            energy = en,
            bfield_type="igrf"
        )
        
        traj.get_trajectory(dt = dt, max_time=max_time)

        assert np.allclose(traj.final_time, expected_times[iexp])

def test_trajectories_unvectorized():
    '''
    Test the final times of the trajectory evaluation in the igrf field for 
    the unvectorized Runge Kutta version
    '''
    
    expected_times= [
        0.004199999999999997,
        0.2999300000001592,
        0.19221000000005145,
        0.20289000000006213,
        0.21144000000007068,
        0.2024600000000617,
        0.19869000000005793,
        0.2169600000000762,
        0.19499000000005423,
        0.23231000000009155,
        0.007719999999999854,
        0.019809999999999363,
        0.19331000000005255
    ]

    dt = 1e-5
    max_time = 1.

    for iexp, initial_variables in enumerate(initial_variable_list):

        (plabel, zenith, azimuth, palt, lat, lng, dalt, rig, en) = initial_variables
        
        traj = Trajectory(
            plabel=plabel,
            zenith_angle = zenith,
            azimuth_angle = azimuth,
            particle_altitude = palt,
            latitude = lat,
            longitude = lng,
            detector_altitude = dalt,
            rigidity = rig,
            energy = en,
            bfield_type="igrf"
        )
        
        traj.get_trajectory(dt = dt, max_time=max_time, use_unvectorized=True)

        assert np.allclose(traj.final_time, expected_times[iexp])

def test_trajectories_dates():
    '''
    Test the final times of the trajectory evaluation in the igrf field for 
    different dates
    '''
    
    expected_times= [
        0.22074000000007998,
        0.22074000000007998,
        0.22074000000007998,
        0.22074000000007998,
        0.22074000000007998,
        0.22074000000007998,
        0.22074000000007998,
        0.22074000000007998,
        0.22074000000007998,
        0.22074000000007998
    ]

    dt = 1e-5
    max_time = 1.

    dates = [1900., 1905.5, 2010., 2020., 2015., 2025., 1917.25, 1918., 2000., 1990.]

    for iexp, date in enumerate(dates):

        (plabel, zenith, azimuth, palt, lat, lng, dalt, rig, en) = ("p+", 90., 0., 100., 0., 0., 0., 50., None)
        
        traj = Trajectory(
            plabel=plabel,
            zenith_angle = zenith,
            azimuth_angle = azimuth,
            particle_altitude = palt,
            latitude = lat,
            longitude = lng,
            detector_altitude = dalt,
            rigidity = rig,
            energy = en,
            bfield_type="igrf",
            date = date
        )
        
        traj.get_trajectory(dt = dt, max_time=max_time)

        assert np.allclose(traj.final_time, expected_times[iexp])

# def test_dipole_sixvec():
    
#     expected_sixvec= [
#         [ 6.37071829e+06,  1.57079633e+00,  1.94647741e-01, -2.60535458e-18,
#         1.90107218e-33,  1.58675306e-17],
#         [6.37122430e+07, 1.57079633e+00, 5.40334589e+00, 1.55532814e-17,
#         1.78893133e-33, 4.08189149e-18],
#         [6.37129272e+07, 7.97888890e-01, 2.13085193e+00, 2.67479369e-17,
#         7.52425459e-19, 1.49054531e-18],
#         [ 6.37125962e+07,  2.58546013e-01, -1.90468255e+00,  1.07138850e-17,
#         -3.55829151e-19, -6.67352818e-20],
#         [6.37136275e+07, 4.61044272e+00, 6.28478688e+10, 3.07004304e-06,
#         3.13435831e-07, 5.09382752e-11],
#         [ 6.37128663e+07,  1.61314440e+00, -2.79477499e+00,  1.06175779e-17,
#         5.72629792e-19, -1.36291384e-18],
#         [ 6.37139359e+07,  9.04914813e-01, -2.18998315e+00,  1.06076960e-17,
#         8.87330699e-19,  1.26800261e-18],
#         [ 6.37132113e+07,  1.73158257e+00, -1.61984851e+00,  1.05096179e-17,
#         1.50355324e-18,  1.48514594e-18],
#         [ 6.37138238e+07,  4.01888935e-01, -1.52667792e+00,  1.07177358e-17,
#         1.87889123e-19, -1.15052418e-19],
#         [ 6.37143905e+07,  1.55293625e+00,  4.89733120e+00,  1.03355488e-17,
#         -1.02218903e-18,  2.65517609e-18],
#         [ 6.37043199e+06,  1.74171451e+00,  3.33859884e+00, -2.55446081e-18,
#         7.01097812e-19,  4.06929779e-19]
#     ]

#     dt = 1e-5
#     max_time = 1.

#     for iexp, initial_variables in enumerate(initial_variable_list):

#         (plabel, zenith, azimuth, palt, lat, lng, dalt, rig) = initial_variables
        
#         traj = Trajectory(
#             plabel=plabel,
#             zenith_angle = zenith,
#             azimuth_angle = azimuth,
#             particle_altitude = palt,
#             latitude = lat,
#             longitude = lng,
#             detector_altitude = dalt,
#             rigidity = rig
#         )
        
#         traj.get_trajectory(dt = dt, max_time=max_time)

#         assert np.allclose(traj.final_sixvector, np.array(expected_sixvec[iexp]))