'''
Compare the values of the final times and sixvector of the trajectory for the dipole model
'''

import os
import sys
import numpy as np
import pytest

from gtracr.trajectory import Trajectory
from gtracr.lib.trajectorypoint import TrajectoryPoint

def test_diptraj_times():

    expected_times= [
        0.004199999999999997,
        0.19221000000005145
    ]

    # in the form (plabel, zenith, azimuth, particle_altitude, latitude, longitude, detector_altitude, rigidity)
    initial_variable_list = [
        ("p+", 90., 90., 100., 0., 0., 0., 30.),
        ("p+", 0., 25., 100., 50., 100., 0., 50.)
    ]

    dt = 1e-5
    max_time = 1.

    for iexp, initial_variables in enumerate(initial_variable_list):
        plabel= initial_variables[0]
        zenith = initial_variables[1]
        azimuth = initial_variables[2]
        particle_altitude = initial_variables[3]
        latitude = initial_variables[4]
        longitude = initial_variables[5]
        detector_altitude = initial_variables[6]
        rigidity = initial_variables[7]
        
        traj = Trajectory(
            plabel=plabel,
            zenith_angle = zenith,
            azimuth_angle = azimuth,
            particle_altitude = particle_altitude,
            latitude = latitude,
            longitude = longitude,
            detector_altitude = detector_altitude,
            rigidity = rigidity
        )
        
        traj.get_trajectory(dt = dt, max_time=max_time)

        assert np.allclose(traj.final_time, expected_times[iexp])


