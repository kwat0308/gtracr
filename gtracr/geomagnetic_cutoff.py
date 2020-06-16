'''
Obtains the geomagnetic cutoff for each zenith and azimuthal angle

Structure will be much similar to test_trajectory.py
'''

import sys, os
import numpy as np

sys.path.append(os.getcwd())
sys.path.append(os.path.join(os.getcwd(), "gtracr"))

from gtracr.trajectory import ParticleTrajectory
from gtracr.add_location import location_dict


if __name__ == "__main__":

    # define a range of zenith and azimuthal angles
    num = 100
    zenith_arr = np.linspace(0., 180., num)
    azimuth_arr = np.linspace(0., 360., num, endpoint=False)

    # create particle trajectory with desired particle and energy
    # energy_list = [0.5, 10, 20, 50, 100, 1000]
    energy_list = [10]

    # locations: kamioka, icecube, uofa
    for locname, loc in list(location_dict.items()):
        for energy in energy_list:
            traj = ParticleTrajectory("p+", energy, loc.latitude, loc.longitude, loc.altitude)

            # now iterate over each zenith and azimuthal angle
            for zenith in zenith_arr:
                for azimuth in azimuth_arr:
                    traj.getTrajectory(zenith, azimuth)

