import sys, os
# import os
import numpy as np
import matplotlib.pyplot as plt
# import pickle
# from scipy.interpolate import griddata
# import matplotlib.tri as tri

# sys.path.append(os.getcwd())
sys.path.append(os.path.join(os.getcwd(), ".."))
sys.path.append(os.path.join(os.getcwd(), "..", "gtracr"))

from trajectory import Trajectory
from add_location import location_dict
# from add_particle import particle_dict


def get_gmcutoffs(loc, rigidity_list, iter_num):
    result_arr1 = []
    for i in range(iter_num):
        # get a random zenith and azimuth angle
        # zenith angles range from 0 to 180
        # azimuth angles range from 0 to 360
        [azimuth, zenith] = np.random.rand(2)
        azimuth *= 360.
        zenith *= 180.

        #         print("Zenith Angle: {0}, Azimuth Angle {1}\n".format(
        #             zenith, azimuth))
        for k, rigidity in enumerate(rigidity_list):
            #             print("Current rigidity: {:.4e}".format(rigidity))
            traj = Trajectory("p+",
                              latitude=loc.latitude,
                              longitude=loc.longitude,
                              detector_altitude=loc.altitude,
                              zenith_angle=zenith,
                              azimuth_angle=azimuth,
                              particle_altitude=100.,
                              rigidity=rigidity)
            traj.get_trajectory(max_step=10000)

            if traj.particle_escaped == True:
                result_arr1.append(
                    (azimuth, zenith, rigidity))  #, rigidity_analytical))
                break

    return result_arr1


if __name__ == "__main__":
    # create particle trajectory with desired particle and energy
    rigidity_list = np.arange(5, 55, 5)

    loc = location_dict["Kamioka"]
    iter_num = 1000
    result_arr = get_gmcutoffs(loc, rigidity_list, iter_num)