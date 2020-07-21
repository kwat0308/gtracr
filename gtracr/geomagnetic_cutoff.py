'''
Obtains the geomagnetic cutoff for each zenith and azimuthal angle

Structure will be much similar to test_trajectory.py
'''

import sys, os
import numpy as np
import matplotlib.pyplot as plt
import pickle
from scipy.interpolate import griddata
# import matplotlib.tri as tri

sys.path.append(os.getcwd())
sys.path.append(os.path.join(os.getcwd(), "gtracr"))

from gtracr.trajectory import Trajectory
from gtracr.add_location import location_dict
from gtracr.add_particle import particle_dict


def export_as_pkl(fpath, ds):
    with open(fpath, "wb") as f:
        pickle.dump(ds, f, protocol=-1)


def plot_heatmap(rigidity_cutoffdict, iter_num, locname, rigidity_list,
                 particle):
    # fig = plt.figure()
    # ax = plt.subplot(111, projection="3d")
    # ax = plt.subplot(111, projection="hammer")
    fig, ax = plt.subplots(figsize=(12, 9))

    # the azimuth and zenith arrays, obtained from the rigidity cutoff dictionary keys
    # that are in the form (azimuth, zenith) for each point
    azimuth_arr, zenith_arr = zip(*list(rigidity_cutoffdict.keys()))
    # reshape rigidity_cutoff arrays into a matrix-like structure
    # also divide by total number of iterations (normalize)
    rigidity_cutoffarr = np.array(list(rigidity_cutoffdict.values()))

    azimuth_grid = np.linspace(
        np.min(azimuth_arr) - 0.5,
        np.max(azimuth_arr) + 0.5, 500)
    zenith_grid = np.linspace(
        np.max(zenith_arr) + 0.5,
        np.min(zenith_arr) - 0.5, 500)

    # triang = tri.Triangulation(azimuth_arr, zenith_arr)
    # interpolator = tri.LinearTriInterpolator(triang, rigidity_cutoffarr)
    # azimuth_matrix, zenith_matrix = np.meshgrid(azimuth_grid, zenith_grid)
    # rigidity_cutoffgrid = interpolator(azimuth_matrix, zenith_matrix)

    rigidity_cutoffgrid = griddata(points=(azimuth_arr, zenith_arr),
                                   values=rigidity_cutoffarr,
                                   xi=(azimuth_grid[None, :],
                                       zenith_grid[:, None]),
                                   method='nearest')

    # cs = ax.scatter(azimuth_arr, zenith_arr, c=rigidity_cutoffarr, s=2.0)  #,
    # vmin=np.min(rigidity_list) - 0.5,
    # vmax=np.max(rigidity_list) + 0.5)

    # im = ax.pcolormesh(azimuth,
    #                    zenith,
    #                    rigidity_cutoffarr,
    #                    cmap="viridis",
    #                    vmin=np.min(rigidity_list) - 0.5,
    #                    vmax=np.max(rigidity_list) + 0.5)

    cs = ax.contour(
        azimuth_grid,
        zenith_grid,
        rigidity_cutoffgrid,
        #  levels=rigidity_list,
        colors="k",
        linewidths=0.5,
        levels=len(rigidity_list))

    cs = ax.contourf(
        azimuth_grid,
        zenith_grid,
        rigidity_cutoffgrid,
        #  levels=rigidity_list,
        cmap="RdBu_r")

    # #
    # ax.clabel(cs, inline=1, fontsize=10)
    cbar = fig.colorbar(cs, ax=ax)
    cbar.ax.set_ylabel("Rigidity [GV/c]")

    # ax.set_xlim([0., 360.])
    # ax.set_ylim([180., 0.])

    ax.set_xlabel("Azimuthal Angle [Degrees]")
    ax.set_ylabel("Zenith Angle [Degrees]")
    ax.set_title("Geomagnetic Rigidity Cutoffs at {0} for {1}".format(
        locname, particle))

    fig.tight_layout()

    # plt.show()
    plt.savefig(os.path.join(
        os.getcwd(), "..", "gtracr_plots",
        "{0}_{1}_cutoffplot.png".format(locname, particle)),
                dpi=800)


if __name__ == "__main__":

    # define a range of zenith and azimuthal angles
    # we flip zenith since thats how Honda's paper plots it...
    # znum = 15
    # aznum = 30
    # num = 20
    # zenith_arr = np.linspace(180., 0., num, endpoint=False)
    # azimuth_arr = np.linspace(0., 360., num, endpoint=False)

    # # create the matrix for zenith and azimuth to aid with plotting
    # zenith_matrix, azimuth_matrix = np.meshgrid(zenith_arr,
    #                                             azimuth_arr,
    #                                             indexing="ij")

    # create particle trajectory with desired particle and energy
    rigidity_list = np.arange(5, 55, 10)
    particle_list = [("p+", particle_dict["p+"])
                     ]  #, ("e-", particle_dict["e-"])]
    location_list = [("Kamioka", location_dict["Kamioka"])]

    iter_num = 5000  # total number of points used for Monte Carlo process
    # # variables used for determining index in which rigidity cutoffs should be located at
    # zenith_stepsize = 180. / num  # stepsize for zenith angle
    # azimuth_stepsize = 360. / num  # stepsize for azimuth angle

    # arrays to append zenith and azimuth angles
    zenith_arr = []
    azimuth_arr = []

    # dictionary to append rigidity cutoff values
    # this is to allow pair of points that are generated randomly
    # to be stored
    rigidity_cutoffdict = {}

    # geomag_cutoffdict = {
    #     "Zenith": zenith_arr,
    #     "Azimuth": azimuth_arr,
    #     "Location": {}
    # }

    # locations: kamioka, icecube
    # for locname, loc in list(location_dict.items()):
    for (locname, loc) in location_list:
        # geomag_cutoffdict["Location"][locname] = {}
        # for pname, particle in list(particle_dict.items()):
        for (pname, particle) in particle_list:
            # geomag_cutoffdict["Location"][locname][pname] = {}
            # for energy in energy_list:
            # contain the rigidity cutoffs
            # this will be a sparse matrix with a lot of zeros
            rigidity_cutoffarr = np.zeros((iter_num, iter_num))
            for i in range(iter_num):
                # get a random zenith and azimuth angle
                # zenith angles range from 0 to 180
                # azimuth angles range from 0 to 360
                [azimuth, zenith] = np.random.rand(2)
                azimuth *= 360.
                zenith *= 180.

                print("Zenith Angle: {0}, Azimuth Angle {1}\n".format(
                    zenith, azimuth))

                binary_list = np.zeros(len(rigidity_list))
                for k, rigidity in enumerate(rigidity_list):
                    print("Current rigidity: {:.4e}".format(rigidity))
                    traj = Trajectory(pname,
                                      latitude=loc.latitude,
                                      longitude=loc.longitude,
                                      detector_altitude=loc.altitude,
                                      zenith_angle=zenith,
                                      azimuth_angle=azimuth,
                                      particle_altitude=100.,
                                      rigidity=rigidity)
                    traj.get_trajectory(max_step=10000)
                    # append the binary result, whether particle escaped or not (in 0 or 1)
                    binary_list[k] = traj.particle_escaped
                    # cutoff_arrs[k][i][j] = traj.particle_escaped
                    # print(traj.particle_escaped)

                    # figure out starting at which rigidity the particle trajectory is invalid
                    # print(binary_list)
                nonzero_indices = np.nonzero(binary_list)[0]
                # print(nonzero_indices)
                # cutoff_index = nonzero_indices[
                #     len(nonzero_indices) -
                #     1] if len(nonzero_indices) != 0 else 0
                cutoff_index = nonzero_indices[0]
                # print(cutoff_index)

                cutoff_rigidity = rigidity_list[cutoff_index]

                print("Cutoff rigidity: {:.5f}".format(cutoff_rigidity))

                # append the values
                # # indices are determined by index = value / stepsize
                # # this should work since the maximum possible size of the
                # # matrix is num x num, and (for example) the max number of zenith
                # # angles that can be appended is
                # zenith_index = zenith / zenith_stepsize
                # azimuth_index = azimuth / azimuth_stepsize
                # rigidity_cutoffarr[zenith_index][
                #     azimuth_index] = cutoff_rigidity

                # store the rgiditiy cutoff inside a dictionary
                if (azimuth, zenith) in list(rigidity_cutoffdict.keys()):
                    rigidity_cutoffdict[(azimuth, zenith)] += cutoff_rigidity
                else:
                    rigidity_cutoffdict[(azimuth, zenith)] = cutoff_rigidity

                # geomag_cutoffdict["Location"][locname][pname][
                # rigidity] = cutoff_arrs[k]

            # print(rigidity_cutoffarr)

            plot_heatmap(rigidity_cutoffdict, iter_num, locname, rigidity_list,
                         pname)

    # fpath = os.path.join(os.getcwd(), "geomagnetic_cutoff.pkl")
    # export_as_pkl(fpath, geomag_cutoffdict)

# def plot_heatmap(zenith, azimuth, rigidity_cutoffs, locname, rigidity_list,
#                  particle):
#     # fig = plt.figure()
#     # ax = plt.subplot(111, projection="hammer")
#     fig, ax = plt.subplots(figsize=(12, 9))

#     cs = ax.contour(
#         azimuth,
#         zenith,
#         rigidity_cutoffs,
#         #  levels=rigidity_list,
#         cmap="viridis",
#         vmin=np.min(rigidity_list) - 0.5,
#         vmax=np.max(rigidity_list) + 0.5)
#     # im = ax.pcolormesh(azimuth,
#     #                    zenith,
#     #                    rigidity_cutoffarr,
#     #                    cmap="viridis",
#     #                    vmin=np.min(rigidity_list) - 0.5,
#     #                    vmax=np.max(rigidity_list) + 0.5)
#     ax.clabel(cs, inline=1, fontsize=10)
#     cbar = fig.colorbar(cs, ax=ax)
#     cbar.ax.set_ylabel("Rigidity [GV]")

#     # ax.set_xlim([0., 360.])
#     # ax.set_ylim([180., 0.])

#     ax.set_xlabel("Azimuthal Angle [Degrees]")
#     ax.set_ylabel("Zenith Angle [Degrees]")
#     ax.set_title("Geomagnetic Rigidity Cutoffs at {0} for {1}".format(
#         locname, particle))

#     fig.tight_layout()

#     # plt.show()
#     plt.savefig(os.path.join(
#         os.getcwd(), "..", "gtracr_plots",
#         "{0}_{1}_cutoffplot.png".format(locname, particle)),
#                 dpi=800)

# if __name__ == "__main__":

#     # define a range of zenith and azimuthal angles
#     # we flip zenith since thats how Honda's paper plots it...
#     # znum = 15
#     # aznum = 30
#     num = 20
#     zenith_arr = np.linspace(180., 0., num, endpoint=False)
#     azimuth_arr = np.linspace(0., 360., num, endpoint=False)

#     # create the matrix for zenith and azimuth to aid with plotting
#     zenith_matrix, azimuth_matrix = np.meshgrid(zenith_arr,
#                                                 azimuth_arr,
#                                                 indexing="ij")

#     # create particle trajectory with desired particle and energy
#     rigidity_list = np.arange(5, 55, 5)
#     particle_list = [("p+", particle_dict["p+"])
#                      ]  #, ("e-", particle_dict["e-"])]
#     location_list = [("Kamioka", location_dict["Kamioka"])]

#     # geomag_cutoffdict = {
#     #     "Zenith": zenith_arr,
#     #     "Azimuth": azimuth_arr,
#     #     "Location": {}
#     # }

#     # locations: kamioka, icecube
#     # for locname, loc in list(location_dict.items()):
#     for (locname, loc) in location_list:
#         # geomag_cutoffdict["Location"][locname] = {}
#         # for pname, particle in list(particle_dict.items()):
#         for (pname, particle) in particle_list:
#             # geomag_cutoffdict["Location"][locname][pname] = {}
#             # for energy in energy_list:
#             rigidity_cutoffarr = np.zeros((num, num))
#             for i, zenith in enumerate(zenith_arr):
#                 # cutoff_arr[i] = np.zeros(num)
#                 for j, azimuth in enumerate(azimuth_arr):
#                     print("Zenith Angle: {0}, Azimuth Angle {1}\n".format(
#                         zenith, azimuth))
#                     # create a new array to append binary values
#                     binary_list = np.zeros(len(rigidity_list))
#                     for k, rigidity in enumerate(rigidity_list):
#                         print("Current rigidity: {:.4e}".format(rigidity))
#                         traj = Trajectory(pname,
#                                           latitude=loc.latitude,
#                                           longitude=loc.longitude,
#                                           detector_altitude=loc.altitude,
#                                           zenith_angle=zenith,
#                                           azimuth_angle=azimuth,
#                                           particle_altitude=100.,
#                                           rigidity=rigidity)
#                         traj.get_trajectory(max_step=10000)
#                         # append the binary result, whether particle escaped or not (in 0 or 1)
#                         binary_list[k] = traj.particle_escaped
#                         # cutoff_arrs[k][i][j] = traj.particle_escaped
#                         # print(traj.particle_escaped)

#                     # figure out starting at which rigidity the particle trajectory is invalid
#                     # print(binary_list)
#                     nonzero_indices = np.nonzero(binary_list)[0]
#                     # print(nonzero_indices)
#                     # cutoff_index = nonzero_indices[
#                     #     len(nonzero_indices) -
#                     #     1] if len(nonzero_indices) != 0 else 0
#                     cutoff_index = nonzero_indices[0]
#                     # print(cutoff_index)

#                     cutoff_rigidity = rigidity_list[cutoff_index]

#                     print("Cutoff rigidity: {:.5f}".format(cutoff_rigidity))

#                     rigidity_cutoffarr[i][j] = cutoff_rigidity

#                 # geomag_cutoffdict["Location"][locname][pname][
#                 # rigidity] = cutoff_arrs[k]

#             # print(rigidity_cutoffarr)

#             plot_heatmap(zenith_matrix, azimuth_matrix, rigidity_cutoffarr,
#                          locname, rigidity_list, pname)

#     # fpath = os.path.join(os.getcwd(), "geomagnetic_cutoff.pkl")
#     # export_as_pkl(fpath, geomag_cutoffdict)