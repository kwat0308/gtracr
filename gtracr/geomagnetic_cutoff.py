'''
Obtains the geomagnetic cutoff for each zenith and azimuthal angle

Structure will be much similar to test_trajectory.py
'''

import sys, os
import numpy as np
import matplotlib.pyplot as plt
import pickle

sys.path.append(os.getcwd())
sys.path.append(os.path.join(os.getcwd(), "gtracr"))

from gtracr.trajectory import Trajectory
from gtracr.add_location import location_dict
from gtracr.add_particle import particle_dict


def export_as_pkl(fpath, ds):
    with open(fpath, "wb") as f:
        pickle.dump(ds, f, protocol=-1)


def plot_heatmap(zenith, azimuth, cutoff_arrs, locname, rigidity_list,
                 particle):
    fig, ax = plt.subplots(figsize=(12, 9))
    # Z, A = np.meshgrid(zenith, azimuth, indexing="ij")
    # alpha_list = np.linspace(0., 1., num=len(cutoff_arrs))
    for i, cutoff in enumerate(cutoff_arrs):
        im = ax.pcolormesh(zenith,
                           azimuth,
                           cutoff,
                           cmap="viridis",
                           alpha=0.1,
                           vmin=np.min(rigidity_list) - 0.5,
                           vmax=np.max(rigidity_list) + 0.5)
    # cf = ax.contourf(zenith, azimuth, cutoff, levels=rigidity_list, cmap="viridis", vmin=np.min(rigidity_list)-0.5, vmax=np.max(rigidity_list)+0.5)
    cbar = fig.colorbar(im, ax=ax)
    cbar.ax.set_ylabel("Rigidity [GV]")

    ax.set_xlabel("Azimuthal Angle [Degrees]")
    ax.set_ylabel("Zenith Angle [Degrees]")
    ax.set_title("Geomagnetic Cutoffs at {0} for {1}".format(
        locname, particle))

    fig.tight_layout()

    # plt.show()
    plt.savefig(os.path.join(
        os.getcwd(), "..", "gtracr_plots",
        "{0}_{1}_cutoffplot.png".format(locname, particle)),
                dpi=800)


if __name__ == "__main__":

    # define a range of zenith and azimuthal angles
    # znum = 15
    # aznum = 30
    num = 20
    zenith_arr = np.linspace(0., 180., num, endpoint=False)
    azimuth_arr = np.linspace(-180., 180., num, endpoint=False)

    zenith_matrix, azimuth_matrix = np.meshgrid(zenith_arr, azimuth_arr, indexing="ij")

    # create particle trajectory with desired particle and energy
    # rigidity_list = [0.5, 10, 20, 50, 100, 1000]
    # energy_list = [5, 30, 50]
    rigidity_list = [5, 10, 30, 50]
    cutoff_arrs = [np.zeros((num, num))] * len(rigidity_list)
    # rigidity_list = [5, 30]

    # rigidity_list = [5]

    particle_list = [("p+", particle_dict["p+"])]

    location_list = [("Kamioka", location_dict["Kamioka"])]

    geomag_cutoffdict = {
        "Zenith": zenith_arr,
        "Azimuth": azimuth_arr,
        "Location": {}
    }

    # locations: kamioka, icecube, uofa
    # for locname, loc in list(location_dict.items()):
    for (locname, loc) in location_list:
        geomag_cutoffdict["Location"][locname] = {}
        # for pname, particle in list(particle_dict.items()):
        for (pname, particle) in particle_list:
            geomag_cutoffdict["Location"][locname][pname] = {}
            # for energy in energy_list:
            for k, rigidity in enumerate(rigidity_list):
                # get rigidity as this is more conventional for geomagnetic cutoffs
                # particle.set_rigidity_from_energy(energy)
                # rigidity = particle.rigidity
                # cutoff_arr = np.zeros((num, num))
                for i, zenith in enumerate(zenith_arr):
                    # cutoff_arr[i] = np.zeros(num)
                    for j, azimuth in enumerate(azimuth_arr):
                        print("Zenith Angle: {0}, Azimuth Angle {1}\n".format(
                            zenith, azimuth))
                        traj = Trajectory(pname,
                                          loc.latitude,
                                          loc.longitude,
                                          loc.altitude,
                                          zenith,
                                          azimuth,
                                          rigidity=rigidity)
                        traj.get_trajectory(max_step=10000, step_size=0.01)
                        cutoff_arrs[k][i][j] = traj.particle_escaped
                        print(traj.particle_escaped)

                geomag_cutoffdict["Location"][locname][pname][
                    rigidity] = cutoff_arrs[k]

            plot_heatmap(zenith_matrix, azimuth_matrix, cutoff_arrs, locname, rigidity_list, pname)

    fpath = os.path.join(os.getcwd(), "geomagnetic_cutoff.pkl")
    export_as_pkl(fpath, geomag_cutoffdict)
