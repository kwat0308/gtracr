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

from gtracr.trajectory import ParticleTrajectory
from gtracr.add_location import location_dict


def export_as_pkl(fpath, ds):
    with open(fpath, "rb") as f:
        pickle.dump(ds, f, protocol=-1)


def plot_heatmap(zenith, azimuth, cutoff, locname, energy, particle):
    fig, ax = plt.subplots(figsize=(8, 10))
    Z, A = np.meshgrid(zenith, azimuth, indexing="ij")
    im = ax.pcolormesh(A, Z, cutoff, cmap="viridis", vmin=0, vmax=1)

    cbar = fig.colorbar(im, ax=ax)
    cbar.ax.set_ylabel("Trajectory Validity")

    ax.set_xlabel("Azimuthal Angle [degree]")
    ax.set_ylabel("Zenith Angle [degree]")
    ax.set_title("Geomagnetic Cutoffs at {0} for {1} at {2}GeV".format(
        locname, particle, energy))

    fig.tight_layout()

    # plt.show()
    plt.savefig(os.path.join(
        os.getcwd(), "plots",
        "{0}_{1}_{2}_cutoffplot.png".format(locname, particle, energy)),
                dpi=800)



if __name__ == "__main__":

    # define a range of zenith and azimuthal angles
    num = 10
    zenith_arr = np.linspace(0., 180., num)
    azimuth_arr = np.linspace(0., 360., num, endpoint=False)

    # create particle trajectory with desired particle and energy
    # energy_list = [0.5, 10, 20, 50, 100, 1000]
    energy_list = [30, 50]

    particle_list = ["p+"]

    geomag_cutoffdict = {
        "Zenith": zenith_arr,
        "Azimuth": azimuth_arr,
        "Location": {}
    }

    # locations: kamioka, icecube, uofa
    for locname, loc in list(location_dict.items()):
        geomag_cutoffdict["Location"][locname] = {}
        for energy in energy_list:
            geomag_cutoffdict["Location"][locname][energy] = {}
            for particle in particle_list:
                geomag_cutoffdict["Location"][locname][energy][particle] = {
                    "Cutoff": np.zeros(num)
                }
                traj = ParticleTrajectory(particle, energy, loc.latitude,
                                          loc.longitude, loc.altitude)
                # print(traj.startLatitude, traj.startLongitude)
                # now iterate over each zenith and azimuthal angle
                data = []
                for i, azimuth in enumerate(azimuth_arr):
                    data.append(np.zeros(num))
                    for j, zenith in enumerate(zenith_arr):
                        print(zenith, azimuth)

                        cutoff = 1  # binary boolean value
                        (startPoint,
                         endPoint) = traj.getTrajectory(zenith, azimuth)

                        print(endPoint)

                        # if particle touches Earth's surface again 
                        # if startPoint.altitude == endPoint.altitude
                        if endPoint.altitude < 0:
                            cutoff = 0
                            # geomag_cutoffdict["Location"][locname][energy]["Cutoff"].append(True)
                        # else:
                        #     cutoff = False
                        # geomag_cutoffdict["Location"][locname][energy]["Cutoff"].append(False)
                        # geomag_cutoffdict["Location"][locname][energy][particle]["Cutoff"].append(cutoff)
                        data[i][j] = cutoff

                plot_heatmap(zenith_arr, azimuth_arr, data, locname, energy,
                             particle)

                geomag_cutoffdict["Location"][locname][energy][particle] = data

    fpath = os.getcwd()
    export_as_pkl(fpath, geomag_cutoffdict)
