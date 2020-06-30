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


def plot_heatmap(zenith, azimuth, cutoff, locname, energy, particle):
    fig, ax = plt.subplots(figsize=(12, 9))
    Z, A = np.meshgrid(zenith, azimuth, indexing="ij")
    im = ax.pcolormesh(A, Z, cutoff, cmap="binary", vmin=0, vmax=1)

    cbar = fig.colorbar(im, ax=ax)
    cbar.ax.set_ylabel("Trajectory Validity")

    ax.set_xlabel("Azimuthal Angle [Degrees]")
    ax.set_ylabel("Zenith Angle [Degrees]")
    ax.set_title("Geomagnetic Cutoffs at {0} for {1} at {2}GV Rigidity".format(
        locname, particle, energy))

    fig.tight_layout()

    # plt.show()
    plt.savefig(os.path.join(
        os.getcwd(), "..", "gtracr_plots",
        "{0}_{1}_{2}_cutoffplot.png".format(locname, particle, energy)),
                dpi=800)

# def get_cutoffs(zenith, azimuth):
#     data = []
#     for i, azimuth in enumerate(azimuth_arr):
#         data.append(np.zeros(num))
#         for j, zenith in enumerate(zenith_arr):
#             print(zenith, azimuth)

#             cutoff = 1  # binary boolean value
#             (startPoint,
#                 endPoint) = traj.getTrajectory(zenith, azimuth)

#             print(endPoint)

#             # if particle touches Earth's surface again 
#             # if startPoint.altitude == endPoint.altitude
#             if endPoint.altitude < 0:
#                 cutoff = 0
#                 # geomag_cutoffdict["Location"][locname][energy]["Cutoff"].append(True)
#             # else:
#             #     cutoff = False
#             # geomag_cutoffdict["Location"][locname][energy]["Cutoff"].append(False)
#             # geomag_cutoffdict["Location"][locname][energy][particle]["Cutoff"].append(cutoff)
#             data[i][j] = cutoff
    
#     return data

                


if __name__ == "__main__":

    # define a range of zenith and azimuthal angles
    # znum = 15
    # aznum = 30
    num = 20
    zenith_arr = np.linspace(0., 180., num, endpoint=False)
    azimuth_arr = np.linspace(-180., 180., num, endpoint=False)

    # create particle trajectory with desired particle and energy
    # rigidity_list = [0.5, 10, 20, 50, 100, 1000]
    # energy_list = [5, 30, 50]
    # rigidity_list = [5, 10, 30, 50]
    # rigidity_list = [5, 30]

    rigidity_list = [5]

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
            for rigidity in rigidity_list:
                # get rigidity as this is more conventional for geomagnetic cutoffs
                # particle.set_rigidity_from_energy(energy)
                # rigidity = particle.rigidity
                cutoff_arr = []
                for i, azimuth in enumerate(azimuth_arr):
                    cutoff_arr.append(np.zeros(num))
                    for j, zenith in enumerate(zenith_arr):
                        print("Zenith Angle: {0}, Azimuth Angle {1}\n".format(zenith, azimuth))
                # geomag_cutoffdict["Location"][locname][pname][rigidity] = {
                #     "Cutoff": np.zeros(num)
                # }
                # traj = ParticleTrajectory(pname, energy, loc.latitude,
                #                           loc.longitude, loc.altitude)

                # cutoff_arr = get_cutoffs(zenith_arr, azimuth_arr)
                        traj = Trajectory(pname, loc.latitude, loc.longitude, loc.altitude, zenith, azimuth, rigidity=rigidity)
                        traj.getTrajectory()
                        cutoff_arr[i][j] = traj.particleEscaped
                        print(traj.particleEscaped)

                plot_heatmap(zenith_arr, azimuth_arr, cutoff_arr, locname, rigidity,
                            pname)

                geomag_cutoffdict["Location"][locname][pname][rigidity] = cutoff_arr

    fpath = os.path.join(os.getcwd(), "geomagnetic_cutoff.pkl")
    export_as_pkl(fpath, geomag_cutoffdict)
