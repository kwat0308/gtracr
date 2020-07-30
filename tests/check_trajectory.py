import os, sys
import numpy as np
import matplotlib.pyplot as plt

sys.path.append(os.getcwd())
sys.path.append(os.path.join(os.getcwd(), "gtracr"))
# sys.path.append(os.path.join(os.getcwd(), "..", "gtracr"))

from gtracr.trajectory import Trajectory
# from gtracr.utils import spherical_to_cartesian
from gtracr.constants import EARTH_RADIUS

if __name__ == "__main__":
    traj1 = Trajectory(
        "p+",
        latitude=0.,
        longitude=0.,
        detector_altitude=0.,
        zenith_angle=90.,
        azimuth_angle=90.,
        particle_altitude=100.,
        rigidity=30.,
    )
    # for i in range(10):
    result = traj1.get_trajectory(max_step=10000,
                                  step_size=1e-5,
                                  get_data=True)

    # result = traj1.get_plotting_variables()

    t_arr = result["t"]
    x_arr = result["x"] / EARTH_RADIUS
    y_arr = result["y"] / EARTH_RADIUS
    z_arr = result["z"] / EARTH_RADIUS

    fig_3d = plt.figure()
    ax_3d = fig_3d.add_subplot(111, projection="3d")
    cm_3d = ax_3d.scatter(x_arr, y_arr, z_arr, c=t_arr, marker='o')
    # ax_3d.scatter(x2, y2, z2, c="r", marker='o')
    cbar_3d = fig_3d.colorbar(cm_3d, ax=ax_3d)
    # plt.colorbar()
    # plt.savefig("test.png")
    ax_3d.set_xlim([-3, 3])
    ax_3d.set_ylim([-3, 3])
    ax_3d.set_zlim([-3, 3])
    ax_3d.set_xlabel(r"x [$R_E$]")
    ax_3d.set_ylabel(r"y [$R_E$]")
    ax_3d.set_zlabel(r"z [$R_E$]")
    cbar_3d.ax.set_ylabel("Time [s]")
    plt.show()

    fig_xy, ax_xy = plt.subplots()
    cm_xy = ax_xy.scatter(x_arr, y_arr, c=t_arr)
    cbar_xy = fig_xy.colorbar(cm_xy, ax=ax_xy)
    ax_xy.set_xlim([-3, 3])
    ax_xy.set_ylim([-3, 3])
    ax_xy.set_xlabel(r"x [$R_E$]")
    ax_xy.set_ylabel(r"y [$R_E$]")
    cbar_xy.ax.set_ylabel("Time [s]")
    plt.show()

    fig_xz, ax_xz = plt.subplots()
    cm_xz = ax_xz.scatter(x_arr, z_arr, c=t_arr)
    cbar_xz = fig_xz.colorbar(cm_xz, ax=ax_xz)
    ax_xz.set_xlim([-3, 3])
    ax_xz.set_ylim([-3, 3])
    ax_xz.set_xlabel(r"x [$R_E$]")
    ax_xz.set_ylabel(r"z [$R_E$]")
    cbar_xz.ax.set_ylabel("Time [s]")
    plt.show()

    fig_yz, ax_yz = plt.subplots()
    cm_yz = ax_yz.scatter(y_arr, x_arr, c=t_arr)
    cbar_yz = fig_yz.colorbar(cm_yz, ax=ax_yz)
    ax_yz.set_xlim([-3, 3])
    ax_yz.set_ylim([-3, 3])
    ax_yz.set_xlabel(r"y [$R_E$]")
    ax_yz.set_ylabel(r"z [$R_E$]")
    cbar_yz.ax.set_ylabel("Time [s]")
    plt.show()
'''
if __name__ == "__main__":
    
    # traj = ParticleTrajectory("e-", startAltitude=100, maxStep=1000)
    # traj.getTrajectory(5)
    # t = traj.results["t"]
    # traj = ParticleTrajectory(
    #     # "p+", 12, startLongitude=89., startLatitude=-63., startAltitude=0.
    #       # "p+", 0.005, startLongitude=137.276599, startLatitude=36.434800, startAltitude=500)
    #       "p+", 30, startAltitude=0., stopAltitude=500)
    # # (startPoint, endPoint) = traj.getTrajectory(zenith=70, azimuth=90)
    # print(startPoint, endPoint)
    # t = traj.results["t"]
    # (x, y, z) = spherical_to_cartesian(traj.results["r"] / EARTH_RADIUS,
    #                                    traj.results["theta"],
    #                                    traj.results["phi"])
    traj1 = ParticleTrajectory("p+", 30, startAltitude=1, stopAltitude=500)
    traj2 = ParticleTrajectory("p+", 30, startAltitude=1, stopAltitude=500)
    data1 = traj1.getTrajectory(zenith=70., azimuth=-90)
    data2 = traj2.getTrajectory(zenith=70., azimuth=90)

    t1 = data1["t"]
    t2 = data2["t"]
    (x1, y1, z1) = spherical_to_cartesian(data1["r"] / EARTH_RADIUS,
                                       data1["theta"],
                                       data1["phi"])

    (x2, y2, z2) = spherical_to_cartesian(data2["r"] / EARTH_RADIUS,
                                       data2["theta"],
                                       data2["phi"])

    fig, ax = plt.subplots()
    ax.scatter(x1, y1, c="k")
    ax.scatter(x2, y2, c="r")
    # plt.colorbar()
    plt.show()

    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    cm = ax.scatter(x1, y1, z1, c="k", marker='o')
    cm = ax.scatter(x2, y2, z2, c="r", marker='o')
    fig.colorbar(cm, ax=ax)
    # plt.colorbar()
    # plt.savefig("test.png")
    plt.show()

    # fig2 = go.Figure(data=[go.Scatter3d(x=x,y=y,z=z,mode="markers",
    #                 marker=dict(size=2.0, color=t, colorscale='Viridis'),
    #                 line=dict(color='darkblue', width=2))])

    # fig2.show()


    ptraj = ParticleTrajectory("p-", 12, startLatitude=0., startLongitude=-60, startAltitude=20, stopAltitude=500)
    results = ptraj.getTrajectory(zenith=0., azimuth=0.)

    t = results["t"]
    (x,y,z) = spherical_to_cartesian(results["r"] / EARTH_RADIUS, results["theta"], results["phi"])

    plt.scatter(x,y,c=t)
    # plt.xlim([-1.5, 1.5])
    # plt.ylim([-1.5, 1.5])
    plt.show()
    # fig = plt.figure()
    # ax = fig.add_subplot(111, projection="3d")
    # cm = ax.scatter(x, y, z, c=t, marker='o')
    # # cm = ax.scatter(x2, y2, z2, c=t, marker='o')
    # fig.colorbar(cm, ax=ax)
    # # plt.colorbar()
    # # plt.savefig("test.png")
    # plt.show()    
'''