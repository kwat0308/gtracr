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
    traj1 = Trajectory("p+",
                       0.,
                       0.,
                       0.,
                       0.,
                       0.,
                       rigidity=10.,
                       escape_altitude=10 * EARTH_RADIUS)
    traj1.get_trajectory(max_step=10000, step_size=1e-7)

    result = traj1.get_plotting_variables()

    t = result["t"]
    x = result["x"] / EARTH_RADIUS
    y = result["y"] / EARTH_RADIUS
    z = result["z"] / EARTH_RADIUS

    # traj2 = Trajectory("p+",
    #                    0.,
    #                    0.,
    #                    20.,
    #                    70.,
    #                    90.,
    #                    rigidity=30.,
    #                    escape_altitude=565)
    # traj2.get_trajectory(max_step=10000, step_size=1e-4)

    # result2 = traj2.get_plotting_variables()

    # x2 = result2["x"] / EARTH_RADIUS
    # y2 = result2["y"] / EARTH_RADIUS
    # z2 = result2["z"] / EARTH_RADIUS

    # fig1, ax1 = plt.subplots()
    # ax1.scatter(x, y, c=t)
    # ax.set_zlabel(r"z [$R_E$]")
    # ax1.scatter(x2, y2, c="r")
    # plt.show()
    fig2, ax2 = plt.subplots()
    ax2.scatter(x, z, c=t)
    # ax2.scatter(x2, z2, c="r")
    ax2.set_xlim([-3, 3])
    ax2.set_ylim([-3, 3])
    # ax.set_zlim([-1, 1])
    ax2.set_xlabel(r"x [$R_E$]")
    ax2.set_ylabel(r"z [$R_E$]")
    plt.show()
    # fig3, ax3 = plt.subplots()
    # ax3.scatter(y, z, c=t)
    # ax3.scatter(y2, z2, c="r")
    # plt.show()

    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    cm = ax.scatter(x, y, z, c=t, marker='o')
    # ax.scatter(x2, y2, z2, c="r", marker='o')
    cbar = fig.colorbar(cm, ax=ax)
    # plt.colorbar()
    # plt.savefig("test.png")
    ax.set_xlim([-1, 1])
    ax.set_ylim([-1, 1])
    ax.set_zlim([-1, 1])
    ax.set_xlabel(r"x [$R_E$]")
    ax.set_ylabel(r"y [$R_E$]")
    ax.set_zlabel(r"z [$R_E$]")
    cbar.ax.set_ylabel("Time [s]")
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