import os, sys
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go

sys.path.append(os.getcwd())
sys.path.append(os.path.join(os.getcwd(), "gtracr"))
# sys.path.append(os.path.join(os.getcwd(), "..", "gtracr"))

from gtracr.trajectory import ParticleTrajectory
from gtracr.utils import spherical_to_cartesian, EARTH_RADIUS

if __name__ == "__main__":
    # traj = ParticleTrajectory("e-", startAltitude=100, maxStep=1000)
    # traj.getTrajectory(5)
    # t = traj.results["t"]
    traj = ParticleTrajectory(
        # "p+", 12, startLongitude=89., startLatitude=-63., startAltitude=0.
          # "p+", 0.005, startLongitude=137.276599, startLatitude=36.434800, startAltitude=500)
          "p+", 30, startAltitude=0., stopAltitude=500)
    # (startPoint, endPoint) = traj.getTrajectory(zenith=70, azimuth=90)
    # print(startPoint, endPoint)
    # t = traj.results["t"]
    # (x, y, z) = spherical_to_cartesian(traj.results["r"] / EARTH_RADIUS,
    #                                    traj.results["theta"],
    #                                    traj.results["phi"])
    data1 = traj.getTrajectory(zenith=70., azimuth=270)
    data2 = traj.getTrajectory(zenith=70., azimuth=90)

    t = data1["t"]
    (x1, y1, z1) = spherical_to_cartesian(data1["r"] / EARTH_RADIUS,
                                       data1["theta"],
                                       data1["phi"])

    (x2, y2, z2) = spherical_to_cartesian(data2["r"] / EARTH_RADIUS,
                                       data2["theta"],
                                       data2["phi"])

    plt.scatter(x1, z1, c=t)
    plt.colorbar()
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