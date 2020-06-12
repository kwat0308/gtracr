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
    # traj = ParticleTrajectory("p+", startAltitude=500, stopAltitude=2, maxStep=1000)
    traj = ParticleTrajectory("p-", startAltitude=500, maxStep=1000)
    traj.getTrajectory(12)#, startLongitude=81, startLatitude=81)#, startZenith=50)#, startAzimuth=10))
    # traj.getTrajectory(1000, startLongitude=81, startLatitude=81, startZenith=50, startAzimuth=10)
    t = traj.results["t"]
    (x,y,z) = spherical_to_cartesian(traj.results["r"] / EARTH_RADIUS, traj.results["theta"], traj.results["phi"])

    plt.scatter(x, y, c=t)
    plt.colorbar()
    plt.show()

    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    cm = ax.scatter(x, y, z, c=t, marker='o')
    fig.colorbar(cm, ax=ax)
    # plt.colorbar()

    plt.show()

    # fig2 = go.Figure(data=[go.Scatter3d(x=x,y=y,z=z,mode="markers",
    #                 marker=dict(size=2.0, color=t, colorscale='Viridis'),
    #                 line=dict(color='darkblue', width=2))])

    # fig2.show()