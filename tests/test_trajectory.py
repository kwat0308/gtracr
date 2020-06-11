import os, sys
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go

sys.path.append(os.getcwd())
sys.path.append(os.path.join(os.getcwd(), "gtracr"))
# sys.path.append(os.path.join(os.getcwd(), "..", "gtracr"))

from gtracr.trajectory import ParticleTrajectory
from gtracr.utils import spherical_to_cartesian

if __name__ == "__main__":
    traj = ParticleTrajectory("p+")

    traj.getTrajectory(12)
    t = traj.results["t"]
    (x,y,z) = spherical_to_cartesian(traj.results["r"], traj.results["theta"], traj.results["phi"])

    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    ax.scatter(x, y, z, c=t, marker='o')

    plt.show()

    fig2 = go.Figure(data=[go.Scatter3d(x=x,y=y,z=z,mode="markers",
                    marker=dict(size=2.0, color=t, colorscale='Viridis'),
                    line=dict(color='darkblue', width=2))])

    fig2.show()