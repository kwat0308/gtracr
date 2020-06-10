import os, sys
import numpy as np
import matplotlib.pyplot as plt

sys.path.append(os.getcwd())
sys.path.append(os.path.join(os.getcwd(), "gtracr"))
# sys.path.append(os.path.join(os.getcwd(), "..", "gtracr"))

from gtracr.trajectory import ParticleTrajectory
from gtracr.utils import spherical_to_cartesian

if __name__ == "__main__":
    traj = ParticleTrajectory("p+", 1, 100)

    (t, r, theta, phi) = traj.getTrajectory(12)

    (x,y,z) = spherical_to_cartesian(r, theta, phi)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    ax.scatter(x, y, z, c=t, marker='o')

    plt.show()