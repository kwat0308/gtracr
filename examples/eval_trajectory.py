'''
A command-line interface to obtain the trajectory plots,
both the projections and the 3-D plot.

This should support some html format with interactive window
like PlotLy in the future. 
'''
import os
import sys
import numpy as np
import argparse

from gtracr.utils import dec_to_dms
from gtracr.lib.constants import EARTH_RADIUS, KG_M_S_PER_GEVC
from gtracr.plotting import plot_3dtraj, plot_2dtraj, plot_trajmomentum
from gtracr.trajectory import Trajectory

CURRENT_DIR = os.path.dirname(os.path.realpath(__file__))
PARENT_DIR = os.path.dirname(CURRENT_DIR)

PLOT_DIR = os.path.join(PARENT_DIR, "..", "gtracr_plots")


def convert_to_cartesian(trajectory_data):

    r_arr = trajectory_data["r"] / EARTH_RADIUS
    theta_arr = trajectory_data["theta"]
    phi_arr = trajectory_data["phi"]

    # convert to cartesian & add to dict
    trajectory_data["x"] = r_arr * np.sin(theta_arr) * np.cos(phi_arr)
    trajectory_data["y"] = r_arr * np.sin(theta_arr) * np.sin(phi_arr)
    trajectory_data["z"] = r_arr * np.cos(theta_arr)


def plot_trajectory(traj_datadict, title, check_3dtraj=False, show_plot=False):
    # # convert to cartesian coordinates
    convert_to_cartesian(traj_datadict)

    # plot the projections
    plot_2dtraj([traj_datadict], plotdir_path=PLOT_DIR)

    # plot the 3-d trajectory with wireframe sphere as the earth
    if check_3dtraj:
        plot_3dtraj([traj_datadict], title_name=title, plotdir_path=PLOT_DIR)


def get_trajectory():
    # set parameters

    # parameters for trajectory
    # particle is assumed to be proton
    q = 1
    plabel = "p+"

    # initial momentum
    p0 = 30.
    rigidity = p0 / np.abs(q)  # convert to rigidity

    # location of detector
    lat = 10.
    lng = 40.
    detector_alt = 0.

    # 3-vector of particle
    zenith = 90.
    azimuth = 0.
    particle_alt = 100.

    # set integration parameters
    dt = 1e-5
    max_time = 1.
    max_step = 10000

    # control variables for the code
    check_pmag = True  # if we want to check the momentum magnitude
    check_3dtraj = True  # if we want to check the 3d trajectory or not
    show_plot = False  # if we want to show the plot on some GUI or not

    # first create plot directory if it doesnt exist
    if not os.path.exists(PLOT_DIR):
        os.mkdir(PLOT_DIR)

    # initialize trajectory
    traj = Trajectory(plabel=plabel,
                      zenith_angle=zenith,
                      azimuth_angle=azimuth,
                      particle_altitude=particle_alt,
                      latitude=lat,
                      longitude=lng,
                      detector_altitude=detector_alt,
                      rigidity=rigidity,
                      bfield_type="igrf")

    # obtain the trajectory result
    traj_datadict = traj.get_trajectory(dt=dt,
                                        max_time=max_time,
                                        get_data=True,
                                        max_step=max_step,
                                        use_python=False,
                                        use_unvectorized=False)

    # convert lat, long in decimal notation to dms
    lat_dms, lng_dms = dec_to_dms(lat, lng)

    title = "Particle Trajectory at {:s}, {:s} with Zenith Angle {:.1f}°, \
           \n Azimuth Angle {:.1f}° and Rigidity R = {:.1f}GV".format(
        lat_dms, lng_dms, zenith, azimuth, rigidity)

    # get momentum only if check_pmag is true
    if check_pmag:
        plot_trajmomentum(traj_datadict, p0, show_plot)

    # plot the trajectory
    plot_trajectory(traj_datadict, title, check_3dtraj, show_plot)


if __name__ == "__main__":
    # should add some argparse thing later on
    get_trajectory()
