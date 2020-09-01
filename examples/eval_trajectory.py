from gtracr.utils import dec_to_dms
from gtracr.lib.constants import EARTH_RADIUS, KG_M_S_PER_GEVC
from gtracr.trajectory import Trajectory
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
import matplotlib.pyplot as plt
# import mpld3
import matplotlib.patches as patches
from mpl_toolkits import mplot3d

CURRENT_DIR = os.path.dirname(os.path.realpath(__file__))
PARENT_DIR = os.path.dirname(CURRENT_DIR)
# sys.path.append(PARENT_DIR)

PLOT_DIR = os.path.join(PARENT_DIR, "..", "gtracr_plots")


def convert_to_cartesian(r_arr, theta_arr, phi_arr):

    # convert to cartesian
    x_arr = r_arr * np.sin(theta_arr) * np.cos(phi_arr)
    y_arr = r_arr * np.sin(theta_arr) * np.sin(phi_arr)
    z_arr = r_arr * np.cos(theta_arr)

    return x_arr, y_arr, z_arr


def plot_momentum(traj_datadict, p0, show_plot=False):
    # check momentum magnitude vs steps since |p| should be
    # constant throughout the trajectory
    p_arr = np.sqrt(traj_datadict["pr"]**2. + traj_datadict["ptheta"]**2. +
                    traj_datadict["pphi"]**2.) / KG_M_S_PER_GEVC
    p_ratio = p_arr / p0

    # figure for momentum vs steps
    # for physics checking
    fig_pmag, ax_pmag = plt.subplots(figsize=(12, 9))

    # momentum ratio vs steps
    ax_pmag.plot(p_ratio, color="b")
    ax_pmag.set_ylim([0.5, 1.5])
    # ax_pmag.set_ylim([-3, 3])
    ax_pmag.set_xlabel(r"Steps")
    ax_pmag.set_ylabel(r"$p/p_0$")
    ax_pmag.set_title(r"Momentum ratio $p/p_0$ per step")

    if show_plot:
        plt.show()
    plt.savefig(os.path.join(PLOT_DIR, "pmag_plot.png"))


def plot_projection(arr1, arr2, t_arr, fig, ax, label1, label2):
    '''
    Plot the projection of the trajectory.

    Parameters
    ----------
    arr1 (numpy array): 
        the array that goes on the x-axis
    arr2 (numpy array): 
        the array that goes on the y-axis
    t_arr (numpy array):
        the time array
    fig, ax:
        matplotlib plt.Figure and plt.Axes.axes objects
    label1, label2 (str):
        the labels that indicate the identity of arr1 and arr2 resp.
    '''
    cm = ax.scatter(arr1, arr2, c=t_arr)
    circ = patches.Circle((0., 0.),
                          1.,
                          alpha=0.8,
                          fc='None',
                          linestyle='-',
                          ec="k",
                          lw=2.0)
    ax.add_patch(circ)
    cbar = fig.colorbar(cm, ax=ax)
    ax.set_xlim([-3, 3])
    ax.set_ylim([-3, 3])
    ax.set_xlabel(r"{:s} [$R_E$]".format(label1))
    ax.set_ylabel(r"{:s} [$R_E$]".format(label2))
    cbar.ax.set_ylabel("Time [s]")
    ax.set_title("Trajectory projected onto {:s}-{:s} plane".format(
        label1, label2))


def plot_trajectory(traj_datadict, title, check_3dtraj=False, show_plot=False):
    # obtain data from dictionary
    t_arr = traj_datadict["t"]
    r_arr = traj_datadict["r"] / EARTH_RADIUS
    theta_arr = traj_datadict["theta"]
    phi_arr = traj_datadict["phi"]

    # perform necessary conversions
    # convert to cartesian coordinates
    x_arr, y_arr, z_arr = convert_to_cartesian(r_arr, theta_arr, phi_arr)
    mag_arr = np.sqrt(x_arr**2. + y_arr**2. + z_arr**2.)

    # figures for the projections
    fig_proj, ax_proj = plt.subplots(ncols=2,
                                     nrows=2,
                                     figsize=(16, 12),
                                     constrained_layout=True)

    # relabel them for each projection type

    ax_xy = ax_proj[0, 0]
    ax_xz = ax_proj[0, 1]
    ax_yz = ax_proj[1, 0]
    ax_mag = ax_proj[1, 1]

    # projection onto xy plane
    plot_projection(x_arr, y_arr, t_arr, fig_proj, ax_xy, "x", "y")
    # projection onto xz plane
    plot_projection(x_arr, z_arr, t_arr, fig_proj, ax_xz, "x", "z")
    # projection onto yz plane
    plot_projection(y_arr, z_arr, t_arr, fig_proj, ax_yz, "y", "z")

    # magnitude vs time
    ax_mag.plot(t_arr, mag_arr)
    ax_mag.set_xlabel("Time [s]")
    ax_mag.set_ylabel(r"$\| \vec{r} \| [$R_E$]$")
    ax_mag.set_title("Magnitude of location throughout trajectory propagation")

    fig_proj.suptitle(title, fontsize=16)
    if show_plot:
        plt.show()
    # html doesnt work with latex labels and suptitle, so its deprecated for now
    # mpld3.save_html(fig_proj,
    #                 os.path.join(PLOT_DIR, "test_trajectory_proj.html"))
    plt.savefig(os.path.join(PLOT_DIR, "test_trajectory_proj.png"))
    plt.clf()

    # plot the 3-d trajectory with wireframe sphere as the earth
    if check_3dtraj:
        # set the figure and axes
        fig_3d = plt.figure(figsize=(12, 9))
        ax_3d = fig_3d.add_subplot(111, projection="3d")

        # set the earth wireframe
        u, v = np.mgrid[0:2 * np.pi:20j, 0:np.pi:
                        10j]  # x in a:b:xj is number of points in [a,b]
        x_sphere = np.sin(u) * np.cos(v)
        y_sphere = np.sin(u) * np.sin(v)
        z_sphere = np.cos(u)

        # plot the sphere
        ax_3d.plot_wireframe(x_sphere, y_sphere, z_sphere, color="b")

        # plot the trajectory
        cm_3d = ax_3d.scatter(x_arr, y_arr, z_arr, c=t_arr, marker='o', s=1.0)

        # labels, colorbars, limits whatnot
        cbar_3d = fig_3d.colorbar(cm_3d, ax=ax_3d)
        ax_3d.set_xlim([-2.5, 2.5])
        ax_3d.set_ylim([-2.5, 2.5])
        ax_3d.set_zlim([-2.5, 2.5])
        ax_3d.set_xlabel(r"x [$R_E$]")
        ax_3d.set_ylabel(r"y [$R_E$]")
        ax_3d.set_zlabel(r"z [$R_E$]")
        cbar_3d.ax.set_ylabel("Time [s]")
        ax_3d.set_title(title)

        # mpld3.save_html(fig_3d, os.path.join(PLOT_DIR, "test_trajectory_3d.html"))
        if show_plot:
            plt.show()
        plt.savefig(os.path.join(PLOT_DIR, "test_trajectory_3d.png"))


def get_trajectory():
    # set parameters

    # parameters for trajectory
    # particle is assumed to be proton
    q = 1
    m = 0.938

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
    dt = 1e-4
    max_time = 1.
    max_step = 10000

    # control variables for the code
    check_pmag = False  # if we want to check the momentum magnitude
    check_3dtraj = False  # if we want to check the 3d trajectory or not
    show_plot = False  # if we want to show the plot on some GUI or not

    # first create plot directory if it doesnt exist
    if not os.path.exists(PLOT_DIR):
        os.mkdir(PLOT_DIR)

    # initialize trajectory
    traj = Trajectory(
        "p+",
        zenith_angle=zenith,
        azimuth_angle=azimuth,
        particle_altitude=particle_alt,
        latitude=lat,
        longitude=lng,
        detector_altitude=detector_alt,
        rigidity=rigidity,
        bfield_type="igrf"
    )

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
        plot_momentum(traj_datadict, p0, show_plot)

    # plot the trajectory
    plot_trajectory(traj_datadict, title, check_3dtraj, show_plot)


if __name__ == "__main__":
    # should add some argparse thing later on
    get_trajectory()
