from gtracr.lib.constants import KG_M_S_PER_GEVC
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import matplotlib.patches as patches
import plotly.graph_objects as go

CURRENT_DIR = os.path.dirname(os.path.realpath(__file__))
ROOT_DIR = os.path.dirname(os.path.dirname(CURRENT_DIR))

# sys.path.append(PARENT_DIR)

PLOT_DIR = os.path.join(ROOT_DIR, "..", "gtracr_plots")


def plot_3dtraj(
    trajectory_data,
    title_name="Particle Trajectory",
    file_name="test_trajectory_3d.html",
    mpl=False,
    plotdir_path=PLOT_DIR
):
    '''
    Plots the trajectory using a 3-dimensional plot using PlotLy, an interactive HTML plotting module. Options are available to plot in usual matplotlib as well.

    Parameters
    ----------

    - trajectory_data : dict(str, np.array)
        The dictionary that stores the information of the trajectory in Cartesian coordinates.
    - title_name : str
        The title name of the plot. Default is setot `"Particle Trajectory"`.
    - file_name : str
        The name of the file that we want to save. This defaults to a filename with a .html extension. *The filename should not be .html if `mpl=True`.*
    - mpl : bool
        Enables plotting with matplotlib instead of PlotLy (default = False).
    - plotdir_path : str
        The path to the directory in which the plots are stored in. Default is set to a directory `gtracr_plots` placed in parallel with the root directory.
    '''

    # unpack dictionary
    x_arr = trajectory_data["x"]
    y_arr = trajectory_data["y"]
    z_arr = trajectory_data["z"]
    t_arr = trajectory_data["t"]

    # set the earth wireframe
    u, v = np.mgrid[0:2 * np.pi:50j, 0:np.pi:
                    50j]  # x in a:b:xj is number of points in [a,b]
    x_sphere = np.sin(u) * np.cos(v)
    y_sphere = np.sin(u) * np.sin(v)
    z_sphere = np.cos(u)

    # plot using matplotlib
    if mpl:
        # set the figure and axes
        fig_3d = plt.figure(figsize=(12, 9))
        ax_3d = fig_3d.add_subplot(111, projection="3d")

        # plot the sphere
        ax_3d.plot_wireframe(x_sphere, y_sphere, z_sphere, color="b")

        # plot the trajectory
        cm_3d = ax_3d.scatter(x_arr, y_arr, z_arr, c=t_arr, marker='o', s=1.0)

        # labels, colorbars, limits whatnot
        cbar_3d = fig_3d.colorbar(cm_3d, ax=ax_3d)
        ax_3d.set_xlim([-3.0, 3.0])
        ax_3d.set_ylim([-3.0, 3.0])
        ax_3d.set_zlim([-3.0, 3.0])
        ax_3d.set_xlabel(r"x [$R_E$]")
        ax_3d.set_ylabel(r"y [$R_E$]")
        ax_3d.set_zlabel(r"z [$R_E$]")
        cbar_3d.ax.set_ylabel("Time [s]")
        ax_3d.set_title(title_name)

        # make file extension to png if it is not png or jpg
        if file_name.find("png") < 0 or file_name.find("jpg") < 0:
            file_name = file_name.split(".")[0] + ".png"

        # mpld3.save_html(fig_3d, os.path.join(PLOT_DIR, "test_trajectory_3d.html"))
        plt.savefig(os.path.join(plotdir_path, file_name))

    # plot using PlotLy
    else:

        # marker setings for trajectory plot
        traj_marker = dict(
            size=4,
            color=t_arr,                # set color to an array/list of desired values
            colorscale='Viridis',   # choose a colorscale
            opacity=0.8,
            colorbar=dict(
                thickness=20,
                title="Time [s]"
            )
        )

        # construct trajectory plot object
        traj_plot = go.Scatter3d(
            x=x_arr,
            y=y_arr,
            z=z_arr,
            mode='markers',
            marker=traj_marker
        )

        # construct wireframe for sphere
        lines = []
        line_marker = dict(color='#0066FF', width=2)
        for i, j, k in zip(x_sphere, y_sphere, z_sphere):
            lines.append(go.Scatter3d(
                x=i, y=j, z=k, mode='lines', line=line_marker))

        # append them all and plot
        data = lines + [traj_plot]
        fig = go.Figure(data=data)

        # additional configurations
        # somehow latex axis labels are still not enabled with plotly?
        fig.update_layout(
            margin=dict(l=0, r=0, b=0, t=0),
            scene_aspectmode='cube',
            scene=dict(
                xaxis=dict(
                    nticks=6,
                    range=[-3.0, 3.0],
                    title=r"x [Re]"
                ),
                yaxis=dict(
                    nticks=6,
                    range=[-3.0, 3.0],
                    title=r"y [Re]"
                ),
                zaxis=dict(
                    nticks=6,
                    range=[-3.0, 3.0],
                    title=r"z [Re]"
                ),
            ),
            showlegend=False
        )

        # make file extension to html if it is not html
        if file_name.find("html") < 0:
            file_name = file_name.split(".")[0] + ".html"

        # write to html
        fig.write_html(os.path.join(plotdir_path, file_name))


def plot_projections(
    trajectory_data,
    title_name="Particle Trajectory",
    file_name="test_trajectory_proj.png",
    mpl=False,
    plotdir_path=PLOT_DIR
):
    '''
    Plots the projections of the trajectory in three different planes, and the time evolution of the magnitude of the trajectory. Only supported with matplotlib for now.

    Parameters
    ----------

    - trajectory_data : dict(str, np.array)
        The dictionary that stores the information of the trajectory in Cartesian coordinates.
    - title_name : str
        The title name of the plot. Default is setot `"Particle Trajectory"`.
    - file_name : str
        The name of the file that we want to save. This defaults to a filename with a .html extension. *The filename should not be .html if `mpl=True`.*
    - mpl : bool
        Enables plotting with matplotlib instead of PlotLy (default = False).
    - plotdir_path : str
        The path to the directory in which the plots are stored in. Default is set to a directory `gtracr_plots` placed in parallel with the root directory.
    '''
    # unpack dictionary
    x_arr = trajectory_data["x"]
    y_arr = trajectory_data["y"]
    z_arr = trajectory_data["z"]
    t_arr = trajectory_data["t"]

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
    ax_mag.set_title("Time Evolution of the Magnitude of the Trajectory")

    fig_proj.suptitle(title_name, fontsize=16)
    # if show_plot:
    #     plt.show()
    # html doesnt work with latex labels and suptitle, so its deprecated for now
    # mpld3.save_html(fig_proj,
    #                 os.path.join(PLOT_DIR, "test_trajectory_proj.html"))
    plt.savefig(os.path.join(PLOT_DIR, "test_trajectory_proj.png"))


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


def plot_momentum(trajectory_data, p0, show_plot=False):
    '''
    Plot the time evolution of the magnitude of the momentum. 
    This is mainly used for debugging purposes as a cross-check 
    to the accuracy of our trajectory simulation, as |p| should
    be constant throughout the trajectory.

    Parameters
    -----------

    - trajectory_data : dict(str, np.array)
        The dictionary that contains the momentum components of the 
        trajectory.
    - p0 : float
        The initial momentum of the trajectory
    - show_plot : bool
        Decide whether to show the plot in a GUI or not (default = False).
    '''
    # check momentum magnitude vs steps since |p| should be
    # constant throughout the trajectory
    p_arr = np.sqrt(trajectory_data["pr"]**2. + trajectory_data["ptheta"]**2. +
                    trajectory_data["pphi"]**2.) / KG_M_S_PER_GEVC
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
