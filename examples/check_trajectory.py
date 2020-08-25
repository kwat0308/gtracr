from gtracr.utils import dec_to_dms
from gtracr.lib.constants import EARTH_RADIUS, KG_M_S_PER_GEVC
from gtracr.trajectory import Trajectory
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
# import mpld3
import matplotlib.patches as patches
from mpl_toolkits import mplot3d

# sys.path.append(os.getcwd())
# sys.path.append(os.path.join(os.getcwd(), "gtracr"))
# sys.path.append(os.path.join(os.getcwd(), "..", "gtracr"))
CURRENT_DIR = os.path.dirname(os.path.realpath(__file__))
PARENT_DIR = os.path.dirname(CURRENT_DIR)
sys.path.append(PARENT_DIR)

PLOT_DIR = os.path.join(PARENT_DIR, "..", "gtracr_plots")

# from gtracr.utils import spherical_to_cartesian

if __name__ == "__main__":
    # set initialization parameters
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

    # initialize trajectory
    traj = Trajectory(
        "p+",
        latitude=lat,
        longitude=lng,
        detector_altitude=detector_alt,
        zenith_angle=zenith,
        azimuth_angle=azimuth,
        particle_altitude=particle_alt,
        rigidity=rigidity,
    )

    # set integration parameters
    dt = 1e-4
    max_time = 1.
    max_step = 10000

    # obtain the trajectory result
    result = traj.get_trajectory(dt=dt,
                                 max_time=max_time,
                                 get_data=True,
                                 max_step=max_step)

    # obtain data from dictionary
    t_arr = result["t"]
    r_arr = result["r"] / EARTH_RADIUS
    theta_arr = result["theta"]
    phi_arr = result["phi"]

    # convert to cartesian
    x_arr = r_arr * np.sin(theta_arr) * np.cos(phi_arr)
    y_arr = r_arr * np.sin(theta_arr) * np.sin(phi_arr)
    z_arr = r_arr * np.cos(theta_arr)

    # check momentum magnitude vs steps since |p| should be
    # constant throughout the trajectory
    p_arr = np.sqrt(result["pr"]**2. + result["ptheta"]**2. +
                    result["pphi"]**2.) / KG_M_S_PER_GEVC
    p_ratio = p_arr / p0

    # convert lat, long in decimal notation to dms
    lat_dms, lng_dms = dec_to_dms(lat, lng)

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
    cm_xy = ax_xy.scatter(x_arr, y_arr, c=t_arr)
    circ_xy = patches.Circle((0., 0.),
                             1.,
                             alpha=0.8,
                             fc='None',
                             linestyle='-',
                             ec="k")
    ax_xy.add_patch(circ_xy)
    cbar_xy = fig_proj.colorbar(cm_xy, ax=ax_xy)
    ax_xy.set_xlim([-3, 3])
    ax_xy.set_ylim([-3, 3])
    ax_xy.set_xlabel(r"x [$R_E$]")
    ax_xy.set_ylabel(r"y [$R_E$]")
    cbar_xy.ax.set_ylabel("Time [s]")
    ax_xy.set_title("Trajectory projected onto x-y plane")

    # projection onto xz plane
    cm_xz = ax_xz.scatter(x_arr, z_arr, c=t_arr)
    circ_xz = patches.Circle((0., 0.),
                             1.,
                             alpha=0.8,
                             fc='None',
                             linestyle='-',
                             ec="k")
    ax_xz.add_patch(circ_xz)
    cbar_xz = fig_proj.colorbar(cm_xz, ax=ax_xz)
    ax_xz.set_xlim([-3, 3])
    ax_xz.set_ylim([-3, 3])
    ax_xz.set_xlabel(r"x [$R_E$]")
    ax_xz.set_ylabel(r"z [$R_E$]")
    cbar_xz.ax.set_ylabel("Time [s]")
    ax_xz.set_title("Trajectory projected onto x-z plane")

    # projection onto yz plane
    cm_yz = ax_yz.scatter(y_arr, x_arr, c=t_arr)
    circ_yz = patches.Circle((0., 0.),
                             1.,
                             alpha=0.8,
                             fc='None',
                             linestyle='-',
                             ec="k")
    ax_yz.add_patch(circ_yz)
    cbar_yz = fig_proj.colorbar(cm_yz, ax=ax_yz)
    ax_yz.set_xlim([-3, 3])
    ax_yz.set_ylim([-3, 3])
    ax_yz.set_xlabel(r"y [$R_E$]")
    ax_yz.set_ylabel(r"z [$R_E$]")
    cbar_yz.ax.set_ylabel("Time [s]")
    ax_yz.set_title("Trajectory projected onto y-z plane")

    # magnitude vs time
    mag_arr = np.sqrt(x_arr**2. + y_arr**2. + z_arr**2.)
    ax_mag.plot(t_arr, mag_arr)
    # ax_mag.set_xlim([-3, 3])
    # ax_mag.set_ylim([-3, 3])
    ax_mag.set_xlabel("Time [s]")
    ax_mag.set_ylabel(r"$\| \vec{r} \|$")
    ax_mag.set_title("Magnitude of location throughout trajectory propagation")

    fig_proj.suptitle(
        "Particle Trajectory at {:s}, {:s} with Zenith Angle {:.1f}°, \
           \n Azimuth Angle {:.1f}° and Rigidity R = {:.1f}GV".format(
            lat_dms, lng_dms, zenith, azimuth, rigidity),
        fontsize=16)
    # plt.show()
    # html doesnt work with latex labels and suptitle, so its deprecated for now
    # mpld3.save_html(fig_proj,
    #                 os.path.join(PLOT_DIR, "test_trajectory_proj.html"))
    plt.savefig(os.path.join(PLOT_DIR, "test_trajectory_proj.png"))

    # plot the 3-d trajectory with wireframe sphere as the earth

    # set the figure and axes
    fig_3d = plt.figure(figsize=(12, 9))
    ax_3d = fig_3d.add_subplot(111, projection="3d")

    # set the earth wireframe
    u, v = np.mgrid[0:2 * np.pi:20j,
                    0:np.pi:10j]  # x in a:b:xj is number of points in [a,b]
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
    ax_3d.set_title(
        "Particle Trajectory at {:s}, {:s} with Zenith Angle {:.1f}°, \
           \n Azimuth Angle {:.1f}° and Rigidity R = {:.1f}GV".format(
            lat_dms, lng_dms, zenith, azimuth, rigidity))

    # mpld3.save_html(fig_3d, os.path.join(PLOT_DIR, "test_trajectory_3d.html"))
    plt.savefig(os.path.join(PLOT_DIR, "test_trajectory_3d.png"))
    '''
    plotly implementation
    # use plotly for the nice user interface
    # import plotly.graph_objects as go
    import plotly.express as px

    df = px.data.iris()
    fig = px.scatter_3d(
        df,
        x=x_arr,
        y=y_arr,
        z=z_arr,
        color=t_arr,
        # size=2.0,
        opacity=0.7,
        labels={
            'x': 'x [$R_E$]',
            'y': 'y [$R_E$]',
            'z': 'z [$R_E$]',
            'color': 'Time [s]'
        })

    fig.update_layout(
        margin=dict(l=0, r=0, b=0, t=0),
        title="Particle Trajectory at {:s}, {:s} with Zenith Angle {:.1f}°, \
             \n Azimuth Angle {:.1f}° and Rigidity R = {:.1f}GV".format(
            lat_dms, lng_dms, zenith, azimuth, rigidity),
        scene=dict(xaxis=dict(
            nticks=6,
            range=[-10, 10],
        ),
                   yaxis=dict(
                       nticks=6,
                       range=[-10, 10],
                   ),
                   zaxis=dict(
                       nticks=6,
                       range=[-10, 10],
                   )),
        coloraxis_colorbar=dict(
            title="Time [s]",
            # thicknessmode="pixels",
            # thickness=50,
            # lenmode="pixels",
            # len=400,
            yanchor="top",
            y=0.9,
            xanchor="right",
            x=0.85,
            # ticks="outside",
            # ticksuffix=" bills",
            # dtick=5
        ))
    # fig.show()
    '''
