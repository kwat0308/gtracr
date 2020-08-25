from gtracr.trajectory import Trajectory
'''
Obtains the geomagnetic cutoff for each zenith and azimuthal angle

Structure will be much similar to test_trajectory.py
'''

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import pickle
from scipy.interpolate import griddata
# import matplotlib.tri as tri

# sys.path.append(os.getcwd())
# sys.path.append(os.path.join(os.getcwd(), "gtracr"))
# add filepath of gtracr to sys.path
CURRENT_DIR = os.path.dirname(os.path.realpath(__file__))
PARENT_DIR = os.path.dirname(CURRENT_DIR)
sys.path.append(PARENT_DIR)

# from gtracr.add_location import location_dict
# from gtracr.add_particle import particle_dict

# path for plots
PLOT_DIR = os.path.join(PARENT_DIR, "..", "gtracr_plots")
# create directory if gtracr_plots dir does not exist
if not os.path.isdir(PLOT_DIR):
    os.mkdir(PLOT_DIR)


def export_as_pkl(fpath, ds):
    with open(fpath, "wb") as f:
        pickle.dump(ds, f, protocol=-1)


def plot_heatmap(data_arr,
                 rigidity_list,
                 locname,
                 pname,
                 ngrid_azimuth=70,
                 ngrid_zenith=70,
                 save_plot=False):
    '''
    Plot the heatmap of the geomagnetic rigidity cutoffs at the specified location and type of particle that the cosmic ray constitutes of. 

    Parameters
    ----------
    - data_arr (list of tuples):
        An array that consists of the azimuthal and zenith components of the particle trajectory and its corresponding rigidity cutoff in the following order: (azimuth, zenith, rigidity cutoff)

    - rigidity_list (array of floats):
        The list of rigidities in which each Monte Carlo iteration had evaluated the rigidity cutoff for. Required to determine colorbar limits and contour levels.

    - locname (str): 
        The name of the detector location.

    - pname (str):
        The label of the particle that constitutes the cosmic ray. 

    - ngrid_azimuth (int, default=70):
        The number of grid points used for the interpolation process for the azimuthal component. The default is set to create a nice plot with linear interpolation.

    - ngrid_zenith (int, default=70):
        The number of grid points used for the interpolation process for the zenith component. The default is set to create a nice plot with linear interpolation.

    - save_plot (bool):
        Decides to choose to save the generated plot or not. If True, saves the plot in a separate parent directory named `gtracr_plots`. If False, only presents the plot in a GUI window. 

    '''

    # get azimuth, zenith, and rigidity cutoff arrays
    azimuth_arr, zenith_arr, rigidity_cutoffarr = zip(*data_arr)

    # interpolate the rigidity cutoffs in a 2-d sense using griddata
    azimuth_grid = np.linspace(np.min(azimuth_arr), np.max(azimuth_arr),
                               ngrid_azimuth)
    zenith_grid = np.linspace(np.max(zenith_arr), np.min(zenith_arr),
                              ngrid_zenith)

    rigidity_cutoffgrid = griddata(points=(azimuth_arr, zenith_arr),
                                   values=rigidity_cutoffarr,
                                   xi=(azimuth_grid[None, :],
                                       zenith_grid[:, None]),
                                   method='linear')

    # plot the contour plot
    # we use imshow to create a mock filled contour plot
    # and plot contour lines over the imshow plot
    fig, ax = plt.subplots(figsize=(12, 9), constrained_layout=True)

    image = ax.imshow(rigidity_cutoffgrid,
                      extent=[-2.5, 362.5, -2.5, 182.5],
                      origin='upper',
                      cmap="RdBu_r",
                      interpolation="bilinear",
                      aspect="auto",
                      vmin=np.min(rigidity_list),
                      vmax=np.max(rigidity_list),
                      alpha=1.)
    ax.axis('image')

    ax.contour(azimuth_grid,
               zenith_grid,
               rigidity_cutoffgrid,
               colors="k",
               linewidths=0.5,
               levels=len(rigidity_list),
               alpha=1.)

    # shrink parameter should change accordingly to
    # figsize (trial and error for now...)
    cbar = fig.colorbar(image, ax=ax, shrink=0.6)
    cbar.ax.set_ylabel("Rigidity [GV]")

    # ylim from 180 to 0 to follow convention in Honda 2002 paper
    ax.set_xlim([0., 360.])
    ax.set_ylim([180., 0.])

    ax.set_xlabel("Azimuthal Angle [Degrees]")
    ax.set_ylabel("Zenith Angle [Degrees]")
    ax.set_title("Geomagnetic Rigidity Cutoffs at {0} for {1}".format(
        locname, pname))

    # only save plot is save_plot is True
    # otherwise show the plot in a GUI window
    if save_plot:
        plt.savefig(
            os.path.join(PLOT_DIR,
                         "{0}_{1}_cutoffplot.png".format(locname, pname)))
    else:
        plt.show()


def plot_scatter(data_arr, locname, pname, save_plot=False):
    '''
    Plot the scatter plot of the geomagnetic rigidity cutoffs at the specified location and type of particle that the cosmic ray constitutes of. 

    Parameters
    ----------
    - data_arr (list of tuples):
        An array that consists of the azimuthal and zenith components of the particle trajectory and its corresponding rigidity cutoff in the following order: (azimuth, zenith, rigidity cutoff)

    - locname (str): 
        The name of the detector location.

    - pname (str):
        The label of the particle that constitutes the cosmic ray. 

    - save_plot (bool):
        Decides to choose to save the generated plot or not. If True, saves the plot in a separate parent directory named `gtracr_plots`. If False, only presents the plot in a GUI window. 

    '''

    # get azimuth, zenith, and rigidity cutoff arrays
    azimuth_arr, zenith_arr, rigidity_cutoffarr = zip(*data_arr)

    fig, ax = plt.subplots()
    sc = ax.scatter(azimuth_arr, zenith_arr, c=rigidity_cutoffarr, s=2.0)
    ax.set_xlabel("Azimuthal Angle [Degrees]")
    ax.set_ylabel("Zenith Angle [Degrees]")
    ax.set_title("Geomagnetic Rigidity Cutoffs at {0} for {1}".format(
        locname, pname))

    cbar = fig.colorbar(sc, ax=ax)
    cbar.ax.set_ylabel("Rigidity [GV/c]")

    ax.set_xlim([0., 360.])
    ax.set_ylim([180., 0.])

    # only save plot is save_plot is True
    # otherwise show the plot in a GUI window
    if save_plot:
        plt.savefig(os.path.join(
            PLOT_DIR, "{0}_{1}_scatterplot.png".format(locname, pname)),
            dpi=800)
    else:
        plt.show()


def evaluate_rcutoff(rigidity_list,
                     detector_coord,
                     particle_label="p+",
                     iter_num=10000):
    '''
    Evaluate the rigidity cutoff value at some provided location
    on Earth for a given cosmic ray particle.

    Parameters
    ----------
    - rigidity_list (list / NumPy array):
        The list of rigidities to evalue the minimum rigidity cutoff value for. 

    - detector_coord (NumPy array, float, size=3):
        The geodesic coordinates of the location of the detector. This must provide the latitude [decimal], longitude [decimal], and altitude [km] of the detector in the given order. 

    - particle_label (str, default="p+") :
        The label of the particle that is acting as the cosmic ray in the simulation. 
        Current valid inputs are: "p+", "p-", "e+", "e-" 

    - iter_num (int, default 10000):
        The number of iterations to perform the Monte Carlo integration with. 

    Returns
    -------
    - data_arr (list of tuple of size 3, size=iter_num):
        An array that contains the tuple of the direction of the particle trajectory and its correlating rigidity cutoff value in the following format: (azimuth, zenith, rigidity_cutoff)

    '''

    # unpack the detector coordinates in latitude, longitude, altitude
    (latitude, longitude, detector_altitude) = detector_coord

    # array that contains the data as a tuple, that is,
    # the azimuth, zenith, and cutoff rigidity
    data_arr = []

    # perform Monte Carlo integration to get cutoff rigidity
    for i in range(iter_num):
        # get a random zenith and azimuth angle
        # zenith angles range from 0 to 180
        # azimuth angles range from 0 to 360
        [azimuth, zenith] = np.random.rand(2)
        azimuth *= 360.
        zenith *= 180.

        # print("Zenith Angle: {0}, Azimuth Angle {1}\n".format(zenith, azimuth))

        # iterate through each rigidity, and break the loop
        # when particle is able to escape earth
        for k, rigidity in enumerate(rigidity_list):
            # print("Current rigidity: {:.4e}".format(rigidity))
            traj = Trajectory(plabel=particle_label,
                              latitude=latitude,
                              longitude=longitude,
                              detector_altitude=detector_altitude,
                              zenith_angle=zenith,
                              azimuth_angle=azimuth,
                              particle_altitude=100.,
                              rigidity=rigidity)
            traj.get_trajectory(max_step=10000)
            # break loop and append direction and current rigidity if particle has escaped
            if traj.particle_escaped == True:
                data_arr.append((azimuth, zenith, rigidity))
                break

        # progress checker
        if i % (iter_num // 10) == 0:
            print("{:d} iterations done.".format(i))

    return data_arr


if __name__ == "__main__":

    # create particle trajectory with desired particle and energy
    rigidity_list = np.arange(5, 55, 5)
    particle_list = [("p+", particle_dict["p+"])
                     ]  # , ("e-", particle_dict["e-"])]
    location_list = [("Kamioka", location_dict["Kamioka"])]

    iter_num = 10000  # total number of points used for Monte Carlo process

    # number of points for azimuth / zenith grid
    ngrid_azimuth = 70
    ngrid_zenith = 70

    # debug mode (bool)
    # should be done in a better fashion like argparse or smth
    debug_mode = True

    # boolean to save plot or not
    save_plot = False if debug_mode else True

    # locations: kamioka, icecube
    # for locname, loc in list(location_dict.items()):
    for (locname, loc) in location_list:
        # print(loc)
        detector_coord = np.array([loc.latitude, loc.longitude, loc.altitude])
        # geomag_cutoffdict["Location"][locname] = {}
        # for pname, particle in list(particle_dict.items()):
        for (pname, particle) in particle_list:
            # geomag_cutoffdict["Location"][locname][pname] = {}
            # evaluate rigidity cutoffs at the location for
            # the specified particle
            data_arr = evaluate_rcutoff(rigidity_list,
                                        detector_coord,
                                        particle_label=pname,
                                        iter_num=iter_num)
            # geomag_cutoffdict["Location"][locname][pname][
            # rigidity] = cutoff_arrs[k]

            # create a debugger / checker as a scatter plot
            # of the dataset
            if debug_mode:
                plot_scatter(data_arr, locname, pname, save_plot=save_plot)

            plot_heatmap(data_arr,
                         rigidity_list,
                         locname=locname,
                         pname=pname,
                         ngrid_azimuth=ngrid_azimuth,
                         ngrid_zenith=ngrid_zenith,
                         save_plot=save_plot)

    # fpath = os.path.join(os.getcwd(), "geomagnetic_cutoff.pkl")
    # export_as_pkl(fpath, geomag_cutoffdict)
'''

Below info will be used in a future release so its just 
stored here for now.

    - min_rigidity (float, default 5):
        The minimum rigidity in GV to evaluate the minimum rigidity cutoff value. 
        Note: For energies ~ particle mass, the particle may be trapped within the Earth's magnetic field lines, and thus such low energies need not be considered. 

    - max_rigidity (float, default 60):
        The maximum rigidity in GV to evaluate the minimum rigidity cutoff value. The rigidity cutoff value diminishes quickly at high rigidities, so rigidity values > 100 GV are usually not necessary.

    - delta_rigidity (float, default 5):
        The spacing 
'''
