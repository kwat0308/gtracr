from gtracr.trajectory import Trajectory
from gtracr.utils import get_locationdict, get_particledict
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
import argparse
from tqdm import tqdm
# import matplotlib.tri as tri

# sys.path.append(os.getcwd())
# sys.path.append(os.path.join(os.getcwd(), "gtracr"))
# add filepath of gtracr to sys.path
CURRENT_DIR = os.path.dirname(os.path.realpath(__file__))
PARENT_DIR = os.path.dirname(CURRENT_DIR)
PLOT_DIR = os.path.join(PARENT_DIR, "..", "gtracr_plots")
# sys.path.append(PARENT_DIR)

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
                 plabel,
                 ngrid_azimuth=70,
                 ngrid_zenith=70,
                 show_plot=False):
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

    - plabel (str):
        The label of the particle that constitutes the cosmic ray. 

    - ngrid_azimuth (int, default=70):
        The number of grid points used for the interpolation process for the azimuthal component. The default is set to create a nice plot with linear interpolation.

    - ngrid_zenith (int, default=70):
        The number of grid points used for the interpolation process for the zenith component. The default is set to create a nice plot with linear interpolation.

    - show_plot (bool):
        Decides to choose to show the generated plot or not. If True, presents the plot in a GUI window. 

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
        locname, plabel))

    # only save plot is save_plot is True
    # otherwise show the plot in a GUI window
    if show_plot:
        plt.show()

    plt.savefig(
        os.path.join(PLOT_DIR,
                     "{0}_{1}_cutoffplot.png".format(locname, plabel)))


def plot_scatter(data_arr, locname, plabel, show_plot=False):
    '''
    Plot the scatter plot of the geomagnetic rigidity cutoffs at the specified location and type of particle that the cosmic ray constitutes of. 

    Parameters
    ----------
    - data_arr (list of tuples):
        An array that consists of the azimuthal and zenith components of the particle trajectory and its corresponding rigidity cutoff in the following order: (azimuth, zenith, rigidity cutoff)

    - locname (str): 
        The name of the detector location.

    - plabel (str):
        The label of the particle that constitutes the cosmic ray. 

    - show_plot (bool):
        Decides to choose to show the generated plot or not. If True, presents the plot in a GUI window.

    '''

    # get azimuth, zenith, and rigidity cutoff arrays
    azimuth_arr, zenith_arr, rigidity_cutoffarr = zip(*data_arr)

    fig, ax = plt.subplots()
    sc = ax.scatter(azimuth_arr, zenith_arr, c=rigidity_cutoffarr, s=2.0)
    ax.set_xlabel("Azimuthal Angle [Degrees]")
    ax.set_ylabel("Zenith Angle [Degrees]")
    ax.set_title("Geomagnetic Rigidity Cutoffs at {0} for {1}".format(
        locname, plabel))

    cbar = fig.colorbar(sc, ax=ax)
    cbar.ax.set_ylabel("Rigidity [GV]")

    ax.set_xlim([0., 360.])
    ax.set_ylim([180., 0.])

    plt.savefig(os.path.join(PLOT_DIR,
                             "{0}_{1}_scatterplot.png".format(locname,
                                                              plabel)),
                dpi=800)

    if args.show_plot:
        plt.show()


def evaluate_rcutoff(rigidity_list,
                     location_name,
                     iter_num,
                     bfield_type,
                     particle_label="p+"):
    '''
    Evaluate the rigidity cutoff value at some provided location
    on Earth for a given cosmic ray particle.

    Parameters
    ----------
    - rigidity_list (list / NumPy array):
        The list of rigidities to evalue the minimum rigidity cutoff value for.  

    - location_name (str, default="Kamioka"):
        The name of the detector / location in which the neutrinos are detected. Must be a valid location in location_dict. Contained in args if eval_all is False.

    - iter_num (int, default 10000):
        The number of iterations to perform the Monte Carlo integration with. Contained in args.
        
    - bfield_type (str, default ="igrf") :
        The type of the magnetic field model to use. Contained in args.
    
    - particle_label (str, default="p+") :
        The label of the particle that is acting as the cosmic ray in the simulation. 
        Current valid inputs are: "p+", "p-", "e+", "e-" 

    Returns
    -------
    - data_arr (list of tuple of size 3, size=iter_num):
        An array that contains the tuple of the direction of the particle trajectory and its correlating rigidity cutoff value in the following format: (azimuth, zenith, rigidity_cutoff)

    '''

    # unpack the detector coordinates in latitude, longitude, altitude
    # (latitude, longitude, detector_altitude) = detector_coord

    # array that contains the data as a tuple, that is,
    # the azimuth, zenith, and cutoff rigidity
    data_arr = []

    # perform Monte Carlo integration to get cutoff rigidity
    for i in tqdm(range(iter_num)):
        # get a random zenith and azimuth angle
        # zenith angles range from 0 to 180
        # azimuth angles range from 0 to 360
        [azimuth, zenith] = np.random.rand(2)
        azimuth *= 360.
        zenith *= 180.

        # print("Zenith Angle: {0}, Azimuth Angle {1}\n".format(zenith, azimuth))

        # iterate through each rigidity, and break the loop
        # when particle is able to escape earth
        for rigidity in rigidity_list:
            # print("Current rigidity: {:.1f}".format(rigidity))
            traj = Trajectory(
                plabel=particle_label,
                #   latitude=latitude,
                #   longitude=longitude,
                #   detector_altitude=detector_altitude,
                location_name=location_name,
                zenith_angle=zenith,
                azimuth_angle=azimuth,
                particle_altitude=100.,
                rigidity=rigidity,
                bfield_type=bfield_type)
            traj.get_trajectory(dt=1e-5, max_time=1.)
            # print("Final time : {:.2e}".format(traj.final_time))
            # break loop and append direction and current rigidity if particle has escaped
            if traj.particle_escaped == True:
                # print("particle escaped")
                data_arr.append((azimuth, zenith, rigidity))
                break

        # print(data_arr)

        # progress checker
        # if i % (iter_num // 10) == 0:
        #     print("{:d} iterations done.".format(i))

    return data_arr


def eval_gmcutoff(args):
    # create particle trajectory with desired particle and energy
    rigidity_list = np.arange(5, 55, 5)
    # print(rigidity_list)
    # particle_dict = get_particledict()
    location_dict = get_locationdict()

    # particle_list = ["p+"]
    plabel = "p+"

    # number of points for azimuth / zenith grid
    ngrid_azimuth = 70
    ngrid_zenith = 70

    # change initial parameters if debug mode is set
    if args.debug_mode:
        args.iter_num = 10
        args.show_plot = True

    if args.eval_all:
        for locname in list(location_dict.keys()):
            # for (locname, loc) in location_list:
            # print(loc)
            # geomag_cutoffdict["Location"][locname] = {}
            # for pname, particle in list(particle_dict.items()):

            # for pname in particle_list:
            # print("Current configuration: {:s}, {:s}".format(
            #     loc.name, particle.name))
            # geomag_cutoffdict["Location"][locname][pname] = {}
            # evaluate rigidity cutoffs at the location for
            # the specified particle
            data_arr = evaluate_rcutoff(rigidity_list,
                                        locname,
                                        iter_num=args.iter_num,
                                        bfield_type=args.bfield_type,
                                        particle_label=plabel)

            # print(data_arr)
            # geomag_cutoffdict["Location"][locname][pname][
            # rigidity] = cutoff_arrs[k]

            # create a debugger / checker as a scatter plot
            # of the dataset
            if args.debug_mode:
                plot_scatter(data_arr,
                             locname,
                             plabel,
                             show_plot=args.show_plot)

            plot_heatmap(data_arr,
                         rigidity_list,
                         locname=locname,
                         plabel=plabel,
                         ngrid_azimuth=ngrid_azimuth,
                         ngrid_zenith=ngrid_zenith,
                         show_plot=args.show_plot)

    else:
        # evaluate rigidity cutoffs at the location for
        # the specified particle
        data_arr = evaluate_rcutoff(rigidity_list,
                                    args.locname,
                                    iter_num=args.iter_num,
                                    bfield_type=args.bfield_type,
                                    particle_label=plabel)

        # print(data_arr)
        # geomag_cutoffdict["Location"][locname][pname][
        # rigidity] = cutoff_arrs[k]

        # create a debugger / checker as a scatter plot
        # of the dataset
        if args.debug_mode:
            plot_scatter(data_arr,
                         args.locname,
                         plabel,
                         show_plot=args.show_plot)

        plot_heatmap(data_arr,
                     rigidity_list,
                     locname=args.locname,
                     plabel=plabel,
                     ngrid_azimuth=ngrid_azimuth,
                     ngrid_zenith=ngrid_zenith,
                     show_plot=args.show_plot)

    # fpath = os.path.join(os.getcwd(), "geomagnetic_cutoff.pkl")
    # export_as_pkl(fpath, geomag_cutoffdict)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=
        'Evaluates the geomagnetic cutoff rigidities of some location for N iterations using a Monte-Carlo scheme, and produces a heatmap for such geomagnetic cutoff rigidities.'
    )
    parser.add_argument(
        '-ln',
        '--locname',
        dest="locname",
        # action="store_const",
        default="Kamioka",
        type=str,
        help="Detector location to evaluate GM cutoffs.")
    parser.add_argument(
        '-n',
        '--iter_num',
        dest="iter_num",
        # action="store_const",
        default=10000,
        type=int,
        help="Number of iterations for Monte-Carlo.")
    parser.add_argument(
        '-bf',
        '--bfield',
        dest="bfield_type",
        # action="store_const",
        default="igrf",
        type=str,
        help="The geomagnetic field model used.")
    parser.add_argument('-a',
                        '--all',
                        dest="eval_all",
                        action="store_true",
                        help="Evaluate GM cutoffs for all locations.")
    parser.add_argument('--show',
                        dest="show_plot",
                        action="store_true",
                        help="Show the plot in an external display.")
    parser.add_argument('-d',
                        '--debug',
                        dest="debug_mode",
                        action="store_true",
                        help="Enable debug mode.")
    parser.add_argument('-c',
                        '--clean',
                        dest="clean_dict",
                        action="store_true",
                        help='Clean the dataset. ')

    args = parser.parse_args()
    eval_gmcutoff(args)
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
