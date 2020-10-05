'''
Evaluates the benchmarks between different versions of the code
'''

import os, sys
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import pickle

from gtracr.trajectory import Trajectory

import time

CURRENT_DIR = os.path.dirname(os.path.realpath(__file__))
PARENT_DIR = os.path.dirname(CURRENT_DIR)

DATA_DIR = os.path.join(PARENT_DIR, "gtracr", "data")
PLOT_DIR = os.path.join(PARENT_DIR, "..", "gtracr_plots")


def read_pkl(fpath):
    with open(fpath, "rb") as f:
        benchmark_data = pickle.load(f)
    return benchmark_data


def write_pkl(fpath, datadict):
    with open(fpath, "wb") as f:
        pickle.dump(datadict, f, protocol=-1)


def get_evaltime(iter_num,
                 initial_variables,
                 bfield_type,
                 use_unvec=False,
                 use_python=False):
    '''
    Evaluates iter_num number of iterations of the same trajectory calculation and return the average evaluation time for those number of iterations

    Parameters
    ----------

    - iter_num : int
        the number of iterations to evaluate for 
    - initial_variables : tuple
        initial parameters to initialize trajectory
    - bfield_type : str
        type of magnetic field model to use
    - use_unvec : bool
        whether to use unvectorized version of Trajectory evaluator or not (default False)
    - use_python : bool
        whether to use Python version of Trajectory evaluator or not (default False)
    '''

    plabel, zenith, azimuth, part_alt, lat, lng, dec_alt, rig = initial_variables

    trajectory = Trajectory(plabel=plabel,
                            latitude=lat,
                            longitude=lng,
                            detector_altitude=dec_alt,
                            zenith_angle=zenith,
                            azimuth_angle=azimuth,
                            particle_altitude=part_alt,
                            rigidity=rig,
                            bfield_type=bfield_type)

    eval_time = 0.  # initialize evaluation time
    dt = 1e-5  # step size in integration
    max_step = 10000  # max steps in integration

    for i in tqdm(range(iter_num)):
        # start counter
        # perf counter for more precise time evaluations
        start_time = time.perf_counter()
        trajectory.get_trajectory(dt=dt,
                                  max_step=max_step,
                                  use_python=use_python,
                                  use_unvectorized=use_unvec)
        stop_time = time.perf_counter()
        eval_time += stop_time - start_time

    # evaluate average
    avg_evaltime = eval_time / iter_num

    return avg_evaltime


def get_evaltime_data(iternum_list):
    '''
    Evaluates the evaluation time for different iteration values, and store them into a dictionary.

    Parameters
    -----------

    - iternum_list : list
        list of maximum iterations to perform
    '''
    # set geographic location parameters

    # initial variables chosen since it has a well-defined trajectory with # steps with dt=1e-5, max_time=0.1
    # total step size for integrator is 2119
    initial_variables = ("p+", 20., -30., 100., 0., 0., 0., 40.)

    # initialize performance benchmark time arrays

    avg_evaltime_dict = {
        "pydip": np.zeros(len(iternum_list)),
        "pyigrf": np.zeros(len(iternum_list)),
        "cppdip_novec": np.zeros(len(iternum_list)),
        "cppigrf_novec": np.zeros(len(iternum_list)),
        "cppdip_vec": np.zeros(len(iternum_list)),
        "cppigrf_vec": np.zeros(len(iternum_list))
    }

    # evaluate average evaluation times for each iteration
    for i, iter_num in enumerate(iternum_list):

        # python, dipole
        avg_evaltime_dict["pydip"][i] = get_evaltime(iter_num,
                                                     initial_variables,
                                                     bfield_type="dipole",
                                                     use_python=True)
        # c++, dipole, scalar form
        avg_evaltime_dict["cppdip_novec"][i] = get_evaltime(
            iter_num, initial_variables, bfield_type="dipole", use_unvec=True)
        # C++, dipole, vector form
        avg_evaltime_dict["cppdip_vec"][i] = get_evaltime(iter_num,
                                                          initial_variables,
                                                          bfield_type="dipole")
        if iter_num <= 100:
            # python, igrf
            avg_evaltime_dict["pyigrf"][i] = get_evaltime(iter_num,
                                                          initial_variables,
                                                          bfield_type="igrf",
                                                          use_python=True)

        # C++, igrf, scalar form
        avg_evaltime_dict["cppigrf_novec"][i] = get_evaltime(
            iter_num, initial_variables, bfield_type="igrf", use_unvec=True)

        # C++, igrf, vector form
        avg_evaltime_dict["cppigrf_vec"][i] = get_evaltime(iter_num,
                                                           initial_variables,
                                                           bfield_type="igrf")

    # return the data
    return avg_evaltime_dict


def plot_benchmarks(benchmark_data, iternum_list):
    '''
    Plot benchmark results for each maximal iteration
    '''

    label_arr = [
        "Python, Dipole (x1e-3)", "Python, IGRF (x1e-4)",
        "C++, Dipole (Scalar)", "C++, IGRF (Scalar)", "C++, Dipole (Vector)",
        "C++, IGRF (Vector)"
    ]

    color_arr = ["b", "m", "g", "c", "r", "y"]

    fig, ax = plt.subplots(figsize=(12, 9), constrained_layout=True)

    for i, (code_type,
            avg_evaltime) in enumerate(list(benchmark_data.items())):
        if code_type.find("pyigrf") != -1:
            avg_evaltime *= 1e-4
        elif code_type.find("pydip") != -1:
            avg_evaltime *= 1e-3
        ax.plot(iternum_list,
                avg_evaltime,
                label=label_arr[i],
                color=color_arr[i],
                marker="o",
                ms=3.0)

    ax.set_xlabel("Number of Iterations", fontsize=14)
    ax.set_ylabel("Average Evaluation Time [s]", fontsize=14)
    ax.set_title("Performance Benchmarks for Trajectory Evaluations",
                 fontsize=16)
    ax.legend()

    plt.savefig(os.path.join(PLOT_DIR, "benchmark_plot.png"))


if __name__ == "__main__":
    # set geographic location parameters

    # initial variables chosen since it has a well-defined trajectory with # steps with dt=1e-5, max_time=0.1
    # total step size for integrator is 2119
    initial_variables = ("p+", 20., 25., 100., 0., 0., 0., 10.)
    plabel, zenith, azimuth, part_alt, lat, lng, dec_alt, rig = initial_variables

    # integration parameters
    iternum_list = [10, 50, 100, 500, 1000, 5000,
                    10000]  # number of iterations
    # iternum_list = [10, 100]

    # check if the benchmark data file is already contained in data/, otherwise run the above
    fpath = os.path.join(DATA_DIR, "benchmark_data.pkl")
    if os.path.exists(fpath):
        benchmark_data = read_pkl(fpath)

    else:
        benchmark_data = get_evaltime_data(iternum_list)
        write_pkl(fpath, benchmark_data)

    # plot results
    plot_benchmarks(benchmark_data, iternum_list)
