# Index

## Introduction

This package is a high-performance 3-D tracking simulation of cosmic ray trajectories for different locations around the globe. The code can simulate cosmic ray trajectories that are detected at a certain altitude from some detector location, be it at Super-K in Kamioka, Japan or IceCube at the South Pole. The code utilizes the [International Geomagnetic Reference Field (IGRF) model](https://www.ngdc.noaa.gov/IAGA/vmod/igrf.html) to represent the main magnetic field from Earth's outer core.

The cosmic ray trajectories are described by a set of coupled ordinary differential equations based on the Lorentz force ($F_L = q(\vec{v} \times \vec{B})$). The equations are expressed in spherical coordinates, the canonical coordinate system for performing calculations with spherical objects (i.e. the Earth).

Since the analytical expression of the coordinate and momentum vector cannot be obtained, the evaluation of the trajectories are performed numerically by a 4th-order Runge-Kutta algorithm. This is considered one of the most simplest yet accurate methods used for numerical integration of all kinds. There are more accurate techniques that can be implemented to improve the accuracy of the results, but with a truncation error of $O(h^5)$, we can consider this algorithm already pretty good. _In future versions,ff we will implement a better algorithm_.

The notion "high-performance" does not come without a meaning. We intend to create our package such that the evaluation of these trajectories can be performed with **1,000,000** particles **simultaneously** at approximately **5 seconds max!** This is achieved by the aid of multi-threading, i.e. parallelization of processes across different threads of CPU / GPU cores. We intend to mainly utilize the thousands of GPU cores to perform such evaluations in a very fast manner. See the performance section for more details and benchmarks. **_This is not implemented yet, please wait for future versions!_**

## Installation

The package can be installed by using `pip`:

```
pip install gtracr
```

This allows module-level usage of the package.

Alternative methods are available to download this package.

The first method is to use `setup.py` with `python`:

```
python setup.py install --user
```

_Note_: The `--user` flag is required for installing the package in this manner.

The lastest (unstable) version can be obtained by cloning the repository from the `master` branch in [GitHub](https://github.com/kwat0308/gtracr).

## Quickstart

### Evaluating a single trajectory

Here we present a simple example that can be run in a Python script (to run in Jupyter cell, set `jupyter=True` in `plot_3dtraj()`):

```
from gtracr.trajectory import Trajectory
# matplotlib not required, just required for plotting purposes
import matplotlib.pyplot as plt

# initialize a cosmic ray trajectory that arrives at the horizon
# from the West at 100km above sea level at the Kamioka site with
# a rigidity of 30 GV:

traj = Trajectory(
        location_name="Kamioka",
        zenith_angle=90.,
        azimuth_angle=90.,
        particle_altitude=100.,
        rigidity=30.
    )

# evaluate the trajectory and get the trajectory data
trajectory_data = traj.get_trajectory(get_data=True)

# plot the results
from gtracr.plotting import plot_3dtraj

plot_3dtraj(
		trajectory_data,
        file_name = "plot.html",
        plot_path="../gtracr_plots"
    )
```

This will return a PlotLy interactive plot that will display on the default web browser, and saves the html for future access.
