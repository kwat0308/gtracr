# Evaluating a Single Trajectory

Evaluating a single trajectory can be done seamlessly with a few steps.

## 1. Initialize a Trajectory object.

The first thing to do is to initalize the Trajectory object, which is a container for the specifications and initial configurations of the cosmic ray trajectory.

To get the relavant module for the initialization of the trajectory:

```
from gtracr.trajectory import Trajectory
```

The required parameters to initialize the trajectory are as follows:

- `location_name (str)` : The name of the location in which we want to the detect the cosmic rays
- `zenith_angle (float)` : the zenith angle relative to the detector location in degrees ([0, 180]). Angles > 90 are for upwards-moving particles from the other end of the Earth.
- `azimuth_angle (float)` : the azimuthal angle relative to the detector location in degrees ([0, 360]). An angle of 0 degrees points to the Geographic South Pole.
- `particle_altitude (float)` : the altitude above sea level in which the cosmic ray arrives at in km
- `rigidity (float)` : the momentum per charge of the cosmic ray at the arrival location in GV

For example, assume we want to simulate a cosmic ray that arrives at the horizon (`zenith_angle=90`) from the West (`azimuth_angle=90`) at 100km above sea level (`particle_altitude=100`) at the Kamioka site (`location_name="Kamioka"`) with a rigidity of 30 GV (`rigidity=30`):

```
traj = Trajectory(
        location_name="Kamioka",
        zenith_angle=90,
        azimuth_angle=90,
        particle_altitude=100,
        rigidity=30
    )
```

There are other optional parameters that one can configure to initialize the Trajectory object.

- `plabel (str)` : the type of particle of the cosmic ray (default=`"p+"`)
- `energy (float)` : the kinetic energy of the cosmic ray, which can be used instead of the rigidity to set the momentum of the cosmic ray (default=`None`). Cannot be used concurrently with the rigidity.
- `bfield_type (str)` : the type of magnetic field model that is being used (default=`igrf`). One can set this to `dipole` to use the dipole approximation of Earth's magnetic field instead.
- `date (str)` : the date (yyyy/mm/dd) in which we want to simulate the cosmic ray trajectory (default is set to the current date). This will change the value of the Gauss coefficients associated with the IGRF model.

## 2. Evaluate the trajectory

The next step is to evaluate the cosmic ray trajectory by using the numerical integrator.

There are a few parameters that can be set to modify the integration:

- `dt` : the step size of the integration, i.e. the time difference between each point in the trajectory (default : 1e-5s)
- `max_time` : the maximum time in which the integration occurs. No trajectory will be evaluated longer than this time (default : 1s).

Additionally, one can choose either to obtain the data of the trajectory in a dictionary format by setting the option `get_data=True`. This dictionary will return the trajectory data of the six-vector in spherical coordinates as well as the time as an array.

The keys of the dictionary are as follows: `["t", "r", "theta", "phi", "pr", "ptheta", "pphi"]`

The evaluation is performed simply by one function call:

```
traj.get_trajectory()  # without trajectory data
trajectory_data = traj.get_trajectory(get_data=True) # with data
```

## Optional steps

The following step are performed one desires to visualize the trajectory as a plot.

### 3. Plot the results

We can now plot the trajectory. We can either plot this using `matplotlib`, a standard Python module for plotting, or `PlotLy`, a HTML interface for interactive plots. The desired method is using `PlotLy`, which allows for plots obtained from both Python scripts and Jupyter cells to produce an interactive plot.

This can be done by the following script:

```
from gtracr.plotting import plot_3dtraj

plot_3dtraj(
		trajectory_data,
        file_name = "plot.png",
		mpl=False,
        plot_path="../gtracr_plots",
        title_name="Particle Trajectory"
    )
```

There are several options available when visualizing the plot:

- `file_name` : the name of the file in which we save the plot to. This needs only to be the relative path within the directory in which the plot is saved in. This should have extensions in `.html` format for PlotLy plots, and some image format (`.png`, `.jpg` etc) for matplotlib plots.
- `mpl` : choose whether to plot the trajectory using the `matplotlib` module with `plt.axes(projection="3d")` (default : False).
- `plot_dir` : the path to the directory in which we store the plots. The default is set to "gtracr_plots", placed in parallel with the parent directory.
- `title_name` : the title of the plot. This should be changed to appropriately describe the trajectory plot.

The plotter also implicitly converts the data (which is in spherical coordinates) into Cartesian coordinates.

_Note_ : To run this in JupyterLab, additional modeles are required to allow PlotLy to work. Check the [PlotLy Installation Guide](https://plotly.com/python/getting-started/) for more details. Running PlotLy on Jupyter Notebook still works without any additional requirements.

## Example

Here we present a simple example that can be run in a Python script or a Jupyter Cell:

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
