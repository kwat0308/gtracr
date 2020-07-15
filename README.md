# gtracr

gTraCR - A **G**PU-based **Tra**cking simulation for **C**osmic **R**ays

This code simulates a cosmic ray from any longitude, latitude, altitude, local zenith angle and local azimuthal angle at any energy (in GeV) or rigidity (in GV) simultaneously with the aid of GPU acceleration.

## Current progress

- Simulates for proton / anti-proton (p+/-), electron / positron (e-/+)
- Uses the ideal dipole approximation for the Earth's magnetic field
- Works at low energies (~50GeV) but doesnt work well at higher energies (>1PeV)
- ~~Assumes a non-relativistic particle~~ This now fully accounts for relativistic charged particles.
- We have established a (preliminary) geomagnetic cutoff evaluator that evaluates the valid zenith and azimuth angles at different locations on Earth

  - Kamioka (the location of Super-K)
  - South Pole (location of IceCube)
  - University of Alberta

- Different locations can easily be set / integrated as long as one knows the geographic latitude and longitude of the location.

## Requirements

- Python 3 and above
- NumPy, matplotlib
- **pip** (if Python version is <=3.4)
- C++11 or above
- A C++ compiler (g++, CLANG++, or Microsoft C++ compiler)
- pybind11

## Installation

By either using **pip** or **python setup.py install** 

<!-- The module containing the relavent classes can be installed in the following ways: <!-- 1\. using **pip**: **pip install .** 2\. using **python**: **python setup.py build_ext**, then **python setup.py install** The modules for the user-defined Matrix class (matrices created from C-arrays) can be imported with **Matrix**. The module for matrices constructed from Boost can be imported with **BoostMatrix**. 3\. If a different C++ compiler is being used, make sure to change the compiler configurations (**os.environ["CC"]**) in the top of the file **setup.py**. 4\. If you do not have **pybind11** installed, install from **pip** by using the following command: **pip install pybind11**. Alternatively, this can be downloaded using **git clone** from [here](https://github.com/pybind/pybind11/tree/stable). 5\. If you do not have **Boost** installed, install it from [here](https://www.boost.org/doc/libs/1_73_0/more/getting_started/windows.html). -->

 ## Running the simulation

A single trajectory can be simulated by performing the following steps:

1. Create a Trajectory object, which contains the following information:

  - Particle label: the label of the particle in which we want to simulate.

    - List of available particles are indicated above.

  - The geodesic coordinates of the detector (i.e. the geographic latitude, longitude, and altitude of the detector).

    - The latitude and longitudes are defined by their decimal degree notation. The latitude is valid in the interval [+90, -90], and the longitude is valid in the interval [-180, 180].

      - For example, New York, defined in GPS coordinates as (40.7128° N, 74.0060° W), is defined in decimal degrees as (40.730610, -73.935242).

    - The altitude, that is, the height above sea level of the detector, is defined in kilometers.

    - The latitude and longitudes are defined in a geographic sense, which indicates that a latitude of 90 degrees would correlate to the North Pole.

  - The arrival direction and altitude of the cosmic ray, i.e. the zenith angle, azimuthal angle, and altitude of the cosmic ray

    - All such values are defined within the detector's reference frame, that is, such angles and altitude defined from the location of the detector.

    - The zenith angle is defined as the angle from the local zenith, where a zenith angle of 0 degrees indicate the local zenith itself, and that of 90 degrees indicate the horizon. The valid angles are in the interval [0, 180].

    - The azimuthal angle is defined where 0 degrees points to the North Pole. Valid angles are in the interval of [0, 360].

    - The altitude, defined in km, is the altitude in which the cosmic ray collides with the atmosphere and creates particle showers.

  - The energy (center-of-mass energy) or the rigidity of the cosmic ray, defined in GeV (for energy) or GV (for rigidity).

    - Only one or the other may be used as valid inputs, otherwise an error will be thrown.
    - The rigidity is defined as such: pc / Ze, where p is the momentum of the particle, Z is the charge of the particles (defined in e), and c is the speed of light in vacuum.

  - Optional parameters:

    - escape_altitude: the altitude in which we declare that the cosmic ray trajectory is allowed. Default to 10 times the Earth's radius. Also defined in terms of kilometers.

2. Run get_trajectory().

  - This will run the simulation by performing the integration to the particle trajectory with the specifications.
  - Optional parameters:

    - max_steps: the maximum number of iterations to make for the integrator. Default 10000.
    - step_size: the step size, i.e. the time step for the integration process. Default 1e-5.

3. Run get_plotting_variables().

  - This will return a dictionary of Cartesian coordinates that contain the information of the particle trajectory. Using this, one can plot the trajectory onto some device of their wishing.

An example of the simulation can be observed by running **test_trajectory.py**. The specifications can be changed manually within the code.

## Geomagnetic Cutoffs

The geomagnetic cutoffs are the cutoffs to the valid particle trajectories at a specified location on Earth. A heatmap of the valid trajectories (with the necessary configurations (particle type, particle energy / rigidity, location, altitude)) can be obtained by running **geomagnetic_cutoff.py**.

- The locations must be set from **add_location.py**, which will create a global database for different locations on Earth. 

<!-- - Performance tests can be run by **python tests/test_from_nparray.py** and **python tests/test_from_datafile.py**, where tests are run by using numpy arrays or constructed data files as valid inputs respectively. The maximum matrix dimension (n) for a n-by-n matrix is set to 1000\. A plot showing the performance as a function of dimension will be displayed. The constructed numpy array / datafile will have random values from 0-1 as elements with the provided shape. Available flags: - **--verbosity, -v**: Set verbosity level (integer from 0-4). Default is 0\. - **--debug, -d**: Activate debugging mode. Sets verbosity to level 4 and presets maximum dimension (n) to 50\. <!-- Performance test can also be run by **python tests/test_matrix.py**. This compares the performance for a small and large matrix. The rows and columns of the smaller matrix, as well as the scaling factor of the smaller vs larger matrix, is set by user input. The constructed numpy array / datafile will have random values from 0-1 as elements with the provided shape. Available flags: - **--verbosity, -v**: Set verbosity level (integer from 0-4). Default is 0\. - **--mode, -m**: Set whether to benchmark performance from a datafile, a numpy array, or both (inputs: datafile, np_array, or both). Default is both. - **--debug, -d**: Activate debugging mode. Sets verbosity to level 4 and presets rowsize = 3, columnsize = 4, scale = 50\. ## Tasks - [x] Implement performance benchmarks in Python with datafiles - [x] Implement performance benchmarks in Python with iterations - [x] Implement constructors that pass numpy arrays by reference - [x] Integrate Boost performance test with numpy array tests - [x] Create performance plot (dimension vs time) - [x] Fix slow performance for C++ matrices - [ ] Add more descriptions for distribution in setup.py - [ ] Implement performance benchmarks in C++ - [ ] Create a Makefile for performance benchmarks in C++ - [ ] Add uneven matrices for plotting -->
