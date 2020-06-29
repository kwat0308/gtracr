# gtracr
gTraCR - A **G**PU-based **Tra**cking simulation for **C**osmic **R**ays

This code simulates a cosmic ray from any longitude, latitude, altitude, local zenith angle and local azimuthal angle at any energy (in GeV) or rigidity (in GV). 

## Current progress
- Simulates for proton / anti-proton, electron / positron
- Uses the ideal dipole approximation for the Earth's magnetic field
- Works at low energies (~50GeV) but doesnt work well at higher energies (>1PeV)
- ~~Assumes a non-relativistic particle~~ This now fully accounts for relativistic charged particles.
- We have established a geomagnetic cutoff evaluator that evaluates the valid zenith and azimuth angles at different locations on Earth
  - Kamioka (the location of Super-K)
  - South Pole (location of IceCube)
  - University of Alberta
- Different locations can easily be set / integrated as long as one knows the geographic latitude and longitude of the location.

## Requirements
- Python 3 and above
- NumPy, matplotlib
<!-- - Python3 or above 
- **pip** (if Python version is <=3.4)
- numpy and matplotlib modules in Python
- C++11 or above
- uBLAS library in Boost
- A C++ compiler (g++, CLANG++, or Microsoft C++ compiler)
- pybind11  -->

## Installation
From the usual method: either using **pip** or **python setup.py install**
<!-- The module containing the relavent classes can be installed in the following ways:
1. using **pip**: **pip install .**
2. using **python**: **python setup.py build_ext**, then **python setup.py install**
The modules for the user-defined Matrix class (matrices created from C-arrays) can be imported with **Matrix**. The module for matrices constructed from Boost can be imported with **BoostMatrix**.
- If a different C++ compiler is being used, make sure to change the compiler configurations (**os.environ["CC"]**) in the top of the file **setup.py**.
- If you do not have **pybind11** installed, install from **pip** by using the following command: **pip install pybind11**. Alternatively, this can be downloaded using **git clone** from [here](https://github.com/pybind/pybind11/tree/stable).
- If you do not have **Boost** installed, install it from [here](https://www.boost.org/doc/libs/1_73_0/more/getting_started/windows.html). -->

## Testing for performance
This can be done by running **test_trajectory.py**. The specifications should be changed manually within the code.


## Geomagnetic Cutoffs
The geomagnetic cutoffs are the cutoffs to the valid particle trajectories at a specified location on Earth. A heatmap of the valid trajectories (with the necessary configurations (particle type, particle energy / rigidity, location, altitude)) can be obtained by running **geomagnetic_cutoff.py**. 
- The locations must be set from **add_location.py**, which will create a global database for different locations on Earth. 
<!-- Performance tests can be run by **python tests/test_from_nparray.py** and **python tests/test_from_datafile.py**, where tests are run by using numpy arrays or constructed data files as valid inputs respectively. The maximum matrix dimension (n) for a n-by-n matrix is set to 1000. A plot showing the performance as a function of dimension will be displayed. The constructed numpy array / datafile will have random values from 0-1 as elements with the provided shape.
Available flags:
- **--verbosity, -v**: Set verbosity level (integer from 0-4). Default is 0. 
- **--debug, -d**: Activate debugging mode. Sets verbosity to level 4 and presets maximum dimension (n) to 50.

<!--Performance test can also be run by **python tests/test_matrix.py**. This compares the performance for a small and large matrix. The rows and columns of the smaller matrix, as well as the scaling factor of the smaller vs larger matrix, is set by user input. The constructed numpy array / datafile will have random values from 0-1 as elements with the provided shape. 
Available flags:
- **--verbosity, -v**: Set verbosity level (integer from 0-4). Default is 0. 
- **--mode, -m**: Set whether to benchmark performance from a datafile, a numpy array, or both (inputs: datafile, np_array, or both). Default is both.
- **--debug, -d**: Activate debugging mode. Sets verbosity to level 4 and presets rowsize = 3, columnsize = 4, scale = 50. -->

<!--## Tasks
 - [x] Implement performance benchmarks in Python with datafiles
- [x] Implement performance benchmarks in Python with iterations
- [x] Implement constructors that pass numpy arrays by reference
- [x] Integrate Boost performance test with numpy array tests
- [x] Create performance plot (dimension vs time)
- [x] Fix slow performance for C++ matrices 
- [ ] Add more descriptions for distribution in setup.py
- [ ] Implement performance benchmarks in C++
- [ ] Create a Makefile for performance benchmarks in C++
- [ ] Add uneven matrices for plotting -->
