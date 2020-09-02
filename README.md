# gtracr - A **G**PU-accelerated **Tra**cking simulation for **C**osmic **R**ays

**gtracr** is a 3-D simulation package that simulates the trajectories of cosmic rays that arrive at a certain location on any location around the globe :earth_americas: . The package uses the [IGRF (International Geomagnetic Reference Field) model](https://www.ngdc.noaa.gov/IAGA/vmod/igrf.html) as the Earth's magnetic field and simulate trajectories using a 4th-order Runge Kutta numerical integration method.

The main components as well as the user interface of the package is written in Python, so using this package is straighforward with minimal steps for evaluation of a trajectory. The core of the package (that is, the evaluation of the geomagnetic field and the numerical integration) is written in C++, however, and as such each trajectory is optimized to perform evaluations at around 1 ms :watch:. 

The code can further produce geomagnetic cutoff rigidities that either validates or invalidates a cosmic ray based on its trajectory, which is a key feature necessary to distinguish between allowed and forbidden trajectories.

*Note*: The current version does ***not*** support GPU parallelization. This will be done in future versions, check the [CHANGELOG](CHANGELOG) for more details.

## Dependencies

- Python 3 and above
- NumPy
- SciPy
- datetime (for obtaining the current date)
- tqdm

All such dependencies will be installed with the package.

### Optional requirements

These packages are required to observe plots and test different trajectory cases:

- `matplotlib, plotly` for plots  
- `pytest, pytest-benchmark` for testing
<!-- 
## Installation

The package can be installed by using **pip**:

```
pip install gtracr
``` -->

## Documentation

Check out the [documentation](docs) for more details on installing the package as well as instructions to start using this package with detailed examples. 

## Copyright and License

This project is under the BSD 3-Clause License. See [LICENSE](LICENSE) for more details.
