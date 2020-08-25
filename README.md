# gtracr - A **G**PU-accelerated **Tra**cking simulation for **C**osmic **R**ays

gtracr is a 3-D simulation package that simulates the trajectories of cosmic rays that arrive at a certain location on any location around the globe. The package uses the [IGRF (International Geomagnetic Reference Field) model](https://www.ngdc.noaa.gov/IAGA/vmod/igrf.html) as the Earth's magnetic field and simulate trajectories using a 4th-order Runge Kutta numerical integration method.

The code can produce geomagnetic cutoff rigidities that either validates or invalidates a cosmic ray based on its trajectory. With the aid of using C++ programs at the core of the computation, the evaluation of each trajectory can be

## Dependencies

- Python 3 and above
- NumPy
- SciPy

### Optional requirements

These packages are required to observe plots and test different trajectory cases:

- matplotlib <!-- - PlotLy --> 
- pytest

## Installation

The package can be installed by using **pip**:

```
pip install gtracr
```

There are other alternatives to the installation, for example by using `setup.py` (make sure to add the `--user` flag when installing in this fashion).

## Documentation

Check out the [documentation](docs) for more details on using the package.

## Copyright and License

This project is under the BSD 3-Clause License. See <LICENSE> for more details.
