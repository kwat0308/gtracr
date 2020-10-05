# Geomagnetic Cutoff Rigidities

One of the key features of this code is to evaluate the geomagnetic cutoff rigidities of a certain location on the globe.

The geomagnetic cutoff rigidity for a given location on Earth is defined as the minimum rigidity in which a cosmic ray can enter Earth's atmosphere at a certain arrival direction. We distinguish those that can enter the Earth's atmosphere as an _allowed_ trajectory, and those which cannot as a _forbidden_ one.

As the cutoff between an allowed and forbidden trajectory is (in most cases) smooth,[^1] we can construct clear cutoff rigidities at each arrival direction of the particle.

[^1]: At low enough rigidities and with lighter particles such smoothness may not apply (check [Smart & Shea](https://ui.adsabs.harvard.edu/link_gateway/2005AdSpR..36.2012S/doi:10.1016/j.asr.2004.09.015) for more details).

## Evaluating the Geomagnetic Cutoff Rigidities

The evaluation of the geomagnetic cutoff rigidities can be done in the following steps.

### 1. Initialize the Geomagnetic Cutoff Rigidity Evaluator

We first have to initialize the container that performs the evaluation of the geomagntic cutoff rigidities and stores the relavent cutoff rigidities.

We first import the module as such:

```
from gtracr.geomagnetic_cutoffs import GMRC
```

The container can then be initialized by providing the name of the location.

For example, if we want to initialize the evaluator to determine the cutoff rigidities for the Kamioka site, we write the following code block:

```
gmrc = GMRC("Kamioka")
```

One can additionally add some optional arguments for the evaluator:

- `iter_num (int)` : number of iterations for Monte Carlo sampling (default:10000)
- `particle_altitude (float)` : altitude in which cosmic ray collides with atmosphere (default:100)
- `bfield_type (str)` : type of magnetic field model to use (default: igrf)

Note that only the name of the location is required to initialize the evaluator.

### 2. Evaluate the cutoff rigidities

We then evaluate the geomagnetic cutoffs by using a Monte-Carlo sampling scheme. This can be done by using the following code:

```
gmrc.evaluate()
```

Some additional configurations can be set as follows:

- `dt` : the step size of the integration, i.e. the time difference between each point in the trajectory (default : 1e-5s)
- `max_time` : the maximum time in which the integration occurs. No trajectory will be evaluated longer than this time (default : 1s).

### 3. Plot the results

We can then plot the results as a heatmap using the in-built heatmap plotting function. In order to plot this first, however, we need to interpolate the 2-D results, which is done by `scipy.interpolate.griddata`.

We can perform this as such:

```
interpd_gmrc_data = gmrc.interpolate_results()
```

We can then plot our results on the heatmap:

```
from gtracr.plotting import plot_gmrc_heatmap

plot_gmrc_heatmap(interpd_gmrc_data,
                        gmrc.rigidity_list,
                        locname=gmrc.location,
                        plabel=gmrc.plabel)
```

### Example

The following example evaluates the geomagnetic cutoff rigidities and returns a 2-D heatmap of
the interpolated results.

```
from gtracr.geomagnetic_cutoffs import GMRC

# initialize geomagnetic rigidity cutoff evaluator at Kamioka
# with 10000 iterations
gmrc = GMRC("Kamioka")

# evaluate with default stepsize and max_time
gmrc.evaluate()

# interpolate results
interpd_gmrc_data = gmrc.interpolate_results()

# create heatmap and save as png
from gtracr.plotting import plot_gmrc_heatmap

plot_gmrc_heatmap(interpd_gmrc_data,
                        gmrc.rigidity_list,
                        locname=gmrc.location,
                        plabel=gmrc.plabel)
```
