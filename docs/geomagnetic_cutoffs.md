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
from gtracr.geomagnetic_cutoffs import GMCutoffEvaluator
```

The container can then be initialized by providing the name of the location.

For example, if we want to initialize the evaluator to determine the cutoff rigidities for the Kamioka site, we write the following code block:

```
gmcutoff_evaluator = GMCutoffEvaluator("Kamioka")
```

One can additionally add some optional arguments for the evaluator:

Note that only the name of the location is required to initialize the evaluator.

### 2. Evaluate the cutoff rigidities

We then evaluate the geomagnetic cutoffs by using a Monte-Carlo sampling scheme. This can be done by using the following code:

```
gmcutoff_evaluator.evaluate()
```

Some additional configurations can be set as follows:

### 3. Plot the results

We can then plot the results as a heatmap using the in-built heatmap plotting function. In order to plot this first, however, we need to interpolate the 2-D results, which is done by `scipy.interpolate.griddata`.

We can perform this as such:

```
gmcutoff_evaluator.interpolate_results()
```

We can then plot our results on the heatmap:

```
from gtracr.lib.plotting import plot_gmcutoffs_heatmap

plot_gmcutoff_heatmap(interpd_gmcutoff_data,
                        gmcutoff_evaluator.rigidity_list,
                        locname=gmcutoff_evaluator.location,
                        plabel=plabel)
```

### Example

Coming soon!
