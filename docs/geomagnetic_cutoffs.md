# Geomagnetic Cutoff Rigidities

The `GMRC` class computes geomagnetic rigidity cutoffs — the minimum
rigidity a cosmic ray must have to reach Earth from a given direction at a
given location.

Cutoffs are determined via Monte Carlo: random (zenith, azimuth) directions
are sampled, and for each direction the code scans rigidities from low to
high until the trajectory escapes Earth's field.

## Two evaluation modes

### `evaluate()` — Python-orchestrated

Supports all field types (`"igrf"`, `"dipole"`, `"table"`). Uses
`ProcessPoolExecutor` for `igrf`/`dipole` and `ThreadPoolExecutor` for
`table` (GIL released in C++).

```python
from gtracr.geomagnetic_cutoffs import GMRC

gmrc = GMRC(location="Kamioka", iter_num=10000, bfield_type="igrf",
            solver="rk4", n_workers=8)
gmrc.evaluate(dt=1e-5, max_time=1.)
```

### `evaluate_batch()` — C++ batch mode (fastest)

The entire MC loop — RNG, coordinate transforms, rigidity scanning, and
threading — runs in a single C++ call via `BatchGMRC`. Requires
`bfield_type="table"`. Achieves ~35k trajectories/second with `solver="rk45"`.

```python
gmrc = GMRC(location="Kamioka", iter_num=10000,
            bfield_type="table", solver="rk45")
gmrc.evaluate_batch(dt=1e-5, max_time=1., base_seed=42)
```

Pass `base_seed` for reproducible results. If omitted, a random seed is
chosen at runtime.

**Safety limit**: the C++ loop retries failed directions up to
`30 × iter_num` total attempts. If this limit is reached before
`iter_num` successful samples are collected, a `UserWarning` is raised and
`bin_results()` / `interpolate_results()` will have sparser coverage.

## Getting results

### `bin_results()` (recommended)

Bins MC samples into a regular (azimuth, zenith) grid. Fast and appropriate
for large sample counts.

```python
az_centres, zen_centres, cutoff_grid = gmrc.bin_results(
    nbins_azimuth=72, nbins_zenith=36
)
# cutoff_grid shape: (36, 72), NaN where no samples fell
```

### `interpolate_results()` (legacy)

Scattered interpolation via `scipy.interpolate.griddata`.

```python
az_grid, zen_grid, cutoff_grid = gmrc.interpolate_results()
```

## Parameters

`location` also accepts a `gtracr.location.Location` object with explicit
latitude, longitude and altitude (km). This is useful when a downstream
atmosphere model supplies the site geometry. Its coordinates take precedence
over the built-in location table, even if its name matches a registered site.
The same object is accepted by `Trajectory(location_name=...)`.

### `GMRC` constructor

| Parameter | Default | Description |
|-----------|---------|-------------|
| `location` | `"Kamioka"` | Predefined detector location. |
| `iter_num` | `10000` | Number of MC samples. |
| `bfield_type` | `"igrf"` | Field model: `"igrf"`, `"dipole"`, or `"table"`. |
| `particle_type` | `"p+"` | Particle label. |
| `date` | today | IGRF date (`"yyyy-mm-dd"`). |
| `min_rigidity` | `5.` (`0.3` for penumbra) | Lower bound of rigidity scan (GV). |
| `max_rigidity` | `55.` (`100.` for penumbra) | Upper bound of rigidity scan (GV). |
| `delta_rigidity` | `1.` | Rigidity step size (GV). |
| `n_workers` | physical CPUs | Number of parallel workers. |
| `solver` | `"rk4"` | Integrator: `"rk4"`, `"boris"`, or `"rk45"`. |
| `atol` | `1e-3` | Absolute tolerance (RK45). |
| `rtol` | `1e-6` | Relative tolerance (RK45). |
| `method` | `"bisection"` | `"bisection"` or `"penumbra"` (see below). |
| `ladder_delta` | `0.005` | Penumbra: relative rung spacing. |
| `coarse`, `guard` | `8`, `2` | Penumbra: adaptive-scan parameters. |
| `store_ladder` | `False` | Penumbra: keep the allowed/forbidden flag of every rung. |

### `evaluate_batch()` parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `dt` | `1e-5` | Time step in seconds. |
| `max_time` | `1.` | Maximum integration time per trajectory (seconds). |
| `base_seed` | `None` | Integer RNG seed for reproducibility. Random if `None`. |
| `zenith`, `azimuth` | `None` | Penumbra only: explicit directions (degrees) instead of `iter_num` random draws. |

## Penumbra-resolving cutoffs (`method="penumbra"`)

The default `method="bisection"` assumes that every rigidity above the cutoff
is allowed and every rigidity below it is forbidden. It returns the first
escaping rigidity that the bisection lands on. Inside the penumbra, where
allowed and forbidden bands alternate, that value is none of the standard
cutoffs. Directions whose cutoff lies above `max_rigidity` are dropped
and redrawn, so maps are biased where cutoffs are high.

`method="penumbra"` scans a geometric rigidity ladder
R_k = R_min (1 + δ)^k, k = 0 … K−1, K = ⌈ln(R_max/R_min)/ln(1+δ)⌉ + 1. Rung k
covers [R_k/√q, R_k √q] with q = 1 + δ. For every direction, `evaluate_batch()`
returns the cutoffs of Cooke et al. (1991):

| key | definition |
|-----|------------|
| `r_upper` (R_U) | upper edge of the highest forbidden rung |
| `r_lower` (R_L) | lower edge of the lowest allowed rung |
| `r_eff` (R_eff) | R_U − Σ allowed rung widths in [R_L, R_U] (= bottom ladder edge + Σ forbidden widths); also stored as `rcutoff` |
| `flag` | 0 resolved; 1 nothing forbidden (cutoff below the ladder, R_L = R_U = R_eff = bottom edge); 2 top rung forbidden (R_U, R_eff are lower bounds) |

Every direction is kept, including flag 1 and 2. R_eff is the cutoff to use for
fluxes: a sharp cut at R_eff gives the same transmitted rigidity range as
the band structure, measured in linear rigidity. The top rung centre may
lie slightly above `max_rigidity`.

```python
gmrc = GMRC(location="Kamioka", bfield_type="table", solver="rk45",
            method="penumbra", min_rigidity=0.3, max_rigidity=100.,
            ladder_delta=0.005)
# fixed grid of directions
zen, az = np.meshgrid(np.arange(0., 181., 5.), np.arange(0., 360., 5.))
gmrc.evaluate_batch(zenith=zen.ravel(), azimuth=az.ravel())
r_eff = gmrc.data_dict["r_eff"]
```

**Adaptive scan.** First every `coarse`-th rung is traced, plus the bottom
rung. Each coarse interval within `guard` intervals of an allowed/forbidden
change is then traced rung by rung. The other intervals take the common
value of their two end rungs. With `coarse=1` every rung is traced and the
result is exact on the ladder. The adaptive scan can miss a band narrower
than one coarse interval if it lies away from every detected change, most
often an isolated allowed band low in the penumbra. When that happens R_L
moves a lot, but R_eff hardly changes. Two IGRF directions at Kamioka
(δ = 0.5 %, 0.3–100 GV) each missed a single allowed rung: R_L moved by 22 % and 12 %,
R_eff by 0.4 %. In a 0.3–100 GV global study with the same scheme, adaptive
and full-ladder R_eff agreed to 0.41 % rms. R_eff changed by 0.35 % rms
between δ = 0.5 % and 0.25 %. Use R_L only from full-ladder scans.

**Transmission functions.** With `store_ladder=True`, `gmrc.ladder` holds the
rung centres and `gmrc.ladder_allowed` (shape (n, K), bool) the
allowed/forbidden flag of every rung and direction. Averaging the flags over
directions gives a transmission function.

**Cost.** About 190 trajectories per direction with the defaults
(0.3–100 GV, δ = 0.5 %, coarse 8, guard 2), against about 7 for bisection on
5–55 GV in 1 GV steps. That is about 30× the wall time at equal sample count. The
full ladder takes 1166 trajectories per direction. Kamioka, `bfield_type="table"`, 2000
random directions, 16 threads:

| solver | bisection | penumbra |
|--------|-----------|----------|
| `rk45` | 0.04 s (13.2k traj) | 1.2 s (382k traj) |
| `rk4`, dt = 1e-5 | 3.9 s (13.2k traj) | 114 s (382k traj) |

**When to use.** Use penumbra mode whenever the cutoff feeds a flux
calculation, a cutoff map at low or mid latitudes, or any direction close
to the horizon. These are the cases where the penumbra is wide. Bisection is
fine for quick looks where an error of a few GV inside the penumbra does not
matter.

**Escape criterion.** Both modes call a trajectory allowed when it reaches
10 R_E within `max_time`. Near the cutoff, escaping trajectories can loop
for longer than the default `max_time = 1 s`. They are then counted as
forbidden, which raises R_L. For example, in the dipole at 60° latitude
the vertical R_L is 0.99 GV at 1 s and 0.87 GV at 10 s, against a Störmer
value of 0.85 GV. Raise `max_time` (cheap with `rk45`) for cutoffs below about
1 GV. Below about 0.5 GV the Störmer length exceeds 10 R_E, so "escaped"
no longer implies "connected to infinity".

## Solver recommendation

For GMRC evaluation, **`solver="rk45"`** with **`bfield_type="table"`** is
the fastest combination. The adaptive stepper takes ~100x fewer steps than
fixed-step RK4 for allowed trajectories, and the tabulated field eliminates
spherical harmonic computation.

## Complete example

```python
import numpy as np
from gtracr.geomagnetic_cutoffs import GMRC
from gtracr.plotting import plot_gmrc_heatmap

gmrc = GMRC(
    location="Kamioka",
    iter_num=100000,
    bfield_type="table",
    solver="rk45",
    min_rigidity=1.,
    max_rigidity=55.,
    delta_rigidity=0.5,
)
gmrc.evaluate_batch(dt=1e-5, max_time=1.)

grids = gmrc.bin_results(nbins_azimuth=72, nbins_zenith=36)
print(f"Mean cutoff: {np.nanmean(grids[2]):.1f} GV")

plot_gmrc_heatmap(grids, gmrc.rigidity_list,
                  locname="Kamioka", plabel="p+", bfield_type="table")
```

![Kamioka cutoff map — 100k events, table+RK45](img/Kamioka_p+_table_cutoffplot.png)
