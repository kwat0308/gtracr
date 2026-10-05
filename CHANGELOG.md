# CHANGELOG

Presents the version history of the package.

### Unreleased

- `GMRC(method="penumbra")`: penumbra-resolving cutoffs in the C++ batch engine
  (`batch_gmrc_penumbra`). Each direction gets the Cooke et al. (1991) R_L, R_U and R_eff
  on a geometric rigidity ladder, via an adaptive scan, and every direction is kept with a flag.
  `evaluate_batch(zenith=..., azimuth=...)` takes explicit direction grids.
  `store_ladder=True` returns the allowed/forbidden ladder for transmission functions.
  The default `method="bisection"` is unchanged.
- Rigidity-aware escape criterion (breaking: changes results below about 10 GV).
  The escape radius is max(escape_radius, k·C_St(p)) with the Störmer length C_St of
  the field model's dipole. The time budget is max(max_time, N·r_esc/(βc)). The
  defaults k = 4 and N = 80 come from a convergence study (≤ 0.5 % in R_eff for dipole
  and Oulu; at the South Pole the fixed criterion was +213 %, now 3 %). The criterion
  applies to all solvers, `Trajectory` and both GMRC modes. `escape_factor` and
  `time_factor` set k and N; `0, 0` restores the fixed 10 R_E / `max_time` criterion.
  The logic lives in `IntegratorControl.hpp` and `gtracr.escape`.
- RK45 error norm on the dimensionless state (breaking: changes RK45 results). Position
  errors are measured in R_E with tolerance atol + rtol·r/R_E and momentum errors in |p|
  with tolerance atol + rtol. The norm is rotation invariant. The defaults change from
  `atol=1e-3, rtol=1e-6` (raw SI components, which left the momentum uncontrolled) to
  `atol=0, rtol=1e-7`, chosen by a convergence study against rtol = 1e-10. The RK4
  default `dt = 1e-5` is unchanged; it is accurate to one 0.5 % ladder rung.
- Penumbra ladder: rung edges are shared, so R_L ≤ R_eff ≤ R_U holds exactly.
- Batch mode with `bfield_type="dipole"` now builds the dipole tracer directly instead of
  reaching it through the table tracer's out-of-range fallback. The results do not change.

### Version 1.0.0

Initial release of the gtracr package.
Major bugs are now fixed and is now ready for initial release.

- Examples are cleaned up and are ready to be used with users.
- Plotter now plots multiple trajectories
- Created a new class for array storage / operator overloading in replace to std::array
- Bug fixes for compatibility with UNIX systems and Windows

### version 0.10.0

Refactoring of the code was performed:

- Evaluation of the geomagnetic cutoffs are now done by class-level objects
- Plotters are separated with the package / evaluation process

### version 0.7.0

The initial pre-release of the beta version.

- The code evaluates the cosmic ray trajectories and the geomagnetic cutoff rigidities.
- It does _not_ support any sort of parallelization (GPU / CPU).
