# CHANGELOG

Presents the version history of the package.

### Version 3.0.0 (unreleased)

Breaking changes:

- **Azimuth convention.** `azimuth_angle` (Trajectory) and the GMRC
  azimuths are now the azimuth of the arrival direction measured from
  geographic north, clockwise (0 = from north, 90 = from east; astropy AltAz,
  MCEq). Before, 0 meant from south (Honda's convention) while the docstring
  claimed a north-based angle. Convert old inputs and maps with
  `(180 - azimuth_old) % 360`; see docs/conventions.md.
- `batch_gmrc_evaluate` returns `(zenith, azimuth, rcutoff, above_max,
  total_trajectories)`; GMRC keeps directions whose cutoff exceeds
  `max_rigidity` (clipped to it, flagged in `data_dict["above_rmax"]`)
  instead of dropping them (batch) or storing 0 (Python). The Python
  rigidity list includes `max_rigidity`.
- IGRF-14 replaces IGRF-13 (`igrf14.json`, `IGRF14.shc`, `IGRF14.COF`);
  `gtracr.bfield.IGRF13` is a deprecated alias of `IGRF14`.

Fixes:

- Exact SI constants: the old 1 GeV/c = 5.36e-19 kg m/s with e = 1.602e-19 C
  made every rigidity label 0.30 % low. Proton mass 0.937272 → 0.938272088.
- RK4: the last stage was evaluated at `y + h*k3` with `k3` already scaled
  by h (not fourth order; |p| drifted by up to 40 % on long forbidden paths).
- RK4/RK45 integrate in Cartesian coordinates: no NaN at the geographic
  poles and no |p| blow-up for paths near the polar axis.
- C++ IGRF epoch selection: dates more than 2.5 years after an epoch were
  interpolated from the following epoch (up to 0.9 % field error).
- CTA-South longitude (#23).
- Docs: conventions page; corrected azimuth and energy docstrings.

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
