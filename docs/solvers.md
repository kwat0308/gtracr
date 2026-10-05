# Solvers

gtracr provides three integration methods for the relativistic Lorentz force
equation, all implemented in C++.

## Comparison

| Solver | `solver=` | B-field evals/step | Step size | Key properties |
|--------|-----------|-------------------|-----------|----------------|
| Frozen-field RK4 | `"rk4"` | 1 | Fixed | Default; O(h^4) with frozen-field O(h^2) error |
| Boris pusher | `"boris"` | 1 | Fixed | Symplectic; ~30% faster than RK4; 2x better momentum conservation |
| Adaptive RK45 | `"rk45"` | ~6 (FSAL) | Adaptive | Dormand-Prince; uses `atol`/`rtol`; ~100x fewer steps for escaping trajectories |

## When to use each solver

### Single trajectory visualization

Use **`"rk4"`** or **`"boris"`** with `bfield_type="igrf"`:

- Fixed time steps produce regularly spaced output points, good for plotting.
- Boris is faster and conserves energy better, but outputs Cartesian
  internally (converted back to spherical).
- RK4 is the safe default.

### GMRC cutoff scanning

Use **`"rk45"`** with `bfield_type="table"`:

- Adaptive stepping takes far fewer steps for escaping trajectories.
- Combined with the tabulated IGRF (~7x faster field evaluation), this gives
  the highest throughput (~35k trajectories/s in batch mode).
- Default tolerances (`atol=0`, `rtol=1e-7`) resolve cutoffs to the 0.5 %
  ladder rung (see below).

### High-precision work

Tighten the RK45 tolerance:

```python
traj = Trajectory(..., solver="rk45", rtol=1e-9)
```

## RK45 error norm and tolerances

The RK45 step error is measured on the dimensionless state. Position errors are
lengths in units of R_E, with tolerance `atol + rtol·r/R_E`. Momentum errors are
in units of |p|, with tolerance `atol + rtol`. The norm is the RMS over the six
components, applied to the error vectors, so it does not depend on the coordinate
system. Before this change the norm was applied to the raw components (r in m,
angles in rad, momenta in kg m/s). There `atol = 1e-3` meant 6 km in angle and
left the momentum (~1e-19 kg m/s) uncontrolled. With the dimensionless norm,
`atol = 1e-3` still means 1e-3 R_E and 1e-3 |p| per step. That is far too loose:
|p| drifts by 1.3 % over a 1 s bound orbit.

Convergence study (IGRF, Kamioka / equator 100°E / Oulu, 6 directions each,
penumbra ladder 0.3–60 GV, δ = 0.5 %, fixed escape criterion, Cartesian
integrator; reference RK45 rtol = 1e-10; wall time for 6 directions × ~300
trajectories on 48 threads):

| solver, tolerance | max \|ΔR_eff\| Kamioka / equator / Oulu | wall (Oulu) |
|-------------------|-------------------------------------------|-------------|
| rk45 `atol=1e-3`, any rtol ≤ 1e-6 (old default) | 1.0 / 0.5 / 7.7 % | 0.2 s |
| rk45 `atol=0`, `rtol=1e-5` | 0 / 0.5 / 3.0 % | 0.3 s |
| rk45 `atol=0`, `rtol=1e-6` | 0 / 0 / 0.5 % | 0.4 s |
| **rk45 `atol=0`, `rtol=1e-7` (default)** | **0 / 0 / 0** | **0.6 s** |
| rk45 `atol=0`, `rtol=1e-8` | 0 / 0 / 0 | 0.9 s |
| rk4 dt = 2e-5 | 0.55 / 0.5 / 0 % | 7.2 s |
| **rk4 dt = 1e-5 (default)** | **0.06 / 0.5 / 0 %** | **14 s** |
| rk4 dt = 5e-6 | 0.86 / 0 / 0 % | 29 s |
| rk4 dt = 2e-6 | 0.43 / 0 / 0 % | 72 s |
| boris dt = 1e-5 | 0.5 / 0 / 0.5 % | 14 s |

One ladder rung is 0.5 %. RK4 at dt = 1e-5 is at the one-rung level, and
refining dt does not improve it monotonically, because single penumbra bands
flip. The default stays at dt = 1e-5. RK45 at the new default is 10–25× faster
than RK4 for the same accuracy. The tolerance table below gives the trade-off.

| `atol` | `rtol` | Use case |
|--------|--------|----------|
| `0` | `1e-6` | Quick maps (≤ 1 rung in R_eff) |
| `0` | `1e-7` | Default: cutoffs converged on a 0.5 % ladder |
| `0` | `1e-10` | Reference solutions |

## Benchmarking solvers

Run the solver comparison script to compare accuracy and performance:

```bash
python examples/eval_solver_comparison.py --n-perf 100
```
