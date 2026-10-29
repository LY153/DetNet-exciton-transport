# Sensitivity analysis version

DetNet with disorder and ensemble options, intended for sensitivity scans
of the efficiency. Run: in MATLAB, `cd` to this folder, then type

```matlab
DeterminedNet
```

## Requirements (same as the original)

- MATLAB R2018b or later
- Parallel Computing Toolbox
- Statistics and Machine Learning Toolbox

## What changed (vs. the original code in `../DetNet`)

- **A2** Time-invariant quantities (`Rpc*dt`, `Rnc.'*dt` and the accumarray
  target base) are precomputed once before the time loop.
- **A3** The 4D overlap tensor Jcoe (size N^4, ~3.4 GB for x=8) is no longer
  built; `Rate_t_LMR` computes the required entries on the fly.
- **A4** The bath rate and pure-dephasing rate are evaluated only at the
  final time RT instead of over the full 1000-point time grid.
- **B1** One ensemble per initial site, run in parallel with `parfor` inside `DetNetcalculation`; `DeterminedNet` now
  loops over the (m, x) combinations serially.
- **C1** The per-step site-population time series are no longer recorded to
  disk. `DetNetcalculation` now returns only the efficiencies (one per
  initial site); `DeterminedNet` saves the per-(m,x) efficiency results to a
  .mat file and CSV tables (summary + per-site details).

New parameters in `DeterminedNet.m` (all optional, defaults reproduce the
original behavior): `R_cr`, `R_fl`, `Jc`, `energy_gradient_total`,
`dipole_mode`/`dipole_sigma`, `energy_disorder`, `rng_seed`,
`initial_state_mode`, `N_r`.

## Notes

- The number of parallel workers is set in `DeterminedNet` (`n_workers`,
  default `min(12, numcores)`). Raise it to your CPU core count for more
  speedup; the trajectories are independent, so it scales well.
- Outputs of `DetNetcalculation` are now `[eta, t1, c]`; the population time
  series outputs (`population_s_cell`, `totalpopulation_cell`) and the
  temporary per-trajectory .mat files have been removed. Efficiency results
  are saved by `DeterminedNet` as
  `eta_results_m<m>_<initial_state_mode>_Nr<N_r>_ED<energy_disorder>.mat`
  plus `*_summary.csv` (one row per (m, x)) and `*_details.csv` (per-site
  eta).
- To disable parallelization entirely (e.g. for debugging), replace the
  `parfor` in `DetNetcalculation` with `for`.
- New parameter `N_r` (default 1): the N-system ensemble is split into N_r
  independent parts. Each part builds its own Hamiltonian realization, so
  the random parts of the Hamiltonian construction (per-layer energy
  disorder and the gaussian_angle dipole directions)
  use independent random numbers.
  If N_r does not divide N, the first N_r-1 parts get floor(N/N_r) systems
  and the remainder is added to the last part. The final result is the
  system-count-weighted average over all parts, i.e. still an average over
  the whole N-system ensemble. N_r = 1 reproduces the original behavior;
  resume is only supported when N_r = 1.
- New parameter `energy_disorder` (default 0, in cm^-1): the nominal energy
  of layer i, `E_i = (i-1)*energy_gradient_total/(m-1) + 14100`, is replaced
  by `E_i + dE_i` with `dE_i` uniform in (-energy_disorder,
  +energy_disorder). One `dE_i` is drawn per layer and shared by all sites
  of that layer, so the intra-layer degeneracy of the model is preserved.
  The reaction-center energy is not randomized.
  With `energy_disorder > 0` the result becomes a disorder average, so N_r
  should be > 1 (e.g. 20-100); each of the N_r parts then samples an
  independent set of layer energies and the final eta is still the average
  over the whole ensemble.
- New parameter `rng_seed` (default `[]`, i.e. random each run): part r of
  an N_r run uses seed (`rng_seed + r - 1`) for its Hamiltonian, and task k
  of the parfor loop uses seed (`rng_base_seed + 1e6 + k`) for its
  Monte-Carlo draws. Setting `rng_seed` to an integer therefore makes the
  whole run bit-for-bit reproducible, which is convenient when comparing
  different `energy_disorder` values.
