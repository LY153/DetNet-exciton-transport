# Improved version

Optimized version of the DetNet simulation. Run: in MATLAB, `cd` to this
folder, then type

```matlab
DeterminedNet
```

## Requirements (same as the original)

- MATLAB R2018b or later
- Parallel Computing Toolbox
- Statistics and Machine Learning Toolbox

## What changed (vs. the original code in `../DetNet`)

- **A1** Site populations are now obtained from the diagonal of V*rho*V'
  with O(N^2) matrix-vector products instead of two full O(N^3) matmuls.
- **A2** Time-invariant quantities (`Rpc*dt`, `Rnc.'*dt`, `|vec_exciton|.^2`
  and the accumarray target base) are precomputed once before the time loop.
- **A3** The 4D overlap tensor Jcoe (size N^4, ~3.4 GB for x=8) is no longer
  built; `Rate_t_LMR` computes the required entries on the fly.
- **A4** The bath rate and pure-dephasing rate are evaluated only at the
  final time RT instead of over the full 1000-point time grid.
- **B1** The independent trajectories (one per initial site) are run in
  parallel with `parfor` inside `DetNetcalculation`; `DeterminedNet` now
  loops over the (m, x) combinations serially.

## Diagnostic scripts and saved results

- `run_cutoff_compare.m` runs `CompareHamiltonians` for m = 5, x = 4..8 and
  saves the results to `CutoffCompare.mat`.
- `CompareCutoffVsJc.m` sweeps the cutoff parameter Jc and plots the
  relative eigenvalue shift and the eigenstate deviation versus x, saving
  the figure to `CutoffCompareResult.fig`.
- The two figures included here come from two runs of
  `CompareCutoffVsJc.m`: `Result_1_5_10.png`/`.pdf` with `Jc_list = [1 5 10]`
  and `Result_1_5_25.png` with `Jc_list = [1 5 25]`. The committed script
  hardcodes `Jc_list = [1 5 10]`; edit that line and rerun to reproduce the
  other figure.

## Notes

- The number of parallel workers is set in `DeterminedNet` (`n_workers`,
  default `min(12, numcores)`). Raise it to your CPU core count for more
  speedup; the trajectories are independent, so it scales well.
- Outputs and interface of `DetNetcalculation` are unchanged.
- To disable parallelization entirely (e.g. for debugging), replace the
  `parfor` in `DetNetcalculation` with `for`.
