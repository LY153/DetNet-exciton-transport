# DetNet — Kinetic Monte-Carlo simulation of exciton transport in multi-layer ring-shaped antenna networks

DetNet simulates excitation-energy transfer in a multi-layer ring-shaped
antenna network and computes the trapping efficiency eta: the fraction of
an initially excited ensemble of N systems that is captured by the
reaction center. The exciton dynamics combine coherent evolution under the
network Hamiltonian with incoherent transitions given by a
modified-Redfield treatment of a 300 K Ohmic bath; the populations are
propagated with a kinetic Monte-Carlo jump process.

## Repository structure

| Folder | Description |
|---|---|
| `DetNet/` | Original version. |
| `Improved_version/` | Optimized version: O(N^2) population readout, no 4D overlap tensor, bath rates evaluated at the final time only, parallel trajectories. Also contains the Hamiltonian/cutoff diagnostic scripts. See its README. |
| `Sensitivity_analysis/` | Version with disorder and ensemble options: coupling cutoff, recombination/fluorescence rates, energy gradient, dipole and site-energy disorder, ensemble splitting, reproducible seeding. See its README. |

Each folder is self-contained: it only calls files within its own folder,
and the three folders are not cross-compatible (the function signatures of
`DetNetcalculation`/`DetNetGeneration` differ between them).

## Requirements

- MATLAB R2018b or later (developed and tested with R2025a)
- Statistics and Machine Learning Toolbox (`binornd`, `mnrnd`)
- Parallel Computing Toolbox (`parpool`, `parfor`)

## How to run

In MATLAB, `cd` into one of the three folders and run

```matlab
DeterminedNet
```

`DeterminedNet.m` is a script; all run parameters (network size, `T_max`,
convergence threshold, and — in `Sensitivity_analysis` — the disorder and
ensemble settings) are set at the top of the file and can be edited before
running.

Outputs:

- `DetNet/` and `Improved_version/`: on-screen efficiency table and eta
  plots. Per-trajectory time series are buffered to temporary `.mat` files
  in the system temp folder (`sim_*.mat`, `state_*.mat`); they can be
  deleted after a run.
- `Sensitivity_analysis/`: efficiency results saved to the current folder
  as `eta_results_*.mat` plus `*_summary.csv` and `*_details.csv`.

A run with the full parameter set (T_max up to 1e7 steps, one trajectory
per initial site) takes hours per (m, x) combination; reduce `T_max` for a
quick test. The diagnostic scripts in `Improved_version/`
(`run_cutoff_compare.m`, `CompareCutoffVsJc.m`) run in seconds to minutes.
