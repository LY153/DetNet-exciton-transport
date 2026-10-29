# DetNet (original version)

Kinetic Monte-Carlo simulation of excitation-energy transfer in a
multi-layer ring-shaped antenna network, with a modified-Redfield bath.
This folder contains the original version of the code.

## Requirements

- MATLAB R2018b or later
- Parallel Computing Toolbox
- Statistics and Machine Learning Toolbox

## How to run

1. In MATLAB, `cd` to this folder (all source files are already in one
   directory, so no path setup is needed).
2. Run:

   ```matlab
   DeterminedNet
   ```

The script prints the per-site efficiency table, plots the site-averaged
efficiency eta versus x and versus the number of sites NH, and reports the
elapsed time. Run parameters (`m_f`, `x_r`, `T_max`, convergence threshold)
can be edited at the top of `DeterminedNet.m`.
