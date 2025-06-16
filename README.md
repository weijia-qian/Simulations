# 202208_wrobel_permutations

## Project Directory

This project runs simulations for the Ripley's K permutations project


## Workflow

1. Run the file `source/simulations/20240429_simulations.R`. This file runs simulations for different scenarios.The code is parallelized. Sources the following files:
  * `source/simulate_ppp.R`: contains functions used to define a simulated point process dataset
  * `source/utils_fast.R`: contains function for calculating the different methods compared, including K inhomogeneous, fperm, perm. For a given simulated dataset, this function returns a matrix of the estimated K values for the different methods, as well as computation times.

Running this file outputs the simulation results to the `output/simulation_results/` folder.


2. Run the Rmarkdown file `20240429_simulation_results.Rmd`. This file will pull together the simulation results in `output` folder and generated some plots. Plots generated are automatically saved to the `output` folder as `.png` files.
