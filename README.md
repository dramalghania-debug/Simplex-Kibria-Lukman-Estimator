# Simplex Kibria-Lukman Estimator (SKLE)

Supplementary code for the manuscript "Kibria-Lukman Estimator for Simplex
Regression Model under Multicollinearity: Theory, Simulation and Applications".

## Files

- `Simulation.R` — Monte Carlo simulation study (Section 4). Computes the SMLE,
  SRRE, and SKLE at both shrinkage parameters k1 and k2.
- `Body fat new.R` — analysis of the Body Fat data (Section 5.1); data from
  Johnson (1996).
- `CI hald data.R` — analysis of the Hald Cement data (Section 5.2); data from
  Montgomery, Peck & Vining (2021).

## Software requirements

- R 4.5.1 (2025-06-13) or later
- VGAM (1.1-12), MASS, parallel

`parallel` is part of base R; VGAM and MASS install from CRAN:

    install.packages(c("VGAM", "MASS"))

**Operating system.** `Simulation.R` uses a FORK cluster, which is available on
macOS and Linux only. The simulations reported in the manuscript were run on
macOS. Windows users will need to substitute a PSOCK cluster in Section 5 of
the script (`type = "PSOCK"`); the script exports all required objects to the
workers, but this configuration was not used for the published results.

## Reproducing other dimensions (p = 4, 12, 16)

The number of explanatory variables is set by `P_VAL` at the top of the script
(line 8). To reproduce the other dimensionalities, change this value and rerun:

    P_VAL <- 8   # Set to 4, 8, 12, or 16

All other factors (sample size, multicollinearity level, dispersion, and link
function) are controlled by the same configuration block and need no
modification.

## Output

Results are written to a folder named `New_SKLE_Simulation_p<P_VAL>` on the
Desktop, as `Simulation_Results_p<P_VAL>.csv`, one row per scenario. The script
appends after each scenario and resumes from the existing file if interrupted —
delete the CSV to start a scenario grid from scratch.

The seed is fixed at `set.seed(2025)`. Each of the 256 scenarios runs 1000
replications; non-converging replications are discarded and counted in the
`Successful_Reps` column.
