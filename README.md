# Simplex Kibria-Lukman Estimator (SKLE)

Code and results accompanying the manuscript
*Kibria-Lukman estimator for the simplex regression model under
multicollinearity: Theory, simulation and applications*.

## Contents

| File | Purpose |
|---|---|
| `simulation.R` | Monte Carlo study. One dimension per run; four runs give 1,024 scenarios. |
| `figures.R` | Produces the figures from the simulation output. |
| `body_fat.R` | Body Fat application. |
| `hald_cement.R` | Hald Cement application. |
| `hald_gof.R` | Goodness-of-fit tests for the transformed Hald response. |
| `Simulation_Results_p4.csv` | Simulation output, p = 4. |
| `Simulation_Results_p8.csv` | Simulation output, p = 8. |
| `Simulation_Results_p12.csv` | Simulation output, p = 12. |
| `Simulation_Results_p16.csv` | Simulation output, p = 16. |
| `mse_tables.R` | Writes the MSE tables for p = 8 at the dispersion levels the manuscript does not print. |
| `tables/` | Those twelve tables, one CSV per link function and dispersion level. |

All scripts read and write in the working directory, so they run from a
clone without editing.

## Reproducing the simulation

The dimension and the number of replications are read from the command line,
one dimension per run:

```
Rscript simulation.R 4  1000
Rscript simulation.R 8  1000
Rscript simulation.R 12 1000
Rscript simulation.R 16 1000
```

Each run crosses four sample sizes, four correlation levels, four dispersion
levels and four link functions at the given dimension, giving 256 scenarios
per run and 1,024 in total. Each scenario is replicated 1,000 times.

Each run writes one `Simulation_Results_p*.csv` to the working directory, and
resumes from the last completed scenario if the file already exists. Set
`SKLE_OUT` to write elsewhere.

The study is seeded, so a rerun reproduces the deposited files. Replications
run in parallel on `detectCores() - 1` cores; the full study takes some hours.

The four CSV files in this repository are the output used in the manuscript,
so the figures and tables can be reproduced without rerunning the simulation.

## Reproducing the figures

```
Rscript figures.R
```

Reads the four CSV files from the working directory and writes 13 figures to
`./figures`. Seven of them appear in the manuscript; the rest cover
combinations of the design factors reported in the manuscript tables and in
`tables/`.
Set `SKLE_DATA` to read the CSV files from elsewhere and `SKLE_FIGS` to write
the figures elsewhere:

```
SKLE_DATA=/path/to/csvs SKLE_FIGS=/path/to/figures Rscript figures.R
```

## MSE tables

```
Rscript mse_tables.R
```

Writes twelve CSV files to `./tables`, one for each link function at
sigma^2 = 0.5, 1 and 1.5, with p = 8. These are the dispersion levels the
manuscript does not print; its Tables 1 to 4 give the same quantities at
sigma^2 = 2, and the files follow that layout. Columns are the sample size, the
correlation, the SMLE, and the SRRE and SKLE at each of the four shrinkage
parameters, to four decimal places.

Reads `Simulation_Results_p8.csv` from the working directory. Set `SKLE_DATA`
to read it from elsewhere and `SKLE_TABLES` to write the tables elsewhere:

```
SKLE_DATA=/path/to/csvs SKLE_TABLES=/path/to/tables Rscript mse_tables.R
```

The twelve CSV files in `tables/` are the output of this script, so they can be
read directly without running it.

## Reproducing the applications

```
Rscript body_fat.R
Rscript hald_cement.R
```

Both print their results, including the estimated coefficients, the estimated
MSE and the shrinkage parameters. `hald_cement.R` also refits each model by
direct numerical maximum likelihood and reports the agreement with the IRLS
fit. The goodness-of-fit tests for the Hald response are in `hald_gof.R`,
below.

## Goodness of fit

```
Rscript hald_gof.R
```

Tests the simplex distribution against the transformed Hald response, with the
mean and dispersion estimated from that response. Reports the
Kolmogorov-Smirnov, Cramer-von Mises and Anderson-Darling statistics with
asymptotic p-values, chi-square statistics on equiprobable bins, and
Kolmogorov-Smirnov and Cramer-von Mises p-values from a parametric bootstrap
that re-estimates the parameters on each of 4,000 samples.

## Data

Neither application needs a data file. The Body Fat data are loaded from the
`mfp` package; the Hald Cement data are given in the script.

## Requirements

R, with `MASS` and `parallel` for the simulation, `ggplot2`, `dplyr`, `tidyr`
and `patchwork` for the figures, `mfp` for the Body Fat application and
`goftest` for the goodness-of-fit tests. `mse_tables.R` uses base R only.
