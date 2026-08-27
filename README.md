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
| `verify_table18.R` | Recomputes every number in Table 18 from the CSV files, printing the column arithmetic used for each. |
| `Simulation_Results_p4.csv` | Simulation output, p = 4. |
| `Simulation_Results_p8.csv` | Simulation output, p = 8. |
| `Simulation_Results_p12.csv` | Simulation output, p = 12. |
| `Simulation_Results_p16.csv` | Simulation output, p = 16. |

All four scripts read and write in the working directory, so they run from a
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
combinations of the design factors that the manuscript reports in tables.
Set `SKLE_DATA` to read the CSV files from elsewhere and `SKLE_FIGS` to write
the figures elsewhere:

```
SKLE_DATA=/path/to/csvs SKLE_FIGS=/path/to/figures Rscript figures.R
```

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

## Verifying Table 18

```
Rscript verify_table18.R
```

Recomputes every number in Table 18 from the four CSV files and prints the
column arithmetic behind each one. Base R only.

## Data

Neither application needs a data file. The Body Fat data are loaded from the
`mfp` package; the Hald Cement data are given in the script.

## Requirements

R, and the packages named in the `library()` calls at the head of each script.
Run under R 4.5.1.
