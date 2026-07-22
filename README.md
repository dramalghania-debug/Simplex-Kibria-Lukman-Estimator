# Simplex Kibria-Lukman Estimator (SKLE)

This repository contains the supplementary code and datasets for the manuscript: 
**"Kibria-Lukman Estimator for Simplex Regression under Multicollinearity: Simulation and Applications"**

## Files Included:
* **`Simulation.R`**: The primary R script used to execute the Monte Carlo simulations, calculate the $k_1$ and $k_2$ shrinkage parameters, and formulate the SKLE matrix.
* **`Body fat new.R`, `CI hald data.R` **: The real-world dataset utilized in Section 5 of the manuscript.
  
* ## Reproducing Other Dimensions (p = 4, 12, 16)

`Simulation.R` is parameterized, not hardcoded — the number of explanatory 
variables is set by the `P_VAL` variable at the top of the script (line 8). 
To reproduce the other dimensionalities reported in the manuscript, change 
this value and rerun the script:

```r
P_VAL <- 8   # Set to 4, 8, 12, or 16 to reproduce the corresponding scenario
```

All other simulation factors (sample size, multicollinearity level, 
dispersion, and link function) are controlled by the same configuration 
block and do not require modification.

## Software Requirements:
The simulations were conducted using the R programming environment. 
