# Simulation Workflow

This repository contains scripts to reproduce the simulation results and generate the final plots.

## Requirements
- **R** (version ≥ 4.0)
- Required R packages (install if not already available):
  ```r
  install.packages(c("glmnet", "CVXR", "stringr", "dplyr", "glmtrans", "ggplot2", "ggpubr", "RColorBrewer"))
  ```

## Workflow

1. **Generate simulation results for Settings 1–3**
   ```bash
   Rscript main_1.R
   ```

2. **Generate simulation results for Settings 4.1 and 4.2**
   ```bash
   Rscript main_2.R
   ```

3. **Aggregate results and produce final plot**
   ```bash
   Rscript summ.R
   ```

## Notes
- `main_1.R` will create and save the simulation results for **Settings 1–3**.
- `main_2.R` will create and save the simulation results for **Settings 4.1 and 4.2**.
- `summ.R` will read the output files from the above steps, aggregate the results, and generate the final plots.
