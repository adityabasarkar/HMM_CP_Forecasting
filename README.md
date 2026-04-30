# HMM CP Forecasting

This project applies conformal prediction methods for hidden Markov models to conflict forecasting data. The maintained workflow preprocesses raw conflict data, runs simulation and fatality forecasting analyses, and writes publication-ready plots.

## Quick Start in RStudio

1. Open RStudio.
2. Select **File > Open Project...** if an `.Rproj` file is available, or **File > Open...** and choose this repository folder.
3. Set the working directory to the repository root, the folder that contains `README.md`:

   ```r
   setwd("[YOUR PATH TO REPOSITORY]")
   getwd()
   ```

4. Install the required R packages:

   ```r
   install.packages(c(
     "arrow",
     "dplyr",
     "forcats",
     "future",
     "future.apply",
     "ggplot2",
     "gridExtra",
     "Rcpp",
     "RcppArmadillo",
     "RcppDist",
     "scales",
     "stringr",
     "tibble",
     "tidyr"
   ))
   ```

5. On Windows, install Rtools if `sourceCpp()` cannot compile the C++ code. The analysis scripts compile `src/cpp_files/hmm_conf_pred_final.cpp` through `Rcpp`.
6. Add the raw input files:

   ```text
   data/raw/markov_data2507.csv
   data/raw/country_list.csv
   ```

7. Run the scripts in the order listed in `workflow.txt`.

## Project Layout

```text
data/raw/                         Raw input data, not tracked by git
outputs/r_objects/                Generated .rds analysis outputs, not tracked by git
outputs/plots/                    Generated PDF plots, not tracked by git
src/cpp_files/                    C++ conformal prediction implementation
src/entry_points/preprocessing/   Data preparation scripts
src/entry_points/simulations_and_analysis/
                                  Analysis scripts that write .rds outputs
src/entry_points/visualization/   Plot scripts that read .rds outputs
```

## Running the Workflow

Run scripts from the repository root. In RStudio, open a script and click **Source**, or run:

```r
source("src/entry_points/preprocessing/preprocess.R")
```

Start with preprocessing. It reads the raw CSV files from `data/raw/` and creates:

```text
outputs/r_objects/preprocessed_data_full.rds
outputs/r_objects/preprocessed_data.rds
```

Then run the simulation and fatality analysis scripts. These scripts may take a while because several of them repeatedly compile or call the C++ conformal prediction routines.

Finally, run the visualization scripts. They expect the corresponding `.rds` files to already exist in `outputs/r_objects/` and write PDFs to `outputs/plots/`.

## Expected Generated Files

The repository intentionally ignores large or reproducible files:

```text
data/raw/
outputs/plots/
outputs/r_objects/
```

If another user clones the project, they must provide the raw CSV files and rerun the workflow to regenerate the ignored outputs.

## Notes for Troubleshooting

- Always run scripts from the repository root. Most paths are relative to that location.
- If `sourceCpp("src/cpp_files/hmm_conf_pred_final.cpp")` fails, confirm that Rtools is installed and available on your `PATH`.
- If a visualization script fails with a missing `.rds` file, run the matching analysis script first.
- If `arrow` fails to install on Windows, restart RStudio and try installing it again from a clean R session.
- Several analysis scripts use random simulation. The maintained scripts call `set.seed()` so results should be reproducible for the same R/package environment.

## Script Naming

The main workflow lives under `src/entry_points/`. Older exploratory scripts are kept in `src/` with descriptive names so they are easier to identify, but they are not required for the standard workflow.
