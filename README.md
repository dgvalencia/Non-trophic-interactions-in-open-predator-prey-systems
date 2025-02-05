# Non-Trophic Interactions in Open Predator-Prey Systems

## Overview
This repository contains code and data associated with the study:

**"The complexity of a ‘simple' predator-prey system: Nontrophic positive interactions generate unsuspected dynamics and dependencies"**

**Authors:** Daniel E. Valencia, Alexandre Génin, Sergio Rojas, and Sergio A. Navarrete

This study explores the effects of non-trophic interactions (NTIs) on predator-prey dynamics in an open system. The repository provides the simulation scripts, dataset, and code used to generate the figures presented in the manuscript.

---

## Repository Contents

### 1. `numerical_simulation.R`
This R script runs numerical simulations of the predator-prey system under different scenarios:
- **Only trophic interactions**
- **Recruitment facilitation** (prey facilitating predator recruitment)
- **Refuge provision** (prey reducing predator mortality)
- **Both NTIs combined**

The script solves the system of ordinary differential equations (ODEs) using the `deSolve` package and outputs the temporal trajectories and steady-state abundances of prey and predators.

### 2. `last_time.csv`
A dataset containing the final equilibrium values of predator and prey abundances from the simulation. This dataset is a smaller version of the complete simulation results and is used to generate the figures.

**Columns include:**
- `Prey` - Abundance of prey
- `Predator` - Abundance of predators
- `interaction_type` - Type of interaction (e.g., "Only trophic", "Recruitment facilitation", "Refuge provision", "Both NTI")
- Other relevant simulation parameters

### 3. `results_figures.R`
This R script generates the main figures from the study. It reads the `last_time.csv` dataset, processes the data, and produces visualizations for:
- Predator recruitment rate vs. prey recruitment
- Predator mortality rate vs. prey recruitment
- Predator abundance vs. prey recruitment
- Predator-prey abundance relationships under different NTI scenarios

The script uses `ggplot2` for data visualization and saves the figures as `.png` files.

---

## Requirements
The R scripts require the following R packages:

- `deSolve`
- `tidyverse`
- `plyr`
- `ggplot2`

To install them, run:
```r
install.packages(c("deSolve", "tidyverse", "plyr", "ggplot2"))
```
---

## Running the Simulation

To run the numerical simulation:

 1. Open`numerical_simulation.R` in RStudio or a terminal.
 2. Set up the simulation parameters if needed.
 3. Run the script.

To generate the figures:

 1. Ensure `last_time.csv` is present in the working directory.
 2. Run `results_figures.R`.
 3. Figures will be saved as `.png` files in the working directory. Note that the file extension could be easily changed by replacing `.png` in the `ggsave()` function.

---

## Citation

If you use this code or data in your work, please cite our study:

---
## Contact

For questions or issues, please contact Daniel E. Valencia at `devalencia12@gmail.com`.
