 Nonparameteric Bounds for Evaluating the Clinical Utility of
Treatment Rules - Code Repository
======================================================================

This repository contains the R code for the simulation studies presented in the paper: "Nonparameteric Bounds for Evaluating the Clinical Utility of
Treatment Rules".

Please note: This repository contains the code to reproduce the simulation results. The code for the real-world data example (LEAP study) discussed in the paper is not included.

-----------------------
## Overview
-----------------------

This project provides the code used to run the simulations comparing two strategies for calculating bounds on the effectiveness of clinical treatment rules:

1.  **Reduction Strategy**
2.  **Conditioning Strategy**

The primary script runs simulations across many randomly generated (but valid) probability distributions to empirically compare the performance and width of the bounds produced by these two methods.

-----------------------
## File Descriptions
-----------------------

- `which bounds are better.R`:
  This is the main simulation script. It runs the core simulations described in the paper, comparing the two bounding strategies.

- `helper_function.R`:
  A utility script that is called by the main simulation. Its function is to generate valid random conditional probability distributions for a given causal graph.

- `create_plots.R`:
  This script takes the output from the simulation runs and generates figure used in the paper.

- `README.txt`:
  This file.

-----------------------
## Requirements
-----------------------

- R
- The `causaloptim` R package.
- The parallel and pbapply, when using parallel computing.

You can install the necessary package from CRAN by running the following command in your R console:

```R
install.packages("causaloptim")
