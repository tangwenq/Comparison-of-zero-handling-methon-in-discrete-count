# **Comparing zero-imputation methods for high-dimensional compositional count data**

This project contains code to reproduce the results of the study on evaluating discrete and continuous zero-replacement strategies for high-dimensional compositional datasets, using both real data (Rabbits.xlsx) and synthetic simulation designs.

The goal of the project is to benchmark a broad set of imputation methods, quantify their reconstruction accuracy under different zero proportions and dimension, analyze computational complexity, and visualize the geometric effects of discretization and scaling.

---

## **Main functions**

All code is written in R and requires packages such as  
`compositions`, `robCompositions`, `zCompositions`, `ggplot2`, `patchwork`,  
`foreach`, `doMPI`, and other tidyverse tools.

The main files are:

---

- `generate_figures_coda_discrete.R`: This script generates all geometric
  illustrations related to discrete compositional data (CoDA). It compares
  the integer lattice representation of count data with the continuous
  simplex geometry and visualizes distortions introduced by sequencing
  depth, closure, scaling, and ceiling quantization. The file produces
  Figures 2–4, including lattice vs simplex comparisons, CLR shifts under
  varying depths, and log-ratio effects of discretization.

- `rabbit_comparison.R`: The main driver script for the core simulation
  comparison. It constructs censored datasets from the Rabbits data,
  applies a full suite of zero-imputation methods, and computes the CED
  and ADCS error metrics (raw and ceiling versions) across dimensions
  **m = 50, 500, 1000** and a range of missingness probabilities. Results
  are stored as error matrices and later visualized in Figures 6 and 7.

- `rabbit_pls_check.R`: This script performs a detailed local diagnostic
  study of PLS-based imputation under a single example setting
  (**m = 200**, **p = 0.2**). It extracts PLS component counts, evaluates
  imputed vs. true values for selected tracked variables, and visualizes
  the behavior of the PLS model under left-censoring. The output forms
  the basis of Figure 10.

- `rabbit_time_recording.R`: A high-performance MPI-parallelized script
  designed to benchmark computational runtime for all imputation methods.
  For each choice of dimension *m* and missingness probability *p*, the
  script generates synthetic censored datasets and records execution time
  for methods including `lmrob`, `PLS`, `lrEM`, `lrDA`, `multLN`,
  `multRepl`, `multKMSS`, `lrSVD`, `GBM`, and ad-hoc replacements. The
  results are exported as CSV files and later visualized in Figure 9.

- `rabbit_visulization.R`: This script compiles and visualizes results
  produced by the main comparison experiment. It generates Figures 6–8,
  including boxplots of method performance (Figure 6), mean error curves
  across missingness probabilities (Figure 7), and dimensional-scaling
  trends under fixed missingness (Figure 8). It also prepares publication-
  quality PDF outputs.

- `simulationdata_comparison.R`: The primary script for the appendix
  simulation based on **zero-free synthetic compositional count data**
  derived from the `microbialdata` example. It replicates the main
  comparison framework at **m = 50** components but without naturally
  occurring zeros, enabling assessment of method behavior under a
  fully observed ground truth. Boxplots of error metrics are generated
  for the appendix figures.

- `rabbit_time_plotting.R` (if separated): This script reads the runtime
  CSV files produced by `rabbit_time_recording.R` and constructs the
  combined runtime figure (Figure 9), showing mean execution times as a
  function of dimension *m* and missingness probability *p*.

- `README.md`: The main documentation file summarizing the project,
  experimental designs, methods compared, figure structure, and file
  functionality.

**Authors** 

Tang, W.



**License**
Project license (MIT).

