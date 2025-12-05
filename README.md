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
  illustrations related to discrete data and compositional data analysis. It compares
  the integer lattice representation of count data with the continuous
  simplex geometry and visualizes distortions introduced by sequencing
  depth, closure, scaling, and ceiling quantization. The file produces
  Figures 2–4, including lattice vs simplex comparisons, CLR shifts under
  varying depths, and log-ratio effects of discretization.

- `rabbit_comparison.R`:  
  Main script conducting the full zero-imputation comparison for the
  Rabbits dataset. It constructs left-censored compositional datasets,
  applies all imputation methods, and computes both CED and ADCS error
  metrics (raw and ceiling variants) across dimensions **m = 50, 500,
  1000** and a range of missingness probabilities.  
  The script also compiles and visualizes the resulting error matrices,
  producing Figures 6–8:  
  • Figure 6 — boxplots of method performance under varying p  
  • Figure 7 — mean error curves across missingness probabilities  
  • Figure 8 — scaling of average error as the dimension m increases
    under fixed missingness levels  
  Publication-ready PDF figures are automatically generated.


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

- `simulationdata_comparison_appendix.R`: The primary script for the appendix
  simulation based on **zero-free synthetic compositional count data**
  derived from the `microbialdata` example. It replicates the main
  comparison framework at **m = 50** components but without naturally
  occurring zeros, enabling assessment of method behavior under a
  fully observed ground truth. Boxplots of error metrics are generated
  for the appendix figures.

- `rabbit_visulization.R`:  
  Visualization script for the Rabbits count dataset. It produces
  Figure 5, including:  
  • Figure 5(a) — a heatmap of the Rabbits count table after sorting
    taxa by mean abundance  
  • Figure 5(b) — a log–log mean–variance plot illustrating
    sparsity and overdispersion in the raw count data  
  These plots provide an overview of structural characteristics of the
  dataset before imputation. Publication-ready PDF files are generated.

- `README.md`: The main documentation file summarizing the project,
  experimental designs, methods compared, figure structure, and file
  functionality.

**Authors** 

Tang, W.



**License**
Project license (MIT).

