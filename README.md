# **Comparing zero-imputation methods for high-dimensional compositional count data**

This project contains code to reproduce the results of the study on evaluating discrete and continuous zero-replacement strategies for high-dimensional compositional datasets, using both real data (Rabbits.xlsx) and synthetic simulation designs.

The goal of the project is to benchmark a broad set of imputation methods, quantify their reconstruction accuracy under different censoring regimes, analyze computational complexity, and visualize the geometric effects of discretization and scaling.

---

## **Main functions**

All code is written in R and requires packages such as  
`compositions`, `robCompositions`, `zCompositions`, `ggplot2`, `patchwork`,  
`foreach`, `doMPI`, and other tidyverse tools.

The main files are:

---

### **`generate_figures_coda_discrete.R`**  
This file generates all figures illustrating the geometric distortions caused by discretization of compositional data.

Figures produced include:  
- *Figure 2*: Integer lattice **N²** vs continuous 3-part simplex  
- *Figure 3(a)*: Simplex closure at different sequencing depths  
- *Figure 3(b)*: CLR coordinates under varying sequencing depths  
- *Figure 4(a)*: Scaling + ceiling quantization (scatter)  
- *Figure 4(b)*: Scaling + quantization effects on log-ratio shifts  

---

### **`rabbit_comparison.R`**  
Main simulation engine for the core comparison experiment.  
It generates CED and ADCS error matrices (raw and ceiling versions) across:  

- Feature dimensions: **M = 50, 500, 1000**  
- Missingness probabilities: **p ∈ [0.05, 0.80]**

These outputs are used for:  
- **Figure 6** – Boxplots of method performance  
- **Figure 7** – Mean error trends across p  

---

### **`rabbit_pls_check.R`**  
A detailed local inspection of **PLS-based zero-imputation**.

- Uses one replicate at **m = 200**, **p = 0.2**  
- Tracks how PLS chooses the number of latent components  
- Compares imputed vs. true values for selected variables  
Produces:  
- **Figure 10** (PLS diagnostic visualizations)

---

### **`rabbit_time_recording.R`**  
Parallel MPI-based runtime benchmark for all imputation methods.

- Generates censored datasets from Rabbits.xlsx  
- Applies: `lmrob`, `PLS`, `multLN`, `multRepl`, `lrDA`, `lrEM`,  
  `multKMSS`, `lrSVD`, `GBM`, and ad-hoc rules  
- Records elapsed times across **500 repetitions**  
- Saves timing results as CSV files for visualization

These results are used to build:  
- **Figure 9** (runtime vs dimension & runtime vs missingness)

---

### **`rabbit_visulization.R`**  
A complete plotting script for the main comparison experiment.

Produces:  
- **Figure 6** – Boxplots of CED/ADCS errors  
- **Figure 7** – Mean error line plots (vs missingness)  
- **Figure 8** – Mean error vs dimension m under fixed p  

---

### **`simulationdata_comparison.R`**  
Appendix simulation using **zero-free synthetic compositional data** constructed from the `microbialdata` example.

- Repeats the main comparison under **m = 50** components  
- Highlights method behavior when no natural zeros exist  
Provides the appendix figures.

---

## **Authors**
Tang, W.

---

## **License**
LGPL ≥ 2.1
