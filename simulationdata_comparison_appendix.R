############################################################
## Appendix simulation: zero-free synthetic compositional
## count data constructed from the 'microbialdata' example.
## This script repeats the main imputation comparison using
## a simulated zero-free dataset. In this appendix example,
## we use m = 50 components as the feature dimension.
############################################################

## -------- Load packages --------
library(ggplot2)
library(dplyr)
library(gtools)
library(gllvm)

library(compositions)
library(robCompositions)
library(zCompositions)
library(vegan)
library(foreach)
library(doMPI)

set.seed(123)

############################################################
## 1. Construct a zero-free compositional count dataset
##    from the 'microbialdata' example (HTS-like counts)
############################################################

data("microbialdata")
Y <- microbialdata$Y
cat("Loaded Y: dimension =", paste(dim(Y), collapse = " × "), "\n")

## Sort OTUs (columns) by sparsity (number of zeros per column)
data.s <- cbind(1:ncol(Y), apply(Y == 0, 2, sum))
ind <- data.s[order(data.s[, 2], decreasing = FALSE), 1]
Y <- Y[, ind]
cat("Sorting completed: new dimension =", paste(dim(Y), collapse = " × "), "\n")
print(apply(Y == 0, 2, sum)[1:10])

## Optionally, select a subset of OTUs; here we keep all columns
Y_top <- Y
Y_template <- Y_top
n <- nrow(Y_template)
p <- ncol(Y_template)

## Draw a zero-free Dirichlet posterior composition per sample
P_post <- matrix(NA, n, p)
for (i in 1:n) {
  temp <- rdirichlet(1, Y_template[i, ] + 0.5)
  temp[temp < 1e-5] <- 1e-5
  P_post[i, ] <- temp / sum(temp)
}

cat("Any zeros in P_post:", any(P_post == 0), "\n")
cat("Row sums of P_post:\n"); print(summary(rowSums(P_post)))

## Choose a common sequencing depth and build an integer count matrix
depth_full <- 50070
Y_full <- round(P_post * depth_full)

cat("Any zeros in Y_full:", any(Y_full == 0), "\n")
cat("Constructed zero-free integer count matrix:",
    paste(dim(Y_full), collapse = " × "), "\n")

Y_nozero <- as.data.frame(Y_full)
data_used <- Y_nozero

############################################################
## 2. Core functions: data generation and performance metrics
############################################################

## Generate data with artificial detection limits and zeros
generate_data <- function(data, m, probs, i) {
  set.seed(i)
  selected_cols <- sample(ncol(data), m)
  data_m <- data[, selected_cols]
  data_0 <- data_m
  DL <- numeric(m)
  for (j in seq(1, m, 2)) {
    DL[j] <- quantile(data_m[, j], probs = probs, na.rm = TRUE)
  }
  for (j in 1:m) {
    if (DL[j] > 0) {
      data_0[data_0[, j] < DL[j], j] <- 0
    }
  }
  return(list(data_m = data_m, data_0 = data_0, DL = DL))
}

## Average deviation in covariance structure (ADCS) in ilr-space
ADCS <- function(original_data, imputed_data) {
  Z      <- ilr(original_data)
  Z_star <- ilr(imputed_data)
  S      <- cov(Z)
  S_star <- cov(Z_star)
  frob_norm <- norm(S - S_star, type = "F")
  frob_norm / (ncol(original_data) - 1)
}

## Simple +1 replacement baseline
add1_method <- function(data_0) {
  tmp <- data_0
  tmp[tmp == 0] <- 1
  tmp
}

## Deterministic DL-based replacement: 0.65 * DL
dl_065_method <- function(data_0, DL) {
  tmp <- data_0
  for (j in seq_along(DL)) {
    tmp[tmp[, j] == 0, j] <- 0.65 * DL[j]
  }
  tmp
}

## Uniform random replacement on (0.1 * DL, DL)
dl_unif_method <- function(data_0, DL, seed) {
  tmp <- data_0
  set.seed(seed)
  for (j in seq_along(DL)) {
    idx <- which(tmp[, j] == 0)
    if (length(idx) > 0) {
      tmp[idx, j] <- runif(length(idx), 0.1 * DL[j], DL[j])
    }
  }
  tmp
}

## Wrapper: run all imputation methods and compute metrics
imp_comparison <- function(data, m, probs, seed) {
  set.seed(seed)
  result <- generate_data(data = data, m = m, probs = probs, i = seed)
  data_m <- result$data_m
  data_0 <- result$data_0
  DL     <- result$DL
  
  method_names <- c(
    "lmrob", "PLS", "mult_lognorm", "mult_repl",
    "lr_da", "lr_em", "mult_KMSS", "lr_SVD",
    "add1", "dl_065", "dl_unif", "GBM"
  )
  
  ADCS_      <- matrix(NA, nrow = 1, ncol = length(method_names))
  CED_       <- matrix(NA, nrow = 1, ncol = length(method_names))
  ADCS_raw_  <- matrix(NA, nrow = 1, ncol = length(method_names))
  CED_raw_   <- matrix(NA, nrow = 1, ncol = length(method_names))
  colnames(ADCS_) <- colnames(CED_) <-
    colnames(ADCS_raw_) <- colnames(CED_raw_) <- method_names
  
  ## Generic imputation + metric computation wrapper
  impute_wrapper <- function(expr) {
    res <- try(expr, silent = TRUE)
    if (inherits(res, "try-error")) return(rep(NA, 4))
    res <- as.matrix(res)
    if (any(!is.finite(res))) return(rep(NA, 4))
    
    ## Raw evaluation (continuous imputed values)
    res[res <= 0] <- 1
    ced_raw  <- try(ced(data_m, res, ni = sum(data_0 == 0)), silent = TRUE)
    adcs_raw <- try(ADCS(data_m, res), silent = TRUE)
    if (inherits(ced_raw, "try-error") || inherits(adcs_raw, "try-error"))
      return(rep(NA, 4))
    
    ## Ceiled evaluation (counts quantized to integers at zero positions)
    res_corrected <- res
    res_corrected[data_0 == 0] <- ceiling(res_corrected[data_0 == 0])
    ced_ceil  <- try(ced(data_m, res_corrected, ni = sum(data_0 == 0)), silent = TRUE)
    adcs_ceil <- try(ADCS(data_m, res_corrected), silent = TRUE)
    if (inherits(ced_ceil, "try-error") || inherits(adcs_ceil, "try-error"))
      return(rep(NA, 4))
    
    c(ced_ceil, adcs_ceil, ced_raw, adcs_raw)
  }
  
  add1_data   <- add1_method(data_0)
  dl_065_data <- dl_065_method(data_0, DL)
  dl_unif_data <- dl_unif_method(data_0, DL, seed = 123)
  
  CED_ADCS <- rbind(
    impute_wrapper(imputeBDLs(data_0, dl = DL, eps = 0.1, method = "lmrob",
                              verbose = FALSE, variation = FALSE, maxit = 50)$x),
    impute_wrapper(imputeBDLs(data_0, dl = DL, eps = 0.1, method = "pls",
                              verbose = FALSE, R = 50, variation = FALSE, maxit = 50)$x),
    impute_wrapper(multLN(data_0, label = 0, dl = DL)),
    impute_wrapper(multRepl(data_0, label = 0, dl = DL)),
    impute_wrapper(lrDA(data_0, label = 0, dl = DL, ini.cov = "multRepl")),
    impute_wrapper(imputeBDLs(data_0, dl = DL, eps = 0.1, method = "lm",
                              verbose = FALSE, variation = FALSE, R = 5, maxit = 50)$x),
    impute_wrapper(multKM(data_0, label = 0, dl = DL)),
    impute_wrapper(lrSVD(data_0, label = 0, dl = DL)),
    impute_wrapper(add1_data),
    impute_wrapper(dl_065_data),
    impute_wrapper(dl_unif_data),
    impute_wrapper(cmultRepl(data_0, output = "p-counts", z.delete = FALSE))
  )
  
  CED_[1, ]      <- CED_ADCS[, 1]
  ADCS_[1, ]     <- CED_ADCS[, 2]
  CED_raw_[1, ]  <- CED_ADCS[, 3]
  ADCS_raw_[1, ] <- CED_ADCS[, 4]
  
  list(
    CED_ceil = CED_,
    ADCS_ceil = ADCS_,
    CED_raw  = CED_raw_,
    ADCS_raw = ADCS_raw_
  )
}

############################################################
## 3. Parallel simulation setup (appendix setting: m = 50)
############################################################

cl <- startMPIcluster()
registerDoMPI(cl)

batch_size  <- 40
rep_k       <- 500
probs_seq   <- seq(0.05, 0.8, by = 0.05)
output_dir  <- "results_Ynozero"
dir.create(output_dir, showWarnings = FALSE)

set.seed(123)
seeds <- sample(1:1e6, rep_k)
m_compare_results <- list()

############################################################
## 4. Main simulation loop (appendix example m = 50)
############################################################

for (m_val in c(50)) {
  for (prob in probs_seq) {
    
    cat("\nRunning for m =", m_val, ", prob =", prob, "\n")
    result_prefix <- paste0("m", m_val, "_prob", prob)
    output_path   <- file.path(output_dir, paste0("result_", result_prefix, ".rds"))
    
    if (file.exists(output_path)) {
      cat("Loading existing partial results from:", output_path, "\n")
      result_list <- readRDS(output_path)
      done_indices <- which(sapply(result_list, Negate(is.null)))
    } else {
      result_list  <- vector("list", rep_k)
      done_indices <- integer(0)
    }
    
    remaining <- setdiff(1:rep_k, done_indices)
    cat("Remaining simulations:", length(remaining), "\n")
    
    batches <- split(remaining, ceiling(seq_along(remaining) / batch_size))
    
    for (batch_indices in batches) {
      new_results <- foreach(
        i = batch_indices,
        .packages = c("compositions", "robCompositions", "zCompositions")
      ) %dopar% {
        res <- try(
          imp_comparison(data = data_used, m = m_val, probs = prob, seed = seeds[i]),
          silent = TRUE
        )
        if (inherits(res, "try-error")) return(NULL)
        res
      }
      
      for (j in seq_along(batch_indices)) {
        idx <- batch_indices[j]
        result_list[[idx]] <- new_results[[j]]
        cat(sprintf("Simulation %d done.\n", idx))
      }
      
      saveRDS(result_list, output_path)
      cat("Partial results saved after batch.\n")
    }
    
    valid_results <- Filter(Negate(is.null), result_list)
    
    if (length(valid_results) > 0) {
      m_compare_results[[paste0("CEDceil_", result_prefix)]] <-
        do.call(rbind, lapply(valid_results, `[[`, "CED_ceil"))
      m_compare_results[[paste0("ADCSceil_", result_prefix)]] <-
        do.call(rbind, lapply(valid_results, `[[`, "ADCS_ceil"))
      m_compare_results[[paste0("CEDraw_", result_prefix)]]  <-
        do.call(rbind, lapply(valid_results, `[[`, "CED_raw"))
      m_compare_results[[paste0("ADCSraw_", result_prefix)]] <-
        do.call(rbind, lapply(valid_results, `[[`, "ADCS_raw"))
    } else {
      cat("No valid results for m =", m_val, "prob =", prob, "\n")
    }
  }
}

############################################################
## 5. Save results and shut down MPI
############################################################

save(m_compare_results, seeds, file = "results_Ynozero.RData")

closeCluster(cl)
mpi.quit()







