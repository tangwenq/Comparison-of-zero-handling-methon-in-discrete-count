
###############################
## Main comparison script
## Discrete CoDA zero-imputation methods
## Example setting: m = 500 variables
###############################

# ====== Load packages ======
library(readxl)
library(ggplot2)
library(compositions)
library(robCompositions)
library(zCompositions)
library(vegan)
library(foreach)
library(doMPI)

# ====== Download rabbit data (Greenacre, 2021) ======
url <- "https://github.com/michaelgreenacre/CODAinPractice/raw/master/Rabbits.xlsx"
temp_file <- tempfile(fileext = ".xlsx")
download.file(url, destfile = temp_file, mode = "wb")

# First column is an index, second is group label → remove both
rabbits_data <- read_excel(temp_file)[, -1]
rabbits_data <- rabbits_data[, -1]

# ====== Core helper functions ======

# Subsample m columns, generate detection limits and apply left-censoring to create zeros
generate_data <- function(data, m, probs, i) {
  set.seed(i)
  selected_cols <- sample(ncol(data), m)
  data_m <- data[, selected_cols]
  data_0 <- data_m
  DL <- numeric(m)
  
  # Detection limit only for odd-indexed columns
  for (j in seq(1, m, 2)) {
    DL[j] <- quantile(data_m[, j], probs = probs, na.rm = TRUE)
  }
  
  # Set values below DL to zero
  for (j in seq_len(m)) {
    if (DL[j] > 0) {
      data_0[data_0[, j] < DL[j], j] <- 0
    }
  }
  list(data_m = data_m, data_0 = data_0, DL = DL)
}

# Average distortion in covariance structure (Aitchison)
ADCS <- function(original_data, imputed_data) {
  Z <- ilr(original_data)
  Z_star <- ilr(imputed_data)
  S <- cov(Z)
  S_star <- cov(Z_star)
  frob_norm <- norm(S - S_star, type = "F")
  frob_norm / (ncol(original_data) - 1)
}

# Simple ad-hoc zero-handling methods
add1_method <- function(data_0) {
  tmp <- data_0
  tmp[tmp == 0] <- 1
  tmp
}

dl_065_method <- function(data_0, DL) {
  tmp <- data_0
  for (j in seq_along(DL)) {
    tmp[tmp[, j] == 0, j] <- 0.65 * DL[j]
  }
  tmp
}

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

# Wrapper to compare multiple imputation methods on one simulated dataset
imp_comparison <- function(data, m, probs, seed) {
  set.seed(seed)
  
  # Generate censored data with zeros
  result <- generate_data(data = data, m = m, probs = probs, i = seed)
  data_m <- result$data_m     # original (no zeros)
  data_0 <- result$data_0     # censored with zeros
  DL <- result$DL             # detection limits
  
  method_names <- c(
    "lmrob", "PLS", "mult_lognorm", "mult_repl",
    "lr_da", "lr_em", "mult_KMSS", "lr_SVD",
    "add1", "dl_065", "dl_unif", "GBM"
  )
  
  ADCS_ <- matrix(NA, nrow = 1, ncol = length(method_names))
  CED_ <- matrix(NA, nrow = 1, ncol = length(method_names))
  ADCS_raw_ <- matrix(NA, nrow = 1, ncol = length(method_names))
  CED_raw_ <- matrix(NA, nrow = 1, ncol = length(method_names))
  colnames(ADCS_) <- colnames(CED_) <- colnames(ADCS_raw_) <- colnames(CED_raw_) <- method_names
  
  # Evaluate one imputed matrix: return CED/ADCS for raw + ceiling-corrected versions
  impute_wrapper <- function(expr) {
    res <- try(expr, silent = TRUE)
    if (inherits(res, "try-error")) return(rep(NA, 4))
    
    res <- as.matrix(res)
    if (any(!is.finite(res))) return(rep(NA, 4))
    
    # Ensure strictly positive entries
    res[res <= 0] <- 1
    
    ced_raw <- try(ced(data_m, res, ni = sum(data_0 == 0)), silent = TRUE)
    adcs_raw <- try(ADCS(data_m, res), silent = TRUE)
    if (inherits(ced_raw, "try-error") || inherits(adcs_raw, "try-error")) return(rep(NA, 4))
    
    # Ceiling to integer counts at zero positions
    res_corrected <- res
    res_corrected[data_0 == 0] <- ceiling(res_corrected[data_0 == 0])
    
    ced_ceil <- try(ced(data_m, res_corrected, ni = sum(data_0 == 0)), silent = TRUE)
    adcs_ceil <- try(ADCS(data_m, res_corrected), silent = TRUE)
    if (inherits(ced_ceil, "try-error") || inherits(adcs_ceil, "try-error")) return(rep(NA, 4))
    
    c(ced_ceil, adcs_ceil, ced_raw, adcs_raw)
  }
  
  # Construct ad-hoc imputations
  add1_data   <- add1_method(data_0)
  dl_065_data <- dl_065_method(data_0, DL)
  dl_unif_data <- dl_unif_method(data_0, DL, seed = 123)
  
  # Run all methods
  CED_ADCS <- rbind(
    impute_wrapper(imputeBDLs(data_0, dl = DL, eps = 0.1,
                              method = "lmrob", variation = FALSE, maxit = 50)$x),
    impute_wrapper(imputeBDLs(data_0, dl = DL, eps = 0.1,
                              method = "pls", variation = FALSE, R = 50, maxit = 50)$x),
    impute_wrapper(multLN(data_0, label = 0, dl = DL)),
    impute_wrapper(multRepl(data_0, label = 0, dl = DL)),
    impute_wrapper(lrDA(data_0, label = 0, dl = DL, ini.cov = "multRepl")),
    impute_wrapper(lrEM(X = data_0, label = 0, dl = DL, ini.cov = "multRepl")),
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
    CED_raw = CED_raw_,
    ADCS_raw = ADCS_raw_
  )
}

# ====== MPI cluster initialization ======
cl <- startMPIcluster()
registerDoMPI(cl)

# ====== Global settings (main comparison; here m = 500 as example) ======
batch_size  <- 40
rep_k       <- 500
probs_seq   <- seq(0.05, 0.8, by = 0.05)
output_dir  <- "results_500compare"
dir.create(output_dir, showWarnings = FALSE)

set.seed(123)
seeds <- sample(1:1e6, rep_k)
m_compare_results <- list()

# ====== Main simulation loop over m and censoring quantiles ======
for (m_val in c(500)) {                 # main example: m = 500
  for (prob in probs_seq) {
    
    cat("\nRunning: m =", m_val, ", prob =", prob, "\n")
    result_prefix <- paste0("m", m_val, "_prob", prob)
    output_path <- file.path(output_dir, paste0("result_", result_prefix, ".rds"))
    
    # Load partial results if available
    if (file.exists(output_path)) {
      result_list  <- readRDS(output_path)
      done_indices <- which(sapply(result_list, Negate(is.null)))
    } else {
      result_list  <- vector("list", rep_k)
      done_indices <- integer(0)
    }
    
    remaining <- setdiff(seq_len(rep_k), done_indices)
    cat("Remaining simulations:", length(remaining), "\n")
    
    batches <- split(remaining, ceiling(seq_along(remaining) / batch_size))
    
    # Parallel evaluation using MPI
    for (batch_indices in batches) {
      new_results <- foreach(
        i = batch_indices,
        .packages = c("compositions", "robCompositions", "zCompositions")
      ) %dopar% {
        res <- try(
          imp_comparison(data = rabbits_data, m = m_val,
                         probs = prob, seed = seeds[i]),
          silent = TRUE
        )
        if (inherits(res, "try-error")) return(NULL)
        res
      }
      
      for (j in seq_along(batch_indices)) {
        idx <- batch_indices[j]
        result_list[[idx]] <- new_results[[j]]
        cat(sprintf("Simulation %d completed.\n", idx))
      }
      
      saveRDS(result_list, output_path)
      cat("Partial results saved.\n")
    }
    
    # Collect valid results across replications
    valid_results <- Filter(Negate(is.null), result_list)
    
    if (length(valid_results) > 0) {
      m_compare_results[[paste0("CEDceil_",  result_prefix)]] <-
        do.call(rbind, lapply(valid_results, `[[`, "CED_ceil"))
      m_compare_results[[paste0("ADCSceil_", result_prefix)]] <-
        do.call(rbind, lapply(valid_results, `[[`, "ADCS_ceil"))
      m_compare_results[[paste0("CEDraw_",   result_prefix)]]  <-
        do.call(rbind, lapply(valid_results, `[[`, "CED_raw"))
      m_compare_results[[paste0("ADCSraw_",  result_prefix)]] <-
        do.call(rbind, lapply(valid_results, `[[`, "ADCS_raw"))
    } else {
      cat("No valid results for m =", m_val, "prob =", prob, "\n")
    }
  }
}

# ====== Save final aggregated results ======
save(m_compare_results, seeds, file = "results500_rabbit.RData")

# ====== Shutdown MPI ======
closeCluster(cl)
mpi.quit()
