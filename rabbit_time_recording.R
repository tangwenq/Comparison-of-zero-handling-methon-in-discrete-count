#######################################################################
#                                                                     #
#            Compositional Data Imputation — Timing Benchmark         #
#                                                                     #
#  This script compares the computation time of multiple imputation   #
#  methods for left-censored compositional data. It generates         #
#  synthetic censoring patterns, applies different replacement        #
#  algorithms, and records the elapsed time across many repetitions   #
#  using MPI-based parallelization.                                   #
#                                                                     #
#  Methods evaluated include:                                         #
#    - lmrob, PLS, multLN, multRepl, lrDA, lrEM                       #
#    - multKMSS, lrSVD, cmultRepl (GBM), and simple ad-hoc methods    #
#                                                                     #
#                                                                     #
#  Requirements:                                                      #
#    R packages: readxl, compositions, robCompositions,               #
#                zCompositions, foreach, doMPI                        #
#    MPI environment required for parallel execution.                 #
#                                                                     #
#######################################################################








# ====== Load packages ======
library(readxl)
library(compositions)
library(robCompositions)
library(zCompositions)
library(foreach)
library(doMPI)

# ====== Download data ======
url <- "https://github.com/michaelgreenacre/CODAinPractice/raw/master/Rabbits.xlsx"
temp_file <- tempfile(fileext = ".xlsx")
download.file(url, destfile = temp_file, mode = "wb")

rabbits_data <- read_excel(temp_file)[, -1]
rabbits_data <- rabbits_data[, -1]

# ====== Core functions ======
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

# ====== Simple imputation methods ======
add1_method <- function(data_0) {
  tmp <- data_0
  tmp[tmp == 0] <- 1
  tmp
}

dl_065_method <- function(data_0, DL) { 
  tmp <- data_0
  for (j in seq_along(DL)) tmp[tmp[, j] == 0, j] <- 0.65 * DL[j]
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

# ====== MPI initialization ======
cl <- startMPIcluster()
registerDoMPI(cl)

# ====== Global settings ======
batch_size  <- 50
rep_k       <- 500
m_values    <- 200
probs_seq   <- seq(0.05, 0.8, by = 0.05)
output_dir  <- "results_time_only"

dir.create(output_dir, showWarnings = FALSE)

set.seed(123)
seeds <- sample(1:1e6, rep_k)

# ====== Main loop ======
for (m_val in m_values) {
  for (prob in probs_seq) {
    cat("\n📌 Running for m =", m_val, ", prob =", prob, "\n")
    
    time_path <- file.path(output_dir, paste0("time_m", m_val, "_prob", prob, ".csv"))
    all_time  <- matrix(NA, nrow = rep_k, ncol = 12)
    
    colnames(all_time) <- c(
      "lmrob", "PLS", "mult_lognorm", "mult_repl", "lr_da", "lr_em",
      "mult_KMSS", "lr_SVD", "add1", "dl_065", "dl_unif", "GBM"
    )
    
    remaining <- 1:rep_k
    batches   <- split(remaining, ceiling(seq_along(remaining) / batch_size))
    
    for (batch_indices in batches) {
      batch_time <- foreach(
        i = batch_indices,
        .packages = c("compositions", "robCompositions", "zCompositions")
      ) %dopar% {
        
        # Generate censored data
        set.seed(seeds[i])
        res_gen <- generate_data(rabbits_data, m_val, prob, seeds[i])
        data_0  <- res_gen$data_0
        DL      <- res_gen$DL
        
        add1_data    <- add1_method(data_0)
        dl_065_data  <- dl_065_method(data_0, DL)
        dl_unif_data <- dl_unif_method(data_0, DL, seed = 123)
        
        impute_wrapper <- function(expr) {
          t0 <- system.time({ res <- try(expr, silent = TRUE) })
          t0["elapsed"]
        }
        
        elapsed_times <- numeric(12)
        names(elapsed_times) <- c(
          "lmrob", "PLS", "mult_lognorm", "mult_repl", "lr_da", "lr_em",
          "mult_KMSS", "lr_SVD", "add1", "dl_065", "dl_unif", "GBM"
        )
        
        # lmrob
        elapsed_times["lmrob"] <- impute_wrapper(
          imputeBDLs(
            data_0, dl = DL, eps = 0.1,
            method = "lmrob", verbose = FALSE,
            variation = FALSE, maxit = 50
          )$x
        )
        
        # PLS, only when m < 600
        if (m_val < 600) {
          elapsed_times["PLS"] <- impute_wrapper(
            imputeBDLs(
              data_0, dl = DL, eps = 0.1,
              method = "pls", verbose = FALSE,
              R = 50, variation = FALSE, maxit = 50
            )$x
          )
        } else {
          elapsed_times["PLS"] <- NA
        }
        
        # Other methods
        elapsed_times["mult_lognorm"] <- impute_wrapper(
          multLN(data_0, label = 0, dl = DL)
        )
        
        elapsed_times["mult_repl"] <- impute_wrapper(
          multRepl(data_0, label = 0, dl = DL)
        )
        
        elapsed_times["lr_da"] <- impute_wrapper(
          lrDA(data_0, label = 0, dl = DL, ini.cov = "multRepl")
        )
        
        elapsed_times["lr_em"] <- impute_wrapper(
          lrEM(X = data_0, label = 0, dl = DL, ini.cov = "multRepl")
        )
        
        elapsed_times["mult_KMSS"] <- impute_wrapper(
          multKM(data_0, label = 0, dl = DL)
        )
        
        elapsed_times["lr_SVD"] <- impute_wrapper(
          lrSVD(data_0, label = 0, dl = DL)
        )
        
        elapsed_times["add1"] <- impute_wrapper(add1_data)
        elapsed_times["dl_065"] <- impute_wrapper(dl_065_data)
        elapsed_times["dl_unif"] <- impute_wrapper(dl_unif_data)
        
        elapsed_times["GBM"] <- impute_wrapper(
          cmultRepl(data_0, output = "p-counts", z.delete = FALSE)
        )
        
        elapsed_times
      }
      
      # Save batch results into matrix
      for (j in seq_along(batch_indices)) {
        all_time[batch_indices[j], ] <- batch_time[[j]]
      }
      
      cat("✅ Batch time recorded.\n")
    }
    
    # Write CSV
    write.csv(all_time, file = time_path, row.names = FALSE)
  }
}

# ====== Close MPI ======
closeCluster(cl)
mpi.quit()


