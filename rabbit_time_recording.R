###########################################################################
#                                                                         #
#   Figure 9 – Compositional Data Imputation: Timing Benchmark            #
#                                                                         #
#   Part 1: Simulation                                                    #
#   ------------------------------------------------------------------    #
#   - Compares computation time of multiple imputation methods for        #
#     left-censored compositional data.                                   #
#   - Generates synthetic censoring patterns, applies different           #
#     replacement algorithms, and records elapsed time over many          #
#     repetitions using MPI-based parallelization.                        #
#   - Output: CSV files in a relative folder, e.g. "results_time_only",   #
#     named as "time_m<m>_prob<p>.csv".                                   #
#                                                                         #
#   Part 2: Visualisation                                                 #
#   ------------------------------------------------------------------    #
#   - Uses the saved CSV files to produce Figure 9 with two panels:       #
#       (a) Runtime vs dimension m at fixed p = 0.5.                      #
#       (b) Runtime vs missingness probability p at fixed m = 200.        #
#   - Plots mean runtime (seconds, log10 scale) for selected methods.     #
#                                                                         #
###########################################################################

############################
# Part 1 – Simulation      #
############################

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
  
  list(data_m = data_m, data_0 = data_0, DL = DL)
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
m_values    <- 200                      # can be a vector if needed
probs_seq   <- seq(0.05, 0.8, by = 0.05)
output_dir  <- "results_time_only"

dir.create(output_dir, showWarnings = FALSE)

set.seed(123)
seeds <- sample(1:1e6, rep_k)

# ====== Main loop ======
for (m_val in m_values) {
  for (prob in probs_seq) {
    cat("\nRunning for m =", m_val, ", prob =", prob, "\n")
    
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
        
        # PLS only when m < 600
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
        
        elapsed_times["add1"]    <- impute_wrapper(add1_data)
        elapsed_times["dl_065"]  <- impute_wrapper(dl_065_data)
        elapsed_times["dl_unif"] <- impute_wrapper(dl_unif_data)
        
        elapsed_times["GBM"] <- impute_wrapper(
          cmultRepl(data_0, output = "p-counts", z.delete = FALSE)
        )
        
        elapsed_times
      }
      
      for (j in seq_along(batch_indices)) {
        all_time[batch_indices[j], ] <- batch_time[[j]]
      }
      
      cat("Batch timing recorded.\n")
    }
    
    write.csv(all_time, file = time_path, row.names = FALSE)
  }
}

# ====== Close MPI ======
closeCluster(cl)
mpi.quit()

############################
# Part 2 – Figure 9 plots  #
############################

library(ggplot2)
library(dplyr)
library(readr)
library(patchwork)
library(stringr)

# ====== Methods and aesthetics ======
methods <- c(
  "lmrob","PLS","mult_lognorm","mult_repl",
  "lr_da","lr_em","mult_KMSS","lr_SVD",
  "add1","dl_065","dl_unif","GBM"
)

plot_methods <- setdiff(
  methods,
  c("lmrob", "dl_065", "dl_unif", "add1", "lr_em", "lr_da")
)

color_values <- c(
  "PLS"          = "#E41A1C",
  "mult_lognorm" = "#377EB8",
  "mult_repl"    = "#4DAF4A",
  "mult_KMSS"    = "#FB9A99",
  "lr_SVD"       = "#A65628",
  "GBM"          = "#66C2A5"
)

shape_values <- c(
  "PLS"          = 16,
  "mult_lognorm" = 17,
  "mult_repl"    = 15,
  "mult_KMSS"    = 8,
  "lr_SVD"       = 18,
  "GBM"          = 2
)

base_theme <- theme_minimal(base_size = 14) +
  theme(
    legend.title   = element_blank(),
    legend.position = "top",
    legend.text    = element_text(size = 13),
    axis.text.x    = element_text(angle = 45, hjust = 1, size = 13),
    axis.text.y    = element_text(size = 13),
    axis.title     = element_text(size = 14),
    strip.text     = element_text(size = 12, face = "bold")
  )

# =========================================================
# (a) Runtime vs dimension m at fixed probability p = 0.5
# =========================================================

m_values_plot <- c(50, 100, 200, 400, 500, 600, 700, 800, 1000)
prob_fixed    <- 0.5
all_time_df   <- data.frame()
time_dir      <- "results_time_only"

for (m in m_values_plot) {
  file_path <- file.path(time_dir, sprintf("time_m%d_prob%.1f.csv", m, prob_fixed))
  if (file.exists(file_path)) {
    df <- read_csv(file_path, show_col_types = FALSE)
    df$m <- m
    all_time_df <- bind_rows(all_time_df, df)
  } else {
    message("File not found: ", file_path)
  }
}

mean_time_df <- all_time_df %>%
  group_by(m) %>%
  summarise(across(all_of(methods), ~ mean(.x, na.rm = TRUE))) %>%
  tidyr::pivot_longer(-m, names_to = "method", values_to = "mean_time") %>%
  filter(method %in% plot_methods) %>%
  mutate(
    mean_time = ifelse(method == "PLS" & m >= 600, NA, mean_time),
    method    = factor(method, levels = names(color_values))
  )

p1 <- ggplot(mean_time_df,
             aes(x = m, y = mean_time,
                 color = method, shape = method, group = method)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  scale_y_log10() +
  scale_color_manual(values = color_values) +
  scale_shape_manual(values = shape_values) +
  labs(x = "Dimension (m)", y = "Mean Runtime (seconds, log10)") +
  base_theme

# =========================================================
# (b) Runtime vs missingness probability at fixed m = 200
# =========================================================

m_val_plot <- 200
probs_plot <- seq(0.05, 0.8, by = 0.05)
time_prob_dir <- "results_time_only"   # same folder as in simulation

time_results <- lapply(probs_plot, function(p) {
  file <- file.path(time_prob_dir, sprintf("time_m%d_prob%g.csv", m_val_plot, p))
  if (!file.exists(file)) {
    warning(paste("Missing file:", file))
    return(NULL)
  }
  dat <- read_csv(file, show_col_types = FALSE)
  colnames(dat) <- methods
  mean_times <- colMeans(dat, na.rm = TRUE)
  data.frame(prob = p, method = names(mean_times), mean_time = mean_times)
})

time_df <- do.call(rbind, time_results) %>%
  filter(method %in% plot_methods) %>%
  mutate(
    method    = factor(method, levels = names(color_values)),
    mean_time = ifelse(method == "PLS" & prob >= 0.6, NA, mean_time)
  )

p2 <- ggplot(time_df,
             aes(x = prob, y = mean_time,
                 color = method, shape = method, group = method)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  scale_y_log10() +
  scale_color_manual(values = color_values) +
  scale_shape_manual(values = shape_values) +
  labs(x = "Missingness Probability", 
       y = "Mean Runtime (seconds, log10)") +
  base_theme

# =========================================================
# Combined Figure 9 layout and export
# =========================================================

p_combined <- p1 / p2 +
  plot_layout(guides = "collect") &
  theme(legend.position = "top")

pdf(file = "Figure9_runtime_combined.pdf",
    width = 9, height = 12, useDingbats = FALSE)
print(p_combined)
dev.off()

p_combined

