
###########################################################################
## Main comparison script and plots
## Discrete CoDA zero-imputation methods
## Example base setting: m = 500 variables
###########################################################################
#                                                                         #
#   This script contains the main simulation comparison and three figure  #
#   generations based on two experimental setups.                         #
#                                                                         #
#   Figures 6 & 7 – Fixed m, increasing p                                 #
#   -------------------------------------------------------------------   #
#   • Setting: fixed dimension m (e.g. M = 50, 500, 1000),                #
#              varying proportion of missingness p.                       #
#   • Figure 6: boxplots of CED / ADCS errors (raw and ceiling versions)  #
#               across methods, probabilities, and M.                     #
#   • Figure 7: mean error line plots (CED / ADCS, raw and ceiling)       #
#               across methods, probabilities, and M.                     #
#   • Input: RData files (m_compare_results) for M = 50, 500, 1000.       #
#                                                                         #
#   Figure 8 – Fixed p, increasing m                                      #
#   -------------------------------------------------------------------   #
#   • Setting: fixed proportions of missingness p = 0.2 and p = 0.5,      #
#              varying dimension m.                                       #
#   • Figure 8: mean error line plots (CED / ADCS, raw and ceiling)       #
#               as a function of m for each method.                       #
#   • Input: RDS files of the form                                        #
#       P=0.2/result_m<m>_prob0.2.rds                                     #
#       P=0.5/result_m<m>_prob0.5.rds                                     #
#     where each file contains a list of simulation runs with matrices    #
#     named e.g. ADCS_raw, ADCS_ceil, CED_raw, CED_ceil.                  #
#                                                                         #
###########################################################################
# ====== Load packages ======
library(readxl)
library(ggplot2)
library(compositions)
library(robCompositions)
library(zCompositions)
library(vegan)
library(foreach)
library(doMPI)
library(dplyr)
library(patchwork)
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




# ============================== #
#   Load simulation result files #
# ============================== #

load(file.path("M=50",   "results50_rabbit.RData"))
res_m50 <- m_compare_results

load(file.path("M=500",  "results500_rabbit.RData"))
res_m500 <- m_compare_results

load(file.path("M=1000", "results1000_rabbit.RData"))
res_m1000 <- m_compare_results

# ============================== #
#     Method sets and metadata   #
# ============================== #

methods_all <- c(
  "lmrob","PLS","mult_lognorm","mult_repl",
  "lr_da","lr_em","mult_KMSS","lr_SVD",
  "add1","dl_065","dl_unif","GBM"
)

fixed_methods <- c(
  "mult_lognorm","mult_repl","mult_KMSS",
  "lr_da","lr_em","PLS","lr_SVD",
  "dl_065","dl_unif","GBM","add1"
)

fixed_methods <- setdiff(fixed_methods, c("lmrob","dl_065"))

probs <- c(0.05, 0.4, 0.6, 0.8)
prob_strs <- sub("\\.?0+$", "", as.character(probs))

color_values <- c(
  "PLS"          = "#E41A1C",
  "mult_lognorm" = "#377EB8",
  "mult_repl"    = "#4DAF4A",
  "lr_da"        = "#984EA3",
  "lr_em"        = "#FF7F00",
  "mult_KMSS"    = "#FB9A99",
  "lr_SVD"       = "#A65628",
  "add1"         = "#F781BF",
  "dl_unif"      = "#999999",
  "GBM"          = "#66C2A5"
)

shape_values <- c(
  "PLS"          = 16,
  "mult_lognorm" = 17,
  "mult_repl"    = 15,
  "lr_da"        = 3,
  "lr_em"        = 7,
  "mult_KMSS"    = 8,
  "lr_SVD"       = 18,
  "add1"         = 4,
  "dl_unif"      = 9,
  "GBM"          = 2
)

# ====================================================== #
#   Figure 6 – Boxplots over p and M (CED / ADCS)        #
# ====================================================== #

extract_results <- function(res_list, M_value, metric_prefix){
  df_list <- list()
  
  for(p in prob_strs){
    ceil_key <- sprintf("%sceil_m%d_prob%s", metric_prefix, M_value, p)
    raw_key  <- sprintf("%sraw_m%d_prob%s",  metric_prefix, M_value, p)

    if(is.null(res_list[[ceil_key]]) || is.null(res_list[[raw_key]])) next

    B <- nrow(res_list[[ceil_key]])
    methods_sel <- fixed_methods
    if(M_value != 50) methods_sel <- setdiff(methods_sel, "lr_em")
    methods_exist <- intersect(methods_sel, colnames(res_list[[ceil_key]]))
    if(length(methods_exist) == 0) next

    df_ceil <- data.frame(
      error  = as.numeric(res_list[[ceil_key]][, methods_exist, drop = FALSE]),
      method = factor(rep(methods_exist, each = B), levels = fixed_methods),
      metric = "Ceil",
      prob   = p,
      M_dim  = factor(paste0("M=", M_value), levels = c("M=50","M=500","M=1000"))
    )

    df_raw <- data.frame(
      error  = as.numeric(res_list[[raw_key]][, methods_exist, drop = FALSE]),
      method = factor(rep(methods_exist, each = B), levels = fixed_methods),
      metric = "Raw",
      prob   = p,
      M_dim  = factor(paste0("M=", M_value), levels = c("M=50","M=500","M=1000"))
    )

    df_list[[paste0("ceil_", p)]] <- df_ceil
    df_list[[paste0("raw_",  p)]] <- df_raw
  }
  bind_rows(df_list)
}

CED_box_df <- bind_rows(
  extract_results(res_m50,   50,  "CED"),
  extract_results(res_m500,  500, "CED"),
  extract_results(res_m1000, 1000,"CED")
)

ADCS_box_df <- bind_rows(
  extract_results(res_m50,   50,  "ADCS"),
  extract_results(res_m500,  500, "ADCS"),
  extract_results(res_m1000, 1000,"ADCS")
)

plot_box <- function(df){
  ggplot(df, aes(x = method, y = error, fill = metric)) +
    geom_boxplot(alpha = 0.45, outlier.size = 0.35,
                 position = position_dodge(width = 0.8)) +
    scale_y_log10() +
    facet_grid(rows = vars(M_dim), cols = vars(prob), scales = "free_y") +
    labs(x = "", y = "Error (log10)") +
    theme_bw(base_size = 16) +
    theme(
      legend.position = "top",
      legend.title    = element_blank(),
      axis.text.x     = element_text(angle = 45, hjust = 1, size = 11),
      strip.text      = element_text(size = 13, face = "bold")
    )
}

p_box_CED  <- plot_box(CED_box_df)
p_box_ADCS <- plot_box(ADCS_box_df)

ggsave("Figure6_CED_boxplot.pdf",  p_box_CED,  width = 12, height = 8)
ggsave("Figure6_ADCS_boxplot.pdf", p_box_ADCS, width = 12, height = 8)

combined_fig6 <- p_box_CED + p_box_ADCS + 
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "top")

ggsave("Figure6_combined.pdf", combined_fig6, width = 16, height = 8)

# ====================================================== #
#   Figure 7 – Mean lines over p and M (CED / ADCS)      #
# ====================================================== #

extract_means <- function(res_list, metric_prefix){
  df_list <- list()
  for(M_value in c(50, 500, 1000)){
    for(p in prob_strs){

      key_ceil <- sprintf("%sceil_m%d_prob%s", metric_prefix, M_value, p)
      key_raw  <- sprintf("%sraw_m%d_prob%s",  metric_prefix, M_value, p)

      if(!(key_ceil %in% names(res_list))) next

      methods_sel <- fixed_methods
      if(M_value != 50) methods_sel <- setdiff(methods_sel, "lr_em")
      methods_exist <- intersect(methods_sel, colnames(res_list[[key_ceil]]))
      if(length(methods_exist) == 0) next

      mat_ceil <- res_list[[key_ceil]][, methods_exist, drop = FALSE]
      mat_raw  <- res_list[[key_raw]][,  methods_exist, drop = FALSE]

      means_ceil <- colMeans(mat_ceil, na.rm = TRUE)
      means_raw  <- colMeans(mat_raw,  na.rm = TRUE)

      df_list[[paste0("ceil_", M_value, p)]] <- data.frame(
        prob     = as.numeric(p),
        M_value  = M_value,
        method   = names(means_ceil),
        mean_val = means_ceil,
        type     = "Ceil"
      )
      df_list[[paste0("raw_",  M_value, p)]] <- data.frame(
        prob     = as.numeric(p),
        M_value  = M_value,
        method   = names(means_raw),
        mean_val = means_raw,
        type     = "Raw"
      )
    }
  }
  bind_rows(df_list)
}

CED_line_df  <- bind_rows(
  extract_means(res_m50,  "CED"),
  extract_means(res_m500, "CED"),
  extract_means(res_m1000,"CED")
)

ADCS_line_df <- bind_rows(
  extract_means(res_m50,  "ADCS"),
  extract_means(res_m500, "ADCS"),
  extract_means(res_m1000,"ADCS")
)

CED_line_df$method  <- factor(CED_line_df$method,  levels = fixed_methods)
ADCS_line_df$method <- factor(ADCS_line_df$method, levels = fixed_methods)

plot_line <- function(df){
  ggplot(df, aes(x = prob, y = mean_val, color = method)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2.5) +
    scale_y_log10() +
    facet_grid(rows = vars(M_value), cols = vars(type), scales = "free_y") +
    labs(x = "Probability of Missingness", y = "Mean Error (log10)") +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "top",
      legend.title    = element_blank(),
      axis.text.x     = element_text(angle = 45, hjust = 1),
      strip.text      = element_text(size = 13, face = "bold")
    )
}

p_line_CED  <- plot_line(CED_line_df)
p_line_ADCS <- plot_line(ADCS_line_df)

ggsave("Figure7_CED_line.pdf",  p_line_CED,  width = 12, height = 8)
ggsave("Figure7_ADCS_line.pdf", p_line_ADCS, width = 12, height = 8)

combined_fig7 <- p_line_CED + p_line_ADCS +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "top")

ggsave("Figure7_combined.pdf", combined_fig7, width = 16, height = 8)

# ====================================================== #
#   Figure 8 – Fixed p, varying m (CED / ADCS)           #
# ====================================================== #

m_values  <- c(50,100,200,300,400,500,600,700,800,900,1000,1400,1800)
probs_dim <- c(0.2, 0.5)

plot_methods_dim <- setdiff(methods_all, c("lmrob","dl_065","lr_da","lr_em"))

extract_means_prob <- function(metric_field, prob_val) {
  prob_str   <- sub("\\.?0+$", "", as.character(prob_val))
  means_list <- list()
  
  for (m in m_values) {
    base_path <- sprintf("P=%.1f", prob_val)
    file_path <- file.path(base_path, sprintf("result_m%d_prob%s.rds", m, prob_str))
    
    if (!file.exists(file_path)) next
    
    tmp <- readRDS(file_path)
    mat_list <- lapply(tmp, function(x) x[[metric_field]])
    mat <- do.call(rbind, mat_list)
    if (is.null(mat) || ncol(mat) == 0) next
    
    if (m > 1200) {
      methods <- setdiff(methods_all, c("PLS","lr_SVD","lr_da","lr_em","lmrob"))
    } else if (m > 500) {
      methods <- setdiff(methods_all, c("PLS","lmrob"))
    } else {
      methods <- methods_all
    }
    methods <- intersect(methods, colnames(mat))
    if (length(methods) == 0) next
    
    mean_val <- colMeans(mat[, methods, drop = FALSE], na.rm = TRUE)
    
    df <- data.frame(
      prob     = factor(prob_val),
      M        = m,
      method   = names(mean_val),
      mean_val = mean_val,
      metric   = metric_field
    )
    means_list[[as.character(m)]] <- df
  }
  
  bind_rows(means_list)
}

get_metric_data <- function(metric_field) {
  bind_rows(lapply(probs_dim, function(p) extract_means_prob(metric_field, p))) %>%
    filter(method %in% plot_methods_dim) %>%
    mutate(method = factor(method, levels = names(color_values)))
}

ADCS_raw_df  <- get_metric_data("ADCS_raw")  %>% mutate(type = "Raw",  group = "ADCS")
ADCS_ceil_df <- get_metric_data("ADCS_ceil") %>% mutate(type = "Ceil", group = "ADCS")
CED_raw_df   <- get_metric_data("CED_raw")   %>% mutate(type = "Raw",  group = "CED")
CED_ceil_df  <- get_metric_data("CED_ceil")  %>% mutate(type = "Ceil", group = "CED")

ADCS_df <- bind_rows(ADCS_raw_df, ADCS_ceil_df)
CED_df  <- bind_rows(CED_raw_df,  CED_ceil_df)

plot_lines_dim <- function(df, x_label) {
  ggplot(df, aes(x = M, y = mean_val, color = method, shape = method)) +
    geom_line(aes(group = method), linewidth = 1) +
    geom_point(size = 3) +
    scale_y_log10() +
    scale_color_manual(values = color_values) +
    scale_shape_manual(values = shape_values) +
    facet_grid(rows = vars(prob), cols = vars(type), scales = "free_y") +
    labs(x = x_label, y = "Mean Value (log10)") +
    theme_minimal(base_size = 14) +
    theme(
      legend.title    = element_blank(),
      legend.position = "top",
      legend.text     = element_text(size = 13),
      axis.text.x     = element_text(angle = 45, hjust = 1, size = 13),
      axis.text.y     = element_text(size = 13),
      axis.title      = element_text(size = 14),
      strip.text      = element_text(size = 12, face = "bold")
    )
}

p_CED_dim  <- plot_lines_dim(CED_df,  x_label = "CED: dimension m")
p_ADCS_dim <- plot_lines_dim(ADCS_df, x_label = "ADCS: dimension m")

combined_fig8 <- p_CED_dim + p_ADCS_dim +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "top")

ggsave(
  "Figure8_CED_ADCS_dimension_m.pdf",
  combined_fig8, width = 14, height = 9, useDingbats = FALSE
)









