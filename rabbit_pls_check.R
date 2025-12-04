###########################################################################
#                                                                         #
#      Local Examination of PLS Imputation Behavior (Example: m=200, p=0.2)
#                                                                         #
#  This script performs a detailed local inspection of PLS-based          #
#  imputation for left-censored compositional data. For a single          #
#  simulation replicate (m=200, p=0.2 here as an example). 
#  Result show in Figure 10
###########################################################################

# ================== Load packages ==================
library(readxl)
library(compositions)
library(robCompositions)
library(zCompositions)
library(moments)
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)

# ================== Download data ==================
url <- "https://github.com/michaelgreenacre/CODAinPractice/raw/master/Rabbits.xlsx"
temp_file <- tempfile(fileext = ".xlsx")
download.file(url, destfile = temp_file, mode = "wb")
rabbits_data <- read_excel(temp_file)[, -1]
rabbits_data <- rabbits_data[, -1]

# ================== Group columns based on means ==================
col_means <- colMeans(rabbits_data, na.rm = TRUE)
sorted_cols <- order(col_means)
sorted_names <- colnames(rabbits_data)[sorted_cols]
n <- length(sorted_cols)

group_size <- ceiling(n / 3)
low_cols    <- sorted_names[1:group_size]
medium_cols <- sorted_names[(group_size + 1):(2 * group_size)]
high_cols   <- sorted_names[(2 * group_size + 1):n]

# ================== Parameter settings ==================
rep_k1 <- 1
m1 <- 200
probs1 <- 0.2
set.seed(123)
seeds1 <- sample(1:1e6, rep_k1)

methods_list <- c(
  "add1", "dl_unif", "PLS", "mult_lognorm",
  "mult_repl", "lr_da", "mult_KMSS", "lr_SVD", "GBM"
)

# ================== Safe wrapper for imputation ==================
safe_impute <- function(expr_fun){
  res <- try(expr_fun(), silent = TRUE)
  if(inherits(res, "try-error") || is.null(res)) return(NULL)
  res
}

# ================== Data generation with structured tracked columns ==================
generate_data <- function(data, m, probs, seed, low_cols, medium_cols, high_cols){
  set.seed(seed)
  
  colnames(data) <- trimws(colnames(data))
  all_cols <- colnames(data)
  
  low_cols    <- trimws(low_cols)
  medium_cols <- trimws(medium_cols)
  high_cols   <- trimws(high_cols)
  
  tracked_low    <- intersect(sample(low_cols, 2), all_cols)
  tracked_medium <- intersect(sample(medium_cols, 2), all_cols)
  tracked_high   <- intersect(sample(high_cols, 2), all_cols)
  
  tracked_cols <- unique(c(tracked_low, tracked_medium, tracked_high))
  
  if(length(tracked_cols) < 6){
    tracked_cols <- c(tracked_cols,
                      sample(setdiff(all_cols, tracked_cols), 6 - length(tracked_cols)))
  }
  
  other_cols <- setdiff(all_cols, tracked_cols)
  num_needed <- max(0, m - length(tracked_cols))
  random_cols <- if(num_needed > 0) sample(other_cols, num_needed) else character(0)
  
  final_cols <- character(m)
  odd_pos <- seq(1, m, by = 2)
  if(length(odd_pos) < 6) stop("m too small")
  
  final_cols[odd_pos[1:6]] <- tracked_cols[1:6]
  remaining_slots <- which(final_cols == "")
  final_cols[remaining_slots] <- head(setdiff(random_cols, tracked_cols), length(remaining_slots))
  
  final_cols <- intersect(final_cols, all_cols)
  
  data_m <- data[, final_cols, drop = FALSE]
  data_0 <- data_m
  DL <- numeric(m)
  
  for(j in seq(1, m, by = 2))
    DL[j] <- quantile(data_m[, j], probs = probs, na.rm = TRUE)
  
  for(j in 1:m){
    col_j <- as.numeric(data_0[[j]])
    if(DL[j] > 0)
      data_0[col_j < DL[j], j] <- 0
  }
  
  list(data_m = data_m, data_0 = data_0, DL = DL,
       tracked_cols = tracked_cols, final_cols = final_cols)
}

# ================== Simple methods ==================
add1_method <- function(data_0){
  tmp <- data_0
  tmp[tmp == 0] <- 1
  tmp
}

dl_unif_method <- function(data_0, DL, seed){
  tmp <- data_0
  set.seed(seed)
  for(j in seq_along(DL)){
    idx <- which(tmp[, j] == 0)
    if(length(idx) > 0)
      tmp[idx, j] <- runif(length(idx), 0.1 * DL[j], DL[j])
  }
  tmp
}

# ================== Single-run imputation ==================
all_results1 <- list()
pls_components <- list()

data_list <- generate_data(
  rabbits_data, m1, probs1, seed = seeds1[1],
  low_cols, medium_cols, high_cols
)

data_0 <- data_list$data_0
DL <- data_list$DL
tracked_cols <- data_list$tracked_cols

data_0 <- as.matrix(sapply(data_0, as.numeric))
colnames(data_0) <- make.unique(colnames(data_0))

sim_results <- lapply(methods_list, function(method_name){
  
  imp_vals_df <- safe_impute(function(){
    if(method_name == "add1") return(add1_method(data_0))
    if(method_name == "dl_unif") return(dl_unif_method(data_0, DL, seed = 123))
    
    if(method_name == "PLS"){
      pls_res <- imputeBDLs(
        data_0, dl = DL, eps = 0.1, method = "pls",
        verbose = FALSE, R = 50, variation = FALSE, maxit = 50
      )
      
      imp_vals_df <- as.data.frame(pls_res$x)
      colnames(imp_vals_df) <- colnames(data_0)
      
      ncomp_df <- data.frame(
        tracked_col = colnames(imp_vals_df),
        nComp = pls_res$nComp,
        stringsAsFactors = FALSE
      )
      
      pls_components[[1]] <<- ncomp_df %>% filter(tracked_col %in% tracked_cols)
      return(imp_vals_df)
    }
    
    if(method_name == "mult_lognorm") return(multLN(data_0, label = 0, dl = DL))
    if(method_name == "mult_repl")   return(multRepl(data_0, label = 0, dl = DL))
    if(method_name == "lr_da")       return(lrDA(data_0, label = 0, dl = DL, ini.cov="multRepl"))
    if(method_name == "mult_KMSS")   return(multKM(data_0, label = 0, dl = DL))
    if(method_name == "lr_SVD")      return(lrSVD(data_0, label = 0, dl = DL))
    if(method_name == "GBM")         return(cmultRepl(data_0, output="p-counts", z.delete=FALSE))
    
    NULL
  })
  
  if(is.null(imp_vals_df)) return(NULL)
  
  cols_to_keep <- intersect(tracked_cols, colnames(imp_vals_df))
  if(length(cols_to_keep) == 0) return(NULL)
  
  res_list <- list()
  for(col in cols_to_keep){
    zero_idx <- which(data_0[, col] == 0)
    if(length(zero_idx) > 0){
      res_list[[col]] <- data.frame(
        sim = 1,
        method = method_name,
        tracked_col = col,
        imp_value = imp_vals_df[zero_idx, col]
      )
    }
  }
  
  if(length(res_list) == 0) return(NULL)
  do.call(rbind, res_list)
})

sim_results <- sim_results[!sapply(sim_results, is.null)]
if(length(sim_results) > 0){
  all_results1[[1]] <- do.call(rbind, sim_results)
}

all_results1_df <- do.call(rbind, all_results1)
pls_components_df <- do.call(rbind, pls_components)

# ================== Log10 boxplot ==================
facet_labels <- pls_components_df %>%
  mutate(label = paste0(tracked_col, " (nComp=", nComp, ")")) %>%
  select(tracked_col, label)

all_results1_df <- all_results1_df %>%
  left_join(facet_labels, by = "tracked_col")

p_box <- ggplot(all_results1_df, aes(x = method, y = imp_value, fill = method)) +
  geom_boxplot() +
  facet_wrap(~ label, scales = "free_y") +
  scale_y_log10() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = paste0("Log10 Imputed Values for Tracked Columns (p=", probs1, " m=", m1, ")"),
    y = "log10(Imputed Value)", x = "Method"
  )

print(p_box)

ggsave(paste0("boxplot_all_methods_p", probs1, "_m", m1, ".pdf"),
       p_box, width = 12, height = 6)

# ================== Scatter plots for 6 tracked columns ==================
tracked_cols_all <- data_list$tracked_cols
plot_list <- list()

for(colname in tracked_cols_all){
  
  subset_df <- all_results1_df %>%
    filter(tracked_col == colname, method %in% c("PLS", "lr_SVD"))
  
  if(nrow(subset_df) == 0) next
  
  zero_idx <- which(data_list$data_0[, colname] == 0)
  true_values <- data_list$data_m[zero_idx, colname]
  
  plot_df <- subset_df %>%
    mutate(true_value = true_values) %>%
    filter(imp_value > 0, true_value > 0)
  
  ncomp_label <- pls_components_df %>%
    filter(tracked_col == colname) %>% pull(nComp) %>% unique()
  
  ncomp_text <- if(length(ncomp_label) > 0)
    paste0(" (PLS nComp=", ncomp_label, ")") else ""
  
  p <- ggplot(plot_df, aes(x = true_value, y = imp_value, color = method)) +
    geom_point(size = 2, alpha = 0.7) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
    scale_x_continuous(trans = "log10") +
    scale_y_continuous(trans = "log10") +
    theme_bw() +
    labs(
      title = paste0("Tracked: ", colname, ncomp_text),
      x = "True Value (log10)",
      y = "Imputed Value (log10)"
    ) +
    theme(legend.position = "top")
  
  plot_list[[colname]] <- p
}

final_plot <- wrap_plots(plot_list, ncol = 3, nrow = 2) +
  plot_annotation(
    title = paste0("p", probs1, "m", m1),
    theme = theme(plot.title = element_text(size = 18, face = "bold"))
  )

print(final_plot)

ggsave(paste0("six_tracked_scatter_p", probs1, "_m", m1, ".pdf"),
       final_plot, width = 14, height = 8)


