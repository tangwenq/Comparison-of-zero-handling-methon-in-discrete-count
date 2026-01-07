############################################################
## Puropose:  Generate all figures related to discrete CoDA,
##            including lattice vs simplex, depth resolution,
##            scaling & quantization effects, log-ratio shifts.
## 
##   - Fig 2  : Integer lattice N^2 vs continuous 3-part simplex
##   - Fig 3(a): Simplex closure at different sequencing depths
##   - Fig 3(b): CLR coordinates at different sequencing depths
##   - Fig 4(a): Scaling + ceiling quantization of counts (scatter)
##   - Fig 4(b): Scaling + ceiling quantization (log-ratio shifts)
############################################################

library(ggplot2)
library(extraDistr)
library(dplyr)
library(patchwork)
library(zCompositions)
set.seed(123)

############################################################
## Global theme, colors and helper functions
############################################################

## Base theme used across all ggplots
title_size       <- 10
axis_title_size  <- 10
axis_text_size   <- 10
strip_title_size <- 13

theme_coda <- theme_bw(base_size = axis_text_size) +
  theme(
    plot.title   = element_text(size = title_size, face = "bold", hjust = 0.5),
    axis.title   = element_text(size = axis_title_size),
    strip.text   = element_text(size = strip_title_size),
    legend.title = element_text(size = axis_title_size),
    legend.text  = element_text(size = axis_text_size)
  )

## Shared color palette (Morandi-like)
col_blue   <- "#5C6D82"
col_green  <- "#7C9885"
col_red    <- "#BC6C7C"

## Simplex projection of (A,B,C) onto an equilateral triangle
simplex_xy <- function(A, B, C) {
  x <- 0.5 * (2 * B + C)
  y <- (sqrt(3) / 2) * C
  data.frame(x = x, y = y)
}


## Triangle vertices for the 3-part simplex
triangle_vertices <- simplex_xy(
  A = c(1, 0, 0, 1),
  B = c(0, 1, 0, 0),
  C = c(0, 0, 1, 0)
)

## General saving helper: saves PDF + PNG
save_figure <- function(plot_obj, filename, width = 8, height = 4, dpi = 300) {
  ggsave(
    paste0(filename, ".pdf"),
    plot   = plot_obj,
    width  = width,
    height = height,
    device = cairo_pdf
  )
  ggsave(
    paste0(filename, ".png"),
    plot   = plot_obj,
    width  = width,
    height = height,
    dpi    = dpi
  )
  message("Saved: ", filename, ".pdf / .png")
}

############################################################
## Fig 2: Integer lattice N^2 vs continuous 3-part simplex
############################################################

## Left panel: 2D integer lattice N^2, 0..10
lattice_2d <- expand.grid(x = 0:10, y = 0:10)

p_A_lattice <- ggplot(lattice_2d, aes(x, y)) +
  geom_point(size = 2, color = col_blue) +
  scale_x_continuous(breaks = 0:10) +
  scale_y_continuous(breaks = 0:10) +
  coord_fixed() +
  theme_coda +
  labs(
    x = expression(x %in% N),
    y = expression(y %in% N)
  )

## Right panel: 3-part simplex with Dirichlet(1,1,1) samples
n_comp        <- 1000
alpha_simplex <- c(1, 1, 1)
P_simplex     <- extraDistr::rdirichlet(n_comp, alpha_simplex)
simp_xy       <- simplex_xy(P_simplex[, 1], P_simplex[, 2], P_simplex[, 3])

p_A_simplex <- ggplot() +
  geom_polygon(
    data = triangle_vertices,
    aes(x, y),
    fill  = NA,
    color = "grey40"
  ) +
  geom_point(
    data = simp_xy,
    aes(x, y),
    alpha = 0.5,
    size  = 1.8,
    color = col_green
  ) +
  coord_equal() +
  theme_coda +
  theme(
    axis.title = element_blank(),
    axis.text  = element_blank(),
    axis.ticks = element_blank()
  )

## Combine into Fig A
fig_A <- p_A_lattice + p_A_simplex + plot_layout(ncol = 2)

save_figure(fig_A, "fig_A_lattice_simplex")

############################################################
## Fig 3: Closure at different sequencing depths (simplex)
############################################################

## Simulate Dirichlet-multinomial counts at three depths
n_samples  <- 500
alpha_true <- c(4, 2, 1)
P_true     <- extraDistr::rdirichlet(n_samples, alpha_true)

depths      <- c(10, 100, 1000)
depth_label <- paste0("Depth = ", depths)

sim_list <- lapply(seq_along(depths), function(i) {
  N   <- depths[i]
  lab <- depth_label[i]
  
  counts <- t(apply(P_true, 1, function(p) {
    rmultinom(1, size = N, prob = p)
  }))
  colnames(counts) <- c("x", "y", "z")
  
  total <- rowSums(counts)
  comp  <- counts / total
  
  xy       <- simplex_xy(comp[, 1], comp[, 2], comp[, 3])
  clr_mat  <- t(apply(comp, 1, clr))
  colnames(clr_mat) <- c("clr_x", "clr_y", "clr_z")
  
  data.frame(
    x     = xy$x,
    y     = xy$y,
    clr_x = clr_mat[, 1],
    clr_y = clr_mat[, 2],
    clr_z = clr_mat[, 3],
    depth = lab
  )
})

sim_data <- do.call(rbind, sim_list)
sim_data$depth <- factor(sim_data$depth, levels = depth_label)

## Colors for depth-specific plots
depth_cols <- c(
  "Depth = 10"   = col_red,
  "Depth = 100"  = col_green,
  "Depth = 1000" = col_blue
)

## Fig C1: simplex closure at different depths (facets)
fig_C1 <- ggplot(sim_data, aes(x, y)) +
  geom_polygon(
    data = triangle_vertices,
    aes(x, y),
    inherit.aes = FALSE,
    fill        = NA,
    color       = "grey40"
  ) +
  geom_point(alpha = 0.6, size = 1.6, aes(color = depth)) +
  scale_color_manual(values = depth_cols) +
  coord_equal() +
  facet_wrap(~ depth, nrow = 1) +
  theme_coda +
  theme(
    axis.title      = element_blank(),
    axis.text       = element_blank(),
    axis.ticks      = element_blank(),
    plot.title      = element_blank(),
    legend.position = "none"
  )

save_figure(fig_C1, "fig_C1_simplex_depth")

############################################################
## Fig 3: CLR coordinates at three sequencing depths
############################################################

## Left: Depth = 10
fig_C2_left <- ggplot(subset(sim_data, depth == "Depth = 10"),
                      aes(x = clr_x, y = clr_y)) +
  geom_point(alpha = 0.7, size = 2, color = depth_cols["Depth = 10"]) +
  theme_coda +
  labs(x = "clr(x)", y = "clr(y)") +
  theme(plot.title = element_blank())

## Middle: Depth = 100
fig_C2_mid <- ggplot(subset(sim_data, depth == "Depth = 100"),
                     aes(x = clr_x, y = clr_y)) +
  geom_point(alpha = 0.7, size = 2, color = depth_cols["Depth = 100"]) +
  theme_coda +
  labs(x = "clr(x)", y = NULL) +
  theme(plot.title = element_blank())

## Right: Depth = 1000
fig_C2_right <- ggplot(subset(sim_data, depth == "Depth = 1000"),
                       aes(x = clr_x, y = clr_y)) +
  geom_point(alpha = 0.7, size = 2, color = depth_cols["Depth = 1000"]) +
  theme_coda +
  labs(x = "clr(x)", y = NULL) +
  theme(plot.title = element_blank())

fig_C2 <- fig_C2_left + fig_C2_mid + fig_C2_right +
  plot_layout(nrow = 1)

save_figure(fig_C2, "fig_C2_clr_depth")

############################################################
## Fig 4: Scaling + ceiling quantization (scatter)
############################################################

## Parameters for S1/S2 simulation
n_samples_dm <- 100
N_total_dm   <- 1000
alpha_dm     <- c(6, 3, 1)

P_dm <- extraDistr::rdirichlet(n_samples_dm, alpha_dm)

count_dm <- t(apply(P_dm, 1, function(p) {
  rmultinom(1, size = N_total_dm, prob = p)
}))
colnames(count_dm) <- c("x1", "x2", "x3")

x1_dm <- count_dm[, "x1"]
x2_dm <- count_dm[, "x2"]

scales      <- c(1, 0.1, 0.01, 0.001)
scale_label <- c("1x", "0.1x", "0.01x", "0.001x")

## Ceiling at four scales
data_list <- lapply(seq_along(scales), function(i) {
  s <- scales[i]
  data.frame(
    x1    = ceiling(x1_dm * s),
    x2    = ceiling(x2_dm * s),
    scale = scale_label[i]
  )
})

all_data <- do.call(rbind, data_list)
all_data$scale <- factor(all_data$scale, levels = rev(scale_label))

## Only 1x used for density contours
high_data <- subset(all_data, scale == "1x")

## Mean log-ratio per scale
lr_stats <- all_data %>%
  group_by(scale) %>%
  summarise(
    mean_x1 = mean(x1),
    mean_x2 = mean(x2),
    lr      = log10(mean_x1 / mean_x2),
    .groups = "drop"
  )

## Legend labels: "ceil(1x) (0.320)" etc.
legend_labels <- setNames(
  paste0("ceil(", lr_stats$scale, ") (", round(lr_stats$lr, 3), ")"),
  lr_stats$scale
)

morandi_scales <- c(
  "0.001x" = "#F4A259",
  "0.01x"  = "#7C9885",
  "0.1x"   = "#BC6C7C",
  "1x"     = "#5C6D82"
)

fig_S1 <- ggplot() +
  geom_abline(
    slope = seq(0.1, 2, length.out = 25),
    intercept = 0,
    alpha = 0.05,
    color  = "grey70"
  ) +
  geom_abline(
    slope = 0.5,
    intercept = 0,
    linetype = "dashed",
    color    = "grey40"
  ) +
  geom_point(
    data = all_data,
    aes(x1, x2, color = scale),
    alpha = 0.75, size = 2
  ) +
  stat_density_2d(
    data = high_data,
    aes(x1, x2),
    color = "grey35",
    size = 0.9
  ) +
  scale_color_manual(
    values = morandi_scales,
    labels = legend_labels
  ) +
  coord_cartesian(
    xlim = c(0, max(high_data$x1) * 1.02),
    ylim = c(0, max(high_data$x2) * 1.02)
  ) +
  labs(x = "x1", y = "x2", color = "scale") +
  theme_bw() +
  theme(
    plot.title = element_blank(),
    legend.position = c(0.14, 0.87),
    legend.background = element_rect(
      fill = alpha("white", 0.85),
      color = "grey80"
    ),
    legend.title = element_text(size = 9),
    legend.text  = element_text(size = 8)
  )

fig_S1

save_figure(fig_S1, "Fig_S1_scaling_ceiling", width = 5.5, height = 5.5)

############################################################
## Fig 4: Scaling + ceiling quantization (log-ratio shifts)
############################################################

valid_idx <- which(x1_dm > 0 & x2_dm > 0)
log_raw   <- log10(x1_dm[valid_idx] / x2_dm[valid_idx])

per_point <- lapply(seq_along(scales), function(i) {
  s   <- scales[i]
  lab <- scale_label[i]
  
  x1_c <- ceiling(x1_dm[valid_idx] * s)
  x2_c <- ceiling(x2_dm[valid_idx] * s)
  
  keep <- (x1_c > 0 & x2_c > 0)
  data.frame(
    scale = lab,
    shift = log10(x1_c[keep] / x2_c[keep]) - log_raw[keep]
  )
})

shift_df <- do.call(rbind, per_point)

fig_S2 <- ggplot(shift_df, aes(scale, shift)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_boxplot(fill = "#002FA7", alpha = 0.8) +
  theme_bw() +
  labs(
    y = expression(Delta * log[10] (x[1] / x[2])),
    x = "Scale"
  )







save_figure(fig_S2, "Fig_S2_logratio_shift",  width = 5.5, height = 5.5)


