###############################
## Rabbit count data: visualization scripts
## - Figure 1: count distribution (binned)
## - Figure 2: heatmap of sorted counts (by mean abundance)
## - Figure 3: log–log mean–variance plot
###############################

## ====== Load packages ======
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape2)
library(scales)

## ====== Download and prepare Rabbit data ======
# Source: Greenacre, Rabbits.xlsx from CODAinPractice GitHub repo
url <- "https://github.com/michaelgreenacre/CODAinPractice/raw/master/Rabbits.xlsx"
temp_file <- tempfile(fileext = ".xlsx")

download.file(url, destfile = temp_file, mode = "wb")

# Drop the first two columns (non-count columns)
rabbits_data <- read_excel(temp_file)[, -1]
rabbits_data <- rabbits_data[, -1]

dim(rabbits_data)
head(rabbits_data)

## =========================================================
## Figure 1: Binned count distribution across all OTUs/samples
## =========================================================

# Reshape to long format: each row is one count value
rabbits_long <- rabbits_data %>%
  pivot_longer(cols = everything(), names_to = "OTU", values_to = "Count")

# Define bins on the original count scale
bins <- c(0, 1, 10, 20, 50, 100, 200, 500, 1000, 5000)

rabbits_long <- rabbits_long %>%
  mutate(bin = cut(Count, breaks = bins, include.lowest = TRUE, right = FALSE))

# Compute percentage of entries in each bin
bin_summary <- rabbits_long %>%
  group_by(bin) %>%
  summarise(percent = n() / nrow(rabbits_long) * 100, .groups = "drop")

# Plot and save Figure 1
p_dist <- ggplot(bin_summary, aes(x = bin, y = percent)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  scale_y_continuous(expand = c(0, 0)) +
  labs(
    x = "Count bin",
    y = "Percentage of entries",
    title = "Rabbit data: empirical count distribution"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave("rabbit_count_distribution.pdf", p_dist, width = 7, height = 5)
ggsave("rabbit_count_distribution.png", p_dist, width = 7, height = 5, dpi = 300)

## =========================================================
## Figure 2: Heatmap of Rabbit counts sorted by mean abundance
## =========================================================

# Sort OTUs (columns) by their mean abundance
col_means <- colMeans(rabbits_data, na.rm = TRUE)
rabbits_sorted <- rabbits_data[, order(col_means)]

# Convert to long format for heatmap
rabbits_long_hm <- melt(as.matrix(rabbits_sorted))
colnames(rabbits_long_hm) <- c("Sample", "Species", "Value")

# Heatmap with custom color scale and low/medium/high annotation
p_hm <- ggplot(rabbits_long_hm, aes(x = Species, y = Sample, fill = Value)) +
  geom_tile() +
  scale_fill_gradientn(
    colors = rev(hcl.colors(8, "Spectral")),
    values = rescale(c(0, 2, 5, 10, 15, 50, 100, 200, 1024)),
    limits = c(0, 1024),
    oob = squish,
    breaks = c(10, 100, 500, 1000),
    labels = c("0–10", "10–100", "100–500", ">500")
  ) +
  theme_minimal(base_size = 18) +
  labs(
    x = NULL,            # remove x-axis title
    y = "Sample",
    fill = "Count range"
  ) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y  = element_blank(),
    axis.title.y = element_text(size = 20, face = "bold"),
    legend.title = element_text(size = 20, face = "bold"),
    legend.text  = element_text(size = 16),
    legend.position = "bottom",
    legend.direction = "horizontal",
    plot.margin = margin(10, 30, 10, 10)
  ) +
  # Vertical dashed lines at tertiles of the sorted species
  geom_vline(
    xintercept = c(ncol(rabbits_sorted) / 3, 2 * ncol(rabbits_sorted) / 3),
    linetype = "dashed",
    color = "black",
    size = 0.8
  ) +
  # Annotate low / medium / high abundance regions
  annotate(
    "text",
    x = ncol(rabbits_sorted) / 6,
    y = -3,
    label = "Low",
    size = 6,
    fontface = "bold"
  ) +
  annotate(
    "text",
    x = ncol(rabbits_sorted) / 2,
    y = -3,
    label = "Medium",
    size = 6,
    fontface = "bold"
  ) +
  annotate(
    "text",
    x = 5 * ncol(rabbits_sorted) / 6,
    y = -3,
    label = "High",
    size = 6,
    fontface = "bold"
  ) +
  guides(
    fill = guide_colorbar(
      barwidth = 40,
      barheight = 2.5,
      title.position = "top",
      title.hjust = 0.5
    )
  )

# Save heatmap as PDF and PNG
ggsave(
  filename = "rabbits_heatmap.pdf",
  plot = p_hm,
  width = 10,
  height = 8,
  device = cairo_pdf,
  useDingbats = FALSE
)

ggsave(
  filename = "rabbits_heatmap.png",
  plot = p_hm,
  width = 10,
  height = 8,
  dpi = 300
)

## =========================================================
## Figure 3: Log–log mean–variance relationship by mean-abundance tertile
## =========================================================

species_stats <- data.frame(
  Species  = colnames(rabbits_data),
  Mean     = colMeans(rabbits_data, na.rm = TRUE),
  Variance = apply(rabbits_data, 2, var, na.rm = TRUE)
)

# Define low / medium / high groups based on tertiles of mean
q <- quantile(species_stats$Mean, probs = c(1/3, 2/3))

species_stats <- species_stats %>%
  mutate(
    Group = case_when(
      Mean <= q[1] ~ "Low",
      Mean <= q[2] ~ "Medium",
      TRUE         ~ "High"
    )
  )

p_mv <- ggplot(species_stats, aes(x = Mean, y = Variance, color = Group)) +
  geom_point(size = 3) +
  scale_color_manual(
    values = c("Low" = "#1b9e77", "Medium" = "#d95f02", "High" = "#7570b3")
  ) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    color = "black",
    size = 1
  ) +  # Poisson reference line
  scale_x_log10() +
  scale_y_log10() +
  labs(
    x = "Mean count",
    y = "Count variance",
    color = "Abundance group"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    legend.position = "bottom",
    legend.title    = element_text(size = 16, face = "bold"),
    legend.text     = element_text(size = 14),
    axis.title      = element_text(size = 18),
    axis.text       = element_text(size = 14)
  )

ggsave(
  filename = "mean_variance_loglog.pdf",
  plot = p_mv,
  width = 8,
  height = 6,
  device = cairo_pdf,
  useDingbats = FALSE
)

ggsave(
  filename = "mean_variance_loglog.png",
  plot = p_mv,
  width = 8,
  height = 6,
  dpi = 300
)
