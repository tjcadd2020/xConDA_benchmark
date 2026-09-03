#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(grid)
})

# Standalone reproduction and redesign of the disease-level TOPSIS correlation plot.
# All 17,750 observations and all 71 evaluated strategies are retained.

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) {
  stop("Run this file with Rscript.")
}

script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
output_dir <- dirname(script_path)
project_dir <- normalizePath(file.path(output_dir, ".."), mustWork = TRUE)
input_file <- file.path(project_dir, "all_method_all_simulation_all_metrics_performance.tsv")

if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

performance <- read.delim(input_file, check.names = FALSE)
required_columns <- c("Disease", "Method", "TOPSIS")
missing_columns <- setdiff(required_columns, colnames(performance))
if (length(missing_columns) > 0L) {
  stop("Missing required columns: ", paste(missing_columns, collapse = ", "))
}

disease_order <- c("ASD", "CRC", "IBD", "RA", "T2D")

# Data-integrity checks for the confirmed full-factorial dataset.
stopifnot(
  nrow(performance) == 17750L,
  dplyr::n_distinct(performance$Method) == 71L,
  all(as.integer(table(factor(performance$Disease, levels = disease_order))) == 3550L)
)

method_means <- performance %>%
  filter(Disease %in% disease_order) %>%
  group_by(Disease, Method) %>%
  summarise(TOPSIS = mean(TOPSIS, na.rm = TRUE), .groups = "drop") %>%
  complete(Disease = disease_order, Method) %>%
  pivot_wider(names_from = Disease, values_from = TOPSIS) %>%
  arrange(Method)

stopifnot(
  nrow(method_means) == 71L,
  !anyNA(method_means[, disease_order])
)

score_matrix <- as.matrix(method_means[, disease_order])
correlation_matrix <- cor(score_matrix, method = "pearson")

pvalue_matrix <- matrix(
  NA_real_,
  nrow = length(disease_order),
  ncol = length(disease_order),
  dimnames = list(disease_order, disease_order)
)

for (i in seq_along(disease_order)) {
  for (j in seq_along(disease_order)) {
    pvalue_matrix[i, j] <- if (i == j) {
      1
    } else {
      cor.test(score_matrix[, i], score_matrix[, j], method = "pearson")$p.value
    }
  }
}

plot_data <- expand.grid(
  row = seq_along(disease_order),
  column = seq_along(disease_order),
  KEEP.OUT.ATTRS = FALSE
) %>%
  mutate(
    row_disease = disease_order[row],
    column_disease = disease_order[column],
    correlation = correlation_matrix[cbind(row, column)],
    p_value = pvalue_matrix[cbind(row, column)]
  )

lower_data <- plot_data %>%
  filter(row > column) %>%
  mutate(
    coefficient_label = sprintf("%.2f", correlation),
    significance_label = case_when(
      p_value <= 0.001 ~ "***",
      p_value <= 0.01 ~ "**",
      p_value <= 0.05 ~ "*",
      TRUE ~ "ns"
    ),
    coefficient_y = row - 0.11,
    significance_y = row + 0.18
  )

diagonal_data <- plot_data %>%
  filter(row == column) %>%
  mutate(label = disease_order[row])

# Light sequential blue family, focused on correlations from 0.50 to 1.00.
# The upper end remains deliberately moderate rather than dark navy.
correlation_palette <- c(
  "#EEF5FA", "#DAEAF4", "#C3DEEE",
  "#A9D0E6", "#8FC0DB", "#73AECF"
)

correlation_plot <- ggplot() +
  geom_tile(
    data = lower_data,
    aes(x = column, y = row, fill = correlation),
    width = 1,
    height = 1,
    colour = "#C8C8C8",
    linewidth = 0.35
  ) +
  geom_tile(
    data = diagonal_data,
    aes(x = column, y = row),
    width = 1,
    height = 1,
    fill = "white",
    colour = "#C8C8C8",
    linewidth = 0.35
  ) +
  geom_text(
    data = lower_data,
    aes(x = column, y = coefficient_y, label = coefficient_label),
    family = "Helvetica",
    fontface = "bold",
    size = 7.5 / ggplot2::.pt,
    colour = "#2F6E9F"
  ) +
  geom_text(
    data = lower_data,
    aes(x = column, y = significance_y, label = significance_label),
    family = "Helvetica",
    fontface = "bold",
    size = 10.0 / ggplot2::.pt,
    colour = "#242424"
  ) +
  geom_text(
    data = diagonal_data,
    aes(x = column, y = row, label = label),
    family = "Helvetica",
    size = 7.5 / ggplot2::.pt,
    colour = "#74599C"
  ) +
  scale_fill_gradientn(
    colours = correlation_palette,
    limits = c(0.50, 1.00),
    breaks = seq(0.50, 1.00, by = 0.10),
    labels = function(x) sprintf("%.2f", x),
    guide = guide_colourbar(
      title = NULL,
      direction = "horizontal",
      barwidth = unit(50.9, "mm"),
      barheight = unit(2.4, "mm"),
      ticks = TRUE,
      frame.colour = "#333333",
      frame.linewidth = 0.4
    )
  ) +
  scale_x_continuous(limits = c(0.5, 5.5), expand = c(0, 0)) +
  scale_y_reverse(limits = c(5.5, 0.5), expand = c(0, 0)) +
  coord_fixed(clip = "off") +
  theme_void(base_family = "Helvetica", base_size = 7) +
  theme(
    legend.position = "bottom",
    legend.justification = "center",
    legend.box.margin = margin(t = 2, r = 0, b = 0, l = 0, unit = "pt"),
    legend.margin = margin(0, 0, 0, 0, unit = "pt"),
    legend.text = element_text(
      family = "Helvetica",
      size = 6.2,
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      colour = "#333333"
    ),
    plot.margin = margin(t = 0.8, r = 0.8, b = 0.5, l = 0.8, unit = "mm"),
    plot.background = element_rect(fill = "white", colour = NA)
  )

width_mm <- 210 / 4
height_mm <- width_mm * 688 / 568

pdf_file <- file.path(output_dir, "disease_topsis_correlation_lower_triangle.pdf")
png_file <- file.path(output_dir, "disease_topsis_correlation_lower_triangle.png")

grDevices::cairo_pdf(
  filename = pdf_file,
  width = width_mm / 25.4,
  height = height_mm / 25.4,
  family = "Helvetica",
  bg = "white"
)
print(correlation_plot)
dev.off()

grDevices::png(
  filename = png_file,
  width = width_mm,
  height = height_mm,
  units = "mm",
  res = 600,
  type = "cairo-png",
  family = "Helvetica",
  bg = "white"
)
print(correlation_plot)
dev.off()

message("Created: ", pdf_file)
message("Created: ", png_file)
