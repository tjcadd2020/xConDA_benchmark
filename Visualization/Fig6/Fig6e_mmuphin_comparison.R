#!/usr/bin/env Rscript

# Fig. 6e | Multi-cohort strategies compared with MMUPHin
#
# Layout: simulation data on the left; real-data consistency on the right.
# This panel keeps the main structure used in the Illustrator-composed figure.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(tidyr)
})

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_directory <- if (length(script_argument) == 1L) {
  dirname(normalizePath(sub("^--file=", "", script_argument)))
} else {
  normalizePath(getwd())
}

input_file <- file.path(script_directory, "meta_plot_2606.RData")
output_directory <- file.path(script_directory, "Fig6e_output")
output_prefix <- file.path(
  output_directory,
  "Fig6e_mmuphin_comparison"
)

if (!file.exists(input_file)) {
  stop("Input environment not found: ", input_file)
}
if (dir.exists(output_directory)) {
  stop("Output directory already exists; refusing to overwrite: ", output_directory)
}

figure_environment <- new.env(parent = globalenv())
load(input_file, envir = figure_environment)

required_objects <- c(
  "evaluation_result",
  "consistency_summary_performance",
  "integrated_consistency_summary_performance"
)
missing_objects <- setdiff(required_objects, ls(figure_environment))
if (length(missing_objects) > 0L) {
  stop("Required object(s) missing: ", paste(missing_objects, collapse = ", "))
}

evaluation_result <- figure_environment$evaluation_result
consistency_result <- figure_environment$consistency_summary_performance
integrated_consistency <- figure_environment$integrated_consistency_summary_performance

strategy_order_top_to_bottom <- c(
  "ANCOM-BC2_TMM",
  "ANCOM-BC2_count",
  "ANCOM-BC2_CSS",
  "edgeR_CSS",
  "MMUPHin"
)
strategy_labels <- c(
  `ANCOM-BC2_TMM` = "ANCOM-BC2_TMM_Mega",
  `ANCOM-BC2_count` = "ANCOM-BC2_count_Mega",
  `ANCOM-BC2_CSS` = "ANCOM-BC2_CSS_Mega",
  edgeR_CSS = "edgeR_CSS_Mega",
  MMUPHin = "MMUPHin"
)

simulation_metric_order <- c(
  "FDR", "Sensitivity", "AUROC", "AUPR", "Macro-F1", "MCC"
)
simulation_metric_labels <- c(
  FDR = "FDR\n(lower)",
  Sensitivity = "Sens.",
  AUROC = "AUROC",
  AUPR = "AUPR",
  `Macro-F1` = "Macro-\nF1",
  MCC = "MCC"
)

performance_required_columns <- c(
  "Method", "FDR", "Sensitivity", "AUC", "AUPR", "Macro.F1", "MCC"
)
consistency_required_columns <- c(
  "Disease", "Strategy", "Mean_Jaccard_index", "Mean_DA_count"
)
integrated_required_columns <- c(
  "Strategy", "Mean_Jaccard_index", "Mean_DA_count"
)

if (!all(performance_required_columns %in% colnames(evaluation_result))) {
  stop("evaluation_result does not contain all required Fig. 5e columns.")
}
if (!all(consistency_required_columns %in% colnames(consistency_result))) {
  stop("consistency_summary_performance does not contain all required Fig. 5e columns.")
}
if (!all(integrated_required_columns %in% colnames(integrated_consistency))) {
  stop("integrated consistency data do not contain all required Fig. 5e columns.")
}

performance_wide <- evaluation_result %>%
  filter(Method %in% strategy_order_top_to_bottom) %>%
  transmute(
    Strategy = Method,
    FDR = FDR,
    Sensitivity = Sensitivity,
    AUROC = AUC,
    AUPR = AUPR,
    `Macro-F1` = Macro.F1,
    MCC = MCC
  )

if (nrow(performance_wide) != length(strategy_order_top_to_bottom) ||
    anyDuplicated(performance_wide$Strategy) > 0L ||
    !setequal(performance_wide$Strategy, strategy_order_top_to_bottom)) {
  stop("Each requested strategy must occur exactly once in evaluation_result.")
}

performance_long <- performance_wide %>%
  pivot_longer(
    cols = all_of(simulation_metric_order),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    Strategy = factor(Strategy, levels = rev(strategy_order_top_to_bottom)),
    Metric = factor(Metric, levels = simulation_metric_order),
    Metric_x = match(as.character(Metric), simulation_metric_order),
    Metric_family = if_else(as.character(Metric) == "FDR", "FDR", "Other simulation metrics")
  )

if (anyNA(performance_long) ||
    any(!is.finite(performance_long$Value)) ||
    any(performance_long$Value < 0 | performance_long$Value > 1)) {
  stop("Simulation metrics must be complete, finite, and within [0, 1].")
}

consistency_points <- consistency_result %>%
  filter(Strategy %in% strategy_order_top_to_bottom) %>%
  select(all_of(consistency_required_columns))

consistency_means <- integrated_consistency %>%
  filter(Strategy %in% strategy_order_top_to_bottom) %>%
  select(all_of(integrated_required_columns))

point_counts <- consistency_points %>% count(Strategy, name = "n_diseases")
if (nrow(point_counts) != length(strategy_order_top_to_bottom) ||
    any(point_counts$n_diseases != 4L)) {
  stop("Fig. 5e expects four disease observations per strategy.")
}
if (nrow(consistency_means) != length(strategy_order_top_to_bottom) ||
    anyDuplicated(consistency_means$Strategy) > 0L) {
  stop("Each requested strategy must occur exactly once in the consistency means.")
}
if (anyNA(consistency_points) || anyNA(consistency_means)) {
  stop("Consistency inputs contain missing values; no rows were silently removed.")
}

consistency_metric_labels <- c(
  Mean_Jaccard_index = "Jaccard",
  Mean_DA_count = "DA count"
)

consistency_points_long <- consistency_points %>%
  pivot_longer(
    cols = c(Mean_Jaccard_index, Mean_DA_count),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    Strategy = factor(Strategy, levels = rev(strategy_order_top_to_bottom)),
    Metric = factor(
      Metric,
      levels = names(consistency_metric_labels),
      labels = unname(consistency_metric_labels)
    )
  )

consistency_means_long <- consistency_means %>%
  pivot_longer(
    cols = c(Mean_Jaccard_index, Mean_DA_count),
    names_to = "Metric",
    values_to = "Mean"
  ) %>%
  mutate(
    Strategy = factor(Strategy, levels = rev(strategy_order_top_to_bottom)),
    Metric = factor(
      Metric,
      levels = names(consistency_metric_labels),
      labels = unname(consistency_metric_labels)
    )
  )

# Low-saturation purple-blue-green palette adapted from the provided reference.
palette <- c(
  purple = "#76539A",
  purple_outline = "#5D397F",
  blue = "#397FBB",
  blue_outline = "#2369A6",
  green = "#7F9F55",
  green_outline = "#5F803C",
  point = "#2D4058",
  separator = "#D7D7D7",
  row_guide = "#E5E5E5"
)

simulation_plot <- ggplot(
  performance_long,
  aes(x = Metric_x, y = Strategy)
) +
  geom_vline(
    xintercept = seq(1.5, length(simulation_metric_order) - 0.5, by = 1),
    colour = palette[["separator"]],
    linetype = "dotted",
    linewidth = 0.28
  ) +
  geom_point(
    aes(size = Value, fill = Metric_family, colour = Metric_family),
    shape = 21,
    stroke = 0.38,
    alpha = 0.96,
    show.legend = FALSE
  ) +
  scale_size_area(limits = c(0, 1), max_size = 4.8) +
  scale_fill_manual(
    values = c(
      FDR = palette[["purple"]],
      `Other simulation metrics` = palette[["blue"]]
    )
  ) +
  scale_colour_manual(
    values = c(
      FDR = palette[["purple_outline"]],
      `Other simulation metrics` = palette[["blue_outline"]]
    )
  ) +
  scale_x_continuous(
    position = "top",
    breaks = seq_along(simulation_metric_order),
    labels = unname(simulation_metric_labels[simulation_metric_order]),
    limits = c(0.52, length(simulation_metric_order) + 0.48),
    expand = expansion(mult = 0)
  ) +
  scale_y_discrete(labels = strategy_labels, drop = FALSE) +
  labs(tag = "e", title = "Simulation data", x = NULL, y = NULL) +
  coord_cartesian(clip = "off") +
  theme_void(base_family = "Arial", base_size = 6.2) +
  theme(
    axis.text.x.top = element_text(
      family = "Arial",
      size = 6.0,
      face = "bold",
      colour = "#202020",
      lineheight = 0.90,
      margin = margin(b = 2.5, unit = "pt")
    ),
    axis.text.y = element_text(
      family = "Arial",
      size = 6.2,
      colour = "#202020",
      margin = margin(r = 3, unit = "pt")
    ),
    plot.title = element_text(
      family = "Arial",
      size = 7.2,
      face = "bold",
      hjust = 0.5,
      colour = "#111111",
      margin = margin(b = 3, unit = "pt")
    ),
    plot.tag = element_text(
      family = "Arial",
      size = 10,
      face = "bold",
      colour = "#111111"
    ),
    plot.tag.position = c(-0.01, 1.08),
    plot.margin = margin(t = 8, r = 3, b = 10, l = 3, unit = "pt")
  )

consistency_plot <- ggplot(
  consistency_means_long,
  aes(x = Mean, y = Strategy)
) +
  geom_col(
    width = 0.58,
    fill = palette[["green"]],
    colour = palette[["green_outline"]],
    linewidth = 0.32
  ) +
  geom_point(
    data = consistency_points_long,
    aes(x = Value, y = Strategy),
    inherit.aes = FALSE,
    position = position_jitter(width = 0, height = 0.072, seed = 2606),
    size = 0.82,
    colour = palette[["point"]],
    alpha = 0.94
  ) +
  facet_wrap(~ Metric, nrow = 1, scales = "free_x") +
  scale_x_continuous(
    breaks = scales::breaks_pretty(n = 3),
    expand = expansion(mult = c(0, 0.035))
  ) +
  scale_y_discrete(drop = FALSE) +
  labs(title = "Real-data consistency", x = NULL, y = NULL) +
  theme_classic(base_family = "Arial", base_size = 6.2) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(
      family = "Arial",
      size = 6.0,
      face = "bold",
      colour = "#202020",
      margin = margin(b = 2.5, unit = "pt")
    ),
    axis.line = element_line(linewidth = 0.32, colour = "black"),
    axis.ticks = element_line(linewidth = 0.28, colour = "black"),
    axis.ticks.length = grid::unit(1.3, "pt"),
    axis.text.x = element_text(
      family = "Arial",
      size = 6.0,
      colour = "#202020",
      margin = margin(t = 1.5, unit = "pt")
    ),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_line(
      colour = palette[["row_guide"]],
      linewidth = 0.24
    ),
    panel.grid.minor = element_blank(),
    panel.spacing.x = grid::unit(4, "pt"),
    plot.title = element_text(
      family = "Arial",
      size = 7.2,
      face = "bold",
      hjust = 0.5,
      colour = "#111111",
      margin = margin(b = 3, unit = "pt")
    ),
    plot.margin = margin(t = 8, r = 3, b = 4, l = 3, unit = "pt")
  )

combined_plot <- (simulation_plot | consistency_plot) +
  plot_layout(widths = c(1.72, 1.00))

# Preserve the requested half-A4 width and reduce v1 height by one half.
fig_width_mm = 105
fig_height_mm <- 62.5
figure_width_in <- fig_width_mm / 25.4
figure_height_in <- fig_height_mm / 25.4

dir.create(output_directory, recursive = FALSE)

output_files <- paste0(output_prefix, c(".svg", ".pdf", ".tiff", ".png"))
if (any(file.exists(output_files))) {
  stop("One or more Fig. 5e v2 outputs already exist; refusing to overwrite.")
}

# Render through an ASCII-safe temporary path for graphics-device compatibility.
temporary_prefix <- file.path(
  tempdir(),
  paste0("Fig6e_mmuphin_comparison_", Sys.getpid())
)
temporary_files <- paste0(temporary_prefix, c(".svg", ".pdf", ".tiff", ".png"))

svglite::svglite(
  paste0(temporary_prefix, ".svg"),
  width = figure_width_in,
  height = figure_height_in,
  bg = "white"
)
print(combined_plot)
dev.off()

grDevices::cairo_pdf(
  paste0(temporary_prefix, ".pdf"),
  width = figure_width_in,
  height = figure_height_in,
  family = "Arial",
  bg = "white"
)
print(combined_plot)
dev.off()

ragg::agg_tiff(
  paste0(temporary_prefix, ".tiff"),
  width = figure_width_in,
  height = figure_height_in,
  units = "in",
  res = 600,
  background = "white",
  scaling = 1
)
print(combined_plot)
dev.off()

ragg::agg_png(
  paste0(temporary_prefix, ".png"),
  width = figure_width_in,
  height = figure_height_in,
  units = "in",
  res = 300,
  background = "white",
  scaling = 1
)
print(combined_plot)
dev.off()

copy_success <- file.copy(temporary_files, output_files, overwrite = FALSE)
if (!all(copy_success)) {
  stop(
    "Failed to copy one or more completed Fig. 5e v2 exports: ",
    paste(basename(output_files[!copy_success]), collapse = ", ")
  )
}

write.table(
  performance_wide,
  file = file.path(output_directory, "Fig6e_source_data_simulation.tsv"),
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)
write.table(
  consistency_points,
  file = file.path(output_directory, "Fig6e_source_data_consistency_by_disease.tsv"),
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)
write.table(
  consistency_means,
  file = file.path(output_directory, "Fig6e_source_data_consistency_means.tsv"),
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)

writeLines(
  c(
    "Fig. 6e: multi-cohort strategies compared with MMUPHin",
    sprintf("Final size: %.1f x %.1f mm", fig_width_mm, fig_height_mm),
    "Layout: simulation data left; real-data consistency right.",
    "Palette: FDR purple; other simulation metrics blue; consistency green.",
    "Simulation bubble area encodes the original 0-1 value; no within-column normalization.",
    "FDR is explicitly labelled lower; other simulation metrics are higher-is-better.",
    "Consistency bars show across-disease means; dark points show four individual diseases.",
    "Typography: Arial; minimum rendered text size: 6 pt.",
    "No rows were excluded after selecting the five requested strategies.",
    sprintf("Simulation observations: %d strategies x %d metrics", nrow(performance_wide), length(simulation_metric_order)),
    sprintf("Consistency observations: %d strategy-disease rows", nrow(consistency_points)),
    "Source environment: meta_plot_2606.RData"
  ),
  con = file.path(output_directory, "README.txt")
)

message("Created Fig. 5e v2 in: ", output_directory)
