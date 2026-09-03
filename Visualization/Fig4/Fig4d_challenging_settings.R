suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
})

# Rebuild the six-panel extreme-condition comparison from the exported source
# data. The original file and plot are left unchanged.
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg) == 1L) {
  dirname(normalizePath(sub("^--file=", "", script_arg)))
} else {
  normalizePath(getwd())
}

input_file <- file.path(script_dir, "worst_case_performance.tsv")
output_prefix <- file.path(script_dir, "worst_case_performance_with_AUROC_A4")

required_columns <- c(
  "Worst_scenario",
  "Extreme_scenario_setting",
  "Method",
  "Accuracy_metric",
  "Accuracy_value"
)

raw_data <- readr::read_tsv(
  input_file,
  show_col_types = FALSE,
  progress = FALSE,
  locale = readr::locale(encoding = "UTF-8")
)

missing_columns <- setdiff(required_columns, colnames(raw_data))
if (length(missing_columns) > 0L) {
  stop("Missing required columns: ", paste(missing_columns, collapse = ", "))
}

strategy_order <- c(
  "ANCOM-BC2_CSS",
  "edgeR_CSS",
  "ANCOM-BC2_count",
  "ANCOM-BC2_TMM",
  "MaAslin2_CPLM_TSS",
  "edgeR_count",
  "edgeR_TMM",
  "MaAslin2_NEGBIN_CSS",
  "MaAslin2_NEGBIN_TMM",
  "MaAslin2_ZINB_count",
  "DESeq2_count",
  "MaAslin2_NEGBIN_count"
)

setting_order <- c(
  "25 Cases + 25 Controls",
  "30 Cases + 570 Controls",
  "Primary effect = 2, Confounding effect = 2",
  "Primary effect = 4, Confounding effect = 5",
  "1 continuous variable + 4 categorical variables",
  "Prevalence:0~0.05"
)

panel_labels <- c(
  "25 Cases + 25 Controls" = "Sample size\n(n = 25 per group)",
  "30 Cases + 570 Controls" = "Case-control ratio\n(30:570)",
  "Primary effect = 2, Confounding effect = 2" = "Phenotype effect\n(effect = 2)",
  "Primary effect = 4, Confounding effect = 5" = "Confounding effect\n(effect = 5)",
  "1 continuous variable + 4 categorical variables" = "Confounder complexity\n(1 continuous + 4 categorical)",
  "Prevalence:0~0.05" = "Taxa prevalence\n(<0.05)"
)

metric_order <- c("AUPR", "AUROC", "MCC", "Macro-F1")
metric_colors <- c(
  "AUPR" = "#264653",
  "AUROC" = "#7A5195",
  "MCC" = "#E76F51",
  "Macro-F1" = "#2A9D8F"
)

plot_data <- raw_data %>%
  filter(
    Extreme_scenario_setting %in% setting_order,
    Accuracy_metric %in% c("AUPR", "AUC", "MCC", "Macro-F1")
  ) %>%
  mutate(
    Method = sub("^EdgeR", "edgeR", Method),
    Method = sub("^Maaslin2", "MaAslin2", Method),
    Accuracy_metric = if_else(Accuracy_metric == "AUC", "AUROC", Accuracy_metric),
    Method = factor(Method, levels = strategy_order),
    Accuracy_metric = factor(Accuracy_metric, levels = metric_order),
    Panel = factor(Extreme_scenario_setting, levels = setting_order)
  )

if (nrow(plot_data) != length(setting_order) * length(strategy_order) * length(metric_order)) {
  stop(
    "Unexpected selected row count: ", nrow(plot_data),
    "; expected ", length(setting_order) * length(strategy_order) * length(metric_order)
  )
}

if (anyNA(plot_data$Method) || anyNA(plot_data$Accuracy_metric) || anyNA(plot_data$Panel)) {
  stop("Unmapped strategy, metric, or panel label detected.")
}

duplicate_cells <- plot_data %>%
  count(Panel, Method, Accuracy_metric, name = "n") %>%
  filter(n != 1L)
if (nrow(duplicate_cells) > 0L) {
  stop("Each panel-strategy-metric combination must occur exactly once.")
}

axis_color <- "#4D4D4D"

p <- ggplot(
  plot_data,
  aes(x = Method, y = Accuracy_value, group = Accuracy_metric, color = Accuracy_metric)
) +
  geom_line(alpha = 0.88, linewidth = 0.48) +
  geom_point(size = 1.05, alpha = 0.95) +
  facet_grid(
    . ~ Panel,
    labeller = labeller(Panel = panel_labels),
    drop = FALSE
  ) +
  scale_color_manual(values = metric_colors, breaks = metric_order, drop = FALSE) +
  scale_y_continuous(
    breaks = seq(0, 1, by = 0.25),
    labels = sprintf("%.2f", seq(0, 1, by = 0.25)),
    expand = expansion(mult = c(0.015, 0.025))
  ) +
  scale_x_discrete(drop = FALSE, expand = expansion(mult = c(0.035, 0.035))) +
  coord_cartesian(ylim = c(0, 1.03), clip = "on") +
  labs(x = NULL, y = NULL, color = NULL) +
  theme_classic(base_size = 7, base_family = "Arial") +
  theme(
    axis.text.x = element_text(
      size = 6,
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      color = axis_color,
      margin = margin(t = 1.5)
    ),
    axis.text.y = element_text(size = 6.5, color = axis_color),
    axis.line = element_line(color = axis_color, linewidth = 0.35),
    axis.ticks = element_line(color = axis_color, linewidth = 0.35),
    axis.ticks.length = grid::unit(1.6, "pt"),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.text = element_text(size = 7, face = "bold", color = "#222222"),
    legend.key.width = grid::unit(12, "pt"),
    legend.key.height = grid::unit(7, "pt"),
    legend.spacing.x = grid::unit(2.5, "pt"),
    legend.margin = margin(t = 1, r = 0, b = 0, l = 0),
    strip.background = element_rect(
      fill = "#F2F2F2",
      color = axis_color,
      linewidth = 0.35
    ),
    strip.text.x = element_text(
      size = 6.8,
      face = "bold",
      color = axis_color,
      lineheight = 0.98,
      margin = margin(t = 3.4, b = 3.4)
    ),
    panel.border = element_rect(color = axis_color, fill = NA, linewidth = 0.35),
    panel.spacing.x = grid::unit(1.2, "pt"),
    plot.margin = margin(t = 3, r = 4, b = 2, l = 3)
  )

width_mm <- 210
height_mm <- 74
width_in <- width_mm / 25.4
height_in <- height_mm / 25.4

grDevices::cairo_pdf(
  filename = paste0(output_prefix, ".pdf"),
  width = width_in,
  height = height_in,
  family = "Arial",
  onefile = TRUE
)
print(p)
grDevices::dev.off()

svglite::svglite(
  filename = paste0(output_prefix, ".svg"),
  width = width_in,
  height = height_in,
  bg = "white"
)
print(p)
grDevices::dev.off()

ragg::agg_png(
  filename = paste0(output_prefix, ".png"),
  width = width_in,
  height = height_in,
  units = "in",
  res = 600,
  background = "white"
)
print(p)
grDevices::dev.off()

message("Created: ", paste0(output_prefix, ".pdf"))
message("Created: ", paste0(output_prefix, ".svg"))
message("Created: ", paste0(output_prefix, ".png"))
