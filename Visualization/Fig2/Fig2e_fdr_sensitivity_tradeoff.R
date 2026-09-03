#!/usr/bin/env Rscript

# Full_factorial trade-off between mean FDR and mean sensitivity.
# One point represents one strategy averaged over all 250 scenarios.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(ragg)
})

required_packages <- c("ggplot2", "dplyr", "ragg")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop("Missing R packages: ", paste(missing_packages, collapse = ", "))
}

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) {
  stop("Run this file with Rscript.")
}
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
output_dir <- dirname(script_path)
project_dir <- normalizePath(file.path(output_dir, ".."), mustWork = TRUE)

input_path <- file.path(
  project_dir,
  "all_method_all_simulation_all_metrics_performance.tsv"
)
rdata_path <- file.path(project_dir, "rank_metrics.RData")

if (!file.exists(input_path)) {
  stop("Missing input: all_method_all_simulation_all_metrics_performance.tsv")
}
if (!file.exists(rdata_path)) {
  stop("Missing input: rank_metrics.RData")
}

# ----- Data validation -----------------------------------------------------

performance <- read.delim(
  input_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

required_columns <- c(
  "Index", "Disease", "Method", "FDR", "Sensitivity"
)
if (!all(required_columns %in% names(performance))) {
  stop("The TSV does not contain all required fields: ", paste(required_columns, collapse = ", "))
}

method_counts <- table(performance$Method)
scenario_counts <- table(interaction(
  performance$Index,
  performance$Disease,
  drop = TRUE
))

stopifnot(
  nrow(performance) == 17750L,
  length(unique(performance$Method)) == 71L,
  length(unique(performance$Index)) == 250L,
  length(unique(performance$Disease)) == 5L,
  length(method_counts) == 71L,
  all(method_counts == 250L),
  length(scenario_counts) == 250L,
  all(scenario_counts == 71L),
  all(is.finite(performance$FDR)),
  all(is.finite(performance$Sensitivity)),
  all(performance$FDR >= 0 & performance$FDR <= 1),
  all(performance$Sensitivity >= 0 & performance$Sensitivity <= 1),
  !any(grepl("spiked", performance$Method, ignore.case = TRUE))
)

# Cross-check the complete TSV against the corresponding object in RData.
data_env <- new.env(parent = baseenv())
loaded_objects <- load(rdata_path, envir = data_env)
if (!"all_simulation_performance_output" %in% loaded_objects) {
  stop("rank_metrics.RData is missing all_simulation_performance_output.")
}

rdata_output <- as.data.frame(
  data_env$all_simulation_performance_output,
  stringsAsFactors = FALSE
)
rdata_output$Method <- as.character(rdata_output$Method)
common_columns <- intersect(names(performance), names(rdata_output))
row_key_file <- do.call(paste, performance[c("Index", "Disease", "Method")])
row_key_rdata <- do.call(paste, rdata_output[c("Index", "Disease", "Method")])
rdata_output <- rdata_output[
  match(row_key_file, row_key_rdata),
  common_columns,
  drop = FALSE
]
rdata_match <- isTRUE(all.equal(
  performance[, common_columns, drop = FALSE],
  rdata_output,
  check.attributes = FALSE,
  tolerance = 1e-12
))
stopifnot(rdata_match)

# ----- Strategy-level means -----------------------------------------------

classify_method <- function(x) {
  case_when(
    grepl("^(edgeR|limma|DESeq2)_", x) ~ "RNA-seq-derived",
    grepl("^(ALDEx2|ANCOM-BC2|Corncob|fastANCOM|ZicoSeq|MaAslin2)_", x) ~
      "Microbiome-tailored",
    TRUE ~ "Classical statistical"
  )
}

strategy_means <- performance %>%
  group_by(Method) %>%
  summarise(
    Mean_FDR = mean(FDR),
    Mean_Sensitivity = mean(Sensitivity),
    Scenario_count = n(),
    .groups = "drop"
  ) %>%
  mutate(
    Mean_1_minus_FDR = 1 - Mean_FDR,
    Method_class = classify_method(Method),
    Method_class = factor(
      Method_class,
      levels = c("Microbiome-tailored", "RNA-seq-derived", "Classical statistical")
    ),
    Is_labelled = Method %in% c("ALDEx2_count", "edgeR_CSS") |
      grepl("^ANCOM-BC2_", Method)
  ) %>%
  arrange(desc(Mean_Sensitivity), Mean_FDR, Method)

stopifnot(
  nrow(strategy_means) == 71L,
  all(strategy_means$Scenario_count == 250L),
  sum(strategy_means$Is_labelled) == 5L
)

labelled_methods <- strategy_means %>% filter(Is_labelled)
expected_labels <- c(
  "ALDEx2_count",
  "ANCOM-BC2_CSS",
  "ANCOM-BC2_TMM",
  "ANCOM-BC2_count",
  "edgeR_CSS"
)
stopifnot(setequal(labelled_methods$Method, expected_labels))

manual_label_positions <- data.frame(
  Method = c(
    "ANCOM-BC2_TMM",
    "ANCOM-BC2_count",
    "ANCOM-BC2_CSS",
    "edgeR_CSS",
    "ALDEx2_count"
  ),
  Label_x = c(0.40, 0.59, 0.54, 0.87, 0.75),
  Label_y = c(0.98, 0.98, 0.86, 0.85, 0.445),
  Segment_x = c(0.42, 0.57, 0.62, 0.82, 0.855),
  Segment_y = c(0.975, 0.975, 0.89, 0.825, 0.445),
  Hjust = c(1, 0, 0.5, 0.5, 0.5),
  stringsAsFactors = FALSE
)
manual_labelled_methods <- labelled_methods %>%
  inner_join(manual_label_positions, by = "Method")
stopifnot(nrow(manual_labelled_methods) == 5L)

write.table(
  strategy_means,
  file = file.path(output_dir, "fdr_sensitivity_strategy_means.tsv"),
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)

# ----- Publication figure -------------------------------------------------

class_palette <- c(
  "Microbiome-tailored" = "#5FA8CC",
  "RNA-seq-derived" = "#978DCA",
  "Classical statistical" = "#8C979D"
)
label_palette <- c(
  "Microbiome-tailored" = "#2F78A6",
  "RNA-seq-derived" = "#746AAF",
  "Classical statistical" = "#667177"
)

font_family <- "Helvetica"
base_font_pt <- 6.5
body_font_pt <- 6.2
label_font_pt <- 6.1
title_font_pt <- 7.5
pt_to_mm <- 1 / ggplot2::.pt

tradeoff_plot <- ggplot(
  strategy_means,
  aes(x = Mean_1_minus_FDR, y = Mean_Sensitivity)
) +
  geom_point(
    aes(fill = Method_class),
    shape = 21,
    size = 1.55,
    stroke = 0.25,
    colour = "#65727A",
    alpha = 0.96,
    show.legend = FALSE
  ) +
  geom_segment(
    data = manual_labelled_methods,
    aes(
      x = Mean_1_minus_FDR,
      y = Mean_Sensitivity,
      xend = Segment_x,
      yend = Segment_y
    ),
    inherit.aes = FALSE,
    colour = "#6F7478",
    linewidth = 0.28,
    show.legend = FALSE
  ) +
  geom_text(
    data = manual_labelled_methods,
    aes(
      x = Label_x,
      y = Label_y,
      label = Method,
      colour = Method_class,
      hjust = Hjust
    ),
    inherit.aes = FALSE,
    family = font_family,
    size = label_font_pt * pt_to_mm,
    fontface = "plain",
    vjust = 0.5,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = class_palette, drop = FALSE, guide = "none") +
  scale_colour_manual(
    values = label_palette,
    drop = FALSE,
    guide = "none"
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.25),
    labels = function(x) sprintf("%.2f", x),
    expand = expansion(mult = c(0.015, 0.015))
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.25),
    labels = function(x) sprintf("%.2f", x),
    expand = expansion(mult = c(0.015, 0.015))
  ) +
  labs(
    title = "The trade-off between FDR and sensitivity",
    x = "Mean (1- FDR) across scenarios",
    y = "Mean sensitivity across scenarios"
  ) +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = base_font_pt, base_family = font_family) +
  theme(
    axis.line = element_line(linewidth = 0.35, colour = "#262626"),
    axis.ticks = element_line(linewidth = 0.35, colour = "#262626"),
    axis.title = element_text(
      size = base_font_pt,
      face = "bold",
      colour = "#262626"
    ),
    axis.text = element_text(
      size = body_font_pt,
      face = "bold",
      colour = "#262626"
    ),
    plot.title = element_text(
      size = title_font_pt,
      face = "bold",
      hjust = 0.5,
      colour = "#202020",
      margin = margin(b = 4, unit = "pt")
    ),
    legend.position = "none",
    panel.grid = element_blank(),
    plot.margin = margin(3, 3, 3, 3, unit = "pt")
  )

save_png_pdf <- function(
  plot,
  stem,
  width_mm = 70,
  height_mm = 70 * 776 / 742,
  dpi = 600
) {
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4
  final_png <- file.path(output_dir, paste0(stem, ".png"))
  final_pdf <- file.path(output_dir, paste0(stem, ".pdf"))

  temp_png <- tempfile(pattern = paste0(stem, "_"), fileext = ".png")
  temp_pdf <- tempfile(pattern = paste0(stem, "_"), fileext = ".pdf")

  ragg::agg_png(
    temp_png,
    width = width_in,
    height = height_in,
    units = "in",
    res = dpi,
    pointsize = 12,
    background = "white"
  )
  print(plot)
  grDevices::dev.off()

  grDevices::cairo_pdf(
    temp_pdf,
    width = width_in,
    height = height_in,
    family = font_family,
    bg = "white",
    onefile = FALSE
  )
  print(plot)
  grDevices::dev.off()

  copied <- c(
    file.copy(temp_png, final_png, overwrite = TRUE),
    file.copy(temp_pdf, final_pdf, overwrite = TRUE)
  )
  if (!all(copied)) {
    stop("Failed to copy PNG/PDF into the output directory.")
  }
  unlink(c(temp_png, temp_pdf))
}

save_png_pdf(
  tradeoff_plot,
  stem = "full_factorial_mean_fdr_sensitivity_tradeoff",
  width_mm = 70,
  height_mm = 70 * 776 / 742,
  dpi = 600
)

validation_lines <- c(
  "Full_factorial mean FDR-sensitivity trade-off validation",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "Input: all_method_all_simulation_all_metrics_performance.tsv",
  paste0("Data rows: ", nrow(performance)),
  paste0("File lines including header: ", nrow(performance) + 1L),
  paste0("Strategies: ", length(method_counts)),
  paste0("Rows per strategy: ", paste(range(method_counts), collapse = "-")),
  paste0("Index-Disease scenarios: ", length(scenario_counts)),
  paste0("Strategies per scenario: ", paste(range(scenario_counts), collapse = "-")),
  paste0("TSV matches rank_metrics.RData output: ", rdata_match),
  paste0("Strategy-level points: ", nrow(strategy_means)),
  paste0("Labelled strategies: ", paste(sort(labelled_methods$Method), collapse = ", ")),
  "Aggregation: arithmetic mean of FDR and Sensitivity across all 250 scenarios for each strategy.",
  "Displayed x-coordinate: 1 - Mean_FDR, equivalent to the mean of (1 - FDR) across scenarios.",
  "Figure size: 70.00 mm wide x 73.21 mm high (reference aspect ratio 742:776).",
  "No rows were sampled or excluded.",
  "Outputs: PNG at 600 dpi and editable-text PDF only."
)
writeLines(validation_lines, file.path(output_dir, "data_validation.txt"))

capture.output(
  sessionInfo(),
  file = file.path(output_dir, "R_sessionInfo.txt")
)

message("Completed FDR-sensitivity trade-off figure in: ", output_dir)
