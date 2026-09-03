#!/usr/bin/env Rscript

# Publication-ready composite rank figure for the Full_factorial simulation.
# This script loads rank_metrics.RData directly into an isolated environment.
# It does not source rank_test.R and cannot enter the spiked0.05 code branch.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(ragg)
  library(scales)
})

required_packages <- c(
  "ggplot2", "dplyr", "tidyr", "patchwork", "ragg", "scales"
)
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
rdata_path <- file.path(project_dir, "rank_metrics.RData")
raw_full_factorial_path <- file.path(
  project_dir, "all_simulation_all_method_all_metrics_8.tsv"
)
disk_rank_path <- file.path(project_dir, "full_factorial_metrics_rank_table.tsv")

if (!file.exists(rdata_path)) {
  stop("Required input is missing: rank_metrics.RData")
}

# ----- Data-integrity gate -------------------------------------------------

data_env <- new.env(parent = baseenv())
loaded_objects <- load(rdata_path, envir = data_env)

required_objects <- c(
  "all_simulation_performance",
  "sub_all_simulation_performance_topsis",
  "topsis_by_method",
  "rank_table"
)
missing_objects <- setdiff(required_objects, loaded_objects)
if (length(missing_objects) > 0L) {
  stop("rank_metrics.RData is missing objects: ", paste(missing_objects, collapse = ", "))
}

if (length(grep("005|spiked", loaded_objects, ignore.case = TRUE)) > 0L) {
  stop("Unexpected spiked0.05-named object detected in rank_metrics.RData.")
}

all_performance <- data_env$all_simulation_performance
topsis_rows <- data_env$sub_all_simulation_performance_topsis
rank_table_rdata <- data_env$rank_table

expected_metrics <- c(
  "AUC", "AUPR", "FDR", "FPR", "Sensitivity", "MCC", "Macro.F1", "nMCC"
)
required_performance_columns <- c(
  "Index", "Disease", "Method", expected_metrics
)
if (!all(required_performance_columns %in% names(all_performance))) {
  stop("Full_factorial performance object does not contain all expected columns.")
}

method_counts <- table(all_performance$Method)
scenario_counts <- table(interaction(
  all_performance$Index, all_performance$Disease, drop = TRUE
))

stopifnot(
  nrow(all_performance) == 17750L,
  length(unique(all_performance$Method)) == 71L,
  length(unique(all_performance$Index)) == 250L,
  length(unique(all_performance$Disease)) == 5L,
  length(method_counts) == 71L,
  all(method_counts == 250L),
  length(scenario_counts) == 250L,
  all(scenario_counts == 71L),
  nrow(topsis_rows) == 17750L,
  !any(grepl("spiked0\\.05", as.character(all_performance$Method), ignore.case = TRUE)),
  !"maaslin2_cplm_LOG_TSS" %in% all_performance$Method
)

# Cross-check the RData object against the Full_factorial source table.
raw_source_match <- NA
raw_source_rows <- NA_integer_
excluded_rows <- NA_integer_
if (file.exists(raw_full_factorial_path)) {
  raw_full_factorial <- read.table(
    raw_full_factorial_path,
    sep = ",",
    header = TRUE,
    row.names = 1,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  excluded_flag <- raw_full_factorial$Method == "maaslin2_cplm_LOG_TSS"
  raw_source_rows <- nrow(raw_full_factorial)
  excluded_rows <- sum(excluded_flag)
  raw_source_kept <- raw_full_factorial[!excluded_flag, , drop = FALSE]
  raw_source_match <- isTRUE(all.equal(
    raw_source_kept,
    all_performance,
    check.attributes = FALSE,
    tolerance = 1e-12
  ))
  stopifnot(
    raw_source_rows == 18000L,
    excluded_rows == 250L,
    nrow(raw_source_kept) == 17750L,
    raw_source_match
  )
}

# Recalculate all eight metric ranks from the 17,750-row object, then verify
# them against the rank_table stored in rank_metrics.RData.
metric_means <- all_performance %>%
  group_by(Method) %>%
  summarise(across(all_of(expected_metrics), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

all_methods_sorted <- sort(unique(metric_means$Method))
rank_check <- data.frame(Method = all_methods_sorted, stringsAsFactors = FALSE)
for (metric_name in expected_metrics) {
  ranked_methods <- if (metric_name %in% c("FDR", "FPR")) {
    metric_means %>% arrange(.data[[metric_name]], Method) %>% pull(Method)
  } else {
    metric_means %>% arrange(desc(.data[[metric_name]]), Method) %>% pull(Method)
  }
  rank_check[[metric_name]] <- match(all_methods_sorted, ranked_methods)
}

rank_table_subset <- rank_table_rdata[, c("Method", expected_metrics), drop = FALSE]
rank_table_subset <- rank_table_subset[
  match(rank_check$Method, rank_table_subset$Method),
  ,
  drop = FALSE
]
rank_rdata_match <- isTRUE(all.equal(
  rank_check,
  rank_table_subset,
  check.attributes = FALSE
))
stopifnot(rank_rdata_match)

# The existing on-disk Full_factorial rank table is also checked. Its
# RankAggreg column is all NA and may differ only by logical/numeric storage.
disk_rank_match <- NA
if (file.exists(disk_rank_path)) {
  disk_rank <- read.table(
    disk_rank_path,
    sep = "\t",
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    quote = "\""
  )
  disk_rank_match <- isTRUE(all.equal(
    disk_rank[
      match(rank_check$Method, disk_rank$Method),
      c("Method", expected_metrics),
      drop = FALSE
    ],
    rank_table_subset,
    check.attributes = FALSE
  ))
  stopifnot(disk_rank_match)
}

# Recalculate the method-level TOPSIS mean and SD from all 17,750 scores and
# verify the result against the current RData environment.
topsis_summary <- topsis_rows %>%
  group_by(Method) %>%
  summarise(
    TOPSIS_score_mean = mean(TOPSIS_score, na.rm = TRUE),
    TOPSIS_score_sd = sd(TOPSIS_score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(TOPSIS_score_mean), Method)

topsis_rdata <- data_env$topsis_by_method %>%
  select(Method, TOPSIS_score_mean, TOPSIS_score_sd) %>%
  mutate(Method = as.character(Method))
topsis_check_sorted <- topsis_summary %>% arrange(Method)
topsis_rdata <- topsis_rdata[
  match(topsis_check_sorted$Method, topsis_rdata$Method),
  ,
  drop = FALSE
]
topsis_rdata_match <- isTRUE(all.equal(
  topsis_check_sorted,
  topsis_rdata,
  check.attributes = FALSE,
  tolerance = 1e-12
))
stopifnot(topsis_rdata_match)

# ----- Explicit field mapping ---------------------------------------------

method_rename <- function(x) {
  x <- gsub("maaslin2_lm_res", "MaAslin2_GLM", x, fixed = TRUE)
  x <- gsub("maaslin2_lm", "MaAslin2_GLM", x, fixed = TRUE)
  x <- gsub("maaslin2", "MaAslin2", x, fixed = TRUE)
  x <- gsub("cplm", "CPLM", x, fixed = TRUE)
  x <- gsub("negbin", "NEGBIN", x, fixed = TRUE)
  x <- gsub("megbin", "NEGBIN", x, fixed = TRUE)
  x <- gsub("zinb", "ZINB", x, fixed = TRUE)
  x <- gsub("LM-fixed", "LFEM", x, fixed = TRUE)
  x <- gsub("lmem", "LMEM", x, fixed = TRUE)
  x <- gsub("ANCOMBC", "ANCOM-BC2", x, fixed = TRUE)
  x
}

classify_method <- function(x) {
  case_when(
    grepl("^(edgeR|limma|DESeq2)_", x) ~ "RNA-seq-derived",
    grepl("^(ALDEx2|ANCOMBC|Corncob|fastANCOM|ZicoSeq|maaslin2)_", x) ~
      "Microbiome-tailored",
    TRUE ~ "Classical statistical"
  )
}

plot_metrics <- c("AUC", "AUPR", "MCC", "Macro.F1", "FDR", "FPR", "Sensitivity")
metric_display <- c(
  AUC = "AUROC",
  AUPR = "AUPR",
  MCC = "MCC",
  `Macro.F1` = "Macro-F1",
  FDR = "FDR",
  FPR = "FPR",
  Sensitivity = "Sensitivity"
)

# A small horizontal gap separates the four TOPSIS components from the three
# operating characteristics while retaining one aligned heatmap panel.
metric_positions <- c(
  AUC = 1,
  AUPR = 2,
  MCC = 3,
  `Macro.F1` = 4,
  FDR = 5.28,
  FPR = 6.28,
  Sensitivity = 7.28
)
metric_group_centres <- c(2.5, mean(metric_positions[c("FDR", "FPR", "Sensitivity")]))
metric_group_labels <- c("TOPSIS\ncomponents", "Operating\ncharacteristics")

method_summary <- topsis_summary %>%
  left_join(rank_table_subset, by = "Method") %>%
  mutate(
    Overall_rank = row_number(),
    Method_label = method_rename(Method),
    Method_class = classify_method(Method),
    Method_class = factor(
      Method_class,
      levels = c("Microbiome-tailored", "RNA-seq-derived", "Classical statistical")
    ),
    TOPSIS_lower = pmax(0, TOPSIS_score_mean - TOPSIS_score_sd),
    TOPSIS_upper = pmin(1, TOPSIS_score_mean + TOPSIS_score_sd)
  )

if (anyDuplicated(method_summary$Method_label)) {
  duplicated_labels <- unique(method_summary$Method_label[
    duplicated(method_summary$Method_label) |
      duplicated(method_summary$Method_label, fromLast = TRUE)
  ])
  stop("Display-name collision after renaming: ", paste(duplicated_labels, collapse = ", "))
}

write.table(
  method_summary,
  file = file.path(output_dir, "full_factorial_method_summary_all71.tsv"),
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)

class_palette <- c(
  "Microbiome-tailored" = "#5FA8CC",
  "RNA-seq-derived" = "#978DCA",
  "Classical statistical" = "#8C979D"
)

class_legend_labels <- c(
  "Microbiome-tailored" = "Microbiome-tailored methods",
  "RNA-seq-derived" = "RNA-seq-derived methods",
  "Classical statistical" = "Classical statistical models"
)

font_family <- "Helvetica"
base_font_pt <- 6.5
body_font_pt <- 6.2
tile_font_pt <- 6.0
title_font_pt <- 7.0
legend_font_pt <- 6.0
pt_to_mm <- 1 / ggplot2::.pt

theme_rank <- function() {
  theme_classic(base_size = base_font_pt, base_family = font_family) +
    theme(
      axis.line = element_line(linewidth = 0.3, colour = "#262626"),
      axis.ticks = element_line(linewidth = 0.3, colour = "#262626"),
      axis.title = element_text(size = base_font_pt, colour = "#262626"),
      axis.text = element_text(size = body_font_pt, colour = "#262626"),
      plot.title = element_text(
        size = title_font_pt, face = "bold", colour = "#202020",
        margin = margin(b = 3, unit = "pt")
      ),
      legend.title = element_text(size = legend_font_pt, face = "bold"),
      legend.text = element_text(size = legend_font_pt),
      legend.box = "vertical",
      legend.box.spacing = unit(2, "pt"),
      legend.key.width = unit(3.5, "mm"),
      legend.key.height = unit(3.5, "mm"),
      legend.margin = margin(t = 2, unit = "pt"),
      panel.grid = element_blank(),
      plot.margin = margin(2, 1.5, 2, 1.5, unit = "pt")
    )
}

panel_blank_theme <- function(title_hjust = 0.5) {
  theme_rank() +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank(),
      plot.title = element_text(
        size = title_font_pt, face = "bold", hjust = title_hjust,
        colour = "#202020", margin = margin(b = 3, unit = "pt")
      )
    )
}

build_rank_figure <- function(top_n) {
  plot_df <- method_summary %>%
    slice_head(n = top_n) %>%
    mutate(
      Method_plot = factor(
        Method_label,
        levels = c(rev(Method_label), ".metric_group_header")
      ),
      cap_ymin = as.numeric(Method_plot) - 0.14,
      cap_ymax = as.numeric(Method_plot) + 0.14
    )

  stopifnot(
    min(plot_df$TOPSIS_lower) >= 0.4,
    max(plot_df$TOPSIS_upper) <= 1
  )

  write.table(
    plot_df %>% select(-Method_plot, -cap_ymin, -cap_ymax),
    file = file.path(output_dir, paste0("top", top_n, "_plot_data.tsv")),
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
  )

  shared_y <- scale_y_discrete(drop = FALSE, expand = expansion(add = c(0.3, 0.75)))

  p_rank <- ggplot(plot_df, aes(y = Method_plot)) +
    geom_text(
      aes(x = 0, label = Overall_rank),
      family = font_family,
      size = body_font_pt * pt_to_mm,
      colour = "#222222",
      hjust = 0.5
    ) +
    shared_y +
    scale_x_continuous(limits = c(-0.5, 0.5), expand = c(0, 0)) +
    labs(title = "Rank") +
    panel_blank_theme(0.5)

  p_method <- ggplot(plot_df, aes(y = Method_plot)) +
    geom_text(
      aes(x = 0, label = Method_label),
      family = font_family,
      size = body_font_pt * pt_to_mm,
      colour = "#222222",
      hjust = 0
    ) +
    shared_y +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    coord_cartesian(clip = "off") +
    labs(title = "Strategy") +
    panel_blank_theme(0)

  p_class <- ggplot(plot_df, aes(y = Method_plot)) +
    geom_point(
      aes(x = 0, colour = Method_class),
      shape = 16,
      size = 2.6
    ) +
    shared_y +
    scale_x_continuous(limits = c(-0.5, 0.5), expand = c(0, 0)) +
    scale_colour_manual(
      values = class_palette,
      labels = class_legend_labels,
      drop = FALSE,
      name = NULL,
      guide = guide_legend(
        nrow = 1,
        byrow = TRUE,
        override.aes = list(shape = 16, size = 3.2)
      )
    ) +
    labs(title = "Type") +
    panel_blank_theme(0.5)

  p_lollipop <- ggplot(plot_df, aes(y = Method_plot)) +
    geom_segment(
      aes(x = TOPSIS_lower, xend = TOPSIS_upper, yend = Method_plot),
      linewidth = 0.48,
      colour = "#7C858C",
      lineend = "round"
    ) +
    geom_segment(
      aes(x = TOPSIS_lower, xend = TOPSIS_lower, y = cap_ymin, yend = cap_ymax),
      linewidth = 0.38,
      colour = "#7C858C"
    ) +
    geom_segment(
      aes(x = TOPSIS_upper, xend = TOPSIS_upper, y = cap_ymin, yend = cap_ymax),
      linewidth = 0.38,
      colour = "#7C858C"
    ) +
    geom_point(
      aes(x = TOPSIS_score_mean),
      shape = 21,
      size = 1.7,
      stroke = 0.35,
      fill = "#111111",
      colour = "#111111"
    ) +
    shared_y +
    scale_x_continuous(
      limits = c(0.4, 1),
      breaks = seq(0.4, 1, by = 0.2),
      labels = label_number(accuracy = 0.1),
      expand = c(0, 0)
    ) +
    labs(title = "TOPSIS Score", x = NULL, y = NULL) +
    theme_rank() +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.y = element_blank(),
      plot.title = element_text(
        size = title_font_pt, face = "bold", hjust = 0.5,
        colour = "#202020", margin = margin(b = 3, unit = "pt")
      )
    )

  heatmap_df <- plot_df %>%
    select(Method_plot, all_of(plot_metrics)) %>%
    pivot_longer(
      cols = all_of(plot_metrics),
      names_to = "Metric",
      values_to = "Metric_rank"
    ) %>%
    mutate(
      Metric_position = unname(metric_positions[Metric]),
      Rank_text_colour = "#3B2B20"
    )

  heatmap_group_df <- data.frame(
    Metric_position = metric_group_centres,
    Method_plot = factor(
      rep(".metric_group_header", length(metric_group_centres)),
      levels = levels(plot_df$Method_plot)
    ),
    Group_label = metric_group_labels
  )

  p_heatmap <- ggplot(heatmap_df, aes(x = Metric_position, y = Method_plot)) +
    geom_tile(
      aes(fill = Metric_rank),
      colour = "white",
      linewidth = 0.22,
      width = 0.98
    ) +
    geom_text(
      aes(label = Metric_rank, colour = Rank_text_colour),
      family = font_family,
      size = tile_font_pt * pt_to_mm,
      fontface = "plain"
    ) +
    geom_text(
      data = heatmap_group_df,
      aes(x = Metric_position, y = Method_plot, label = Group_label),
      inherit.aes = FALSE,
      family = font_family,
      size = 6.0 * pt_to_mm,
      lineheight = 0.92,
      colour = "#4A4038"
    ) +
    shared_y +
    scale_x_continuous(
      breaks = unname(metric_positions[plot_metrics]),
      labels = unname(metric_display[plot_metrics]),
      limits = c(0.5, 7.78),
      expand = c(0, 0)
    ) +
    scale_fill_gradient(
      low = "#C6A476",
      high = "#F4F0EA",
      limits = c(1, 71),
      breaks = c(1, 20, 40, 60, 71),
      name = "Rank (1 = best)",
      guide = guide_colorbar(
        direction = "horizontal",
        title.position = "top",
        title.hjust = 0.5,
        barwidth = unit(34, "mm"),
        barheight = unit(2.2, "mm"),
        ticks = FALSE,
        frame.colour = "#555555",
        frame.linewidth = 0.25
      )
    ) +
    scale_colour_identity(guide = "none") +
    labs(title = "Metric-specific Rank", x = NULL, y = NULL) +
    theme_minimal(base_size = base_font_pt, base_family = font_family) +
    theme(
      axis.title = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.text.x = element_text(
        size = body_font_pt,
        angle = 45,
        hjust = 1,
        vjust = 1,
        colour = "#222222",
        margin = margin(t = 2, unit = "pt")
      ),
      panel.grid = element_blank(),
      plot.title = element_text(
        size = title_font_pt, face = "bold", hjust = 0.5,
        colour = "#202020", margin = margin(b = 3, unit = "pt")
      ),
      legend.title = element_text(size = legend_font_pt, face = "bold"),
      legend.text = element_text(size = legend_font_pt),
      legend.box = "vertical",
      legend.margin = margin(t = 2, unit = "pt"),
      plot.margin = margin(2, 1.5, 2, 1.5, unit = "pt")
    )

  p_rank + p_method + p_class + p_lollipop + p_heatmap +
    plot_layout(
      widths = c(0.055, 0.255, 0.055, 0.285, 0.35),
      guides = "collect"
    ) &
    theme(
      text = element_text(family = font_family),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "vertical",
      legend.justification = "center"
    )
}

save_figure_bundle <- function(plot, stem, width_mm, height_mm) {
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4

  pdf_file <- file.path(output_dir, paste0(stem, ".pdf"))
  png_file <- file.path(output_dir, paste0(stem, ".png"))

  # The local R runtime starts in the C locale and its graphics devices cannot
  # open a path containing non-ASCII characters. Draw to an ASCII-only
  # temporary path, then copy the unchanged files into this isolated folder.
  temp_pdf <- tempfile(pattern = paste0(stem, "_"), fileext = ".pdf")
  temp_png <- tempfile(pattern = paste0(stem, "_"), fileext = ".png")

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

  ragg::agg_png(
    temp_png,
    width = width_in,
    height = height_in,
    units = "in",
    res = 600,
    pointsize = 12,
    background = "white"
  )
  print(plot)
  grDevices::dev.off()

  copied <- c(
    file.copy(temp_pdf, pdf_file, overwrite = TRUE),
    file.copy(temp_png, png_file, overwrite = TRUE)
  )
  if (!all(copied)) {
    stop("Failed to copy one or more rendered files into the output folder.")
  }
  unlink(c(temp_pdf, temp_png))
  unlink(file.path(output_dir, paste0(stem, c(".svg", ".tiff"))))
}

figure_top20 <- build_rank_figure(20L)
figure_top30 <- build_rank_figure(30L)

# Width is exactly two-thirds of A4 (210 mm). Heights retain at least 6 pt text.
save_figure_bundle(
  figure_top20,
  "full_factorial_rank_top20",
  width_mm = 140,
  height_mm = 118
)
save_figure_bundle(
  figure_top30,
  "full_factorial_rank_top30",
  width_mm = 140,
  height_mm = 160
)

validation_lines <- c(
  "Full_factorial rank figure data validation",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("Input RData: rank_metrics.RData"),
  paste0("Loaded object count: ", length(loaded_objects)),
  paste0("Spiked/005-named objects: ", length(grep("005|spiked", loaded_objects, ignore.case = TRUE))),
  paste0("Full_factorial rows: ", nrow(all_performance)),
  paste0("Methods: ", length(method_counts)),
  paste0("Rows per method: ", paste(range(method_counts), collapse = "-")),
  paste0("Index values: ", length(unique(all_performance$Index))),
  paste0("Disease values: ", length(unique(all_performance$Disease))),
  paste0("Index-Disease scenarios: ", length(scenario_counts)),
  paste0("Methods per scenario: ", paste(range(scenario_counts), collapse = "-")),
  paste0("Raw Full_factorial rows before exclusion: ", raw_source_rows),
  paste0("Excluded maaslin2_cplm_LOG_TSS rows: ", excluded_rows),
  paste0("Raw source after exclusion matches RData: ", raw_source_match),
  paste0("Recomputed metric ranks match RData rank_table: ", rank_rdata_match),
  paste0("Existing Full_factorial rank table matches RData values: ", disk_rank_match),
  paste0("Recomputed TOPSIS mean/SD match RData topsis_by_method: ", topsis_rdata_match),
  "TOPSIS definition: within-scenario min-max normalization of AUPR, AUROC, MCC, and Macro-F1 with equal weights.",
  "Variability: SD across 250 Index-Disease scenarios for each method.",
  "Heatmap: method rank per metric; FDR/FPR ascending, all other metrics descending.",
  "Heatmap order: AUROC, AUPR, MCC, Macro-F1 | FDR, FPR, Sensitivity.",
  "Heatmap groups: TOPSIS components | Operating characteristics.",
  "TOPSIS display axis: 0.4 to 1.0; all Top 20 and Top 30 mean-SD intervals remain inside the axis.",
  "Outputs: editable-text PDF and 600 dpi PNG only.",
  "No rows were sampled or omitted after loading the verified 17,750-row Full_factorial object.",
  "The plotting script does not source or execute rank_test.R."
)
writeLines(validation_lines, file.path(output_dir, "data_validation.txt"))

capture.output(
  sessionInfo(),
  file = file.path(output_dir, "R_sessionInfo.txt")
)

message("Completed Full_factorial Top 20 and Top 30 figure bundles in: ", output_dir)
