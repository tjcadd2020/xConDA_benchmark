#!/usr/bin/env Rscript

# Fig. 3c: oral real-data recovery breadth versus non-prior calls.
# Layout, legend, relative point sizes and axes follow the paired Fig. 3d style.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(grid)
  library(svglite)
  library(ragg)
})

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
script_dir <- if (length(file_arg) == 1L) {
  dirname(normalizePath(sub("^--file=", "", file_arg)))
} else {
  getwd()
}

old_working_dir <- getwd()
setwd(script_dir)
on.exit(setwd(old_working_dir), add = TRUE)

input_file <- "matched_result.RData"
output_prefix <- "Fig3c"

if (!file.exists(input_file)) {
  stop("Missing input file: ", input_file)
}

plot_env <- new.env(parent = globalenv())
load(input_file, envir = plot_env)

if (!exists("df_plot", envir = plot_env, inherits = FALSE)) {
  stop("df_plot was not found in matched_result.RData")
}
if (!exists("pal", envir = plot_env, inherits = FALSE)) {
  stop("pal was not found in matched_result.RData")
}

plot_data <- as.data.frame(plot_env$df_plot, stringsAsFactors = FALSE)
required_columns <- c(
  "Method", "Method_group", "True_positive", "False_positive"
)
missing_columns <- setdiff(required_columns, names(plot_data))
if (length(missing_columns) > 0L) {
  stop("Missing required columns: ", paste(missing_columns, collapse = ", "))
}

if (nrow(plot_data) != 59L) {
  stop("Expected 59 plotted strategies; found ", nrow(plot_data))
}
if (anyDuplicated(as.character(plot_data$Method))) {
  stop("Method names are not unique")
}
if (any(!is.finite(plot_data$True_positive)) ||
    any(!is.finite(plot_data$False_positive))) {
  stop("Plotting coordinates contain missing or non-finite values")
}

# Each method family has its own independent colour. Cross-panel colour
# harmonisation is intentionally not imposed for this figure.
method_levels <- names(plot_env$pal)
method_palette <- plot_env$pal
missing_colours <- setdiff(
  unique(as.character(plot_data$Method_group)),
  names(method_palette)
)
if (length(missing_colours) > 0L) {
  stop("Missing method colours: ", paste(missing_colours, collapse = ", "))
}

# Separate exact coordinate duplicates deterministically while preserving their
# integer source coordinates. This matches the paired figure point layout.
plot_data <- plot_data %>%
  mutate(
    source_row = row_number(),
    Method = as.character(Method),
    Method_group = factor(as.character(Method_group), levels = method_levels)
  ) %>%
  group_by(True_positive, False_positive) %>%
  arrange(Method, .by_group = TRUE) %>%
  mutate(
    duplicate_n = n(),
    duplicate_index = row_number(),
    display_radius = case_when(
      duplicate_n == 1L ~ 0,
      duplicate_n <= 3L ~ 0.20,
      duplicate_n <= 5L ~ 0.28,
      TRUE ~ 0.38
    ),
    display_angle = if_else(
      duplicate_n == 1L,
      0,
      2 * pi * (duplicate_index - 1) / duplicate_n + pi / 8
    ),
    display_x = True_positive + display_radius * cos(display_angle),
    display_y = False_positive + display_radius * sin(display_angle)
  ) %>%
  ungroup() %>%
  arrange(source_row)

# Fig. 3c annotations retained from the preceding oral-data version.
annotation_spec <- data.frame(
  method_key = c(
    "ZicoSeq_TSS",
    "MaAslin2_GLM_CLR",
    "Wilcoxon_CSS",
    "MaAslin2_GLM_LOG_TSS",
    "fastANCOM_count"
  ),
  label_x = rep(14.2, 5),
  label_y = c(23.2, 21.7, 20.2, 18.7, 17.2),
  stringsAsFactors = FALSE
)

match_index <- match(annotation_spec$method_key, plot_data$Method)
if (anyNA(match_index)) {
  stop(
    "Requested annotation method(s) not found: ",
    paste(annotation_spec$method_key[is.na(match_index)], collapse = ", ")
  )
}

annotation_data <- bind_cols(
  plot_data[match_index, , drop = FALSE],
  annotation_spec[, c("label_x", "label_y"), drop = FALSE]
)

font_family <- "Arial"

# A4 width is 210 mm; Fig. 3c is exactly half-width. Height preserves the
# supplied 850:500 reference aspect ratio.
width_mm <- 105
height_mm <- width_mm * 500 / 850
width_in <- width_mm / 25.4
height_in <- height_mm / 25.4

fig3c <- ggplot(plot_data, aes(x = display_x, y = display_y)) +
  geom_point(
    aes(fill = Method_group),
    shape = 21,
    colour = "white",
    stroke = 0.22,
    size = 1.75,
    alpha = 0.66
  ) +
  geom_segment(
    data = annotation_data,
    aes(
      x = label_x + 0.10,
      y = label_y,
      xend = display_x,
      yend = display_y
    ),
    inherit.aes = FALSE,
    colour = "#8A8A8A",
    linewidth = 0.34,
    alpha = 0.82,
    show.legend = FALSE
  ) +
  geom_point(
    data = annotation_data,
    aes(x = display_x, y = display_y, fill = Method_group),
    inherit.aes = FALSE,
    shape = 21,
    colour = "#252525",
    stroke = 0.58,
    size = 2.85,
    alpha = 1,
    show.legend = FALSE
  ) +
  geom_label(
    data = annotation_data,
    aes(
      x = label_x,
      y = label_y,
      label = Method,
      colour = Method_group
    ),
    inherit.aes = FALSE,
    hjust = 1,
    family = font_family,
    fontface = "plain",
    size = 2.15,
    fill = "white",
    label.size = 0.28,
    label.padding = unit(0.07, "lines"),
    label.r = unit(0, "lines"),
    show.legend = FALSE
  ) +
  scale_fill_manual(
    values = method_palette,
    breaks = method_levels,
    labels = method_levels,
    drop = FALSE,
    name = NULL
  ) +
  scale_colour_manual(values = method_palette, guide = "none") +
  scale_x_continuous(
    breaks = c(0, 5, 10, 15, 20),
    expand = expansion(mult = 0)
  ) +
  scale_y_continuous(
    breaks = c(0, 5, 10, 15, 20, 25),
    expand = expansion(mult = 0)
  ) +
  coord_cartesian(
    xlim = c(-0.6, 24.0),
    ylim = c(-0.6, 25.0),
    clip = "on"
  ) +
  labs(
    x = "Recovery breadth",
    y = "Non-prior calls"
  ) +
  guides(
    fill = guide_legend(
      ncol = 2,
      byrow = FALSE,
      override.aes = list(
        shape = 21,
        size = 1.55,
        alpha = 0.90,
        colour = "white",
        stroke = 0.22
      )
    )
  ) +
  theme_classic(base_size = 6.5, base_family = font_family) +
  theme(
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.45
    ),
    axis.line = element_blank(),
    axis.ticks = element_line(colour = "#4D4D4D", linewidth = 0.35),
    axis.ticks.length = unit(1.1, "mm"),
    axis.text = element_text(
      colour = "#4D4D4D",
      size = 6.5,
      face = "bold"
    ),
    axis.title = element_text(
      colour = "black",
      size = 8,
      face = "bold"
    ),
    axis.title.x = element_text(margin = margin(t = 2.2, unit = "mm")),
    axis.title.y = element_text(margin = margin(r = 2.2, unit = "mm")),
    legend.position = "right",
    legend.justification = "center",
    legend.direction = "vertical",
    legend.text = element_text(
      colour = "black",
      size = 6,
      face = "bold",
      margin = margin(l = 0.6, unit = "mm")
    ),
    legend.key = element_blank(),
    legend.key.width = unit(3.2, "mm"),
    legend.key.height = unit(4.6, "mm"),
    legend.spacing.x = unit(1.2, "mm"),
    legend.box.spacing = unit(2.0, "mm"),
    plot.margin = margin(1.2, 1.2, 1.2, 1.2, unit = "mm")
  )

# Temporary SVG verifies that all text remains editable. It is deleted after
# QA so the delivered figure formats remain PDF and PNG.
svg_temp <- tempfile(pattern = "Fig3c_QA_", fileext = ".svg")
svglite::svglite(
  svg_temp,
  width = width_in,
  height = height_in,
  bg = "white"
)
print(fig3c)
dev.off()

svg_text <- paste(readLines(svg_temp, warn = FALSE), collapse = "\n")
required_text <- c(
  "Recovery breadth", "Non-prior calls",
  annotation_spec$method_key, method_levels
)
missing_text <- required_text[
  !vapply(required_text, grepl, logical(1), x = svg_text, fixed = TRUE)
]
if (length(missing_text) > 0L) {
  stop("Editable-text QA is missing: ", paste(missing_text, collapse = ", "))
}
unlink(svg_temp)

grDevices::cairo_pdf(
  paste0(output_prefix, ".pdf"),
  width = width_in,
  height = height_in,
  family = font_family,
  bg = "white"
)
print(fig3c)
dev.off()

png_temp <- tempfile(pattern = "Fig3c_", fileext = ".png")
ragg::agg_png(
  png_temp,
  width = width_in,
  height = height_in,
  units = "in",
  res = 600,
  background = "white"
)
print(fig3c)
dev.off()
if (!file.copy(png_temp, paste0(output_prefix, ".png"), overwrite = TRUE)) {
  stop("Failed to copy the PNG output into the project directory")
}
unlink(png_temp)

message("Created Fig3c at ", width_mm, " x ", round(height_mm, 2), " mm")
message(
  "Data check: n = ", nrow(plot_data),
  "; x range = ", paste(range(plot_data$True_positive), collapse = " to "),
  "; y range = ", paste(range(plot_data$False_positive), collapse = " to ")
)
message("Highlighted methods:")
for (i in seq_len(nrow(annotation_data))) {
  message(
    "  ", annotation_data$Method[i],
    ": (", annotation_data$True_positive[i],
    ", ", annotation_data$False_positive[i], ")"
  )
}
