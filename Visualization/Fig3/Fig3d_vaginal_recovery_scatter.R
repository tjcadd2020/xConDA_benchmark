#!/usr/bin/env Rscript

# Recreate Fig. 3d from the data embedded in the original ggplot object.
# The reference image is used for layout and annotation style only.

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

# Relative paths avoid graphics-device encoding problems when the project path
# contains non-ASCII characters under a C locale.
input_file <- "enrichment.RData"
output_prefix <- "Fig3d"

if (!file.exists(input_file)) {
  stop("Missing input file: ", input_file)
}

plot_env <- new.env(parent = emptyenv())
load(input_file, envir = plot_env)

# Use the exact data stored in the original plot whenever available.
if (exists("dot_plot_tp_fp", envir = plot_env, inherits = FALSE)) {
  plot_data <- as.data.frame(plot_env$dot_plot_tp_fp$data)
} else if (exists("df_plot", envir = plot_env, inherits = FALSE)) {
  plot_data <- as.data.frame(plot_env$df_plot)
} else {
  stop("Neither dot_plot_tp_fp nor df_plot was found in enrichment.RData")
}

required_columns <- c(
  "Method", "Method_group", "True_positive", "False_positive"
)
missing_columns <- setdiff(required_columns, names(plot_data))
if (length(missing_columns) > 0L) {
  stop("Missing required columns: ", paste(missing_columns, collapse = ", "))
}

if (nrow(plot_data) != 71L) {
  stop("Expected 71 rows from the original plot; found ", nrow(plot_data))
}
if (any(!is.finite(plot_data$True_positive)) ||
    any(!is.finite(plot_data$False_positive))) {
  stop("The original plotting coordinates contain missing or non-finite values")
}

method_levels <- if (exists("pal", envir = plot_env, inherits = FALSE)) {
  names(plot_env$pal)
} else {
  unique(as.character(plot_data$Method_group))
}

method_palette <- if (exists("pal", envir = plot_env, inherits = FALSE)) {
  plot_env$pal
} else {
  stop("The original method palette was not found in enrichment.RData")
}

missing_colours <- setdiff(unique(as.character(plot_data$Method_group)), names(method_palette))
if (length(missing_colours) > 0L) {
  stop("Missing method colours: ", paste(missing_colours, collapse = ", "))
}

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

# User-requested annotations. Matching is case-insensitive because the stored
# MaAslin2 spelling differs only in capitalization from the display label.
annotation_spec <- data.frame(
  method_key = tolower(c(
    "MaAsLin2_ZINB_count",
    "LFEM_LOG",
    "limma_LOG",
    "LMEM_LOG"
  )),
  display_label = c(
    "MaAsLin2_ZINB_count",
    "LFEM_LOG",
    "limma_LOG",
    "LMEM_LOG"
  ),
  label_x = c(12.2, 12.2, 12.2, 12.2),
  label_y = c(15.8, 31.6, 29.9, 28.2),
  stringsAsFactors = FALSE
)

match_index <- match(annotation_spec$method_key, tolower(plot_data$Method))
if (anyNA(match_index)) {
  stop(
    "Requested annotation method(s) not found: ",
    paste(annotation_spec$display_label[is.na(match_index)], collapse = ", ")
  )
}

annotation_data <- bind_cols(
  plot_data[match_index, , drop = FALSE],
  annotation_spec[, c("display_label", "label_x", "label_y"), drop = FALSE]
)

legend_labels <- setNames(
  sub("^MaAslin2_", "MaAsLin2_", method_levels),
  method_levels
)

font_family <- "Arial"
width_mm <- 100
height_mm <- 58
width_in <- width_mm / 25.4
height_in <- height_mm / 25.4

fig3d <- ggplot(
  plot_data,
  aes(x = display_x, y = display_y)
) +
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
      yend = display_y,
      colour = Method_group
    ),
    inherit.aes = FALSE,
    linewidth = 0.34,
    alpha = 0.78,
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
      label = display_label,
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
    labels = legend_labels,
    drop = FALSE,
    name = NULL
  ) +
  scale_colour_manual(
    values = method_palette,
    guide = "none"
  ) +
  scale_x_continuous(
    breaks = c(0, 5, 10, 15),
    expand = expansion(mult = 0)
  ) +
  scale_y_continuous(
    breaks = c(0, 10, 20, 30),
    expand = expansion(mult = 0)
  ) +
  coord_cartesian(
    xlim = c(-0.6, 19.2),
    ylim = c(-0.6, 33.0),
    clip = "on"
  ) +
  labs(
    x = "Prior-hit count",
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

source_data <- plot_data %>%
  transmute(
    source_row,
    Method,
    Method_group = as.character(Method_group),
    Prior_hit_count = True_positive,
    Non_prior_calls = False_positive,
    display_x,
    display_y,
    highlighted = tolower(Method) %in% annotation_spec$method_key
  )

write.csv(
  source_data,
  paste0(output_prefix, "_source_data.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

svglite::svglite(
  paste0(output_prefix, ".svg"),
  width = width_in,
  height = height_in,
  bg = "white"
)
print(fig3d)
dev.off()

grDevices::cairo_pdf(
  paste0(output_prefix, ".pdf"),
  width = width_in,
  height = height_in,
  family = font_family,
  bg = "white"
)
print(fig3d)
dev.off()

png_temp <- tempfile(pattern = "Fig3d_", fileext = ".png")
ragg::agg_png(
  png_temp,
  width = width_in,
  height = height_in,
  units = "in",
  res = 600,
  background = "white"
)
print(fig3d)
dev.off()
if (!file.copy(png_temp, paste0(output_prefix, ".png"), overwrite = TRUE)) {
  stop("Failed to copy the PNG output into the project directory")
}
unlink(png_temp)

tiff_temp <- tempfile(pattern = "Fig3d_", fileext = ".tiff")
ragg::agg_tiff(
  tiff_temp,
  width = width_in,
  height = height_in,
  units = "in",
  res = 600,
  background = "white"
)
print(fig3d)
dev.off()
if (!file.copy(tiff_temp, paste0(output_prefix, ".tiff"), overwrite = TRUE)) {
  stop("Failed to copy the TIFF output into the project directory")
}
unlink(tiff_temp)

message("Created Fig3d at ", width_mm, " x ", height_mm, " mm")
message(
  "Data check: n = ", nrow(plot_data),
  "; x range = ", paste(range(plot_data$True_positive), collapse = " to "),
  "; y range = ", paste(range(plot_data$False_positive), collapse = " to ")
)
message("Highlighted methods:")
for (i in seq_len(nrow(annotation_data))) {
  message(
    "  ", annotation_data$display_label[i],
    ": (", annotation_data$True_positive[i],
    ", ", annotation_data$False_positive[i], ")"
  )
}
