# Add left-aligned Overall Score labels to the saved funky heatmap.
# The source data and the original plotting script are read only.

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(funkyheatmap)
  library(svglite)
})

data_file <- "summary_plot_2608.RData"
output_prefix <- "summary_performance_with_scores"

if (!file.exists(data_file)) {
  stop("Cannot find ", data_file)
}

load(data_file)

required_objects <- c(
  "merged_all_performance", "column_info", "column_groups",
  "palettes", "legends"
)
missing_objects <- setdiff(required_objects, ls())
if (length(missing_objects) > 0) {
  stop(
    "Missing required object(s) in ", data_file, ": ",
    paste(missing_objects, collapse = ", ")
  )
}

# At the 200 mm target width, give the two long text columns enough physical
# space to retain the 6 pt font floor without overlaps.
column_info$width[column_info$id == "id"] <- 8.0
column_info$width[column_info$id == "DAA_Method_type"] <- 8.5

p <- funky_heatmap(
  merged_all_performance,
  column_info = column_info,
  column_groups = column_groups,
  palettes = palettes,
  legends = legends,
  scale_column = TRUE,
  add_abc = FALSE,
  position_args = position_arguments(
    col_space = 0.1,
    col_bigspace = 0.6,
    col_annot_angle = 90,
    col_annot_offset = 1.5,
    expand_ymax = 2.5
  )
)

if (!inherits(p, "patchwork")) {
  stop("The rebuilt figure is not the expected patchwork object.")
}

main_plot <- p[[1]]

# Reproduce the label organization in the reference layout:
# - the first two column labels remain horizontal above the rows;
# - the Overall Score column label is represented by its purple group header;
# - all metric labels are vertical below the rows;
# - group titles use the requested wording and line breaks.
text_layer_index <- which(vapply(
  main_plot$layers,
  function(layer) {
    layer_data <- layer$data
    is.data.frame(layer_data) &&
      all(c(
        "row_id", "column_id", "label_value", "x", "y",
        "angle", "hjust", "vjust", "size", "fontface"
      ) %in% names(layer_data)) &&
      any(layer_data$label_value == "Strategy Information", na.rm = TRUE) &&
      any(layer_data$label_value == "Strategy", na.rm = TRUE)
  },
  logical(1)
))

if (length(text_layer_index) != 1) {
  stop("Could not uniquely locate the figure annotation layer.")
}

text_data <- main_plot$layers[[text_layer_index]]$data

group_label_map <- c(
  "Strategy Information" = "Strategy Information",
  "Overall" = "Overall\nScore",
  "Simulation-based Accuracy" = "Simulation-based Accuracy",
  "Reference-based Recovery" = "Empirical\nReference\nRecovery",
  "Consistency" = "Consis-\ntency",
  "Stability Under Factor Perturbation" = "Perturbation Stability"
)

for (old_label in names(group_label_map)) {
  text_data$label_value[
    !is.na(text_data$label_value) & text_data$label_value == old_label
  ] <- unname(group_label_map[[old_label]])
}

group_label_rows <-
  is.na(text_data$row_id) &
  is.na(text_data$column_id) &
  !is.na(text_data$label_value) &
  abs(text_data$y - 2) < .Machine$double.eps^0.5

text_data$size[group_label_rows] <- 2.15
text_data$lineheight[group_label_rows] <- 0.78

column_label_rows <-
  is.na(text_data$row_id) &
  is.na(text_data$column_id) &
  !is.na(text_data$label_value) &
  abs(text_data$y) < .Machine$double.eps^0.5

original_column_labels <- text_data$label_value[column_label_rows]
column_width_lookup <- setNames(column_info$width, column_info$name)
label_widths <- unname(column_width_lookup[original_column_labels])

if (anyNA(label_widths)) {
  stop("Could not match every displayed column label to column_info.")
}

# funkyheatmap stores the annotation x position at the right edge of each
# column. Subtracting half of its width centers every label under its column.
text_data$x[column_label_rows] <-
  text_data$x[column_label_rows] - label_widths / 2

top_label_rows <- column_label_rows &
  text_data$label_value %in% c("Strategy", "DAAType")
score_label_row <- column_label_rows &
  text_data$label_value == "Overall Score"
metric_label_rows <- column_label_rows &
  !top_label_rows &
  !score_label_row

text_data$label_value[
  top_label_rows & text_data$label_value == "DAAType"
] <- "DAA Strategy Type"

text_data$y[top_label_rows] <- 0.35
text_data$angle[top_label_rows] <- 0
text_data$hjust[top_label_rows] <- 0.5
text_data$vjust[top_label_rows] <- 0.5
text_data$fontface[top_label_rows] <- "bold"
text_data$size[top_label_rows] <- 2.4

text_data$label_value[score_label_row] <- ""

metric_label_map <- c(
  "FPR" = "FPR",
  "FDR" = "FDR",
  "Sensitivity" = "Sensitivity",
  "AUROC" = "AUROC",
  "AUPR" = "AUPR",
  "MCC" = "MCC",
  "Macro-F1" = "Macro-F1",
  "TOPSIS" = "TOPSIS",
  "Dental plaque scenario" = "Oral",
  "Vaginal scenario" = "Vaginal",
  "Crohn Disease" = "CD",
  "Mean Jaccard Index" = "Mean Jaccard Index",
  "Mean DA count" = "Mean DA count",
  "Sample size" = "Sample size",
  "Case-control ratio" = "Case-control ratio",
  "Primary effect size" = "Phenotype effect size",
  "Confounding effect size" = "Confounding effect size",
  "Confounder complexity" = "Confounder complexity",
  "Feature prevalence" = "Feature prevalence"
)

for (old_label in names(metric_label_map)) {
  text_data$label_value[
    metric_label_rows & text_data$label_value == old_label
  ] <- unname(metric_label_map[[old_label]])
}

text_data$y[metric_label_rows] <- -14.00
text_data$angle[metric_label_rows] <- 90
text_data$hjust[metric_label_rows] <- 1
text_data$vjust[metric_label_rows] <- 0.5
text_data$fontface[metric_label_rows] <- "bold"
text_data$size[metric_label_rows] <- 2.15

method_text_rows <-
  !is.na(text_data$row_id) & text_data$column_id == "id"
type_text_rows <-
  !is.na(text_data$row_id) & text_data$column_id == "DAA_Method_type"

text_data$size[method_text_rows] <- 2.2
text_data$fontface[method_text_rows] <- "bold"
text_data$size[type_text_rows] <- 2.15
text_data$fontface[type_text_rows] <- "italic"

main_plot$layers[[text_layer_index]]$data <- text_data
main_plot$layers[[text_layer_index]]$aes_params$family <- "Arial"

# Reserve vertical space for the labels below the last method row.
blank_layer_index <- which(vapply(
  main_plot$layers,
  function(layer) {
    layer_data <- layer$data
    inherits(layer$geom, "GeomBlank") &&
      is.data.frame(layer_data) &&
      all(c("x", "y") %in% names(layer_data))
  },
  logical(1)
))

if (length(blank_layer_index) != 1) {
  stop("Could not uniquely locate the plot-boundary layer.")
}

blank_data <- main_plot$layers[[blank_layer_index]]$data
blank_data$y[which.min(blank_data$y)] <- -21.5
main_plot$layers[[blank_layer_index]]$data <- blank_data

# Locate the existing Overall Score bar geometry rather than hard-coding row
# positions. label_value contains the unmodified TOPSIS_score used by the bars.
bar_layer_index <- which(vapply(
  main_plot$layers,
  function(layer) {
    layer_data <- layer$data
    is.data.frame(layer_data) &&
      all(c(
        "column_id", "row_id", "xmin", "ymin", "ymax", "label_value"
      ) %in% names(layer_data)) &&
      any(layer_data$column_id == "TOPSIS_score", na.rm = TRUE)
  },
  logical(1)
))

if (length(bar_layer_index) != 1) {
  stop("Could not uniquely locate the Overall Score bar layer.")
}

bar_data <- main_plot$layers[[bar_layer_index]]$data

group_box_rows <-
  is.na(bar_data$row_id) &
  is.na(bar_data$column_id) &
  !is.na(bar_data$ymin) &
  !is.na(bar_data$ymax)

bar_data$ymin[group_box_rows] <- 0.8
bar_data$ymax[group_box_rows] <- 3.2
main_plot$layers[[bar_layer_index]]$data <- bar_data

score_labels <- bar_data[
  !is.na(bar_data$column_id) & bar_data$column_id == "TOPSIS_score",
  c("row_id", "xmin", "ymin", "ymax", "label_value")
]

if (nrow(score_labels) != nrow(merged_all_performance)) {
  stop("The number of score labels does not match the number of methods.")
}

# Verify that the displayed values are identical to the saved data values.
saved_scores <- setNames(
  merged_all_performance$TOPSIS_score,
  merged_all_performance$id
)
matched_scores <- unname(saved_scores[score_labels$row_id])
if (
  anyNA(matched_scores) ||
    !isTRUE(all.equal(
      score_labels$label_value,
      matched_scores,
      tolerance = .Machine$double.eps^0.5,
      check.attributes = FALSE
    ))
) {
  stop("Overall Score values in the plot do not match the saved data.")
}

score_labels$x <- score_labels$xmin + 0.12
score_labels$y <- (score_labels$ymin + score_labels$ymax) / 2
score_labels$label <- formatC(
  score_labels$label_value,
  format = "f",
  digits = 3
)

main_plot_with_scores <- main_plot +
  geom_text(
    data = score_labels,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 0.5,
    colour = "#666666",
    family = "Arial",
    fontface = "bold",
    size = 2.15
  )

p_with_scores <- p
p_with_scores[[1]] <- main_plot_with_scores

# Keep the legend readable at the narrower final width.
legend_plot <- p_with_scores[[3]]
legend_text_layer_index <- which(vapply(
  legend_plot$layers,
  function(layer) inherits(layer$geom, "GeomText"),
  logical(1)
))

if (length(legend_text_layer_index) != 1) {
  stop("Could not uniquely locate the legend text layer.")
}

legend_text_data <-
  legend_plot$layers[[legend_text_layer_index]]$data
legend_text_data$size <- ifelse(
  legend_text_data$label_value == "Scaled score",
  2.4,
  2.15
)
legend_plot$layers[[legend_text_layer_index]]$data <- legend_text_data
legend_plot$layers[[legend_text_layer_index]]$aes_params$family <- "Arial"
p_with_scores[[3]] <- legend_plot

# A4 is 210 mm wide. Use a slightly shorter 200 mm canvas and preserve the
# established panel aspect ratio.
output_width_mm <- 200
output_height_mm <- output_width_mm * 608 / 1053
output_width_in <- output_width_mm / 25.4
output_height_in <- output_height_mm / 25.4

svglite::svglite(
  paste0(output_prefix, ".svg"),
  width = output_width_in,
  height = output_height_in,
  bg = "white"
)
print(p_with_scores)
grDevices::dev.off()

grDevices::cairo_pdf(
  paste0(output_prefix, ".pdf"),
  width = output_width_in,
  height = output_height_in,
  family = "Arial",
  bg = "white"
)
print(p_with_scores)
grDevices::dev.off()

grDevices::png(
  paste0(output_prefix, ".png"),
  width = output_width_in,
  height = output_height_in,
  units = "in",
  res = 300,
  type = "cairo-png",
  bg = "white"
)
print(p_with_scores)
grDevices::dev.off()

message(
  "Created: ",
  paste(
    paste0(output_prefix, c(".svg", ".pdf", ".png")),
    collapse = ", "
  )
)
