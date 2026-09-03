#!/usr/bin/env Rscript

# Fig. 2a: performance distributions before and after confounder adjustment
#
# Figure contract
# Core conclusion: show how the full distribution of each performance metric
# changes from unadjusted to adjusted analysis.
# Archetype: quantitative comparison.
# Evidence: paired half-density shapes share a baseline for each metric;
# unadjusted values are above the baseline, adjusted values are below it; exact
# paired Wilcoxon P values are placed directly to the right of each row.
# Backend/output: R only; 70 x 65 mm PDF and 600 dpi PNG. This is the final
# physical size for a subfigure placed at about one third of an A4 page width.

suppressPackageStartupMessages({
  library(ggplot2)
})

required_packages <- c("ragg")
missing_packages <- required_packages[!vapply(
  required_packages,
  requireNamespace,
  logical(1),
  quietly = TRUE
)]
if (length(missing_packages) > 0L) {
  stop(
    "Missing required package(s): ",
    paste(missing_packages, collapse = ", ")
  )
}

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg) == 1L) {
  dirname(normalizePath(sub("^--file=", "", script_arg)))
} else {
  normalizePath(getwd())
}

args <- commandArgs(trailingOnly = TRUE)
rdata_file <- if (length(args) >= 1L) {
  args[[1L]]
} else {
  file.path(script_dir, "comparison_analysis.RData")
}
output_dir <- if (length(args) >= 2L) {
  args[[2L]]
} else {
  script_dir
}

rdata_file <- normalizePath(rdata_file, mustWork = TRUE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(output_dir, mustWork = TRUE)

analysis_env <- new.env(parent = emptyenv())
load(rdata_file, envir = analysis_env)

required_objects <- c("method_cc_performance", "method_ncc_peformance")
missing_objects <- setdiff(required_objects, ls(analysis_env))
if (length(missing_objects) > 0L) {
  stop(
    "The RData file is missing: ",
    paste(missing_objects, collapse = ", ")
  )
}

adjusted_data <- analysis_env$method_cc_performance
unadjusted_data <- analysis_env$method_ncc_peformance

if (!("Method" %in% names(adjusted_data)) ||
    !("Method" %in% names(unadjusted_data))) {
  stop("Column 'Method' is required to apply the method exclusion.")
}

# This exclusion is inherited from the current upstream analysis.
excluded_methods <- c("maaslin2_cplm_LOG_TSS")
adjusted_data <- adjusted_data[
  !(adjusted_data$Method %in% excluded_methods),
  ,
  drop = FALSE
]
unadjusted_data <- unadjusted_data[
  !(unadjusted_data$Method %in% excluded_methods),
  ,
  drop = FALSE
]

if (nrow(adjusted_data) != nrow(unadjusted_data)) {
  stop("Adjusted and unadjusted data have different row counts after exclusion.")
}
if (!identical(rownames(adjusted_data), rownames(unadjusted_data))) {
  stop("Row names differ; align paired rows before calculating distributions.")
}

metrics <- c("FDR", "FPR", "Sensitivity", "AUPR", "AUROC")
missing_columns <- setdiff(
  metrics,
  intersect(names(adjusted_data), names(unadjusted_data))
)
if (length(missing_columns) > 0L) {
  stop(
    "Required metric columns are missing: ",
    paste(missing_columns, collapse = ", ")
  )
}

font_family <- "Helvetica"
group_colours <- c(
  "Unadjusted" = "#CBC7EA",
  "Adjusted" = "#73C7C4"
)
neutral_dark <- "#24282C"
neutral_mid <- "#646B70"
baseline_colour <- "#D9DEE2"
export_dpi <- 600
combined_width_mm <- 70
combined_height_mm <- 65
density_half_height <- 0.23

format_p_plain <- function(p_value, digits = 3L) {
  if (!is.finite(p_value) || p_value < 0) {
    stop("P value must be finite and non-negative.")
  }
  if (p_value == 0) {
    return("P < 2.23e-308")
  }
  if (p_value < 0.01) {
    return(paste0("P = ", formatC(p_value, format = "e", digits = 2L)))
  }
  paste0("P = ", formatC(p_value, format = "f", digits = digits))
}

format_p_plotmath <- function(p_value, digits = 3L) {
  if (!is.finite(p_value) || p_value < 0) {
    stop("P value must be finite and non-negative.")
  }
  if (p_value == 0) {
    return("bolditalic(P) < bold(2.23 %*% 10^{-308})")
  }
  if (p_value < 0.01) {
    exponent <- floor(log10(p_value))
    coefficient <- p_value / (10^exponent)
    return(paste0(
      "bolditalic(P) == bold(",
      formatC(coefficient, digits = digits, format = "fg", flag = "#"),
      " %*% 10^{",
      exponent,
      "})"
    ))
  }
  paste0(
    "bolditalic(P) == bold(",
    formatC(p_value, digits = digits, format = "fg", flag = "#"),
    ")"
  )
}

paired_source <- list()
test_rows <- list()
density_profiles <- list()

for (metric in metrics) {
  adjusted <- adjusted_data[[metric]]
  unadjusted <- unadjusted_data[[metric]]
  complete <- complete.cases(adjusted, unadjusted)

  adjusted <- adjusted[complete]
  unadjusted <- unadjusted[complete]
  pair_id <- rownames(adjusted_data)[complete]
  n_pairs <- length(pair_id)

  if (n_pairs < 2L) {
    stop("Metric '", metric, "' requires at least two complete pairs.")
  }
  if (any(!is.finite(c(adjusted, unadjusted)))) {
    stop("Metric '", metric, "' contains non-finite values.")
  }
  if (any(c(adjusted, unadjusted) < 0 | c(adjusted, unadjusted) > 1)) {
    stop("Metric '", metric, "' contains values outside the expected 0-1 range.")
  }

  p_value <- stats::wilcox.test(
    x = adjusted,
    y = unadjusted,
    paired = TRUE
  )$p.value

  paired_source[[metric]] <- data.frame(
    Pair_ID = pair_id,
    Metric = metric,
    Unadjusted = unadjusted,
    Adjusted = adjusted,
    stringsAsFactors = FALSE
  )

  test_rows[[metric]] <- data.frame(
    Metric = metric,
    N_pairs = n_pairs,
    N_missing_pairs = sum(!complete),
    Wilcoxon_P_value = p_value,
    P_label = format_p_plain(p_value),
    P_plot_label = format_p_plotmath(p_value),
    stringsAsFactors = FALSE
  )

  pooled_values <- c(unadjusted, adjusted)
  shared_bandwidth <- stats::bw.nrd0(pooled_values)
  if (!is.finite(shared_bandwidth) || shared_bandwidth <= 0) {
    shared_bandwidth <- max(stats::sd(pooled_values), 0.01) * 0.9
  }

  density_unadjusted <- stats::density(
    unadjusted,
    bw = shared_bandwidth,
    from = 0,
    to = 1,
    cut = 0,
    n = 512
  )
  density_adjusted <- stats::density(
    adjusted,
    bw = shared_bandwidth,
    from = 0,
    to = 1,
    cut = 0,
    n = 512
  )

  shared_peak <- max(density_unadjusted$y, density_adjusted$y)
  if (!is.finite(shared_peak) || shared_peak <= 0) {
    stop("Could not estimate a finite density for metric '", metric, "'.")
  }

  density_profiles[[paste(metric, "Unadjusted", sep = "_")]] <- data.frame(
    Metric = metric,
    Group = "Unadjusted",
    X = density_unadjusted$x,
    Density = density_unadjusted$y,
    Scaled_density = density_unadjusted$y / shared_peak,
    Bandwidth = shared_bandwidth,
    stringsAsFactors = FALSE
  )
  density_profiles[[paste(metric, "Adjusted", sep = "_")]] <- data.frame(
    Metric = metric,
    Group = "Adjusted",
    X = density_adjusted$x,
    Density = density_adjusted$y,
    Scaled_density = density_adjusted$y / shared_peak,
    Bandwidth = shared_bandwidth,
    stringsAsFactors = FALSE
  )
}

paired_source <- do.call(rbind, paired_source)
rownames(paired_source) <- NULL
test_table <- do.call(rbind, test_rows)
rownames(test_table) <- NULL
test_table <- test_table[match(metrics, test_table$Metric), , drop = FALSE]
density_table <- do.call(rbind, density_profiles)
rownames(density_table) <- NULL

build_density_figure <- function(
  selected_metrics,
  add_panel_tag = FALSE,
  show_title = TRUE,
  show_legend = TRUE
) {
  selected_metrics <- intersect(metrics, selected_metrics)
  if (length(selected_metrics) == 0L) {
    stop("At least one valid metric is required.")
  }

  layout_data <- data.frame(
    Metric = selected_metrics,
    Y = rev(seq_along(selected_metrics)),
    stringsAsFactors = FALSE
  )

  selected_profiles <- density_table[
    density_table$Metric %in% selected_metrics,
    ,
    drop = FALSE
  ]
  selected_profiles$Y <- layout_data$Y[
    match(selected_profiles$Metric, layout_data$Metric)
  ]

  polygon_data <- do.call(rbind, lapply(
    split(selected_profiles, interaction(
      selected_profiles$Metric,
      selected_profiles$Group,
      drop = TRUE
    )),
    function(profile) {
      profile <- profile[order(profile$X), , drop = FALSE]
      direction <- if (profile$Group[[1L]] == "Unadjusted") 1 else -1
      outer_y <- profile$Y +
        direction * density_half_height * profile$Scaled_density
      data.frame(
        Metric = profile$Metric[[1L]],
        Group = profile$Group[[1L]],
        Polygon_ID = paste(profile$Metric[[1L]], profile$Group[[1L]], sep = "_"),
        X = c(profile$X, rev(profile$X)),
        Y = c(rep(profile$Y[[1L]], nrow(profile)), rev(outer_y)),
        stringsAsFactors = FALSE
      )
    }
  ))
  polygon_data$Group <- factor(
    polygon_data$Group,
    levels = c("Unadjusted", "Adjusted")
  )

  selected_tests <- test_table[
    match(selected_metrics, test_table$Metric),
    ,
    drop = FALSE
  ]
  selected_tests$Y <- layout_data$Y

  top_y <- length(selected_metrics) + 1.32
  title_y <- length(selected_metrics) + 1.22
  legend_y <- length(selected_metrics) + 0.69
  y_limits <- c(0.50, top_y)

  legend_rectangles <- data.frame(
    Group = factor(
      c("Unadjusted", "Adjusted"),
      levels = c("Unadjusted", "Adjusted")
    ),
    Xmin = c(0.05, 0.87),
    Xmax = c(0.12, 0.94),
    Ymin = legend_y - 0.078,
    Ymax = legend_y + 0.078
  )
  legend_labels <- data.frame(
    X = c(0.15, 0.97),
    Y = legend_y,
    Label = c("Unadjusted", "Adjusted")
  )

  plot <- ggplot() +
    geom_segment(
      data = layout_data,
      aes(x = 0, xend = 1, y = Y, yend = Y),
      inherit.aes = FALSE,
      linewidth = 0.30,
      colour = baseline_colour
    ) +
    geom_polygon(
      data = polygon_data,
      aes(x = X, y = Y, group = Polygon_ID, fill = Group),
      colour = NA,
      alpha = 0.98
    ) +
    geom_text(
      data = selected_tests,
      aes(x = 1.025, y = Y, label = P_plot_label),
      inherit.aes = FALSE,
      hjust = 0,
      vjust = 0.5,
      family = font_family,
      size = 6 / ggplot2::.pt,
      parse = TRUE,
      colour = neutral_dark
    ) +
    scale_fill_manual(values = group_colours, guide = "none") +
    scale_x_continuous(
      breaks = seq(0, 1, by = 0.2),
      labels = function(x) formatC(x, format = "f", digits = 1L),
      expand = expansion(mult = c(0, 0))
    ) +
    scale_y_continuous(
      breaks = layout_data$Y,
      labels = layout_data$Metric,
      expand = expansion(mult = c(0, 0))
    ) +
    labs(x = "Mean of Metrics Value", y = NULL) +
    coord_cartesian(xlim = c(0, 1), ylim = y_limits, clip = "off") +
    theme_classic(base_family = font_family, base_size = 6.5) +
    theme(
      axis.line = element_line(colour = neutral_dark, linewidth = 0.38),
      axis.line.y = element_blank(),
      axis.ticks = element_line(colour = neutral_dark, linewidth = 0.35),
      axis.ticks.y = element_blank(),
      axis.text.x = element_text(
        size = 6,
        face = "bold",
        colour = neutral_dark
      ),
      axis.text.y = element_text(
        size = 6.5,
        face = "bold",
        colour = neutral_dark,
        margin = margin(r = 1.5, unit = "mm")
      ),
      axis.title.x = element_text(
        size = 7.5,
        face = "bold",
        margin = margin(t = 1.8, unit = "mm")
      ),
      axis.title.y = element_blank(),
      panel.grid = element_blank(),
      plot.margin = margin(0.5, 17.0, 0.5, 0.5, unit = "mm")
    )

  if (show_title) {
    plot <- plot +
      annotate(
        "text",
        x = 0.5,
        y = title_y,
        label = "Effect of confounder adjustment",
        hjust = 0.5,
        vjust = 1,
        family = font_family,
        fontface = "bold",
        size = 8 / ggplot2::.pt,
        colour = neutral_dark
      )
  }

  if (show_legend) {
    plot <- plot +
      geom_rect(
        data = legend_rectangles,
        aes(xmin = Xmin, xmax = Xmax, ymin = Ymin, ymax = Ymax, fill = Group),
        inherit.aes = FALSE,
        colour = NA
      ) +
      geom_text(
        data = legend_labels,
        aes(x = X, y = Y, label = Label),
        inherit.aes = FALSE,
        hjust = 0,
        vjust = 0.5,
        family = font_family,
        size = 6.5 / ggplot2::.pt,
        colour = neutral_dark
      )
  }

  if (add_panel_tag) {
    plot <- plot +
      annotate(
        "text",
        x = -0.28,
        y = title_y,
        label = "a",
        hjust = 0,
        vjust = 1,
        family = font_family,
        fontface = "bold",
        size = 8.5 / ggplot2::.pt,
        colour = neutral_dark
      )
  }

  plot
}

combined_plot <- build_density_figure(
  selected_metrics = metrics,
  add_panel_tag = TRUE,
  show_title = TRUE,
  show_legend = TRUE
)

copy_export_file <- function(temporary_file, destination_file) {
  if (!file.exists(temporary_file) || file.info(temporary_file)$size <= 0) {
    stop("Export device did not create a valid temporary file.")
  }
  copied <- file.copy(temporary_file, destination_file, overwrite = TRUE)
  unlink(temporary_file)
  if (!copied) {
    stop("Could not copy the exported file to: ", destination_file)
  }
}

# Some graphics devices cannot open non-ASCII paths under a C locale. Export
# through an ASCII-safe temporary path, then copy the completed file back.
save_png <- function(plot, filename, width_mm, height_mm) {
  temporary_file <- tempfile(fileext = ".png")
  ggsave(
    filename = temporary_file,
    plot = plot,
    device = ragg::agg_png,
    width = width_mm,
    height = height_mm,
    units = "mm",
    dpi = export_dpi,
    bg = "white"
  )
  copy_export_file(temporary_file, filename)
}

save_pdf <- function(plot, filename, width_mm, height_mm) {
  temporary_file <- tempfile(fileext = ".pdf")
  grDevices::pdf(
    file = temporary_file,
    width = width_mm / 25.4,
    height = height_mm / 25.4,
    family = font_family,
    bg = "white",
    useDingbats = FALSE
  )
  print(plot)
  grDevices::dev.off()
  copy_export_file(temporary_file, filename)
}

combined_prefix <- file.path(
  output_dir,
  "Fig2a_confounder_adjustment_performance_70x65mm"
)
save_pdf(
  combined_plot,
  paste0(combined_prefix, ".pdf"),
  combined_width_mm,
  combined_height_mm
)
save_png(
  combined_plot,
  paste0(combined_prefix, ".png"),
  combined_width_mm,
  combined_height_mm
)

message("Fig. 2a files saved to: ", output_dir)
print(test_table[, c(
  "Metric",
  "N_pairs",
  "N_missing_pairs",
  "Wilcoxon_P_value"
)], row.names = FALSE)
