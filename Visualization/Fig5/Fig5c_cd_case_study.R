#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

# Figure contract:
# - Left: Active CD differential taxa for three requested strategies, split into
#   paired-sibling reference taxa and other differential taxa.
# - Right: the original confounding-signal pie-chart logic for those strategies.
# - DAA calls and paired-sibling reference taxa use adjusted P (FDR) < 0.10.
# - Confounding-signal taxa use non-relative FDR < 0.10 and paired-sibling
#   FDR >= 0.10.

fig_width_mm = 210
fig_height_mm = 58
font_family <- "Arial"
daa_threshold <- 0.10
reference_threshold <- 0.10
daa_threshold_label <- sprintf("Adjusted p-value (FDR) < %g", daa_threshold)
reference_threshold_label <- sprintf(
  "Paired-sibling adjusted p-value (FDR) < %g",
  reference_threshold
)

reference_file <- "paired_sibling_reference_taxa_FDR_lt_0.1.tsv"
confounding_detail_file <- "nonRelative_not_sibling_vs_all_DAA_detail.tsv"
output_base <- "ActiveCD_three_strategy_compact_figure"
overlap_summary_file <- "ActiveCD_three_strategy_reference_overlap_counts.tsv"
overlap_detail_file <- "ActiveCD_three_strategy_reference_overlap_detail.tsv"
pie_counts_file <- "ActiveCD_three_strategy_confounding_signal_pie_counts.tsv"

periods <- "CD-A"
strategies <- c("ANCOM-BC2_CSS", "ANCOM-BC2_count", "MaAsLin2_CPLM_TSS")

daa_files <- list(
  "CD-A" = c(
    "ANCOM-BC2_CSS" = "CD-A_CSS_ANCOM-BC2_DAA_result.tsv",
    "ANCOM-BC2_count" = "CD-A_Count_ANCOM-BC2_DAA_result.tsv",
    "MaAsLin2_CPLM_TSS" = "CD-A_TSS_MaAsLin2-CPLM_DAA_result.tsv"
  )
)

expected_strategy_labels <- list(
  "ANCOM-BC2_CSS" = c(Method = "ANCOM-BC2", Normalization = "CSS"),
  "ANCOM-BC2_count" = c(Method = "ANCOM-BC2", Normalization = "count"),
  "MaAsLin2_CPLM_TSS" = c(Method = "MaAsLin2_CPLM", Normalization = "TSS")
)
expected_covariates <- list("CD-A" = c("age", "BMI", "smoke"))

bar_colors <- c(
  "Paired-sibling reference" = "#527E71",
  "Other differential taxa" = "#A7B8C2"
)

slice_colors <- c(
  "Not differential" = "#E4A11B",
  "Differential" = "#56A6D1"
)

read_tsv <- function(path) {
  if (!file.exists(path)) {
    stop(sprintf("Input file not found: %s", path))
  }
  read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
}

assert_columns <- function(dat, required, path) {
  missing <- setdiff(required, names(dat))
  if (length(missing) > 0) {
    stop(sprintf(
      "Missing required column(s) in %s: %s",
      path,
      paste(missing, collapse = ", ")
    ))
  }
}

read_header_fields <- function(path) {
  if (!file.exists(path)) {
    stop(sprintf("Input file not found: %s", path))
  }
  strsplit(readLines(path, n = 1L, warn = FALSE), "\t", fixed = TRUE)[[1]]
}

read_first_field <- function(path) {
  lines <- readLines(path, warn = FALSE)
  if (length(lines) < 2L) {
    stop(sprintf("Input table has no data rows: %s", path))
  }
  sub("\t.*$", "", lines[-1L])
}

pretty_covariates <- function(x) {
  replacements <- c(age = "age", sex = "sex", BMI = "BMI", smoke = "smoking")
  out <- unname(ifelse(x %in% names(replacements), replacements[x], x))
  if (length(out) == 0L) {
    return("none")
  }
  if (length(out) == 1L) {
    return(out)
  }
  paste0(paste(out[-length(out)], collapse = ", "), ", and ", out[length(out)])
}

extract_covariates <- function(dat) {
  covariates <- sub("^lfc_", "", grep("^lfc_", names(dat), value = TRUE))
  covariates <- setdiff(covariates, c("(Intercept)", "Group", "Group1"))
  covariates
}

# Map full taxonomic lineages in the supplementary tables to the exact
# compact feature names used in the ANCOM-BC2 outputs. The mapping follows the
# preserved feature order in the source abundance table and both analysis inputs.
full_taxonomy <- read_first_field("genus.tsv")
active_taxa <- read_header_fields("genus_active.tsv")[-1L]
remission_taxa <- read_header_fields("genus_remisson.tsv")[-1L]
if (length(full_taxonomy) != length(active_taxa) ||
    length(full_taxonomy) != length(remission_taxa) ||
    !identical(active_taxa, remission_taxa)) {
  stop("The full-lineage-to-ANCOM feature mapping is not one-to-one across inputs.")
}
if (anyDuplicated(full_taxonomy) || anyDuplicated(active_taxa)) {
  stop("Duplicate taxonomy keys prevent an unambiguous reference mapping.")
}
taxonomy_map <- setNames(active_taxa, full_taxonomy)

reference_dat <- read_tsv(reference_file)
assert_columns(
  reference_dat,
  c(
    "Period", "SourceFile", "Sheet", "SourceRow", "Genus_full",
    "FDR_CD_vs_paired_sibling", "Trend_CD_vs_paired_sibling", "Reference_rule"
  ),
  reference_file
)
reference_dat$FDR_CD_vs_paired_sibling <- suppressWarnings(
  as.numeric(reference_dat$FDR_CD_vs_paired_sibling)
)
if (any(!is.finite(reference_dat$FDR_CD_vs_paired_sibling)) ||
    any(reference_dat$FDR_CD_vs_paired_sibling >= reference_threshold)) {
  stop(sprintf("Reference input must contain only taxa with FDR < %g.", reference_threshold))
}
reference_dat$Reference_taxon <- unname(taxonomy_map[reference_dat$Genus_full])
if (anyNA(reference_dat$Reference_taxon) || any(reference_dat$Reference_taxon == "")) {
  stop("At least one paired-sibling reference taxon could not be mapped to the ANCOM inputs.")
}
reference_key <- paste(reference_dat$Period, reference_dat$Reference_taxon, sep = "\r")
if (anyDuplicated(reference_key)) {
  stop("Paired-sibling reference taxa are not unique within period.")
}

overlap_summary <- list()
overlap_detail <- list()
covariates_by_group <- list()

for (period in periods) {
  reference_period <- reference_dat[reference_dat$Period == period, , drop = FALSE]
  reference_set <- reference_period$Reference_taxon

  for (strategy in strategies) {
    path <- daa_files[[period]][[strategy]]
    daa <- read_tsv(path)
    assert_columns(
      daa,
      c("Taxon", "Adjusted p-value", "Method", "Normalization"),
      path
    )
    adjusted_p <- suppressWarnings(as.numeric(daa[["Adjusted p-value"]]))
    differential <- !is.na(adjusted_p) & adjusted_p < daa_threshold
    reference_flag <- daa$Taxon %in% reference_set
    overlap_flag <- differential & reference_flag
    covariates <- extract_covariates(daa)
    if (length(covariates) == 0L) {
      covariates <- expected_covariates[[period]]
    }

    labels <- expected_strategy_labels[[strategy]]
    if (any(daa$Method != labels[["Method"]]) ||
        any(daa$Normalization != labels[["Normalization"]])) {
      stop(sprintf("Unexpected method or normalization label in %s.", path))
    }
    if (!setequal(covariates, expected_covariates[[period]])) {
      stop(sprintf("Unexpected covariates in %s.", path))
    }

    covariates_by_group[[paste(period, strategy, sep = "\r")]] <- covariates
    overlap_count <- sum(overlap_flag)
    total_differential <- sum(differential)

    overlap_summary[[length(overlap_summary) + 1L]] <- data.frame(
      Period = period,
      Strategy = strategy,
      DAA_source_file = path,
      DAA_differential_taxa_count = total_differential,
      Paired_sibling_reference_set_size = length(reference_set),
      Reference_overlap_count = overlap_count,
      Other_differential_taxa_count = total_differential - overlap_count,
      Reference_overlap_rate_among_DAA_differential = if (
        total_differential > 0
      ) overlap_count / total_differential else NA_real_,
      DAA_diff_rule = daa_threshold_label,
      Reference_rule = reference_threshold_label,
      Adjusted_covariates = paste(covariates, collapse = "; "),
      stringsAsFactors = FALSE
    )

    detail <- data.frame(
      Period = period,
      Strategy = strategy,
      DAA_source_file = path,
      Taxon = daa$Taxon,
      DAA_adjusted_p = adjusted_p,
      DAA_differential = differential,
      Paired_sibling_reference = reference_flag,
      Reference_overlap = overlap_flag,
      DAA_diff_rule = daa_threshold_label,
      Reference_rule = reference_threshold_label,
      Adjusted_covariates = paste(covariates, collapse = "; "),
      stringsAsFactors = FALSE
    )
    reference_fdr_map <- setNames(
      reference_period$FDR_CD_vs_paired_sibling,
      reference_period$Reference_taxon
    )
    detail$Reference_FDR_CD_vs_paired_sibling <- unname(reference_fdr_map[detail$Taxon])
    overlap_detail[[length(overlap_detail) + 1L]] <- detail
  }
}

overlap_summary <- do.call(rbind, overlap_summary)
overlap_detail <- do.call(rbind, overlap_detail)

for (period in periods) {
  period_covariates <- lapply(
    strategies,
    function(strategy) covariates_by_group[[paste(period, strategy, sep = "\r")]]
  )
  covariate_keys <- vapply(
    period_covariates,
    function(x) paste(sort(x), collapse = ";"),
    character(1)
  )
  if (length(unique(covariate_keys)) != 1L) {
    stop(sprintf("DAA outputs use different covariates within %s.", period))
  }
}

write.table(
  overlap_summary,
  file = overlap_summary_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = ""
)
write.table(
  overlap_detail,
  file = overlap_detail_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = ""
)

# Recalculate the right-side pies for the three requested strategies using
# the FDR < 0.10 confounding-signal definition regenerated from Tables S6/S7.
confounding_detail <- read_tsv(confounding_detail_file)
assert_columns(
  confounding_detail,
  c("Period", "Method", "Normalization", "Reference_taxon", "DAA_adjusted_p"),
  confounding_detail_file
)

pie_rows <- list()
for (period in periods) {
  for (strategy in strategies) {
    labels <- expected_strategy_labels[[strategy]]
    rows <- confounding_detail[
      confounding_detail$Period == period &
        confounding_detail$Method == labels[["Method"]] &
        confounding_detail$Normalization == labels[["Normalization"]],
      ,
      drop = FALSE
    ]
    if (nrow(rows) == 0L || anyDuplicated(rows$Reference_taxon)) {
      stop(sprintf("Invalid confounding-signal detail rows for %s / %s.", period, strategy))
    }
    adjusted_p <- suppressWarnings(as.numeric(rows$DAA_adjusted_p))
    differential <- !is.na(adjusted_p) & adjusted_p < daa_threshold
    total <- nrow(rows)
    differential_count <- sum(differential)

    pie_rows[[length(pie_rows) + 1L]] <- data.frame(
      Period = period,
      Strategy = strategy,
      Confounding_signal_taxa = total,
      Not_diff_or_filtered = total - differential_count,
      Differential = differential_count,
      Not_diff_or_filtered_rate = (total - differential_count) / total,
      DAA_diff_rule = paste0(
        daa_threshold_label,
        "; reference taxa absent from DAA output are treated as non-diff/filtered"
      ),
      stringsAsFactors = FALSE
    )
  }
}
pie_dat <- do.call(rbind, pie_rows)
write.table(
  pie_dat,
  file = pie_counts_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

overlap_summary$Strategy <- factor(overlap_summary$Strategy, levels = strategies)
pie_dat$Strategy <- factor(pie_dat$Strategy, levels = strategies)

theme_set(
  theme_classic(base_size = 6.5, base_family = font_family) +
    theme(
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      legend.title = element_blank(),
      legend.text = element_text(size = 6, colour = "#333333"),
      plot.title = element_text(size = 8.5, face = "bold", colour = "#222222"),
      plot.margin = margin(1.2, 2.0, 1.2, 1.2, unit = "mm")
    )
)

make_bar <- function(period, panel_title) {
  dat <- overlap_summary[overlap_summary$Period == period, , drop = FALSE]
  dat$Y <- length(strategies) - as.numeric(dat$Strategy) + 1
  dat$Reference_label <- sprintf("%d", dat$Reference_overlap_count)
  dat$Other_label <- sprintf(
    "/%d",
    dat$DAA_differential_taxa_count - dat$Reference_overlap_count
  )

  segment_dat <- rbind(
    data.frame(
      Strategy = dat$Strategy,
      Y = dat$Y,
      Segment = "Paired-sibling reference",
      Start = 0,
      End = dat$Reference_overlap_count
    ),
    data.frame(
      Strategy = dat$Strategy,
      Y = dat$Y,
      Segment = "Other differential taxa",
      Start = dat$Reference_overlap_count,
      End = dat$DAA_differential_taxa_count
    )
  )
  segment_dat$Segment <- factor(
    segment_dat$Segment,
    levels = c("Paired-sibling reference", "Other differential taxa")
  )

  values <- sort(dat$DAA_differential_taxa_count)
  low_value <- values[1]
  high_value <- values[length(values)]
  # A broken axis would distort the visual share of the stacked segments.
  # Keep a continuous scale so segment length encodes the reported proportion.
  # With the three current values (7, 37, and 48), a break would hide the
  # middle strategy and distort its stacked composition, so use a continuous axis.
  use_break <- FALSE
  y_lower <- 0.4
  y_upper <- length(strategies) + 0.6

  if (!use_break) {
    x_limit <- max(5, ceiling((high_value + max(2, high_value * 0.18)) / 5) * 5)
    segment_dat$Display_start <- segment_dat$Start
    segment_dat$Display_end <- segment_dat$End
    dat$Display_value <- dat$DAA_differential_taxa_count
    dat$Display_label_x <- dat$Display_value + x_limit * 0.06

    return(
      ggplot() +
        geom_rect(
          data = segment_dat,
          aes(
            xmin = Display_start,
            xmax = Display_end,
            ymin = Y - 0.28,
            ymax = Y + 0.28,
            fill = Segment
          ),
          colour = NA
        ) +
        geom_hline(
          yintercept = c(y_lower, y_upper),
          linewidth = 0.4,
          colour = "#333333"
        ) +
        geom_vline(
          xintercept = c(0, x_limit),
          linewidth = 0.4,
          colour = "#333333"
        ) +
        geom_text(
          data = dat,
          aes(x = Display_label_x, y = Y, label = Reference_label),
          hjust = 1,
          size = 2.15,
          family = font_family,
          colour = unname(bar_colors[["Paired-sibling reference"]])
        ) +
        geom_text(
          data = dat,
          aes(x = Display_label_x, y = Y, label = Other_label),
          hjust = 0,
          size = 2.15,
          family = font_family,
          colour = "#333333"
        ) +
        scale_fill_manual(values = bar_colors, guide = "none") +
        scale_x_continuous(
          limits = c(0, x_limit),
          breaks = pretty(c(0, x_limit), n = 5),
          expand = expansion(mult = c(0, 0))
        ) +
        scale_y_continuous(
          limits = c(y_lower, y_upper),
          breaks = seq_along(strategies),
          labels = rev(strategies),
          expand = expansion(mult = c(0, 0))
        ) +
        labs(title = panel_title, x = "Number of differential taxa", y = NULL) +
        coord_cartesian(clip = "off") +
        bar_theme()
    )
  }

  break_start <- max(5, ceiling((low_value + 1) / 5) * 5)
  if (low_value <= 2) {
    break_start <- 6
  }
  break_end <- floor((high_value - low_value) / 5) * 5
  x_limit <- ceiling((high_value + max(3, high_value * 0.08)) / 5) * 5
  if (break_end <= break_start || low_value >= break_start || high_value <= break_end) {
    stop(sprintf("Could not construct a valid broken x-axis for %s.", period))
  }

  right_display_span <- break_start
  gap_width <- break_start * 0.05
  gap_end <- break_start + gap_width
  right_scale <- right_display_span / (x_limit - break_end)
  display_limit <- gap_end + right_display_span
  transform_count <- function(x) {
    ifelse(
      x <= break_start,
      x,
      ifelse(
        x >= break_end,
        gap_end + (x - break_end) * right_scale,
        break_start + gap_width / 2
      )
    )
  }

  segment_dat$Display_start <- transform_count(segment_dat$Start)
  segment_dat$Display_end <- transform_count(segment_dat$End)
  dat$Display_value <- transform_count(dat$DAA_differential_taxa_count)
  dat$Display_label_x <- dat$Display_value + display_limit * 0.06

  left_ticks <- if (break_start == 6) c(0, 5) else seq(0, break_start - 5, by = 5)
  right_original_ticks <- seq(break_end + 5, x_limit, by = 5)
  tick_positions <- c(left_ticks, transform_count(right_original_ticks))
  tick_labels <- c(left_ticks, right_original_ticks)
  slash_dx <- max(0.12, gap_width * 0.28)
  slash_dy <- 0.10

  ggplot() +
    geom_rect(
      data = segment_dat,
      aes(
        xmin = Display_start,
        xmax = Display_end,
        ymin = Y - 0.28,
        ymax = Y + 0.28,
        fill = Segment
      ),
      colour = NA
    ) +
    geom_hline(
      yintercept = c(y_lower, y_upper),
      linewidth = 0.4,
      colour = "#333333"
    ) +
    geom_vline(
      xintercept = c(0, display_limit),
      linewidth = 0.4,
      colour = "#333333"
    ) +
    annotate(
      "rect",
      xmin = break_start,
      xmax = gap_end,
      ymin = -Inf,
      ymax = Inf,
      fill = "white",
      colour = NA
    ) +
    annotate(
      "segment",
      x = c(break_start - slash_dx, gap_end - slash_dx),
      xend = c(break_start + slash_dx, gap_end + slash_dx),
      y = y_lower - slash_dy,
      yend = y_lower + slash_dy,
      linewidth = 0.45,
      colour = "#333333"
    ) +
    annotate(
      "segment",
      x = c(break_start - slash_dx, gap_end - slash_dx),
      xend = c(break_start + slash_dx, gap_end + slash_dx),
      y = y_upper - slash_dy,
      yend = y_upper + slash_dy,
      linewidth = 0.45,
      colour = "#333333"
    ) +
    geom_text(
      data = dat,
      aes(x = Display_label_x, y = Y, label = Reference_label),
      hjust = 1,
      size = 2.15,
      family = font_family,
      colour = unname(bar_colors[["Paired-sibling reference"]])
    ) +
    geom_text(
      data = dat,
      aes(x = Display_label_x, y = Y, label = Other_label),
      hjust = 0,
      size = 2.15,
      family = font_family,
      colour = "#333333"
    ) +
    scale_fill_manual(values = bar_colors, guide = "none") +
    scale_x_continuous(
      limits = c(0, display_limit),
      breaks = tick_positions,
      labels = tick_labels,
      expand = expansion(mult = c(0, 0))
    ) +
    scale_y_continuous(
      breaks = seq_along(strategies),
      labels = rev(strategies),
      expand = expansion(mult = c(0, 0))
    ) +
    labs(title = panel_title, x = "Number of differential taxa", y = NULL) +
    coord_cartesian(ylim = c(y_lower, y_upper), clip = "off") +
    bar_theme()
}

bar_theme <- function() {
  theme(
    axis.line = element_blank(),
    axis.ticks = element_line(linewidth = 0.35, colour = "#333333"),
    axis.ticks.length = grid::unit(1.0, "mm"),
    axis.text.x = element_text(
      size = 6.2,
      colour = "#333333",
      margin = margin(t = 0.55, unit = "mm")
    ),
    axis.text.y = element_text(size = 6.5, face = "bold", colour = "#333333"),
    axis.title.x = element_text(
      size = 6.5,
      face = "bold",
      margin = margin(t = 1.2, unit = "mm")
    ),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    plot.title = element_text(
      size = 8.5,
      face = "bold",
      hjust = 0.5,
      margin = margin(b = 0.8, unit = "mm")
    ),
    plot.margin = margin(1.0, 8.0, 0.7, 1.0, unit = "mm")
  )
}

make_pie <- function(period, strategy) {
  row <- pie_dat[
    pie_dat$Period == period & pie_dat$Strategy == strategy,
    ,
    drop = FALSE
  ]
  dat <- data.frame(
    Status = factor(
      c("Not differential", "Differential"),
      levels = c("Not differential", "Differential")
    ),
    Count = c(row$Not_diff_or_filtered, row$Differential),
    stringsAsFactors = FALSE
  )
  dat$Fraction <- dat$Count / sum(dat$Count)
  dat$Label <- ifelse(
    dat$Status == "Not differential",
    ifelse(
      abs(dat$Fraction * 100 - 100) < 0.05,
      "100%",
      sprintf("%.1f%%", dat$Fraction * 100)
    ),
    ""
  )
  dat$Label_position <- sum(dat$Count) - dat$Count / 2

  ggplot(dat, aes(x = 1, y = Count, fill = Status)) +
    geom_col(width = 1, colour = "white", linewidth = 0.45) +
    geom_text(
      aes(y = Label_position, label = Label),
      family = font_family,
      fontface = "bold",
      size = 3.0,
      colour = "white"
    ) +
    coord_polar(theta = "y", start = 0, clip = "off") +
    scale_fill_manual(values = slice_colors, drop = FALSE) +
    labs(title = strategy) +
    guides(fill = "none") +
    theme_void(base_size = 6.5, base_family = font_family) +
    theme(
      plot.title = element_text(
        size = 7,
        face = "bold",
        hjust = 0.5,
        colour = "#2B2B2B",
        margin = margin(b = 0.3, unit = "mm")
      ),
      legend.position = "none",
      plot.margin = margin(0.3, 0.4, 0.2, 0.4, unit = "mm")
    )
}

bar_panel <- make_bar("CD-A", NULL)
bar_panel_compact <- (
  plot_spacer() / bar_panel / plot_spacer()
) + plot_layout(heights = c(0.08, 0.84, 0.08))
pie_panel <- (
  make_pie("CD-A", "ANCOM-BC2_CSS") |
    make_pie("CD-A", "ANCOM-BC2_count") |
    make_pie("CD-A", "MaAsLin2_CPLM_TSS")
) + plot_layout(ncol = 3)
main_panel <- (bar_panel_compact | pie_panel) + plot_layout(widths = c(0.96, 1.74))

bar_legend_dat <- data.frame(
  Segment = factor(names(bar_colors), levels = names(bar_colors)),
  Label = c("Paired-sibling reference taxa", "Other differential taxa"),
  x = c(0.36, 0.62),
  y = 0.70
)
pie_legend_dat <- data.frame(
  Status = factor(names(slice_colors), levels = names(slice_colors)),
  Label = c(
    "Correctly identified as non-differential",
    "Identified as differential"
  ),
  x = c(0.42, 0.70),
  y = 0.25
)

legend_plot <- ggplot() +
  annotate(
    "text", x = 0.08, y = 0.70,
    label = "Differential taxa in bars:",
    hjust = 0, size = 2.15, family = font_family, fontface = "bold",
    colour = "#333333"
  ) +
  geom_point(
    data = bar_legend_dat,
    aes(x = x, y = y, fill = Segment),
    shape = 22, size = 2.9, stroke = 0, show.legend = FALSE
  ) +
  geom_text(
    data = bar_legend_dat,
    aes(x = x + 0.014, y = y, label = Label),
    hjust = 0, size = 2.15, family = font_family, colour = "#333333"
  ) +
  annotate(
    "text", x = 0.08, y = 0.25,
    label = "Proportion of confounding taxa in pies:",
    hjust = 0, size = 2.15, family = font_family, fontface = "bold",
    colour = "#333333"
  ) +
  geom_point(
    data = pie_legend_dat,
    aes(x = x, y = y, fill = Status),
    shape = 22, size = 2.9, stroke = 0, show.legend = FALSE
  ) +
  geom_text(
    data = pie_legend_dat,
    aes(x = x + 0.014, y = y, label = Label),
    hjust = 0, size = 2.15, family = font_family, colour = "#333333"
  ) +
  scale_fill_manual(values = c(bar_colors, slice_colors)) +
  scale_x_continuous(limits = c(0, 1), expand = expansion(mult = 0)) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
  theme_void(base_size = 6.5, base_family = font_family) +
  theme(plot.margin = margin(-0.2, 0, 0, 0, unit = "mm"))

final_plot <- (main_panel / legend_plot) +
  plot_layout(heights = c(1, 0.18))

width_in <- fig_width_mm / 25.4
height_in <- fig_height_mm / 25.4

svglite::svglite(
  paste0(output_base, ".svg"),
  width = width_in,
  height = height_in,
  bg = "white"
)
print(final_plot)
grDevices::dev.off()

grDevices::cairo_pdf(
  paste0(output_base, ".pdf"),
  width = width_in,
  height = height_in,
  family = font_family,
  bg = "white"
)
print(final_plot)
grDevices::dev.off()

render_raster <- function(extension, dpi) {
  temp_output <- tempfile(fileext = paste0(".", extension))
  final_output <- paste0(output_base, ".", extension)
  if (extension == "png") {
    ragg::agg_png(
      temp_output,
      width = width_in,
      height = height_in,
      units = "in",
      res = dpi,
      background = "white"
    )
  } else if (extension == "tiff") {
    ragg::agg_tiff(
      temp_output,
      width = width_in,
      height = height_in,
      units = "in",
      res = dpi,
      compression = "lzw",
      background = "white"
    )
  } else {
    stop(sprintf("Unsupported raster extension: %s", extension))
  }
  print(final_plot)
  grDevices::dev.off()
  copied <- file.copy(temp_output, final_output, overwrite = TRUE)
  unlink(temp_output)
  if (!copied) {
    stop(sprintf("Could not copy rendered raster to %s", final_output))
  }
}

render_raster("png", dpi = 300)
render_raster("tiff", dpi = 600)

print(overlap_summary[, c(
  "Period", "Strategy", "DAA_differential_taxa_count",
  "Reference_overlap_count", "Reference_overlap_rate_among_DAA_differential",
  "Adjusted_covariates"
)])
print(pie_dat[, c(
  "Period", "Strategy", "Confounding_signal_taxa",
  "Not_diff_or_filtered", "Differential", "Not_diff_or_filtered_rate"
)])
message(sprintf(
  "Saved %s.[pdf|svg|png|tiff] at %.0f x %.0f mm.",
  output_base,
  fig_width_mm,
  fig_height_mm
))
