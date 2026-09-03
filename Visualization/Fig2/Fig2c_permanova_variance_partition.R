#!/usr/bin/env Rscript

# Approximate code for Fig. 2c: PERMANOVA variance partitioning.
# Expected columns: group, factor, R2, p_value.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg) == 1L) {
  dirname(normalizePath(sub("^--file=", "", script_arg)))
} else {
  normalizePath(getwd())
}

args <- commandArgs(trailingOnly = TRUE)
input_file <- if (length(args) >= 1L) args[[1L]] else file.path(script_dir, "Fig2c_permanova_results.tsv")
output_file <- if (length(args) >= 2L) args[[2L]] else file.path(script_dir, "Fig2c_permanova_variance_partition.pdf")

dat <- read.delim(input_file, check.names = FALSE, stringsAsFactors = FALSE)
required <- c("group", "factor", "R2", "p_value")
if (length(setdiff(required, names(dat))) > 0L) {
  stop("Input requires columns: ", paste(required, collapse = ", "))
}

group_order <- c("Scenario-level factor", "Strategy component")
group_colors <- c(
  "Scenario-level factor" = "#73A4CA",
  "Strategy component" = "#BB777A"
)

plot_data <- dat %>%
  mutate(
    group = factor(group, levels = group_order),
    factor = factor(factor, levels = rev(unique(factor))),
    significance = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ ""
    )
  )

p <- ggplot(plot_data, aes(x = R2, y = factor, fill = group)) +
  geom_col(width = 0.72) +
  geom_text(aes(label = significance), hjust = 1.2, size = 3.2, fontface = "bold") +
  facet_grid(group ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(values = group_colors, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(x = expression(R^2), y = NULL) +
  theme_classic(base_family = "Arial", base_size = 9) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.border = element_rect(fill = NA, colour = "black"),
    panel.spacing.y = grid::unit(3, "mm")
  )

ggsave(output_file, p, width = 90, height = 110, units = "mm", device = cairo_pdf)
