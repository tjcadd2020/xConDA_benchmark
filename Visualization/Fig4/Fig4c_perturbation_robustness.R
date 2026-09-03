#!/usr/bin/env Rscript

# Approximate code for Fig. 4c: integrated robustness across perturbation factors.
# Expected columns: strategy, factor, topsis.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg) == 1L) dirname(normalizePath(sub("^--file=", "", script_arg))) else normalizePath(getwd())
args <- commandArgs(trailingOnly = TRUE)
input_file <- if (length(args) >= 1L) args[[1L]] else file.path(script_dir, "Fig4c_perturbation_topsis.tsv")
output_file <- if (length(args) >= 2L) args[[2L]] else file.path(script_dir, "Fig4c_perturbation_robustness.pdf")

dat <- read.delim(input_file, check.names = FALSE, stringsAsFactors = FALSE)
required <- c("strategy", "factor", "topsis")
if (length(setdiff(required, names(dat))) > 0L) stop("Input requires columns: ", paste(required, collapse = ", "))

factor_order <- c(
  "Sample size", "Case-control ratio", "Phenotype effect size",
  "Confounding effect size", "Confounder complexity", "Feature prevalence"
)
strategy_order <- dat %>%
  group_by(strategy) %>%
  summarise(mean_topsis = mean(topsis, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(mean_topsis)) %>%
  pull(strategy)

plot_data <- dat %>%
  mutate(
    factor = factor(factor, levels = factor_order),
    strategy = factor(strategy, levels = strategy_order)
  )

p <- ggplot(plot_data, aes(x = strategy, y = topsis)) +
  geom_col(fill = "#6B2737", width = 0.72) +
  facet_wrap(~ factor, ncol = 2) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.03))) +
  labs(x = NULL, y = "TOPSIS score") +
  theme_classic(base_family = "Arial", base_size = 7) +
  theme(
    strip.background = element_rect(fill = "#F2F2F2", colour = NA),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5.5),
    panel.border = element_rect(fill = NA, colour = "black"),
    axis.line = element_blank()
  )

ggsave(output_file, p, width = 120, height = 120, units = "mm", device = cairo_pdf)
