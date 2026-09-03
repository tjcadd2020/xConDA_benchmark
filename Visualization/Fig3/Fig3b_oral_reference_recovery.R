#!/usr/bin/env Rscript

# Approximate code for Fig. 3b: oral reference-pattern recovery.
# Expected columns: strategy, direction, recovery_breadth, enrichment_p, rres.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg) == 1L) dirname(normalizePath(sub("^--file=", "", script_arg))) else normalizePath(getwd())
args <- commandArgs(trailingOnly = TRUE)
input_file <- if (length(args) >= 1L) args[[1L]] else file.path(script_dir, "Fig3b_oral_reference_recovery.tsv")
output_file <- if (length(args) >= 2L) args[[2L]] else file.path(script_dir, "Fig3b_oral_reference_recovery.pdf")

dat <- read.delim(input_file, check.names = FALSE, stringsAsFactors = FALSE)
required <- c("strategy", "direction", "recovery_breadth", "enrichment_p", "rres")
if (length(setdiff(required, names(dat))) > 0L) stop("Input requires columns: ", paste(required, collapse = ", "))

strategy_order <- dat %>%
  group_by(strategy) %>%
  summarise(rres = first(rres), .groups = "drop") %>%
  arrange(desc(rres)) %>%
  pull(strategy)

plot_data <- dat %>%
  mutate(
    strategy = factor(strategy, levels = strategy_order),
    signed_breadth = if_else(direction == "Subgingival", -recovery_breadth, recovery_breadth),
    significance = case_when(
      enrichment_p < 0.001 ~ "***",
      enrichment_p < 0.01 ~ "**",
      enrichment_p < 0.05 ~ "*",
      enrichment_p < 0.1 ~ "#",
      TRUE ~ ""
    )
  )

direction_colors <- c("Supragingival" = "#C54434", "Subgingival" = "#3486BF")

p_bar <- ggplot(plot_data, aes(x = strategy, y = signed_breadth, fill = direction)) +
  geom_col(width = 0.78) +
  geom_text(
    aes(label = significance, y = signed_breadth + sign(signed_breadth) * 0.6),
    size = 2.5,
    show.legend = FALSE
  ) +
  geom_hline(yintercept = 0, linewidth = 0.45) +
  scale_fill_manual(values = direction_colors) +
  scale_y_continuous(labels = abs, expand = expansion(mult = c(0.08, 0.08))) +
  labs(x = NULL, y = "Recovery breadth", fill = NULL) +
  theme_classic(base_family = "Arial", base_size = 7) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "bottom")

rres_data <- plot_data %>% distinct(strategy, rres)
p_rres <- ggplot(rres_data, aes(x = strategy, y = "RRES", fill = rres)) +
  geom_tile() +
  scale_fill_gradient(low = "#F5F2F0", high = "#4A3B32", name = "RRES") +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_family = "Arial", base_size = 7) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 55, hjust = 1),
    axis.text.y = element_text(face = "bold"),
    legend.position = "bottom"
  )

figure <- p_bar / p_rres + plot_layout(heights = c(5, 1), guides = "collect") & theme(legend.position = "bottom")
ggsave(output_file, figure, width = 180, height = 95, units = "mm", device = cairo_pdf)
