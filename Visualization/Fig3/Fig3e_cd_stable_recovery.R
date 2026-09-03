#!/usr/bin/env Rscript

# Approximate code for Fig. 3e: stable recovery of Crohn's disease reference taxa.
# Expected columns: strategy, taxon, recovered. recovered is 0 or 1.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg) == 1L) dirname(normalizePath(sub("^--file=", "", script_arg))) else normalizePath(getwd())
args <- commandArgs(trailingOnly = TRUE)
input_file <- if (length(args) >= 1L) args[[1L]] else file.path(script_dir, "Fig3e_cd_stable_recovery.tsv")
output_file <- if (length(args) >= 2L) args[[2L]] else file.path(script_dir, "Fig3e_cd_stable_recovery.pdf")

dat <- read.delim(input_file, check.names = FALSE, stringsAsFactors = FALSE)
required <- c("strategy", "taxon", "recovered")
if (length(setdiff(required, names(dat))) > 0L) stop("Input requires columns: ", paste(required, collapse = ", "))

scores <- dat %>%
  group_by(strategy) %>%
  summarise(stable_recovery_score = sum(recovered), .groups = "drop") %>%
  arrange(desc(stable_recovery_score), strategy)

strategy_order <- scores$strategy
plot_data <- dat %>%
  mutate(
    strategy = factor(strategy, levels = strategy_order),
    taxon = factor(taxon, levels = rev(unique(taxon)))
  )
scores$strategy <- factor(scores$strategy, levels = strategy_order)

p_score <- ggplot(scores, aes(x = strategy, y = stable_recovery_score)) +
  geom_col(fill = "#8BAFCC", width = 0.78) +
  scale_y_continuous(breaks = 0:4, limits = c(0, 4.2), expand = c(0, 0)) +
  labs(x = NULL, y = "Stable\nrecovery score") +
  theme_classic(base_family = "Arial", base_size = 7) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

p_matrix <- ggplot(plot_data, aes(x = strategy, y = taxon, fill = interaction(taxon, recovered))) +
  geom_tile(colour = "white", linewidth = 0.25) +
  geom_text(aes(label = recovered), size = 2.1, colour = "#4D4D4D") +
  scale_fill_manual(
    values = c(
      "B. adolescentis.0" = "white", "B. adolescentis.1" = "#B7CEE2",
      "F. prausnitzii.0" = "white", "F. prausnitzii.1" = "#B7CEE2",
      "E. coli.0" = "white", "E. coli.1" = "#F2C18E",
      "R. gnavus.0" = "white", "R. gnavus.1" = "#F2C18E"
    ),
    guide = "none"
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_family = "Arial", base_size = 7) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 55, hjust = 1), axis.ticks = element_blank())

figure <- p_score / p_matrix + plot_layout(heights = c(1.3, 1))
ggsave(output_file, figure, width = 180, height = 75, units = "mm", device = cairo_pdf)
