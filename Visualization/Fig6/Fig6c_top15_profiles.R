#!/usr/bin/env Rscript

# Generate the Fig. 6c Top-15 backup and run an independent audit against:
# 1. the rank values visible in the supplied reference figure;
# 2. ranks recalculated from the unranked simulation metrics;
# 3. the two-dimensional TOPSIS score recalculated from its two inputs; and
# 4. cross-disease means recalculated from all four disease-level observations.

script_argument <- grep(
  "^--file=",
  commandArgs(trailingOnly = FALSE),
  value = TRUE
)
script_directory <- if (length(script_argument) == 1L) {
  dirname(normalizePath(sub("^--file=", "", script_argument)))
} else {
  normalizePath(getwd())
}

Sys.setenv(FIG6C_TOP_N = "15")
source(
  file.path(script_directory, "Fig6c_core.R"),
  local = environment()
)

rank_columns <- c(
  "FPR",
  "FDR",
  "Sensitivity",
  "AUC",
  "AUPR",
  "Macro.F1",
  "MCC",
  "TOPSIS_score"
)

# Values transcribed directly from the first 15 visible rows of the supplied
# reference figure. Column order follows the heatmap from left to right.
reference_rank_matrix <- data.frame(
  Method = c(
    "ANCOM-BC2_TMM",
    "ANCOM-BC2_count",
    "ANCOM-BC2_CSS",
    "edgeR_CSS",
    "ANCOM-BC2_TMM_Fisher",
    "ANCOM-BC2_count_Fisher",
    "ANCOM-BC2_CSS_Fisher",
    "ANCOM-BC2_TMM_Stouffer",
    "ANCOM-BC2_count_Stouffer",
    "ANCOM-BC2_CSS_Stouffer",
    "edgeR_TMM",
    "MaAslin2_NEGBIN_TMM_PM",
    "MaAslin2_NEGBIN_TMM_EB",
    "edgeR_CSS_Stouffer",
    "MaAslin2_NEGBIN_TMM_REML"
  ),
  FPR = c(37, 36, 1, 7, 60, 58, 9, 53, 52, 5, 47, 23, 23, 6, 22),
  FDR = c(35, 34, 1, 7, 60, 59, 8, 51, 49, 5, 40, 22, 22, 6, 24),
  Sensitivity = c(2, 2, 19, 6, 14, 18, 44, 14, 14, 47, 4, 31, 31, 62, 37),
  AUC = c(3, 2, 1, 4, 46, 47, 45, 42, 43, 41, 5, 23, 23, 49, 27),
  AUPR = c(3, 2, 1, 4, 24, 25, 29, 17, 16, 20, 5, 45, 46, 31, 48),
  Macro.F1 = c(35, 34, 1, 2, 60, 59, 8, 51, 49, 6, 39, 19, 19, 22, 23),
  MCC = c(35, 34, 1, 2, 59, 57, 8, 48, 46, 6, 38, 20, 20, 17, 23),
  TOPSIS_score = c(17, 16, 1, 2, 48, 47, 18, 43, 41, 15, 33, 26, 27, 22, 32),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

observed_rank_matrix <- top_rank_matrix %>%
  mutate(Method = as.character(Method)) %>%
  arrange(match(Method, top_methods$Method))

if (!identical(reference_rank_matrix$Method, observed_rank_matrix$Method)) {
  stop("Top-15 method order does not match the supplied reference figure.")
}

reference_crosscheck <- reference_rank_matrix
all_reference_cells_match <- rep(TRUE, nrow(reference_crosscheck))
for (rank_column in rank_columns) {
  observed_column <- paste0(rank_column, "_observed")
  reference_column <- paste0(rank_column, "_reference")
  match_column <- paste0(rank_column, "_match")
  reference_crosscheck[[reference_column]] <- reference_rank_matrix[[rank_column]]
  reference_crosscheck[[observed_column]] <- observed_rank_matrix[[rank_column]]
  reference_crosscheck[[match_column]] <-
    reference_crosscheck[[reference_column]] ==
    reference_crosscheck[[observed_column]]
  all_reference_cells_match <-
    all_reference_cells_match & reference_crosscheck[[match_column]]
  reference_crosscheck[[rank_column]] <- NULL
}
reference_crosscheck$All_rank_cells_match <- all_reference_cells_match

if (!all(reference_crosscheck$All_rank_cells_match)) {
  stop("At least one Top-15 rank cell differs from the supplied reference figure.")
}

recalculate_ranks <- function(metrics) {
  output <- metrics
  for (metric in setdiff(rank_columns, c("FDR", "FPR"))) {
    output[[metric]] <- rank(
      -metrics[[metric]],
      ties.method = "min",
      na.last = "keep"
    )
  }
  for (metric in c("FDR", "FPR")) {
    output[[metric]] <- rank(
      metrics[[metric]],
      ties.method = "min",
      na.last = "keep"
    )
  }
  output
}

simulation_metrics <- as.data.frame(
  figure_environment$evaluation_res_heatmap,
  stringsAsFactors = FALSE
)
simulation_metrics$Method <- as.character(simulation_metrics$Method)
recalculated_ranks <- recalculate_ranks(
  simulation_metrics[, c("Method", rank_columns), drop = FALSE]
)
stored_ranks <- as.data.frame(
  figure_environment$evaluation_res_heatmap_rank,
  stringsAsFactors = FALSE
)
stored_ranks$Method <- as.character(stored_ranks$Method)
stored_ranks <- stored_ranks[, c("Method", rank_columns), drop = FALSE]
recalculated_ranks <- recalculated_ranks[
  match(stored_ranks$Method, recalculated_ranks$Method),
  ,
  drop = FALSE
]

rank_recalculation_match <-
  identical(stored_ranks$Method, recalculated_ranks$Method) &&
  all(
    as.matrix(stored_ranks[, rank_columns]) ==
      as.matrix(recalculated_ranks[, rank_columns])
  )
if (!rank_recalculation_match) {
  stop("Stored simulation ranks differ from ranks recalculated from raw metrics.")
}

recalculate_topsis <- function(data, metric_columns) {
  metric_matrix <- as.matrix(data[, metric_columns, drop = FALSE])
  column_minimum <- apply(metric_matrix, 2, min)
  column_maximum <- apply(metric_matrix, 2, max)
  column_range <- column_maximum - column_minimum
  if (any(column_range == 0)) {
    stop("Cannot recalculate TOPSIS with a constant metric column.")
  }
  normalized <- sweep(metric_matrix, 2, column_minimum, "-")
  normalized <- sweep(normalized, 2, column_range, "/")
  weighted <- normalized / length(metric_columns)
  positive_ideal <- apply(weighted, 2, max)
  negative_ideal <- apply(weighted, 2, min)
  positive_distance <- sqrt(rowSums(
    (weighted - matrix(
      positive_ideal,
      nrow = nrow(weighted),
      ncol = ncol(weighted),
      byrow = TRUE
    ))^2
  ))
  negative_distance <- sqrt(rowSums(
    (weighted - matrix(
      negative_ideal,
      nrow = nrow(weighted),
      ncol = ncol(weighted),
      byrow = TRUE
    ))^2
  ))
  negative_distance / (positive_distance + negative_distance)
}

topsis_input <- integrated_ranking[, c(
  "Method",
  "simulation_TOPSIS",
  "Mean_Jaccard_index",
  "TOPSIS_score"
)]
topsis_recalculated <- recalculate_topsis(
  topsis_input,
  c("simulation_TOPSIS", "Mean_Jaccard_index")
)
topsis_max_absolute_difference <- max(
  abs(topsis_input$TOPSIS_score - topsis_recalculated)
)
topsis_recalculation_match <- topsis_max_absolute_difference < 1e-12
if (!topsis_recalculation_match) {
  stop("Stored integrated TOPSIS differs from the independently recalculated score.")
}

recalculated_consistency <- per_disease_consistency %>%
  mutate(Strategy = as.character(Strategy)) %>%
  group_by(Strategy) %>%
  summarise(
    Mean_Jaccard_index_recalculated = mean(Mean_Jaccard_index),
    Mean_DA_count_recalculated = mean(Mean_DA_count),
    N_disease_recalculated = n(),
    .groups = "drop"
  )
stored_consistency <- integrated_consistency %>%
  mutate(Strategy = as.character(Strategy)) %>%
  select(Strategy, Mean_Jaccard_index, Mean_DA_count, N_disease)
consistency_crosscheck <- stored_consistency %>%
  inner_join(recalculated_consistency, by = "Strategy") %>%
  mutate(
    Jaccard_absolute_difference = abs(
      Mean_Jaccard_index - Mean_Jaccard_index_recalculated
    ),
    DA_count_absolute_difference = abs(
      Mean_DA_count - Mean_DA_count_recalculated
    ),
    N_disease_match = N_disease == N_disease_recalculated,
    All_values_match =
      Jaccard_absolute_difference < 1e-12 &
      DA_count_absolute_difference < 1e-12
  )
consistency_recalculation_match <- all(consistency_crosscheck$All_values_match)
n_disease_metadata_match <- all(consistency_crosscheck$N_disease_match)
if (!consistency_recalculation_match) {
  stop("Stored cross-disease means differ from means recalculated across diseases.")
}

reference_crosscheck_file <- file.path(
  output_directory,
  "Fig6c_top15_reference_figure_crosscheck.tsv"
)
write.table(
  reference_crosscheck,
  file = reference_crosscheck_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

consistency_crosscheck_file <- file.path(
  output_directory,
  "Fig6c_top15_consistency_recalculation_crosscheck.tsv"
)
write.table(
  consistency_crosscheck %>%
    filter(Strategy %in% top_methods$Method) %>%
    arrange(match(Strategy, top_methods$Method)),
  file = consistency_crosscheck_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

audit_summary <- c(
  "Fig. 6c Top-15 independent data audit",
  paste0(
    "Reference figure method order: PASS (",
    nrow(reference_rank_matrix),
    "/",
    nrow(reference_rank_matrix),
    " rows)"
  ),
  paste0(
    "Reference figure rank cells: PASS (",
    nrow(reference_rank_matrix) * length(rank_columns),
    "/",
    nrow(reference_rank_matrix) * length(rank_columns),
    " cells)"
  ),
  paste0(
    "Simulation ranks recalculated from raw metrics: ",
    ifelse(rank_recalculation_match, "PASS", "FAIL")
  ),
  paste0(
    "Integrated TOPSIS independently recalculated: ",
    ifelse(topsis_recalculation_match, "PASS", "FAIL")
  ),
  paste0(
    "Maximum absolute TOPSIS difference: ",
    format(topsis_max_absolute_difference, scientific = TRUE)
  ),
  paste0(
    "Cross-disease means independently recalculated: ",
    ifelse(consistency_recalculation_match, "PASS", "FAIL")
  ),
  paste0(
    "Maximum absolute Jaccard-mean difference: ",
    format(max(consistency_crosscheck$Jaccard_absolute_difference), scientific = TRUE)
  ),
  paste0(
    "Maximum absolute DA-count-mean difference: ",
    format(max(consistency_crosscheck$DA_count_absolute_difference), scientific = TRUE)
  ),
  paste0(
    "N_disease metadata: ",
    ifelse(n_disease_metadata_match, "PASS", "MISMATCH")
  ),
  paste(
    "N_disease is stored as 1 while four disease rows are present per strategy;",
    "this metadata field is not used in the figure or TOPSIS calculation."
  ),
  "Strategy types read directly from Integrate_info$Integrate_type.",
  "No values were inferred from method-name suffixes."
)
writeLines(
  audit_summary,
  con = file.path(output_directory, "Fig6c_top15_data_audit.txt")
)

qa_file <- file.path(output_directory, "Fig6c_top15_QA.txt")
write(
  c("", audit_summary),
  file = qa_file,
  append = TRUE
)

message(
  "Top-15 plotted-data audit: PASS; N_disease metadata mismatch recorded in audit."
)
