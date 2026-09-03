
#########展示方法名字
strategy_name <- c("ANCOM-BC2_count","ANCOM-BC2_CSS","ANCOM-BC2_TMM","DESeq2_count","edgeR_count","edgeR_CSS","edgeR_TMM","MaAslin2_CPLM_TSS","MaAslin2_NEGBIN_CSS","MaAslin2_NEGBIN_TMM","MaAslin2_NEGBIN_count","MaAslin2_ZINB_count")
#方法名转化函数
  method_rename <- function(x){
    x <- gsub(x,pattern = 'maaslin2_lm_res',replacement = 'MaAslin2_GLM',ignore.case = FALSE)
    x <- gsub(x,pattern = 'maaslin2_lm',replacement = 'MaAslin2_GLM',ignore.case = FALSE)
    x <- gsub(x,pattern = 'maaslin2',replacement = 'MaAslin2',ignore.case = FALSE)
    x <- gsub(x,pattern = 'Maaslin2',replacement = 'MaAslin2',ignore.case = FALSE)
    x <- gsub(x,pattern = 'cplm',replacement = 'CPLM',ignore.case = FALSE)
    x <- gsub(x,pattern = 'negbin',replacement = 'NEGBIN',ignore.case = FALSE)
    x <- gsub(x,pattern = 'megbin',replacement = 'NEGBIN',ignore.case = FALSE)
    x <- gsub(x,pattern = 'zinb',replacement = 'ZINB',ignore.case = FALSE)
    x <- gsub(x,pattern = 'lfem',replacement = 'LFEM',ignore.case = FALSE)
    x <- gsub(x,pattern = 'LM_inter',replacement = 'LM-inter',ignore.case = FALSE)
    x <- gsub(x,pattern = 'lmem',replacement = 'LMEM',ignore.case = FALSE)
    x <- gsub(x,pattern = 'ANCOMBC',replacement = 'ANCOM-BC2',ignore.case = FALSE)
    x <- gsub(x,pattern = 'None',replacement = 'count',ignore.case = FALSE)
    return(x)
  }
###导入画图数据
perc5_accuracy_simulation_performance <- read.table(file = 'all_method_all_simulation_all_metrics_performance_5perc.tsv',sep = '\t',header = T, check.names = F, stringsAsFactors = F)
library(dplyr)
perc5_accuracy_simulation_performance <- perc5_accuracy_simulation_performance[perc5_accuracy_simulation_performance$Method %in% strategy_name,]
perc5_accuracy_simulation_performance_mean <- perc5_accuracy_simulation_performance %>%
  group_by(Method) %>%
  summarise(
    across(
      c(AUROC, AUPR, FDR, FPR, Sensitivity, MCC, `Macro-F1`, TOPSIS),
      ~ mean(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  )
accuracy_simulation_performance_mean
plot_rank_comparison <- function(
    data1,
    data2,
    id_col = "Method",
    score_col = "TOPSIS",
    condition_names = c("20% differential features",
                        "5% differential features"),
    higher_score_better = TRUE
) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("请先安装 ggplot2：install.packages('ggplot2')")
  }

  if (!all(c(id_col, score_col) %in% names(data1)) ||
      !all(c(id_col, score_col) %in% names(data2))) {
    stop("两个数据框均需包含：", id_col, " 和 ", score_col)
  }

  if (anyDuplicated(data1[[id_col]]) > 0 ||
      anyDuplicated(data2[[id_col]]) > 0) {
    stop("每个方法必须只有一行。")
  }

  common_id <- data1[[id_col]][data1[[id_col]] %in% data2[[id_col]]]

  d1 <- data1[
    match(common_id, data1[[id_col]]),
    c(id_col, score_col),
    drop = FALSE
  ]

  d2 <- data2[
    match(common_id, data2[[id_col]]),
    c(id_col, score_col),
    drop = FALSE
  ]

  score1 <- d1[[score_col]]
  score2 <- d2[[score_col]]

  if (higher_score_better) {
    rank1 <- rank(-score1, ties.method = "average")
    rank2 <- rank(-score2, ties.method = "average")
  } else {
    rank1 <- rank(score1, ties.method = "average")
    rank2 <- rank(score2, ties.method = "average")
  }

  plot_data <- data.frame(
    Method = common_id,
    Rank_1 = rank1,
    Rank_2 = rank2,
    Rank_change = rank2 - rank1
  )

  ## 按20%条件下的排名排列方法
  plot_data <- plot_data[order(plot_data$Rank_1), ]

  plot_data$Method <- factor(
    plot_data$Method,
    levels = rev(plot_data$Method)
  )

  long_data <- rbind(
    data.frame(
      Method = plot_data$Method,
      Condition = condition_names[1],
      Rank = plot_data$Rank_1
    ),
    data.frame(
      Method = plot_data$Method,
      Condition = condition_names[2],
      Rank = plot_data$Rank_2
    )
  )

  long_data$Condition <- factor(
    long_data$Condition,
    levels = condition_names
  )

  p <- ggplot2::ggplot(plot_data) +
    ggplot2::geom_segment(
      ggplot2::aes(
        x = Rank_1,
        xend = Rank_2,
        y = Method,
        yend = Method
      ),
      linewidth = 0.7,
      alpha = 0.7
    ) +
    ggplot2::geom_point(
      data = long_data,
      ggplot2::aes(
        x = Rank,
        y = Method,
        shape = Condition
      ),
      size = 3
    ) +
    ggplot2::scale_x_continuous(
      breaks = seq_len(nrow(plot_data))
    ) +
    ggplot2::labs(
      x = "TOPSIS rank",
      y = NULL,
      shape = NULL
    ) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      legend.position = "top",
      axis.text.y = ggplot2::element_text(size = 10)
    )

  print(p)

  invisible(
    list(
      result = plot_data,
      plot = p
    )
  )
}
rank_plot_res <- plot_rank_comparison(
  data1 = accuracy_simulation_performance_mean,
  data2 = perc5_accuracy_simulation_performance_mean,
  id_col = "Method",
  score_col = "TOPSIS"
)

rank_plot_res$result
library(dplyr)

topsis_sum_rank <- perc5_accuracy_simulation_performance_mean %>%
  select(Method, TOPSIS_perc5 = TOPSIS) %>%
  inner_join(
    accuracy_simulation_performance_mean %>%
      select(Method, TOPSIS_primary = TOPSIS),
    by = "Method"
  ) %>%
  mutate(
    TOPSIS_sum = TOPSIS_perc5 + TOPSIS_primary
  ) %>%
  arrange(desc(TOPSIS_sum))

topsis_sum_rank
##########
accuracy_simulation_performance <- read.table(file = 'all_method_all_simulation_all_metrics_performance.tsv',sep = '\t',header = T, check.names = F, stringsAsFactors = F)
ref_CD_performance <- read.table(file = 'CD_ref_performance.tsv',sep = '\t',header = T, check.names = F, stringsAsFactors = F)
ref_dental_performance <- read.table(file = 'ref_enrichment_dental_performance.tsv',sep = '\t',header = T, check.names = F, stringsAsFactors = F)
ref_vaginal_performance <- read.table(file = 'ref_enrichment_vaginal_performance.tsv',sep = '\t',header = T, check.names = F, stringsAsFactors = F)
robustness_stability_performance <- read.table(file = 'Perturbation_AUMTC_performance.tsv',sep = '\t',header = T, check.names = F, stringsAsFactors = F)
robustness_consistency_performance <- read.table(file = 'Consistency_Summary_performance.tsv',sep = '\t',header = T, check.names = F, stringsAsFactors = F)
###过滤其他策略
accuracy_simulation_performance <- accuracy_simulation_performance[accuracy_simulation_performance$Method %in% strategy_name,]
ref_CD_performance <- ref_CD_performance[ref_CD_performance$Method %in% strategy_name,]
ref_dental_performance <- ref_dental_performance[ref_dental_performance$Method %in% strategy_name,]
ref_vaginal_performance <- ref_vaginal_performance[ref_vaginal_performance$Method %in% strategy_name,]
robustness_consistency_performance <- robustness_consistency_performance[robustness_consistency_performance$Method %in% strategy_name,]
###提取展示的列
library(dplyr)
accuracy_simulation_performance_mean <- accuracy_simulation_performance %>%
  group_by(Method) %>%
  summarise(
    across(c(AUROC, AUPR, FDR, FPR, Sensitivity, MCC, `Macro-F1`, TOPSIS),
           ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  )
ref_CD_performance <- ref_CD_performance[,c(5,6)]
ref_dental_performance <- ref_dental_performance[,c(1,3)]
ref_vaginal_performance <- ref_vaginal_performance[,c(1,3)]
colnames(ref_vaginal_performance)[2] <- 'RRES_score_v'
robustness_consistency_performance <- robustness_consistency_performance[,c(1,3,6)]
robustness_stability_performance
topsis_minmax <- function(df,
                          metric_cols = c("MCC", "AUPR", "AUC"),
                          group_cols  = NULL,
                          weights = NULL) {

  # 基本检查
  stopifnot(all(metric_cols %in% colnames(df)))

  if (!is.null(group_cols) && length(group_cols) > 0) {
    stopifnot(all(group_cols %in% colnames(df)))
  }

  # 如果没给权重，就等权
  if (is.null(weights)) {
    weights <- rep(1 / length(metric_cols), length(metric_cols))
  } else {
    stopifnot(length(weights) == length(metric_cols))
    weights <- weights / sum(weights)  # 归一化权重
  }

  # 单个数据块的 TOPSIS 计算函数
  compute_topsis <- function(sub_df) {
    mat <- as.matrix(sub_df[, metric_cols, drop = FALSE])

    ## -------- 1. min-max 标准化 --------
    col_min <- apply(mat, 2, min, na.rm = TRUE)
    col_max <- apply(mat, 2, max, na.rm = TRUE)
    rng     <- col_max - col_min

    # 避免 max == min 导致除 0
    rng[rng == 0] <- 1

    norm_mat <- sweep(mat, 2, col_min, `-`)
    norm_mat <- sweep(norm_mat, 2, rng, `/`)

    ## -------- 2. 乘以权重 --------
    w_mat <- sweep(norm_mat, 2, weights, `*`)

    ## -------- 3. 正负理想解（默认都是越大越好）--------
    ideal_pos <- apply(w_mat, 2, max, na.rm = TRUE)
    ideal_neg <- apply(w_mat, 2, min, na.rm = TRUE)

    ## -------- 4. 到理想解的欧氏距离 --------
    d_pos <- sqrt(rowSums(
      (w_mat - matrix(ideal_pos,
                      nrow = nrow(w_mat),
                      ncol = ncol(w_mat),
                      byrow = TRUE))^2
    ))

    d_neg <- sqrt(rowSums(
      (w_mat - matrix(ideal_neg,
                      nrow = nrow(w_mat),
                      ncol = ncol(w_mat),
                      byrow = TRUE))^2
    ))

    ## -------- 5. TOPSIS 得分 --------
    denom <- d_pos + d_neg
    topsis_score <- ifelse(denom == 0, 0.5, d_neg / denom)

    sub_df$TOPSIS_score <- topsis_score
    sub_df
  }

  # 有分组：组内 TOPSIS
  if (!is.null(group_cols) && length(group_cols) > 0) {
    res <- df %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) %>%
      dplyr::group_modify(~ compute_topsis(.x)) %>%
      dplyr::ungroup()
  } else {
    # 无分组：全局 TOPSIS
    res <- compute_topsis(df)
  }

  res
}

topsis_from_long <- function(df_long,
                             perturb_col = "Perturbation_factor",
                             metric_col = "Accuracy_metric",
                             method_col = "Method",
                             value_col = "AUMTC",
                             metric_cols = c("AUC","AUPR", "MCC", "Macro-F1"),
                             weights = NULL) {

  # 先筛选要用于 TOPSIS 的指标
  df_wide <- df_long %>%
    dplyr::filter(.data[[metric_col]] %in% metric_cols) %>%
    dplyr::select(
      all_of(c(perturb_col, method_col, metric_col, value_col))
    ) %>%
    tidyr::pivot_wider(
      names_from = all_of(metric_col),
      values_from = all_of(value_col)
    )

  # 如果有方法缺少任一指标，这里先去掉，避免 TOPSIS 出问题
  df_wide_complete <- df_wide %>%
    dplyr::filter(dplyr::if_all(dplyr::all_of(metric_cols), ~ !is.na(.x)))

  # 组内计算 TOPSIS（每个 Perturbation_factor 内部比较 Method）
  res <- topsis_minmax(
    df = df_wide_complete,
    metric_cols = metric_cols,
    group_cols = perturb_col,
    weights = weights
  )

  # 只返回你要的结果
  res_out <- res %>%
    dplyr::select(
      all_of(perturb_col),
      all_of(method_col),
      TOPSIS_score
    ) %>%
    dplyr::arrange(.data[[perturb_col]], dplyr::desc(TOPSIS_score))

  return(res_out)
}

stability_topsis_performance <- topsis_from_long(robustness_stability_performance)
library(dplyr)
library(tidyr)

stability_topsis_performance <- stability_topsis_performance %>%
  pivot_wider(
    names_from = Perturbation_factor,
    values_from = TOPSIS_score
  )
stability_topsis_performance <- as.data.frame(stability_topsis_performance)
library(dplyr)
library(purrr)

info_list <- list(ref_CD_performance,ref_dental_performance,ref_vaginal_performance,accuracy_simulation_performance_mean, robustness_consistency_performance, stability_topsis_performance)

merged_all_performance <- reduce(info_list, full_join, by = "Method")
merged_all_performance$DAA_Method_type <-
  c(
    'RNA-seq-derived methods',
    'Microbiome-tailored methods',
    'Microbiome-tailored methods',
    'RNA-seq-derived methods',
    'Microbiome-tailored methods',
    'Microbiome-tailored methods',
    'Microbiome-tailored methods',
    'Microbiome-tailored methods',
    'Microbiome-tailored methods',
    'RNA-seq-derived methods',
    'Microbiome-tailored methods',
    'RNA-seq-derived methods'
  )
merged_all_performance <-
  merged_all_performance[, c(
    'Method',
    'DAA_Method_type',
    'FPR',
    'FDR',
    'Sensitivity',
    'AUROC',
    'AUPR',
    'MCC',
    'Macro-F1',
    'TOPSIS',
    'RRES_score',
    'RRES_score_v',
    'SRS',
    'Mean Jaccard Index',
    'Mean DA Count Across Datasets',
    'Sample size',
    'Case-control ratio',
    'Primary effect size',
    'Confounding effect size',
    'Confounder complexity',
    'Feature Prevalence'
  )]
merged_all_performance <- topsis_minmax(
  merged_all_performance,
  metric_cols = c(
    "TOPSIS",
    "RRES_score",
    "RRES_score_v",
    "SRS",
    'Mean Jaccard Index',
    'Sample size',
    'Case-control ratio',
    'Primary effect size',
    'Confounding effect size',
    'Confounder complexity',
    'Feature Prevalence'
  ),
  group_cols  = NULL,
  weights = NULL
)
merged_all_performance <-
  merged_all_performance[, c(
    'Method',
    'DAA_Method_type',
    'TOPSIS_score',
    'FPR',
    'FDR',
    'Sensitivity',
    'AUROC',
    'AUPR',
    'MCC',
    'Macro-F1',
    'TOPSIS',
    'RRES_score',
    'RRES_score_v',
    'SRS',
    'Mean Jaccard Index',
    'Mean DA Count Across Datasets',
    'Sample size',
    'Case-control ratio',
    'Primary effect size',
    'Confounding effect size',
    'Confounder complexity',
    'Feature Prevalence'
  )]
#######画图数据准备完毕
colnames(merged_all_performance)[1] <- 'id'
rownames(merged_all_performance) <-merged_all_performance$id
merged_all_performance <- merged_all_performance[merged_all_performance$id[order(merged_all_performance$TOPSIS_score,decreasing = T)],]

#######准备Funkyheatmap图
#######绘图颜色准备
library(funkyheatmap)

#=============================
# 1. palette：短名字 + 直接单色
#=============================
palettes <- list(
  sim_fill     = "#DDEBF7",
  ref_fill     = "#FCE4D6",
  cons_fill    = "#FFF2CC",
  stab_fill    = "#E2EFDA",
  overall_fill = "#EADBFA",

  strat_bg   = "#595959",
  sim_bg     = "#5B9BD5",
  ref_bg     = "#ED7D31",
  cons_bg    = "#FFC000",
  stab_bg    = "#70AD47",
  overall_bg = "#7030A0"
)

#=============================
# 2. column_info
#=============================
column_info <- data.frame(
  id = colnames(merged_all_performance),

  group = c(
    rep("Strategy Information", 2),
    "Overall",
    rep("Simulation-based Accuracy", 8),
    rep("Reference-based Recovery", 3),
    rep("Consistency", 2),
    rep("Stability Under Factor Perturbation", 6)
  ),

  width = c(
    5.5, 6.7, 2.5,
    rep(1, 19)
  ),

  name = c(
    "Strategy",
    "DAAType",
    "Overall Score",
    "FPR",
    "FDR",
    "Sensitivity",
    "AUROC",
    "AUPR",
    "MCC",
    "Macro-F1",
    "TOPSIS",
    "Dental plaque scenario",
    "Vaginal scenario",
    "Crohn Disease",
    "Mean Jaccard Index",
    "Mean DA count",
    "Sample size",
    "Case-control ratio",
    "Primary effect size",
    "Confounding effect size",
    "Confounder complexity",
    "Feature prevalence"
  ),

  geom = c(
    "text",
    "text",
    "bar",
    rep("funkyrect", 19)
  ),

  palette = c(
    NA, NA,
    "overall_fill",
    rep("sim_fill", 8),
    rep("ref_fill", 3),
    rep("cons_fill", 2),
    rep("stab_fill", 6)
  ),

  stringsAsFactors = FALSE
)

column_info$legend <- FALSE

column_info$draw_outline <- ifelse(column_info$geom == "bar", FALSE, NA)

# 强制转 character，避免任何隐式 factor 问题
column_info$group   <- as.character(column_info$group)
column_info$palette <- as.character(column_info$palette)

#=============================
# 3. column_groups
#=============================
column_groups <- data.frame(
  level1 = c(
    "Strategy Information",
    "Simulation-based Accuracy",
    "Reference-based Recovery",
    "Consistency",
    "Stability Under Factor Perturbation",
    "Overall"
  ),
  group = c(
    "Strategy Information",
    "Simulation-based Accuracy",
    "Reference-based Recovery",
    "Consistency",
    "Stability Under Factor Perturbation",
    "Overall"
  ),
  palette = c(
    "strat_bg",
    "sim_bg",
    "ref_bg",
    "cons_bg",
    "stab_bg",
    "overall_bg"
  ),
  stringsAsFactors = FALSE
)

column_groups$group   <- as.character(column_groups$group)
column_groups$palette <- as.character(column_groups$palette)

#=============================
# 4. 画图
#=============================
legends <- list(
  list(
    palette = "sim_fill",
    geom = "funkyrect",
    title = "Scaled score",
    enabled = TRUE,
    labels = c("0", "", "0.2", "", "0.4", "", "0.6", "", "0.8", "", "1")
  ),

  list(palette = "ref_fill", enabled = FALSE),
  list(palette = "cons_fill", enabled = FALSE),
  list(palette = "stab_fill", enabled = FALSE),
  list(palette = "overall_fill", enabled = FALSE),

  list(palette = "strat_bg", enabled = FALSE),
  list(palette = "sim_bg", enabled = FALSE),
  list(palette = "ref_bg", enabled = FALSE),
  list(palette = "cons_bg", enabled = FALSE),
  list(palette = "stab_bg", enabled = FALSE),
  list(palette = "overall_bg", enabled = FALSE)
)

p <- funky_heatmap(
  merged_all_performance,
  column_info   = column_info,
  column_groups = column_groups,
  palettes      = palettes,
  legends       = legends,
  scale_column  = TRUE,
  add_abc       = FALSE,
  position_args = position_arguments(
    col_space        = 0.1,
    col_bigspace     = 0.6,
    col_annot_angle  = 90,
    col_annot_offset = 1.5,
    expand_ymax      = 2.5
  )
)


#####导出TOPSIS stability评分
write.table(stability_topsis_performance,file='Stability_performance_TOPSIS.tsv',sep = '\t',row.names = F,col.names = T,quote = F )