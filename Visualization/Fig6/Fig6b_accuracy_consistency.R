####导入meta评估结果
evaluation_result <- read.table('metrics_evaluation_result_0.05_mean.tsv',sep = '\t',header = T)
MCC_evaluation_result <- read.table('MCC_evaluation_result_0.05_mean.tsv',sep = '\t',header = T)
evaluation_result <- cbind(evaluation_result, MCC_evaluation_result[,3:4])
evaluation_result$Method <- gsub(evaluation_result$Method,pattern = 'no_meta',replacement = 'mega')
evaluation_result$Method <- gsub(evaluation_result$Method,pattern = 'rma.uni',replacement = '')
method_rename <- function(x){
  x <- gsub(x,pattern = 'maaslin2_lm_res',replacement = 'MaAslin2_GLM',ignore.case = FALSE)
  x <- gsub(x,pattern = 'maaslin2_lm',replacement = 'MaAslin2_GLM',ignore.case = FALSE)
  x <- gsub(x,pattern = 'Maaslin2',replacement = 'MaAslin2',ignore.case = FALSE)
  x <- gsub(x,pattern = 'cplm',replacement = 'CPLM',ignore.case = FALSE)
  x <- gsub(x,pattern = 'negbin',replacement = 'NEGBIN',ignore.case = FALSE)
  x <- gsub(x,pattern = 'megbin',replacement = 'NEGBIN',ignore.case = FALSE)
  x <- gsub(x,pattern = 'zinb',replacement = 'ZINB',ignore.case = FALSE)
  x <- gsub(x,pattern = 'LM-fixed',replacement = 'LFEM',ignore.case = FALSE)
  x <- gsub(x,pattern = 'lmem',replacement = 'LMEM',ignore.case = FALSE)
  x <- gsub(x,pattern = 'ANCOMBC',replacement = 'ANCOM-BC2',ignore.case = FALSE)
  x <- gsub(x,pattern = 'None',replacement = 'count',ignore.case = FALSE)
  return(x)
}
evaluation_result$Method <- method_rename(evaluation_result$Method)
evaluation_result$Method <- gsub(evaluation_result$Method,pattern = '__',replacement = '_')
evaluation_result$Method <- gsub(evaluation_result$Method,pattern = 'None',replacement = 'count')
evaluation_result$Method <- gsub(evaluation_result$Method,pattern = 'MMUPHin_meta',replacement = 'MMUPHin')
evaluation_result$Method <- gsub(evaluation_result$Method,pattern = '_mega',replacement = '')
####估计综合排名
#导入函数
library(RankAggreg)
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
rank_methods_by_metrics <- function(df,
                                    method_col = "Method",
                                    metric_cols = c("MCC", "AUPR", "Macro.F1", "AUC")) {

  # 检查方法列
  if (!method_col %in% colnames(df)) {
    stop(sprintf("'%s' not found in df.", method_col))
  }

  # 如果没有 Macro.F1，但有 F1，则自动用 F1 替代
  metric_cols_used <- metric_cols
  metric_cols_used[metric_cols_used == "Macro.F1" & !("Macro.F1" %in% colnames(df)) & ("F1" %in% colnames(df))] <- "F1"

  # 检查指标列
  missing_cols <- setdiff(metric_cols_used, colnames(df))
  if (length(missing_cols) > 0) {
    stop("These metric columns are missing in df: ", paste(missing_cols, collapse = ", "))
  }

  # 对每个指标分别排序，返回 Method 名称
  rank_list <- lapply(metric_cols_used, function(metric) {
    df[[method_col]][order(df[[metric]], decreasing = TRUE, na.last = TRUE)]
  })

  # 结果转成数据框
  result <- as.data.frame(rank_list, stringsAsFactors = FALSE)

  # 列名保留你想要展示的名字
  colnames(result) <- metric_cols

  return(result)
}
###计算TOPSIS
evaluation_result_without_mmuphin <- evaluation_result[-which(evaluation_result$Method == 'MMUPHin'),]
evaluation_result_without_mmuphin <- evaluation_result_without_mmuphin[,-c(14,15)]
performance_topsis <- topsis_minmax(evaluation_result_without_mmuphin[, c('AUPR','AUC','MCC','Macro.F1','Method'), drop = FALSE],
                                                                metric_cols = c('AUPR','AUC','MCC','Macro.F1'),
                                                                group_cols  = NULL)
###rankaggregation暂时不运行
# rank_df <- rank_methods_by_metrics(evaluation_result)
# performance_rankaggreg <- RankAggreg(
#   x = t(rank_df),
#   k = 64,
#   method = 'CE',
#   distance = 'Spearman',
#   maxIter = 1000,seed = 43
# )
# ###组合信息
# performance_rankaggreg_df <- data.frame(
#   Method = as.character(performance_rankaggreg$top.list),
#   Rank = seq_along(performance_rankaggreg$top.list),
#   stringsAsFactors = FALSE
# )
# evaluation_result<- merge(evaluation_result,performance_rankaggreg_df,by='Method')
# colnames(evaluation_result)[14] <- 'RankAggregation'

evaluation_result_without_mmuphin<- merge(evaluation_result_without_mmuphin,performance_topsis[,c(5,6)],by="Method")

###对指标转化排名
evaluation_res_heatmap <- evaluation_result_without_mmuphin[,c("Method","TOPSIS_score","AUC","AUPR","MCC","Macro.F1","Sensitivity","FDR","FPR")]
rownames(evaluation_res_heatmap) <- evaluation_res_heatmap$Method
rank_df_by_columns <- function(df,
                               method_col = "Method",
                               ascending_cols = c("FDR", "FPR",'RankAggregation'),
                               ties_method = "min") {

  # 1. 基本检查
  if (!method_col %in% colnames(df)) {
    stop(sprintf("'%s' 不在数据框列名中。", method_col))
  }

  metric_cols <- setdiff(colnames(df), method_col)
  asc_cols <- intersect(ascending_cols, metric_cols)
  desc_cols <- setdiff(metric_cols, asc_cols)

  # 2. 复制一个结果表
  out <- df

  # 3. 降序排名：值越大 rank 越靠前
  out[desc_cols] <- lapply(df[desc_cols], function(x) {
    rank(-x, ties.method = ties_method, na.last = "keep")
  })

  # 4. 升序排名：值越小 rank 越靠前
  out[asc_cols] <- lapply(df[asc_cols], function(x) {
    rank(x, ties.method = ties_method, na.last = "keep")
  })

  return(out)
}
evaluation_res_heatmap_rank <- rank_df_by_columns(evaluation_res_heatmap)

###最好的方法排序
###导入consistency结果
consistency_summary_performance <- read.table(file = "summary_multi_cohort_consistency_performance.tsv",sep = '\t',check.names = F,header = T)
consistency_summary_performance$Strategy <- gsub(consistency_summary_performance$Strategy,pattern = '_no_meta',replacement = '')
consistency_summary_performance$Strategy <- gsub(consistency_summary_performance$Strategy,pattern = '_rma.uni',replacement = '')
consistency_summary_performance$Single_strategy <- method_rename(consistency_summary_performance$Single_strategy)
consistency_summary_performance$Strategy <- method_rename(consistency_summary_performance$Strategy)
consistency_summary_performance$Meta_strategy <- gsub(consistency_summary_performance$Meta_strategy,pattern = 'no_meta',replacement = 'Mega')
consistency_summary_performance$Meta_strategy <- gsub(consistency_summary_performance$Meta_strategy,pattern = 'rma.uni_',replacement = '')
library(dplyr)
library(ggplot2)
library(grid)

library(dplyr)
library(ggplot2)

#-----------------------------
# 0. NA处理（只处理数值列）
#-----------------------------
consistency_summary_performance <- consistency_summary_performance %>%
  mutate(across(where(is.numeric), ~ tidyr::replace_na(.x, 0)))

#-----------------------------
# 1. 综合场景表现：同一 Strategy 跨 Disease 取均值
#-----------------------------
integrated_consistency_summary_performance <- consistency_summary_performance %>%
  group_by(Single_strategy, Meta_strategy, Strategy) %>%
  summarise(
    Mean_Jaccard_index = mean(Mean_Jaccard_index, na.rm = TRUE),
    Mean_DA_count      = mean(Mean_DA_count, na.rm = TRUE),
    N_disease          = sum(!is.na(Mean_Jaccard_index)),
    .groups = "drop"
  )
integrated_consistency_summary_performance_without_mmuphin <- integrated_consistency_summary_performance[-which(integrated_consistency_summary_performance$Strategy == 'MMUPHin'),]
#-----------------------------
# 2. Strategy 排序
#-----------------------------
######
evaluation_res_heatmap$Method <- method_rename(evaluation_res_heatmap$Method)
all_strategies_simulation_score <- evaluation_res_heatmap[,c('Method','TOPSIS_score')]
colnames(all_strategies_simulation_score)[2] <-  'simulation_TOPSIS'
all_strategies_consistency_score <- integrated_consistency_summary_performance_without_mmuphin[,c('Strategy','Mean_Jaccard_index')]
colnames(all_strategies_consistency_score)[1] <- 'Method'
all_strategies_2dim_performance <- merge(all_strategies_simulation_score,all_strategies_consistency_score,by = 'Method')
compreh_ranking_2dim <- topsis_minmax(all_strategies_2dim_performance,
                                      metric_cols = c('simulation_TOPSIS','Mean_Jaccard_index'),
                                      group_cols  = NULL)
strategy_order <- compreh_ranking_2dim$Method[order(compreh_ranking_2dim$TOPSIS_score,decreasing = F)]

library(dplyr)
library(ggplot2)

# Disease 顺序
disease_order <- c("adenoma", "CRC", "IBD", "T2D")
disease_levels <- intersect(disease_order, unique(consistency_summary_performance$Disease))

#-----------------------------
# 3. 作图数据
#-----------------------------
consistency_summary_performance_plot <- consistency_summary_performance %>%
  mutate(
    Strategy = factor(Strategy, levels = strategy_order),
    Disease  = factor(Disease, levels = disease_levels)
  ) %>%
  arrange(Strategy, Disease)

integrated_consistency_summary_performance_plot <- integrated_consistency_summary_performance %>%
  mutate(
    Strategy = factor(Strategy, levels = strategy_order)
  ) %>%
  arrange(Strategy)

compreh_ranking_2dim_plot <- compreh_ranking_2dim %>%
  select(Method, TOPSIS_score) %>%
  mutate(Method = factor(Method, levels = strategy_order)) %>%
  arrange(Method)

evaluation_res_heatmap_rank$Method <- method_rename(evaluation_res_heatmap_rank$Method)
evaluation_res_heatmap_rank$Method <- factor(evaluation_res_heatmap_rank$Method,levels = strategy_order)
library(reshape2)
evaluation_res_heatmap_rank_long <- melt(evaluation_res_heatmap_rank,id.vars = "Method")
######删除RankAggregation
evaluation_res_heatmap_rank_long <- evaluation_res_heatmap_rank_long[which(evaluation_res_heatmap_rank_long$variable != 'RankAggregation'),]
evaluation_res_heatmap_rank_long$variable <- factor(evaluation_res_heatmap_rank_long$variable,level = c('FPR','FDR','Sensitivity','AUC','AUPR','Macro.F1','MCC','TOPSIS_score'))


#-----------------------------
# 绘制热图
#-----------------------------
library(dplyr)
library(ggplot2)

# 1. 按你指定的顺序设置列顺序
var_order <- c(
  "FPR",
  "FDR",
  "Sensitivity",
  "AUC",
  "AUPR",
  "Macro.F1",
  "MCC",
  "TOPSIS_score"
)

# 2. 设置显示标签
var_labels <- c(
  FPR = "FPR",
  FDR = "FDR",
  Sensitivity = "Sensitivity",
  AUC = "AUROC",
  AUPR = "AUPR",
  Macro.F1 = "Macro-F1",
  MCC = "MCC",
  TOPSIS_score = "TOPSIS"
)

# 3. 手动设置每一列的位置
# 这里让 TOPSIS 作为最后一列，并和前面的 MCC 稍微隔开
x_pos <- c(
  FPR          = 0.5,
  FDR          = 1.5,
  Sensitivity  = 2.5,
  AUC          = 3.5,
  AUPR         = 4.5,
  Macro.F1     = 5.5,
  MCC          = 6.5,
  TOPSIS_score = 7.9
)

# 4. 构造绘图数据
plot_df <- evaluation_res_heatmap_rank_long %>%
  filter(variable %in% var_order) %>%
  mutate(
    variable = factor(variable, levels = var_order),
    x = x_pos[as.character(variable)]
  )

# 5. 作图
ggplot(plot_df, aes(x = x, y = Method, fill = value)) +
  geom_tile(width = 0.9, height = 0.9) +
  geom_text(aes(label = value), color = "#2F2F2F", size = 2.5) +
  scale_fill_gradient(
    name = "Rank",
    low = "#DEEBF7",
    high = "#084594"
  ) +
  scale_x_continuous(
    breaks = x_pos,
    labels = var_labels[var_order]
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 8),
    panel.grid = element_blank()
  ) +
  labs(title = NULL, x = NULL, y = NULL)
#-----------------------------
# 4. 气泡图
#-----------------------------
bubble_plot <- ggplot(
  consistency_summary_performance_plot,
  aes(
    x = Mean_Jaccard_index,
    y = Strategy,
    size = Mean_DA_count,
    color = Disease
  )
) +
  geom_point(alpha = 0.7) +
  scale_x_continuous(
    limits = c(-0.05, 0.6),
    expand = expansion(mult = c(0.01, 0.03))
  ) +
  scale_size_continuous(
    range = c(2, 7),
    name = "Mean DA count across replicates"
  ) +
  scale_color_manual(
    values = c(
      "CRC" = "#A6CEE3",
      "T2D" = "#FB9A99",
      "adenoma" = "#B2DF8A",
      "IBD" = "#FDBF6F"
    )[disease_levels]
  ) +
  labs(
    x = "Mean Jaccard index across replicates",
    y = NULL,
    color = "Disease"
  ) +
  theme_bw() +
  theme(
    panel.grid.major.y = element_line(color = "grey92", linewidth = 0.35),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 9, color = "black"),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.title.x = element_text(size = 11, color = "black"),
    legend.position = "right"
  )

print(bubble_plot)
######这个bubble plot图暂时不画
#-----------------------------
# 5. 综合表现图：barplot
#-----------------------------
library(dplyr)
library(ggplot2)

#=============================
# 1. 通用函数：准备绘图数据
#=============================
prepare_barplot_data <- function(consistency_df,
                                 integrated_df,
                                 metric_col) {

  # 保留 integrated_df 中已有的 Strategy 顺序
  strategy_levels <- if (is.factor(integrated_df$Strategy)) {
    levels(integrated_df$Strategy)
  } else {
    unique(as.character(integrated_df$Strategy))
  }

  # consistency_df: 用于绘制 bar 上的数据点
  point_df <- consistency_df %>%
    filter(Strategy %in% integrated_df$Strategy) %>%
    transmute(
      Strategy = factor(Strategy, levels = strategy_levels),
      value = .data[[metric_col]]
    )

  # integrated_df: 用于提供 bar 高度
  bar_df <- integrated_df %>%
    filter(Strategy %in% point_df$Strategy) %>%
    transmute(
      Strategy = factor(Strategy, levels = strategy_levels),
      bar_mean = .data[[metric_col]]
    )

  list(bar_df = bar_df, point_df = point_df)
}

#=============================
# 2. 通用函数：画图
#=============================
plot_bar_with_points <- function(consistency_df,
                                 integrated_df,
                                 metric_col,
                                 xlab,
                                 title = NULL,
                                 bar_color) {

  plot_list <- prepare_barplot_data(
    consistency_df = consistency_df,
    integrated_df = integrated_df,
    metric_col = metric_col
  )

  bar_df   <- plot_list$bar_df
  point_df <- plot_list$point_df

  ggplot(bar_df, aes(x = bar_mean, y = Strategy)) +
    geom_col(width = 0.72, fill = bar_color, color = "black", linewidth = 0.3) +
    geom_point(
      data = point_df,
      aes(x = value, y = Strategy),
      inherit.aes = FALSE,
      position = position_jitter(width = 0, height = 0.12, seed = 123),
      size = 0.8,
      alpha = 0.8,
      color = "#1B263B"
    ) +
    labs(
      x = xlab,
      y = NULL,
      title = NULL
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.y = element_text(size = 10),
      axis.text.x = element_text(size = 10),
      axis.title.x = element_text(size = 11),
      panel.grid = element_blank()
    )
}

#=============================
# 3. 绘制 Mean JI barplot
#=============================
p_mean_JI <- plot_bar_with_points(
  consistency_df = consistency_summary_performance_plot,
  integrated_df  = integrated_consistency_summary_performance_plot,
  metric_col     = "Mean_Jaccard_index",
  xlab           = "Mean Jaccard index",
  title          = NULL,
  bar_color      = "#2A9D8F"
)

#=============================
# 4. 绘制 Mean DA count barplot
#=============================
p_mean_DA <- plot_bar_with_points(
  consistency_df = consistency_summary_performance_plot,
  integrated_df  = integrated_consistency_summary_performance_plot,
  metric_col     = "Mean_DA_count",
  xlab           = "Mean DA count",
  title          = NULL,
  bar_color      = "#588157"
)

# 分别显示
p_mean_JI
p_mean_DA

#######barplot
ggplot(compreh_ranking_2dim_plot, aes(y = Method, x = TOPSIS_score)) +
  geom_col(width = 0.8, fill = "#1F4E79") +
  theme_bw() +
  labs(x = NULL, y = "TOPSIS score") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )
####比较基于效应值的Meta分析是否优于仅仅基于P值的Meta分析？是
Integrate_info <- data.frame(Integrate_type = consistency_summary_performance$Meta_strategy,Strategy = consistency_summary_performance$Strategy,Single_strategy = consistency_summary_performance$Single_strategy)
Integrate_info <- Integrate_info[!duplicated(Integrate_info), ]
colnames(Integrate_info)[2] <- 'Method'
evaluation_res_heatmap_rank$Method <- method_rename(evaluation_res_heatmap_rank$Method)
evaluation_rank_comprehensive_info <- merge(evaluation_res_heatmap_rank,Integrate_info,by = 'Method')
library(dplyr)
library(tidyr)
library(purrr)

metric_cols <- c("TOPSIS_score", "RankAggregation", "AUC", "AUPR",
                 "MCC", "Macro.F1", "Sensitivity", "FDR", "FPR")

evaluation_rank_comprehensive_info_long <- evaluation_rank_comprehensive_info %>%
  select(Single_strategy, Integrate_type, all_of(metric_cols)) %>%
  pivot_longer(
    cols = all_of(metric_cols),
    names_to = "Metric",
    values_to = "Value"
  )
paired_group_wilcox <- function(data,
                                better_types = c("PM", "EB", "REML"),
                                worse_types  = c("Stouffer", "Fisher"),
                                lower_is_better = TRUE) {

  tmp <- data %>%
    filter(Integrate_type %in% c(better_types, worse_types)) %>%
    mutate(
      Group = case_when(
        Integrate_type %in% better_types ~ "better",
        Integrate_type %in% worse_types  ~ "worse"
      )
    ) %>%
    group_by(Single_strategy, Group) %>%
    summarise(Value = mean(Value, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = Group, values_from = Value) %>%
    drop_na(better, worse)

  if (nrow(tmp) < 2) {
    return(tibble(
      n = nrow(tmp),
      median_better = NA_real_,
      median_worse = NA_real_,
      median_advantage = NA_real_,
      mean_advantage = NA_real_,
      p.value = NA_real_
    ))
  }

  x <- tmp$better
  y <- tmp$worse

  advantage <- if (lower_is_better) y - x else x - y

  test_res <- wilcox.test(
    x = if (lower_is_better) y else x,
    y = if (lower_is_better) x else y,
    paired = TRUE,
    alternative = "greater",
    exact = FALSE
  )

  tibble(
    n = nrow(tmp),
    median_better = median(x),
    median_worse = median(y),
    median_advantage = median(advantage),
    mean_advantage = mean(advantage),
    p.value = test_res$p.value
  )
}
group_comparison_res <- evaluation_rank_comprehensive_info_long %>%
  group_by(Metric) %>%
  group_modify(~ paired_group_wilcox(.x, lower_is_better = TRUE)) %>%
  ungroup() %>%
  mutate(p.adj.BH = p.adjust(p.value, method = "BH")) %>%
  arrange(Metric)

######source data 此处导出
#多队列DAA策略说明表
strategy_info_output <- evaluation_rank_comprehensive_info[,c(1,11,12)]
colnames(strategy_info_output) <- c('Multi cohort DAA Strategy','Integrating type','Single cohort DAA Strategy')
write.table(strategy_info_output,file='Multi_cohort_Strategy_info.tsv',sep = '\t',col.names = T,row.names = F,quote = F)

#相同single策略内不同Meta分析方法比较
write.table(group_comparison_res,file='meta-analysis_comparision_within_same_singleDAA.tsv',sep = '\t',col.names = T,row.names = F,quote = F)
#多队列DAA模拟数据结果
evaluation_res_heatmap_output <- evaluation_res_heatmap
colnames(evaluation_res_heatmap_output) <- c('Multi cohort DAA Strategy','TOPSIS Score','RankAggregation','AUROC','AUPR','MCC','Macro-F1','Sensitivity','FDR','FPR')
write.table(evaluation_res_heatmap_output,file = 'simulation_evaluation_performance.tsv',sep = '\t',col.names = T,row.names = F,quote = F)
#多队列DAA真实数据结果
consistency_summary_performance_output <- consistency_summary_performance[,c(1,4,5,6)]
colnames(consistency_summary_performance_output) <- c('Disease','Multi cohort DAA Strategy','Mean Jaccard Index Across Replicates','Mean DA Count Across Replicates')
write.table(consistency_summary_performance_output,file = 'consistency_evaluation_performance.tsv',sep = '\t',col.names = T,row.names = F,quote = F)
integrated_consistency_summary_performance_output <- integrated_consistency_summary_performance[,c(3,4,5)]
colnames(integrated_consistency_summary_performance_output) <- c('Multi cohort DAA Strategy','Mean Jaccard Index Across Disease','Mean DA Count Across Diseases')
write.table(integrated_consistency_summary_performance_output,file = 'consistency_evaluation_performance(integrate_disease).tsv',sep = '\t',col.names = T,row.names = F,quote = F)

#####MMUPHin
library(dplyr)
library(tidyr)
library(ggplot2)

# 1. 指定要展示的方法和指标顺序
method_order <- c(
  "ANCOM-BC2_TMM",
  "ANCOM-BC2_count",
  "ANCOM-BC2_CSS",
  "edgeR_CSS",
  "MMUPHin"
)

metric_order <- c(
  "FDR",
  "Sensitivity",
  "AUROC",
  "AUPR",
  "Macro.F1",
  "MCC"
)

# 如果你的列名是 AUC 而不是 AUROC，这里统一改名
plot_df <- evaluation_result %>%
  rename(
    AUROC = AUC
  ) %>%
  filter(Method %in% method_order) %>%
  mutate(
    Method = factor(Method, levels = rev(method_order))
  ) %>%
  select(Method, all_of(metric_order)) %>%
  pivot_longer(
    cols = all_of(metric_order),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    Metric = factor(Metric, levels = metric_order)
  )

p <- ggplot(plot_df, aes(x = Method, y = Value)) +
  geom_col(width = 0.65, fill = "#D58C55") +
  coord_flip() +
  facet_wrap(
    ~ Metric,
    nrow = 2,
    ncol = 3
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "grey90", color = "grey70"),
    strip.text = element_text(size = 11, face = "bold"),

    # 去掉背景网格线
    panel.grid = element_blank(),

    axis.text.y = element_text(size = 9),
    axis.text.x = element_text(size = 9),
    panel.spacing.x = unit(0.6, "lines"),
    panel.spacing.y = unit(0.8, "lines")
  )

p
#####
library(dplyr)

strategy_order <- c(
  "MMUPHin",
  "edgeR_CSS",
  "ANCOM-BC2_CSS",
  "ANCOM-BC2_count",
  "ANCOM-BC2_TMM"
)

consistency_summary_performance_plot_mmuphin <- consistency_summary_performance_plot %>%
  filter(Strategy %in% strategy_order) %>%
  mutate(
    Strategy = factor(Strategy, levels = strategy_order)
  )
integrated_consistency_summary_performance_plot_mmuphin <- integrated_consistency_summary_performance_plot %>%
  filter(Strategy %in% strategy_order) %>%
  mutate(
    Strategy = factor(Strategy, levels = strategy_order)
  )
p_mean_JI <- plot_bar_with_points(
  consistency_df = consistency_summary_performance_plot_mmuphin,
  integrated_df  = integrated_consistency_summary_performance_plot_mmuphin,
  metric_col     = "Mean_Jaccard_index",
  xlab           = "Mean Jaccard index",
  title          = NULL,
  bar_color      = "#2A9D8F"
)
p_mean_DA <- plot_bar_with_points(
  consistency_df = consistency_summary_performance_plot_mmuphin,
  integrated_df  = integrated_consistency_summary_performance_plot_mmuphin,
  metric_col     = "Mean_DA_count",
  xlab           = "Mean DA count",
  title          = NULL,
  bar_color      = "#588157"
)

######导出4strategy指标数据
compreh_ranking_2dim_4stra <- compreh_ranking_2dim[compreh_ranking_2dim$Method %in% c('ANCOM-BC2_CSS','edgeR_CSS','ANCOM-BC2_count','ANCOM-BC2_TMM'),]
evaluation_res_heatmap_4stra <- evaluation_res_heatmap[evaluation_res_heatmap$Method %in% c('ANCOM-BC2_CSS','edgeR_CSS','ANCOM-BC2_count','ANCOM-BC2_TMM'),-2]
integrated_consistency_summary_performance_4stra <- integrated_consistency_summary_performance_output[integrated_consistency_summary_performance_output$`Multi cohort DAA Strategy` %in% c('ANCOM-BC2_CSS','edgeR_CSS','ANCOM-BC2_count','ANCOM-BC2_TMM'),]
colnames(integrated_consistency_summary_performance_4stra)[1] <- 'Method'
integrated_consistency_summary_performance_4stra <- integrated_consistency_summary_performance_4stra[,c(1,3)]
summary_information_4stra <- merge(evaluation_res_heatmap_4stra,compreh_ranking_2dim_4stra,by='Method')
summary_information_4stra <- merge(summary_information_4stra,integrated_consistency_summary_performance_4stra,by = 'Method')
write.table(summary_information_4stra,file = 'Summary_performance.tsv',sep = '\t',row.names = F,col.names = T,quote = F)


#### Fig. 6b：63种策略的模拟表现-一致性二维图
# 数据来自 meta_plot_2606.RData 中的 all_strategies_2dim_performance。
# 三类整合框架严格按 Integrate_info 中的 Integrate_type 定义：
# Mega；REML/PM/EB（效应值 Meta）；Stouffer/Fisher（P 值 Meta）。
library(dplyr)
library(ggplot2)
library(ggrepel)

landscape_required_objects <- c(
  "all_strategies_2dim_performance",
  "Integrate_info",
  "integrated_consistency_summary_performance"
)
landscape_missing_objects <- setdiff(landscape_required_objects, ls())
if (length(landscape_missing_objects) > 0) {
  load("meta_plot_2606.RData")
}
landscape_missing_objects <- setdiff(landscape_required_objects, ls())
if (length(landscape_missing_objects) > 0) {
  stop(
    "Missing required object(s): ",
    paste(landscape_missing_objects, collapse = ", ")
  )
}

landscape_required_columns <- c(
  "Method",
  "simulation_TOPSIS",
  "Mean_Jaccard_index"
)
if (!all(landscape_required_columns %in% colnames(all_strategies_2dim_performance))) {
  stop(
    "all_strategies_2dim_performance must contain: ",
    paste(landscape_required_columns, collapse = ", ")
  )
}

landscape_df <- all_strategies_2dim_performance %>%
  select(all_of(landscape_required_columns))

if (anyDuplicated(landscape_df$Method) > 0) {
  stop("Method must be unique in all_strategies_2dim_performance.")
}
if (any(!complete.cases(landscape_df))) {
  stop("The landscape input contains missing values; no rows were silently removed.")
}

landscape_framework_key <- Integrate_info %>%
  filter(Method %in% landscape_df$Method) %>%
  distinct(Method, Integrate_type) %>%
  mutate(
    Integration_framework = case_when(
      Integrate_type == "Mega" ~ "Mega-analysis",
      Integrate_type %in% c("REML", "PM", "EB") ~ "Effect-size meta",
      Integrate_type %in% c("Stouffer", "Fisher") ~ "P-value meta",
      TRUE ~ NA_character_
    )
  ) %>%
  select(Method, Integration_framework)

landscape_df <- landscape_df %>%
  left_join(landscape_framework_key, by = "Method")

landscape_da_count_key <- integrated_consistency_summary_performance %>%
  distinct(Strategy, Mean_DA_count) %>%
  rename(
    Method = Strategy,
    Mean_detected_DA_taxa = Mean_DA_count
  )

landscape_df <- landscape_df %>%
  left_join(landscape_da_count_key, by = "Method")

if (any(is.na(landscape_df$Integration_framework))) {
  stop(
    "Unclassified strategy/strategies: ",
    paste(landscape_df$Method[is.na(landscape_df$Integration_framework)], collapse = ", ")
  )
}
if (any(is.na(landscape_df$Mean_detected_DA_taxa))) {
  stop(
    "Missing mean detected DA taxa for strategy/strategies: ",
    paste(landscape_df$Method[is.na(landscape_df$Mean_detected_DA_taxa)], collapse = ", ")
  )
}

landscape_framework_levels <- c(
  "Mega-analysis",
  "Effect-size meta",
  "P-value meta"
)
landscape_df$Integration_framework <- factor(
  landscape_df$Integration_framework,
  levels = landscape_framework_levels
)

# 标注现有分析中单独汇总的4个代表性 Mega 策略。
landscape_label_methods <- c(
  "ANCOM-BC2_CSS",
  "edgeR_CSS",
  "ANCOM-BC2_count",
  "ANCOM-BC2_TMM"
)
landscape_label_df <- landscape_df %>%
  filter(Method %in% landscape_label_methods)

landscape_colours <- c(
  "Mega-analysis" = "#168C8C",
  "Effect-size meta" = "#3F78B5",
  "P-value meta" = "#8064A2"
)
landscape_shapes <- c(
  "Mega-analysis" = 21,
  "Effect-size meta" = 24,
  "P-value meta" = 22
)

accuracy_consistency_landscape_plot <- ggplot(
  landscape_df,
  aes(
    x = simulation_TOPSIS,
    y = Mean_Jaccard_index,
    fill = Integration_framework,
    shape = Integration_framework,
    size = Mean_detected_DA_taxa
  )
) +
  geom_point(
    stroke = 0.45,
    colour = "#303030",
    alpha = 0.90
  ) +
  geom_text_repel(
    data = landscape_label_df,
    aes(label = Method),
    family = "Arial",
    size = 2.25,
    fontface = "plain",
    colour = "#202020",
    box.padding = 0.35,
    point.padding = 0.25,
    min.segment.length = 0,
    segment.colour = "#505050",
    segment.size = 0.30,
    direction = "both",
    seed = 2606,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  scale_fill_manual(
    values = landscape_colours,
    breaks = landscape_framework_levels,
    drop = FALSE
  ) +
  scale_shape_manual(
    values = landscape_shapes,
    breaks = landscape_framework_levels,
    drop = FALSE
  ) +
  scale_size_continuous(
    name = "Mean detected\nDA taxa",
    range = c(1.7, 4.2),
    breaks = c(20, 50, 80)
  ) +
  scale_x_continuous(
    name = "Simulation accuracy (TOPSIS score)",
    limits = c(0.19, 1.06),
    breaks = seq(0.2, 1.0, by = 0.2),
    expand = expansion(mult = 0)
  ) +
  scale_y_continuous(
    name = "Cross-cohort consistency\n(Mean Jaccard Index)",
    limits = c(0, 0.36),
    breaks = seq(0, 0.35, by = 0.05),
    expand = expansion(mult = 0)
  ) +
  guides(
    fill = guide_legend(
      title = "Integration framework",
      order = 1,
      override.aes = list(
        shape = unname(landscape_shapes[landscape_framework_levels]),
        size = 2.8,
        alpha = 1
      )
    ),
    shape = "none",
    size = guide_legend(
      title = "Mean detected\nDA taxa",
      order = 2,
      override.aes = list(
        shape = 21,
        fill = "white",
        colour = "#303030",
        alpha = 1
      )
    )
  ) +
  theme_classic(base_size = 7.2, base_family = "Arial") +
  theme(
    plot.title = element_blank(),
    axis.line = element_line(linewidth = 0.42, colour = "black"),
    axis.ticks = element_line(linewidth = 0.35, colour = "black"),
    axis.ticks.length = unit(1.4, "mm"),
    axis.text = element_text(size = 7.2, colour = "black"),
    axis.title = element_text(size = 8.4, face = "bold", colour = "black"),
    axis.title.x = element_text(margin = margin(t = 3.5)),
    axis.title.y = element_text(margin = margin(r = 3.5)),
    legend.position = "right",
    legend.title = element_text(size = 7.2, face = "bold", colour = "black"),
    legend.text = element_text(size = 7.0, colour = "black"),
    legend.key.height = unit(4.7, "mm"),
    legend.key.width = unit(5.0, "mm"),
    legend.spacing.x = unit(1.0, "mm"),
    plot.margin = margin(t = 4, r = 4, b = 4, l = 4, unit = "mm")
  )

# A4宽度略小：190 mm；适当压缩高度以用于 Fig. 6b 横向面板。
landscape_width_mm <- 190
landscape_height_mm <- 70
landscape_width_in <- landscape_width_mm / 25.4
landscape_height_in <- landscape_height_mm / 25.4
landscape_output_prefix <- "Fig6b_accuracy_consistency_landscape_63_strategies"

svglite::svglite(
  paste0(landscape_output_prefix, ".svg"),
  width = landscape_width_in,
  height = landscape_height_in,
  bg = "white"
)
print(accuracy_consistency_landscape_plot)
dev.off()

grDevices::cairo_pdf(
  paste0(landscape_output_prefix, ".pdf"),
  width = landscape_width_in,
  height = landscape_height_in,
  family = "Arial",
  bg = "white"
)
print(accuracy_consistency_landscape_plot)
dev.off()

ragg::agg_tiff(
  file.path(tempdir(), paste0(landscape_output_prefix, ".tiff")),
  width = landscape_width_in,
  height = landscape_height_in,
  units = "in",
  res = 600,
  background = "white",
  scaling = 1
)
print(accuracy_consistency_landscape_plot)
dev.off()
if (!file.copy(
  file.path(tempdir(), paste0(landscape_output_prefix, ".tiff")),
  paste0(landscape_output_prefix, ".tiff"),
  overwrite = TRUE
)) {
  stop("Failed to copy the TIFF export into the working directory.")
}

ragg::agg_png(
  file.path(tempdir(), paste0(landscape_output_prefix, ".png")),
  width = landscape_width_in,
  height = landscape_height_in,
  units = "in",
  res = 300,
  background = "white",
  scaling = 1
)
print(accuracy_consistency_landscape_plot)
dev.off()
if (!file.copy(
  file.path(tempdir(), paste0(landscape_output_prefix, ".png")),
  paste0(landscape_output_prefix, ".png"),
  overwrite = TRUE
)) {
  stop("Failed to copy the PNG export into the working directory.")
}

landscape_framework_counts <- as.data.frame(table(landscape_df$Integration_framework))
colnames(landscape_framework_counts) <- c("Integration_framework", "N_strategies")
print(landscape_framework_counts)
