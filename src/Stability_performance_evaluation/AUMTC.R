######稳定性指标AUMTC构建
library(pracma)
AUMTC <- function(y, levels) {
  # 检查输入
  if (!is.numeric(y)) {
    stop("y 必须是数值向量")
  }
  if (length(y) != length(levels)) {
    stop("y 和 levels 的长度必须一致")
  }
  
  # 生成横坐标 (0 到 1 之间等间距)
  n <- length(levels)
  x <- seq(0, 1, length.out = n)
  
  # 计算曲线下面积
  auc_value <- pracma::trapz(x, y)
  
  return(auc_value)
}

#######AUMTC指标计算
idx_spans <- data.frame(
  name  = c(
    "sample_size",
    "less_case_sample",
    "less_case_sample",
    "less_control_sample",
    "less_control_sample",
    "Feature_effect_strength",
    "Confounding_strength",
    "Number_categorical_variable",
    "Number_continuous_variable",
    "Prevalence"
  ),
  start = c(1, 8, 6, 6, 12, 16, 23, 28, 32, 36),
  end   = c(7, 11, 6, 6, 15, 22, 27, 31, 35, 39),
  stringsAsFactors = FALSE
)
build_index_from_spans <- function(idx_spans) {
  stopifnot(all(c("name", "start", "end") %in% names(idx_spans)))
  # 按 name 拆分
  sp <- split(idx_spans, idx_spans$name)
  lapply(sp, function(d) {
    # 把每一段 start:end 展开，再合并
    rows <- unlist(mapply(seq, d$start, d$end, SIMPLIFY = TRUE))
    sort(unique(rows))  # 排序去重（可选）
  })
}
idx_list <- build_index_from_spans(idx_spans)
compute_AUMTC_by_index <- function(df, y_col, level_col, idx) {
  # 基本校验
  if (!is.data.frame(df)) stop("df 必须是 data.frame")
  if (!all(c(y_col, level_col) %in% names(df))) {
    stop("df 中找不到指定的 y_col 或 level_col")
  }
  if (!is.list(idx) || is.null(names(idx))) {
    stop("idx 必须是带名字的列表：每个元素是一段行号")
  }
  
  out <- lapply(names(idx), function(k) {
    rows <- idx[[k]]
    y_vals <- df[[y_col]][rows]
    lvls   <- df[[level_col]][rows]
    
    # 计算
    auc_val <- AUMTC(y_vals, lvls)
    
    data.frame(
      name      = k,
      AUMTC     = auc_val,
      stringsAsFactors = FALSE
    )
  })
  
  do.call(rbind, out)
}
# 用一个空列表收集结果
robustness_list <- list()
for(colname_method_metrics in colnames(all_sens_performance)[-1]){
  tmp <- compute_AUMTC_by_index(all_sens_performance,colname_method_metrics,'id',idx_list)
  tmp$metrics_method <- colname_method_metrics   # 新增一列，标注方法名
  robustness_list[[colname_method_metrics]] <- tmp
}
# 行拼接
robustness_all_metrics <- do.call(rbind, robustness_list)