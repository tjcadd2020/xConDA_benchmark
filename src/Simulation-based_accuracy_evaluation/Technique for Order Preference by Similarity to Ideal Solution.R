topsis_by_group_minmax <- function(df,
                                   metric_cols = c("MCC", "AUPR","AUC"),
                                   group_cols  = c("Index", "Disease"),
                                   weights = NULL) {
  
  # 如果没给权重，就等权
  if (is.null(weights)) {
    weights <- rep(1 / length(metric_cols), length(metric_cols))
  } else {
    stopifnot(length(weights) == length(metric_cols))
    weights <- weights / sum(weights)  # 归一化权重
  }
  
  df %>%
    group_by(across(all_of(group_cols))) %>%
    group_modify(~ {
      # 当前组的决策矩阵
      mat <- as.matrix(.x[, metric_cols])
      
      ## -------- 1. 组内 min-max 标准化 --------
      # 对每个指标： (x - min) / (max - min)
      col_min <- apply(mat, 2, min, na.rm = TRUE)
      col_max <- apply(mat, 2, max, na.rm = TRUE)
      range   <- col_max - col_min
      
      # 避免 max == min 导致除 0：这种情况下该列全常数，不提供区分信息
      range[range == 0] <- 1
      
      norm_mat <- sweep(mat, 2, col_min, `-`)
      norm_mat <- sweep(norm_mat, 2, range, `/`)
      
      ## -------- 2. 乘以权重 --------
      w_mat <- sweep(norm_mat, 2, weights, `*`)
      
      ## -------- 3. 正负理想解（越大越好型指标）--------
      ideal_pos <- apply(w_mat, 2, max)  # 正理想解
      ideal_neg <- apply(w_mat, 2, min)  # 负理想解
      
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
      
      ## -------- 5. TOPSIS 得分（相对接近度）--------
      denom <- d_pos + d_neg
      # 极端情况下 denom 可能为 0，这时设为 0.5（完全无法区分）
      topsis_score <- ifelse(denom == 0, 0.5, d_neg / denom)
      
      .x$TOPSIS_score <- topsis_score
      .x
    }) %>%
    ungroup()
}