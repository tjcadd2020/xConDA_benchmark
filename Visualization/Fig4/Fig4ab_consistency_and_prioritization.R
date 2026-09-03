#########导入要被去掉的方法
method_rename <- function(x){
  x <- gsub(x,pattern = 'maaslin2_lm_res',replacement = 'MaAslin2_GLM',ignore.case = FALSE)
  x <- gsub(x,pattern = 'maaslin2_lm',replacement = 'MaAslin2_GLM',ignore.case = FALSE)
  x <- gsub(x,pattern = 'maaslin2',replacement = 'MaAslin2',ignore.case = FALSE)
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
###############################
####导入rank信息
CD_rank <- readRDS("CD_rank.rds")
simulation_rank <- readRDS("simulation_rank.rds")
simulation_rank2 <- readRDS("simulation_rank2.rds")
vaginal_rank <- readRDS("vaginal_rank.rds")
dental_rank <- readRDS("dental_rank.rds")
CD_rank <- CD_rank[-which(CD_rank$Method == 'MaAslin2_CPLM_LOG_TSS'),]
vaginal_rank <- vaginal_rank[-which(vaginal_rank$Method == 'MaAslin2_CPLM_LOG_TSS'),]
dental_rank <- dental_rank[-which(dental_rank$Method == 'maaslin2_cplm_LOG_TSS'),]
simulation_rank <- simulation_rank[-which(simulation_rank$Method == 'MaAslin2_CPLM_LOG_TSS'),]

####
CD_filter <- CD_rank$Method[CD_rank$Sum<=quantile(CD_rank$Sum, probs = 0.25)]
vaginal_filter <- vaginal_rank$Method[vaginal_rank$score<=quantile(vaginal_rank$score, probs = 0.25)]
dental_filter <- dental_rank$Method[dental_rank$score<=quantile(dental_rank$score, probs = 0.25)]
dental_filter  <- c(dental_filter,c('LM-inter_CLR','LMEM_CLR','LFEM_CLR'))
dental_filter <- setdiff(dental_filter,'LM_CLR')
CD_filter <- method_rename(CD_filter)
dental_filter<- method_rename(dental_filter)
bottom25_atleast2 <- Reduce(union, list(
  intersect(CD_filter, dental_filter),
  intersect(CD_filter, vaginal_filter),
  intersect(dental_filter, vaginal_filter)
))
bottom25_atleast2 <- sort(unique(bottom25_atleast2))
simulation_rank <- simulation_rank[order(-simulation_rank$TOPSIS_score_mean),]
simulation_rank$Method <- method_rename(simulation_rank$Method)
simulation_rank2 <- simulation_rank2[order(-simulation_rank2$TOPSIS_score_mean),]
simulation_rank2$Method <- method_rename(simulation_rank2$Method)
####模拟数据直接取中位数筛选
simulation_filter <- simulation_rank$Method[simulation_rank$TOPSIS_score_mean<=quantile(simulation_rank$TOPSIS_score_mean, probs = 0.5)]
simulation_filter2 <- simulation_rank2$Method[simulation_rank2$TOPSIS_score_mean<=quantile(simulation_rank2$TOPSIS_score_mean, probs = 0.5)]
filter_method<-intersect(simulation_filter,bottom25_atleast2)
filter_method<-intersect(simulation_filter2,filter_method)
###输出模拟数据不错但真实数据中差
simulation_selected <- simulation_rank$Method[simulation_rank$TOPSIS_score_mean>quantile(simulation_rank$TOPSIS_score_mean, probs = 0.5)]
intersect(simulation_selected,bottom25_atleast2)

rownames(simulation_rank) <- simulation_rank$Method
rownames(simulation_rank2) <- simulation_rank2$Method
cor.test(simulation_rank[simulation_rank$Method,]$TOPSIS_score_mean,simulation_rank2[simulation_rank$Method,]$TOPSIS_score_mean,method = 'spearman')
#S = 13917, p-value = 6.528e-15
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#  rho
#0.7666563
library(dplyr)
library(ggplot2)

# 按 Method 对齐两个数据集
cor_data <- simulation_rank %>%
  select(Method, TOPSIS_20 = TOPSIS_score_mean) %>%
  inner_join(
    simulation_rank2 %>%
      select(Method, TOPSIS_5 = TOPSIS_score_mean),
    by = "Method"
  )

# Spearman相关分析
cor_res <- cor.test(
  cor_data$TOPSIS_20,
  cor_data$TOPSIS_5,
  method = "spearman",
  exact = FALSE
)

rho <- unname(cor_res$estimate)
p_value <- cor_res$p.value

# 绘制相关性散点图
ggplot(cor_data, aes(x = TOPSIS_20, y = TOPSIS_5)) +
  geom_point(size = 2.5, alpha = 0.8) +
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    linewidth = 0.8
  ) +
  annotate(
    "text",
    x = Inf,
    y = -Inf,
    hjust = 1.1,
    vjust = -0.8,
    label = paste0(
      "Spearman \u03c1 = ", round(rho, 2),
      "\nP ", ifelse(
        p_value < 0.001,
        "< 0.001",
        paste0("= ", signif(p_value, 3))
      )
    ),
    size = 4.5
  ) +
  labs(
    x = "Mean TOPSIS score under 20% differential features",
    y = "Mean TOPSIS score under 5% differential features"
  ) +
  theme_classic(base_size = 13)

#####
#[1] "VTwins_TMM"          "MaAslin2_NEGBIN_TMM"
bottom25_atleast2_rank_in_simulation <- as.numeric(rownames(simulation_rank)[simulation_rank$Method %in% bottom25_atleast2])
#########结果导入
subset_res <- read.table(file = './subset_evaluation_result_all.tsv',header = T)
cross_res <- read.table(file = './cross_evaluation_result_all.tsv',header = T)
##########
library(dplyr)

summary_subset_res <- subset_res %>%
  group_by(method) %>%
  summarise(
    n = n(),

    # Jaccard
    jaccard_mean = mean(Jaccard_index, na.rm = TRUE),
    jaccard_median = median(Jaccard_index, na.rm = TRUE),
    jaccard_var  = var(Jaccard_index,  na.rm = TRUE),

    # subset1 count
    subset1_mean = mean(mean_subset1_count, na.rm = TRUE),
    subset1_median = median(mean_subset1_count, na.rm = TRUE),
    subset1_var  = var(mean_subset1_count,  na.rm = TRUE),

    # subset2 count
    subset2_mean = mean(mean_subset2_count, na.rm = TRUE),
    subset2_median = median(mean_subset2_count, na.rm = TRUE),
    subset2_var  = var(mean_subset2_count,  na.rm = TRUE),

    # 两列 count 一起（把两列拼起来后再算均值/方差）
    count_mean_all = mean(c(mean_subset1_count, mean_subset2_count), na.rm = TRUE),
    count_median_all = median(c(mean_subset1_count, mean_subset2_count), na.rm = TRUE),
    count_var_all  = var(c(mean_subset1_count, mean_subset2_count),  na.rm = TRUE),

    .groups = "drop"
  )
summary_subset_res$method <- method_rename(summary_subset_res$method)
summary_subset_res <- summary_subset_res[!(summary_subset_res$method %in% filter_method),]

######subset画图
########
#整体排名图
# 按 Jaccard median 排序；Mean 与 Median 分列展示，避免点位重叠
library(dplyr)
library(ggplot2)
library(patchwork)

benchmark_blue <- "#6F91BB"  # 参考图中 benchmark-prioritized 区域的背景色
mean_blue <- "#0D2F62"
other_text <- "#17353D"
strategy_background <- "#F6F8F8"
priority_text <- "#FFFFFF"
benchmark_header_text <- "#222A31"

# 最新参考图包含 60 个策略；原始数据中的该策略不进入本图
reference_strategy_exclusions <- "MaAslin2_CPLM_TSS_LOG"

summary_subset_res_j <- summary_subset_res %>%
  mutate(
    jaccard_sd = sqrt(jaccard_var),
    method = as.character(method)
  ) %>%
  filter(!method %in% reference_strategy_exclusions)

# 在绘图前计算 benchmark-prioritized 策略，避免依赖工作区中已有对象
summary_subset_res2 <- summary_subset_res_j %>%
  mutate(
    ji_lb = pmax(0, jaccard_mean - jaccard_sd),
    ji_robust = pmin(jaccard_median, ji_lb)
  )

min_count_med <- 5
min_sub_med <- 5
min_ji_med <- 0

kept <- summary_subset_res2 %>%
  filter(
    count_median_all >= min_count_med,
    subset1_median >= min_sub_med,
    subset2_median >= min_sub_med,
    jaccard_median >= min_ji_med
  ) %>%
  arrange(desc(ji_robust), desc(count_median_all), jaccard_var)

kept_methods <- as.character(kept$method)

figure_method_order <- summary_subset_res_j %>%
  mutate(selected = !is.na(jaccard_median) & jaccard_median > 0) %>%
  arrange(
    desc(selected),
    desc(jaccard_median),
    desc(jaccard_mean)
  ) %>%
  pull(method) %>%
  unique()

method_levels <- rev(figure_method_order)
summary_subset_res_j <- summary_subset_res_j %>%
  mutate(
    method = factor(method, levels = method_levels),
    benchmark_prioritized = as.character(method) %in% kept_methods
  )

plot_df <- summary_subset_res_j %>%
  transmute(
    method,
    jaccard_median,
    median_count_all = count_median_all,
    benchmark_prioritized
  )

median_count_limit <- max(
  75,
  ceiling(max(plot_df$median_count_all, na.rm = TRUE) / 25) * 25
)

mean_count_limit <- max(
  75,
  ceiling(max(summary_subset_res_j$count_mean_all, na.rm = TRUE) / 25) * 25
)

theme_consistency <- theme_classic(base_size = 7.5, base_family = "Arial") +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_line(linewidth = 0.35, colour = "black"),
    axis.title = element_text(size = 7.5, face = "bold", colour = "black"),
    axis.text = element_text(size = 6.5, colour = "black"),
    plot.title = element_text(size = 8.5, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 6.2, face = "bold", hjust = 0.5),
    legend.text = element_text(size = 6.2),
    legend.key.height = grid::unit(8, "pt"),
    legend.key.width = grid::unit(13, "pt"),
    legend.spacing.x = grid::unit(1.5, "pt"),
    legend.margin = margin(t = 2, r = 0, b = 0, l = 0),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.45),
    plot.margin = margin(t = 4, r = 3, b = 4, l = 3)
  )

p_mean <- ggplot(
  summary_subset_res_j,
  aes(x = jaccard_mean, y = method)
) +
  geom_errorbar(
    aes(
      xmin = pmax(jaccard_mean - jaccard_sd, 0),
      xmax = jaccard_mean + jaccard_sd
    ),
    width = 0.45,
    linewidth = 0.38,
    colour = "#D8DDE1"
  ) +
  geom_point(
    aes(size = count_mean_all),
    shape = 18,
    colour = mean_blue,
    alpha = 0.95
  ) +
  scale_size_continuous(
    limits = c(0, mean_count_limit),
    breaks = c(0, 25, 50, 75),
    range = c(0.45, 4.0),
    name = "Mean count of differential\ntaxa (across datasets)"
  ) +
  scale_x_continuous(
    limits = c(0, 0.3),
    breaks = seq(0, 0.3, by = 0.1),
    expand = expansion(mult = c(0.015, 0.015))
  ) +
  scale_y_discrete(limits = method_levels, expand = expansion(add = 0.55)) +
  guides(
    size = guide_legend(
      title.position = "top",
      title.hjust = 0.5,
      label.position = "bottom",
      nrow = 1,
      override.aes = list(shape = 18, colour = mean_blue, alpha = 1)
    )
  ) +
  labs(title = NULL, x = "Jaccard index (Mean \u00B1 SD)", y = NULL) +
  theme_consistency +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "bottom",
    legend.justification = "center",
    legend.box = "vertical",
    legend.box.just = "center"
  )

p_median <- ggplot(
  plot_df,
  aes(x = jaccard_median, y = method, size = median_count_all)
) +
  geom_point(colour = benchmark_blue, alpha = 0.78) +
  scale_size_continuous(
    limits = c(0, median_count_limit),
    breaks = c(0, 25, 50, 75),
    range = c(0.45, 4.0),
    name = "Median count of differential\ntaxa (across datasets)"
  ) +
  scale_x_continuous(
    limits = c(0, 0.3),
    breaks = seq(0, 0.3, by = 0.1),
    expand = expansion(mult = c(0.025, 0.025))
  ) +
  scale_y_discrete(limits = method_levels, expand = expansion(add = 0.55)) +
  guides(
    size = guide_legend(
      title.position = "top",
      title.hjust = 0.5,
      label.position = "bottom",
      nrow = 1
    )
  ) +
  labs(title = NULL, x = "Jaccard index (Median)", y = NULL) +
  theme_consistency +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "bottom",
    legend.justification = "center"
  )

# 独立方法标签区：恢复参考图的分组背景与文字配色
label_df <- summary_subset_res_j %>%
  transmute(method, benchmark_prioritized)

priority_y <- as.numeric(label_df$method[label_df$benchmark_prioritized])
benchmark_midpoint <- mean(priority_y)
priority_ymin <- min(priority_y) - 0.5
priority_ymax <- max(priority_y) + 0.55

p_method_labels <- ggplot(label_df, aes(y = method)) +
  annotate(
    "rect",
    xmin = -Inf,
    xmax = Inf,
    ymin = 0.45,
    ymax = priority_ymin,
    fill = strategy_background,
    colour = NA
  ) +
  annotate(
    "rect",
    xmin = -Inf,
    xmax = Inf,
    ymin = priority_ymin,
    ymax = priority_ymax,
    fill = benchmark_blue,
    colour = NA
  ) +
  geom_text(
    aes(x = 0, label = method, colour = benchmark_prioritized),
    hjust = 0,
    family = "Arial",
    fontface = "bold",
    size = 6.2 / 2.845276
  ) +
  annotate(
    "text",
    x = 1.04,
    y = benchmark_midpoint,
    label = "Benchmark-\nprioritized\nstrategies",
    hjust = 0.5,
    vjust = 0.5,
    family = "Arial",
    fontface = "bold",
    size = 6.4 / 2.845276,
    colour = benchmark_header_text
  ) +
  scale_colour_manual(
    values = c(`FALSE` = other_text, `TRUE` = priority_text),
    guide = "none"
  ) +
  scale_x_continuous(limits = c(0, 1.35), expand = expansion(mult = c(0, 0))) +
  scale_y_discrete(limits = method_levels, expand = expansion(add = 0.55)) +
  coord_cartesian(clip = "off") +
  theme_void(base_family = "Arial") +
  theme(
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = strategy_background, colour = NA),
    plot.margin = margin(t = 4, r = 2, b = 4, l = 2)
  )

# Mean:Median = 1.5:1；不添加 a/b/c/d 子图标签
p_overlay <- p_mean + p_median + p_method_labels +
  plot_layout(widths = c(1.50, 1.00, 1.75), guides = "keep")

p_overlay

# 宽度为 A4 的 5/8；高度沿用参考图（626 x 1162）的纵横比
figure_width_in <- (595.276 * 5 / 8) / 72
reference_aspect_ratio <- 1162 / 626
figure_height_in <- figure_width_in * reference_aspect_ratio

svglite::svglite(
  "consistency_mean_median_split.svg",
  width = figure_width_in,
  height = figure_height_in
)
print(p_overlay)
dev.off()

grDevices::cairo_pdf(
  "consistency_mean_median_split.pdf",
  width = figure_width_in,
  height = figure_height_in,
  family = "Arial"
)
print(p_overlay)
dev.off()

preview_png_tmp <- tempfile(fileext = ".png")
ragg::agg_png(
  preview_png_tmp,
  width = figure_width_in,
  height = figure_height_in,
  units = "in",
  res = 300,
  background = "white"
)
print(p_overlay)
dev.off()
file.copy(
  preview_png_tmp,
  "consistency_mean_median_split.png",
  overwrite = TRUE
)
unlink(preview_png_tmp)

submission_tiff_tmp <- tempfile(fileext = ".tiff")
ragg::agg_tiff(
  submission_tiff_tmp,
  width = figure_width_in,
  height = figure_height_in,
  units = "in",
  res = 600,
  background = "white",
  compression = "lzw"
)
print(p_overlay)
dev.off()
file.copy(
  submission_tiff_tmp,
  "consistency_mean_median_split.tiff",
  overwrite = TRUE
)
unlink(submission_tiff_tmp)

# --- 策略筛选图：Median JI > 0 且 Median detected DA count > 5 ---
# 与上图使用同一个最新60策略集合；严格使用用户指定的“>”阈值。
median_ji_filter_threshold <- 0
median_da_filter_threshold <- 5

strategy_filter_df <- summary_subset_res_j %>%
  mutate(
    method = as.character(method),
    selected = jaccard_median > median_ji_filter_threshold &
      count_median_all > median_da_filter_threshold
  )

if (nrow(strategy_filter_df) != 60) {
  stop("Strategy filter figure expects the latest 60-strategy reference set.")
}
if (any(!is.finite(strategy_filter_df$jaccard_median)) ||
    any(!is.finite(strategy_filter_df$count_median_all))) {
  stop("The strategy filter figure contains missing or non-finite values.")
}
if (any(strategy_filter_df$count_median_all < 0)) {
  stop("Median detected DA counts must be non-negative.")
}

selected_strategy_n <- sum(strategy_filter_df$selected)
other_strategy_n <- nrow(strategy_filter_df) - selected_strategy_n

filter_selected_blue <- "#2A78D4"
filter_other_grey <- "#7F7F7F"
filter_threshold_grey <- "#737373"
filter_region_blue <- "#EAF2F8"

# 单一绘图面板的断轴形式：对四个有数据的Y轴区间做分段线性映射。
# 断轴放在9–14、21–39和46–91这三个无观测值区间，避开要显示的刻度数字。
filter_y_segments <- data.frame(
  raw_min = c(0, 14, 39, 91),
  raw_max = c(9, 21, 46, 100),
  display_span = 6.25 * c(1.00, 0.45, 0.50, 0.40)
)
filter_break_gap <- 0.42
filter_y_segments$display_min <- c(
  0,
  cumsum(
    head(filter_y_segments$display_span, -1) + filter_break_gap
  )
)
filter_y_segments$display_max <-
  filter_y_segments$display_min + filter_y_segments$display_span

filter_y_map <- function(y) {
  mapped_y <- rep(NA_real_, length(y))
  for (i in seq_len(nrow(filter_y_segments))) {
    in_segment <- is.finite(y) &
      y >= filter_y_segments$raw_min[i] &
      y <= filter_y_segments$raw_max[i]
    mapped_y[in_segment] <- filter_y_segments$display_min[i] +
      (y[in_segment] - filter_y_segments$raw_min[i]) /
      (filter_y_segments$raw_max[i] - filter_y_segments$raw_min[i]) *
      filter_y_segments$display_span[i]
  }
  mapped_y
}

filter_omitted_ranges <- list(c(9, 14), c(21, 39), c(46, 91))
for (omitted_range in filter_omitted_ranges) {
  if (any(
    strategy_filter_df$count_median_all > omitted_range[1] &
      strategy_filter_df$count_median_all < omitted_range[2]
  )) {
    stop(
      sprintf(
        "The requested %g–%g y-axis break would hide one or more strategies.",
        omitted_range[1], omitted_range[2]
      )
    )
  }
}

# 点仅绘制一次；断轴映射不改变原始数据或筛选结果。
strategy_filter_df <- strategy_filter_df %>%
  mutate(y_plot = filter_y_map(count_median_all))
if (any(!is.finite(strategy_filter_df$y_plot))) {
  stop("One or more strategies fall outside the visible broken-axis intervals.")
}

filter_y_tick_raw <- c(0, 5, 15, 20, 40, 45, 95, 100)
filter_y_tick_plot <- filter_y_map(filter_y_tick_raw)
filter_break_positions <-
  (head(filter_y_segments$display_max, -1) +
     tail(filter_y_segments$display_min, -1)) / 2
filter_y_plot_max <- max(filter_y_segments$display_max)
filter_x_limits <- c(-0.01, 0.225)

# 左右边框的断轴白色遮罩和双斜线标记，保持上一版的跨轴形式。
filter_break_mask_df <- do.call(
  rbind,
  lapply(filter_x_limits, function(axis_x) {
    data.frame(
      xmin = axis_x - 0.0030,
      xmax = axis_x + 0.0030,
      ymin = filter_break_positions - 0.34,
      ymax = filter_break_positions + 0.34
    )
  })
)
filter_break_slash_df <- do.call(
  rbind,
  lapply(filter_x_limits, function(axis_x) {
    do.call(
      rbind,
      lapply(filter_break_positions, function(break_y) {
        data.frame(
          x = axis_x - 0.0027,
          xend = axis_x + 0.0027,
          y = break_y + c(-0.22, 0.05),
          yend = break_y + c(-0.04, 0.23)
        )
      })
    )
  })
)

filter_selected_label <- paste0("Selected strategies (n = ", selected_strategy_n, ")")
filter_other_label <- paste0("Other strategies (n = ", other_strategy_n, ")")
filter_group_levels <- c(filter_selected_label, filter_other_label)
strategy_filter_df <- strategy_filter_df %>%
  mutate(
    display_group = factor(
      if_else(selected, filter_selected_label, filter_other_label),
      levels = filter_group_levels
    )
  )

p_filter_base <- ggplot(
  strategy_filter_df,
  aes(x = jaccard_median, y = y_plot)
) +
  annotate(
    "rect",
    xmin = median_ji_filter_threshold,
    xmax = Inf,
    ymin = filter_y_map(median_da_filter_threshold),
    ymax = filter_y_plot_max,
    fill = filter_region_blue,
    alpha = 0.62
  ) +
  geom_vline(
    xintercept = median_ji_filter_threshold,
    colour = filter_threshold_grey,
    linewidth = 0.45,
    linetype = "dashed"
  ) +
  geom_hline(
    yintercept = filter_y_map(median_da_filter_threshold),
    colour = filter_threshold_grey,
    linewidth = 0.55,
    linetype = "dashed"
  ) +
  geom_point(
    aes(fill = display_group),
    shape = 21,
    size = 2.0,
    stroke = 0.35,
    colour = "white",
    alpha = 1
  ) +
  scale_fill_manual(
    name = NULL,
    values = setNames(
      c(filter_selected_blue, filter_other_grey),
      filter_group_levels
    )
  ) +
  annotate(
    "point", x = 0.082, y = 5.00,
    shape = 21, size = 2.0, stroke = 0.35,
    colour = "white", fill = filter_selected_blue
  ) +
  annotate(
    "text", x = 0.097, y = 5.00,
    label = filter_selected_label,
    hjust = 0, vjust = 0.5,
    family = "Arial", fontface = "bold", size = 6.2 / 2.845276
  ) +
  annotate(
    "point", x = 0.082, y = 4.25,
    shape = 21, size = 2.0, stroke = 0.35,
    colour = "white", fill = filter_other_grey, alpha = 1
  ) +
  annotate(
    "text", x = 0.097, y = 4.25,
    label = filter_other_label,
    hjust = 0, vjust = 0.5,
    family = "Arial", fontface = "bold", size = 6.2 / 2.845276
  ) +
  annotate(
    "text", x = 0.075, y = 2.65,
    label = "Selection criteria:",
    hjust = 0, vjust = 0.5,
    family = "Arial", fontface = "bold", size = 6.2 / 2.845276
  ) +
  annotate(
    "segment", x = 0.080, xend = 0.100, y = 1.65, yend = 1.65,
    colour = filter_threshold_grey, linewidth = 0.45, linetype = "dashed"
  ) +
  annotate(
    "text", x = 0.107, y = 1.65,
    label = "Median JI > 0",
    hjust = 0, vjust = 0.5,
    family = "Arial", size = 6.1 / 2.845276
  ) +
  annotate(
    "segment", x = 0.080, xend = 0.100, y = 0.70, yend = 0.70,
    colour = filter_threshold_grey, linewidth = 0.45, linetype = "dashed"
  ) +
  annotate(
    "text", x = 0.107, y = 0.70,
    label = "Median detected DA count > 5",
    hjust = 0, vjust = 0.5,
    family = "Arial", size = 6.1 / 2.845276
  ) +
  scale_x_continuous(
    breaks = c(0, 0.05, 0.10, 0.15, 0.20),
    labels = c("0", "0.05", "0.10", "0.15", "0.20"),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_continuous(
    limits = c(0, filter_y_plot_max),
    breaks = filter_y_tick_plot,
    labels = as.character(filter_y_tick_raw),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    title = "Prioritization of strategies\n(dual-criterion filter)",
    x = "Median Jaccard index (JI)",
    y = "Median detected DA count"
  ) +
  guides(fill = "none") +
  theme_classic(base_family = "Arial", base_size = 6.5) +
  theme(
    panel.border = element_blank(),
    axis.line = element_blank(),
    axis.title = element_text(face = "bold", size = 7.2, colour = "black"),
    axis.text = element_text(size = 6.2, colour = "black"),
    axis.ticks = element_line(colour = "black", linewidth = 0.35),
    axis.ticks.length = grid::unit(1.4, "mm"),
    plot.title = element_text(
      family = "Arial", face = "bold", size = 8,
      hjust = 0, colour = "black", margin = margin(b = 4)
    ),
    axis.text.y.right = element_blank(),
    axis.ticks.y.right = element_blank(),
    axis.title.y.right = element_blank(),
    legend.position = "none",
    plot.background = element_rect(fill = "white", colour = NA),
    plot.margin = margin(t = 4, r = 4, b = 3, l = 3)
  )

p_strategy_filter <- p_filter_base +
  annotate(
    "rect",
    xmin = filter_x_limits[1], xmax = filter_x_limits[2],
    ymin = 0, ymax = filter_y_plot_max,
    fill = NA, colour = "black", linewidth = 0.55
  ) +
  geom_rect(
    data = filter_break_mask_df,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = "white", colour = NA
  ) +
  geom_segment(
    data = filter_break_slash_df,
    aes(x = x, xend = xend, y = y, yend = yend),
    inherit.aes = FALSE,
    colour = "black", linewidth = 0.55, lineend = "butt"
  ) +
  coord_cartesian(
    xlim = filter_x_limits,
    ylim = c(0, filter_y_plot_max),
    clip = "off",
    expand = FALSE
  )

p_strategy_filter

filter_figure_width_mm <- 210 * 3 / 8  # 78.75 mm
filter_figure_width_in <- filter_figure_width_mm / 25.4
filter_figure_height_mm <- 90 * 6 / 7  # 在上一版基础上整体高度降低1/7（约77.14 mm）
filter_figure_height_in <- filter_figure_height_mm / 25.4

# 本轮仅导出PDF和PNG供预览确认。
grDevices::cairo_pdf(
  "strategy_dual_criterion_filter.pdf",
  width = filter_figure_width_in,
  height = filter_figure_height_in,
  family = "Arial"
)
print(p_strategy_filter)
dev.off()

filter_preview_png_tmp <- tempfile(fileext = ".png")
ragg::agg_png(
  filter_preview_png_tmp,
  width = filter_figure_width_in,
  height = filter_figure_height_in,
  units = "in",
  res = 300,
  background = "white"
)
print(p_strategy_filter)
dev.off()
file.copy(
  filter_preview_png_tmp,
  "strategy_dual_criterion_filter.png",
  overwrite = TRUE
)
unlink(filter_preview_png_tmp)

write.table(
  strategy_filter_df %>%
    transmute(
      Method = method,
      Median_Jaccard_index = jaccard_median,
      Median_detected_DA_count = count_median_all,
      Selected = selected
    ) %>%
    arrange(desc(Selected), desc(Median_Jaccard_index), desc(Median_detected_DA_count)),
  file = "strategy_dual_criterion_filter_data.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
# --- 策略筛选图结束 ---

########
library(dplyr)

summary_subset_res2 <- summary_subset_res %>%
  mutate(
    ji_sd = sqrt(jaccard_var),
    ji_lb = pmax(0, jaccard_mean - ji_sd),
    ji_robust = pmin(jaccard_median, ji_lb)
  )


min_count_med <- 5
min_sub_med   <- 5
min_ji_med    <- 0

kept <- summary_subset_res2 %>%
  filter(
    count_median_all >= min_count_med,
    subset1_median >= min_sub_med,
    subset2_median >= min_sub_med,
    jaccard_median >= min_ji_med
  ) %>%
  arrange(desc(ji_robust), desc(count_median_all), jaccard_var)

kept_methods <- kept$method
kept_methods
###########对数据集的可重复性溯源，是否存在某些数据集任何方法可重复性都很高
library(dplyr)

subset_res_dataset <- subset_res %>%
  mutate(
    dataset = paste(cohort, case, sep = " | "),
    mean_count = (mean_subset1_count + mean_subset2_count) / 2
  )
subset_res_dataset$method <- method_rename(subset_res_dataset$method)
subset_res_dataset <- subset_res_dataset[!(subset_res_dataset$method %in% kept_methods),]
subset_res_dataset <- subset_res_dataset[!(subset_res_dataset$method %in% bottom25_atleast2),]
# 数据集层面：跨方法汇总
subset_res_dataset_sum <- subset_res_dataset %>%
  group_by(dataset, cohort, case) %>%
  summarise(
    n_methods = n(),
    ji_mean   = mean(Jaccard_index, na.rm = TRUE),
    ji_median = median(Jaccard_index, na.rm = TRUE),
    ji_sd     = sd(Jaccard_index, na.rm = TRUE),
    prop_ji_ge_005 = mean(Jaccard_index >= 0.01, na.rm = TRUE),
    count_mean  = mean(mean_count, na.rm = TRUE),
    count_median= median(mean_count, na.rm = TRUE),
    prop_both_nonzero = mean(mean_subset1_count > 0 & mean_subset2_count > 0, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(ji_median), desc(prop_ji_ge_005), desc(count_median))
########################
#提取表现较好/不好的数据集
#count_median大于5，ji_median大于0
library(dplyr)
good_datasets <- subset_res_dataset_sum %>%
  filter(
    !is.na(count_median), !is.na(ji_median),
    count_median > 5,
    ji_median > 0
  ) %>%
  select(dataset, cohort, case)
bad_datasets <- subset_res_dataset_sum %>%
  filter(
    is.na(count_median) | is.na(ji_median) |
      count_median <= 5 |
      ji_median <= 0
  ) %>%
  select(dataset, cohort, case)

###########
subset_res_dataset_good <- subset_res_dataset[subset_res_dataset$dataset %in% good_datasets$dataset,]
subset_res_dataset_bad <- subset_res_dataset[subset_res_dataset$dataset %in% bad_datasets$dataset,]
summary_subset_res_dataset_good <- subset_res_dataset_good %>%
  group_by(method) %>%
  summarise(
    n = n(),

    # Jaccard
    jaccard_mean = mean(Jaccard_index, na.rm = TRUE),
    jaccard_median = median(Jaccard_index, na.rm = TRUE),
    jaccard_var  = var(Jaccard_index,  na.rm = TRUE),

    # subset1 count
    subset1_mean = mean(mean_subset1_count, na.rm = TRUE),
    subset1_median = median(mean_subset1_count, na.rm = TRUE),
    subset1_var  = var(mean_subset1_count,  na.rm = TRUE),

    # subset2 count
    subset2_mean = mean(mean_subset2_count, na.rm = TRUE),
    subset2_median = median(mean_subset2_count, na.rm = TRUE),
    subset2_var  = var(mean_subset2_count,  na.rm = TRUE),

    # 两列 count 一起（把两列拼起来后再算均值/方差）
    count_mean_all = mean(c(mean_subset1_count, mean_subset2_count), na.rm = TRUE),
    count_median_all = median(c(mean_subset1_count, mean_subset2_count), na.rm = TRUE),
    count_var_all  = var(c(mean_subset1_count, mean_subset2_count),  na.rm = TRUE),

    .groups = "drop"
  )
summary_subset_res_dataset_bad <- subset_res_dataset_bad %>%
  group_by(method) %>%
  summarise(
    n = n(),

    # Jaccard
    jaccard_mean = mean(Jaccard_index, na.rm = TRUE),
    jaccard_median = median(Jaccard_index, na.rm = TRUE),
    jaccard_var  = var(Jaccard_index,  na.rm = TRUE),

    # subset1 count
    subset1_mean = mean(mean_subset1_count, na.rm = TRUE),
    subset1_median = median(mean_subset1_count, na.rm = TRUE),
    subset1_var  = var(mean_subset1_count,  na.rm = TRUE),

    # subset2 count
    subset2_mean = mean(mean_subset2_count, na.rm = TRUE),
    subset2_median = median(mean_subset2_count, na.rm = TRUE),
    subset2_var  = var(mean_subset2_count,  na.rm = TRUE),

    # 两列 count 一起（把两列拼起来后再算均值/方差）
    count_mean_all = mean(c(mean_subset1_count, mean_subset2_count), na.rm = TRUE),
    count_median_all = median(c(mean_subset1_count, mean_subset2_count), na.rm = TRUE),
    count_var_all  = var(c(mean_subset1_count, mean_subset2_count),  na.rm = TRUE),

    .groups = "drop"
  )
####bubble_plot
bubble_plot_order <- summary_subset_res_dataset_good$method[order(summary_subset_res_dataset_good$jaccard_median,decreasing = FALSE)]

library(dplyr)
library(ggplot2)
plot_df_good <- summary_subset_res_dataset_good %>%
  select(method, jaccard_median, count_median_all) %>%
  rename(median_count_all = count_median_all)
plot_df_good$method <- factor(plot_df_good$method,levels = bubble_plot_order)


plot_df_bad <- summary_subset_res_dataset_bad %>%
  select(method, jaccard_median, count_median_all) %>%
  rename(median_count_all = count_median_all)
plot_df_bad$method <- factor(plot_df_bad$method,levels = bubble_plot_order)
library(ggplot2)

# 1) 计算全局 size 范围（跨 good + bad）
all_size <- c(plot_df_bad$median_count_all, plot_df_good$median_count_all)
size_lim <- range(all_size, na.rm = TRUE)

# 可选：统一 legend 的刻度（也可以自己指定更“好看”的几个值）
size_brk <- pretty(size_lim, n = 4)

# 2) 定义一个共同的 size scale（两张图复用）
size_scale_shared <- scale_size_continuous(
  limits = size_lim,          # 关键：统一映射范围
  breaks = size_brk,          # 关键：统一 legend 刻度
  range  = c(0.6, 5),         # 你的气泡显示范围（视觉大小）
)

# 3) 在两张图里都加同一个 scale
p_bubble_bad <- ggplot(plot_df_bad,
                       aes(x = jaccard_median, y = method, size = median_count_all)
) +
  geom_point(alpha = 0.7, color = "#E9B99B") +
  size_scale_shared +
  scale_x_continuous(limits = c(0, 0.5)) +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "top") +
  labs(y = NULL, x = "Median Jaccard index", size = "Median count (all)")

p_bubble_good <- ggplot(plot_df_good,
                        aes(x = jaccard_median, y = method, size = median_count_all)
) +
  geom_point(alpha = 0.7, color = "#E9B99B") +
  size_scale_shared +
  scale_x_continuous(limits = c(0, 0.5)) +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "top") +
  labs(y = NULL, x = "Median Jaccard index", size = "Median count (all)")

############################################################################
#########################***以下代码不运行***###############################
############################################################################
# #高重复数据集大多伴随着高检出，画图验证
# library(ggplot2)
# library(ggrepel)
#
# top_label <- subset_res_dataset_sum %>% slice_max(ji_median, n = 10)
#
# ggplot(subset_res_dataset_sum, aes(x = count_median, y = ji_median)) +
#   geom_point(alpha = 0.7,color = "#4D7EAB") +
#   geom_text_repel(data = top_label, aes(label = dataset), size = 3, max.overlaps = 30) +
#   labs(x = "Median detected DA count", y = "Median Jaccard index") +
#   theme_classic()
# ggplot(subset_res_dataset_sum, aes(x = count_mean, y = ji_mean)) +
#   geom_point(alpha = 0.7,color = "#4D7EAB") +
#   geom_text_repel(data = top_label, aes(label = dataset), size = 3, max.overlaps = 30) +
#   labs(x = "Mean detected DA count", y = "Mean Jaccard index") +
#   theme_classic()
# #数据集难度分层
# #选择前三个作为高重复数据集
# high_repro_cohort <- subset_res_dataset_sum$cohort[1:3]
# high_repro_case <- subset_res_dataset_sum$case[1:3]
# high_repro_subset_res <- subset_res[which((subset_res$cohort %in% high_repro_cohort)&(subset_res$case %in% unique(high_repro_case))),]
# summary_high_repro_subset_res <- high_repro_subset_res %>%
#   group_by(method) %>%
#   summarise(
#     n = n(),
#
#     # Jaccard
#     jaccard_mean = mean(Jaccard_index, na.rm = TRUE),
#     jaccard_median = median(Jaccard_index, na.rm = TRUE),
#     jaccard_var  = var(Jaccard_index,  na.rm = TRUE),
#
#     # subset1 count
#     subset1_mean = mean(mean_subset1_count, na.rm = TRUE),
#     subset1_median = median(mean_subset1_count, na.rm = TRUE),
#     subset1_var  = var(mean_subset1_count,  na.rm = TRUE),
#
#     # subset2 count
#     subset2_mean = mean(mean_subset2_count, na.rm = TRUE),
#     subset2_median = median(mean_subset2_count, na.rm = TRUE),
#     subset2_var  = var(mean_subset2_count,  na.rm = TRUE),
#
#     # 两列 count 一起（把两列拼起来后再算均值/方差）
#     count_mean_all = mean(c(mean_subset1_count, mean_subset2_count), na.rm = TRUE),
#     count_median_all = median(c(mean_subset1_count, mean_subset2_count), na.rm = TRUE),
#     count_var_all  = var(c(mean_subset1_count, mean_subset2_count),  na.rm = TRUE),
#
#     .groups = "drop"
#   )
# summary_high_repro_subset_res$method <- method_rename(summary_high_repro_subset_res$method)
# summary_high_repro_subset_res <- summary_high_repro_subset_res[!(summary_high_repro_subset_res$method %in% filter_method),]
# #####低重复数据集
# low_repro_subset_res <- subset_res[which(!((subset_res$cohort %in% high_repro_cohort)&(subset_res$case %in% unique(high_repro_case)))),]
# summary_low_repro_subset_res <- low_repro_subset_res %>%
#   group_by(method) %>%
#   summarise(
#     n = n(),
#
#     # Jaccard
#     jaccard_mean = mean(Jaccard_index, na.rm = TRUE),
#     jaccard_median = median(Jaccard_index, na.rm = TRUE),
#     jaccard_var  = var(Jaccard_index,  na.rm = TRUE),
#
#     # subset1 count
#     subset1_mean = mean(mean_subset1_count, na.rm = TRUE),
#     subset1_median = median(mean_subset1_count, na.rm = TRUE),
#     subset1_var  = var(mean_subset1_count,  na.rm = TRUE),
#
#     # subset2 count
#     subset2_mean = mean(mean_subset2_count, na.rm = TRUE),
#     subset2_median = median(mean_subset2_count, na.rm = TRUE),
#     subset2_var  = var(mean_subset2_count,  na.rm = TRUE),
#
#     # 两列 count 一起（把两列拼起来后再算均值/方差）
#     count_mean_all = mean(c(mean_subset1_count, mean_subset2_count), na.rm = TRUE),
#     count_median_all = median(c(mean_subset1_count, mean_subset2_count), na.rm = TRUE),
#     count_var_all  = var(c(mean_subset1_count, mean_subset2_count),  na.rm = TRUE),
#
#     .groups = "drop"
#   )
# summary_low_repro_subset_res$method <- method_rename(summary_low_repro_subset_res$method)
# summary_low_repro_subset_res <- summary_low_repro_subset_res[!(summary_low_repro_subset_res$method %in% filter_method),]
# ########相关性分析
# all_rank_info <- summary_subset_res[,c(1,3,4)]
# order <- all_rank_info$method
# high_rank_info <-  summary_high_repro_subset_res[,c(1,3,4)]
# low_rank_info <- summary_low_repro_subset_res[,c(1,3,4)]
# rownames(high_rank_info) <- high_rank_info$method
# rownames(low_rank_info) <- low_rank_info$method
# high_rank_info <- high_rank_info[order,]
# low_rank_info <- low_rank_info[order,]
# cor.test(high_rank_info$jaccard_mean,all_rank_info$jaccard_mean)
# cor.test(high_rank_info$jaccard_median,all_rank_info$jaccard_median)
# cor.test(low_rank_info$jaccard_mean,all_rank_info$jaccard_mean)
# cor.test(low_rank_info$jaccard_median,all_rank_info$jaccard_median)
# cor.test(high_rank_info$jaccard_mean,low_rank_info$jaccard_mean)
# cor.test(high_rank_info$jaccard_median,low_rank_info$jaccard_median)
# #####high low排名相关性示意图
# library(ggplot2)
# library(ggpubr)
# library(ggpmisc)
# df1 <- data.frame(high_mean = high_rank_info$jaccard_mean,
#                   low_mean =  low_rank_info$jaccard_mean)
# df2 <- data.frame(high_median = high_rank_info$jaccard_median,
#                   low_median =  low_rank_info$jaccard_median)
# theme_set(ggpubr::theme_pubr()+
#             theme(legend.position = "top"))
# ggplot(df1, aes(x = high_mean, y = low_mean)) + geom_point(color='#3B99B1')+
#   geom_smooth(method = "lm", color = "#0C784F", fill = "#FEFEE3") +
#   stat_cor(method = "pearson",label.x = 0.1, label.y = 0.7) + labs(x = 'Mean Jaccard index in High-reproductivity Datasets',y = 'Median Jaccard index in Low-reproductivity Datasets')+
#   coord_cartesian(ylim = c(0, 0.1))
# ggplot(df2, aes(x = high_median, y = low_median)) + geom_point(color='#3B99B1')+
#   geom_smooth(method = "lm", color = "#0C784F", fill = "#FEFEE3") +
#   stat_cor(method = "pearson",label.x = 0.1, label.y = 0.7) + labs(x = 'Median Jaccard index in High-reproductivity Datasets',y = 'Median Jaccard index in Low-reproductivity Datasets')+
#   coord_cartesian(ylim = c(0, 0.1))
# #分层展示jaccard_median和mean
# sens_method_order <- low_rank_info$method[order(low_rank_info$jaccard_median,decreasing = TRUE)]
# high_rank_info$method <- factor(high_rank_info$method,levels = rev(sens_method_order))
# low_rank_info$method <- factor(high_rank_info$method,levels = rev(sens_method_order))
# library(dplyr)
# library(tidyr)
# library(ggplot2)
# plot_df <- high_rank_info %>%
#   select(method, jaccard_mean, jaccard_median) %>%
#   pivot_longer(
#     cols = c(jaccard_mean, jaccard_median),
#     names_to = "metric",
#     values_to = "value"
#   ) %>%
#   mutate(
#     metric = factor(metric, levels = c("jaccard_median", "jaccard_mean"))
#   )
#
# p_high_ji_barplot <- ggplot(plot_df, aes(x = method, y = value, fill = metric)) +
#   geom_col(position = position_dodge(width = 0.8), width = 0.75) +
#   theme_classic(base_size = 12) + coord_flip() +
#   scale_fill_manual(
#     values = c(
#       jaccard_mean   = "#6A6599",
#       jaccard_median   = "#79AF97"
#     )
#   ) +
#   theme(
#     axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
#     legend.position = "top"
#   ) +
#   labs(x = NULL, y = "Jaccard index", fill = NULL)
#
# p_high_ji_barplot
#
# plot_df <- low_rank_info %>%
#   select(method, jaccard_mean, jaccard_median) %>%
#   pivot_longer(
#     cols = c(jaccard_mean, jaccard_median),
#     names_to = "metric",
#     values_to = "value"
#   ) %>%
#   mutate(
#     metric = factor(metric, levels = c("jaccard_median", "jaccard_mean"))
#   )
#
# p_low_ji_barplot <- ggplot(plot_df, aes(x = method, y = value, fill = metric)) +
#   geom_col(position = position_dodge(width = 0.8), width = 0.75) +
#   theme_classic(base_size = 12) + coord_flip() +
#   scale_fill_manual(
#     values = c(
#       jaccard_mean   = "#6A6599",
#       jaccard_median   = "#79AF97"
#     )
#   ) +
#   theme(
#     axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
#     legend.position = "top"
#   ) +
#   labs(x = NULL, y = "Jaccard index", fill = NULL)
#
# p_low_ji_barplot
########
########
########
#####cross图
library(dplyr)
library(tidyr)

cross_res <- cross_res %>%
  # 保险起见：先把目标列转成 numeric（如果读进来是字符/因子）
  mutate(
    Jaccard_index     = as.numeric(Jaccard_index),
    mean_pairwise_ji  = as.numeric(mean_pairwise_ji)
  ) %>%
  # 只对数值列：NA -> 0
  mutate(across(where(is.numeric), ~replace_na(.x, 0)))

summary_cross_res <- cross_res %>%
  group_by(method) %>%
  summarise(
    n = n(),
    Jaccard_index_mean    = mean(Jaccard_index),
    Jaccard_index_var     = ifelse(n() > 1, var(Jaccard_index), 0),
    mean_pairwise_ji_mean = mean(mean_pairwise_ji),
    mean_pairwise_ji_var  = ifelse(n() > 1, var(mean_pairwise_ji), 0),
    .groups = "drop"
  )

summary_cross_res$method <-method_rename(summary_cross_res$method)
summary_cross_res <- summary_cross_res[!(summary_cross_res$method %in% filter_set2),]

###case_summary_res
case_cross_res <- cross_res
case_cross_res$method <-method_rename(case_cross_res$method)
case_cross_res <- case_cross_res[!(case_cross_res$method %in% filter_set2),]
summary_case_cross_res <- case_cross_res %>%
  group_by(case) %>%
  summarise(
    n = n(),
    Jaccard_index_mean    = mean(Jaccard_index),
    Jaccard_index_var     = ifelse(n() > 1, var(Jaccard_index), 0),
    mean_pairwise_ji_mean = mean(mean_pairwise_ji),
    mean_pairwise_ji_var  = ifelse(n() > 1, var(mean_pairwise_ji), 0),
    .groups = "drop"
  )
summary_case_cross_res_j <- summary_case_cross_res %>%
  mutate(
    jaccard_sd = sqrt(Jaccard_index_var),
    case = as.character(case)
  ) %>%
  arrange(desc(Jaccard_index_mean))
p_case_cross_jaccard <- ggplot(summary_case_cross_res_j, aes(y = Jaccard_index_mean, x = case)) +
  geom_point(size = 3, color = "#C7A76C") +
  geom_errorbar(
    aes(
      ymin = pmax(Jaccard_index_mean - jaccard_sd, 0),
      ymax = Jaccard_index_mean + jaccard_sd
    ),
    width = 0.3,color = '#A3A3A3'
  ) +
  coord_flip() +
  theme_classic(base_size = 12) +
  labs(x = NULL, y = "Jaccard index (mean ± SD)")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

p_case_cross_jaccard

#cross画图
#global_ji
library(dplyr)
library(ggplot2)
summary_cross_res_j <- summary_cross_res %>%
  mutate(
    jaccard_sd = sqrt(Jaccard_index_var),
    method = as.character(method)
  ) %>%
  arrange(desc(Jaccard_index_mean)) %>%
  mutate(method = factor(method, levels = method))

p_cross_jaccard <- ggplot(summary_cross_res_j, aes(x = Jaccard_index_mean, y = method)) +
  geom_point(size = 3, color = "#C7A76C") +
  geom_errorbar(
    aes(
      xmin = pmax(Jaccard_index_mean - jaccard_sd, 0),
      xmax = Jaccard_index_mean + jaccard_sd
    ),
    width = 0.3,color = '#A3A3A3'
  ) +
  coord_flip() +
  theme_classic(base_size = 12) +
  labs(x = NULL, y = "Jaccard index (mean ± SD)")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

p_cross_jaccard

summary_cross_res_mj <- summary_cross_res %>%
  mutate(
    jaccard_sd = sqrt(mean_pairwise_ji_var),
    method = as.character(method)
  ) %>%
  arrange(desc(mean_pairwise_ji_mean)) %>%
  mutate(method = factor(method, levels = method))

p_cross_mjaccard <- ggplot(summary_cross_res_mj, aes(x = mean_pairwise_ji_mean, y = method)) +
  geom_point(size = 3, color = "#C7A76C") +
  geom_errorbar(
    aes(
      xmin = pmax(mean_pairwise_ji_mean - jaccard_sd, 0),
      xmax = mean_pairwise_ji_mean + jaccard_sd
    ),
    width = 0.3,color = '#A3A3A3'
  ) +
  coord_flip() +
  theme_classic(base_size = 12) +
  labs(x = NULL, y = "Pairwise Jaccard index (mean ± SD)")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

p_cross_mjaccard
################################该图用作summary
#########cross_res采用funkyheatmap展示
library(tidyr)
cross_res <- spread(cross_res,key = method,value = Jaccard_index)
cross_res <- cross_res[,-1]
cross_res <- data.frame(cross_res)
library(funkyheatmap)
library(dplyr, warn.conflicts = FALSE)
library(tibble, warn.conflicts = FALSE)
library(paletteer)
palettes <- tribble(
  ~ palette,
  ~ colours,
  'ANCOMBC_count',
  paletteer_c("grDevices::Peach", 30),
  'ANCOMBC_CSS',
  paletteer_c("grDevices::Purp", 30),
  'ANCOMBC_TMM',
  paletteer_c("grDevices::Reds 2", 30),
  "edgeR_CSS",
  rev(paletteer_c("ggthemes::Blue", 30)),
  "fastANCOM_TMM",
  rev(paletteer_c("ggthemes::Classic Area-Brown", 30)),
  "fastANCOM_CSS",
  rev(paletteer_c("ggthemes::Purple", 30)),
  "lfem_CSS",
  paletteer_c("ggthemes::Orange Light", 30),
  "limma_CSS",
  paletteer_c("grDevices::YlOrBr", 30),
  "lmem_CSS",
  rev(paletteer_c("ggthemes::Blue-Teal", 30)),
  "maaslin2_CSS_LOG",
  paletteer_c("grDevices::BrwnYl", 30),
  "maaslin2_TMM_LOG",
  rev(paletteer_c("ggthemes::Classic Area Green", 30))
)

colnames(cross_res)[1] <- 'id'
cross_res <- cross_res[,c('id',order)]
# 提取列名
original_colnames <- colnames(cross_res)[-1]
show_names <- original_colnames
show_names[8:9] <- c('maaslin2_lm_CSS_LOG','maaslin2_lm_TMM_LOG')

cross_res[,2:11] <- apply(cross_res[,2:11],2,as.numeric)
cross_res[,2:11] <- t(apply(cross_res[,2:11],1,scale_minmax))
# 重新设置列名顺序
column_info <- data.frame(id = colnames(cross_res),group = 'Case',width=c(3.5,rep(1.5,10)),
                          name = c('',show_names),
                          geom = c('text',rep('bar',ncol(cross_res)-1)),palette  = c(NA,original_colnames))

funky_heatmap(cross_res,column_info = column_info, palettes = palettes,scale_column = F)
########
subset_res$index <- paste(subset_res$cohort,subset_res$case,sep = '_')
subset_res_wide<-spread(subset_res[,c('Jaccard_index','method','index')],key = index,value = Jaccard_index)
row.names(subset_res_wide) <- subset_res_wide$method
subset_res_wide<-t(subset_res_wide[,-1])
subset_res_wide<-cbind(row.names(subset_res_wide),subset_res_wide)
colnames(subset_res_wide)[1] <- 'id'
subset_res_wide <- subset_res_wide[,c('id',order)]
# 提取列名
original_colnames <- colnames(subset_res_wide)[-1]
show_names <- original_colnames
show_names[8:9] <- c('maaslin2_lm_CSS_LOG','maaslin2_lm_TMM_LOG')
subset_res_wide[is.na(subset_res_wide)] <- 0
subset_res_wide <- data.frame(subset_res_wide)
subset_res_wide[,2:11] <- apply(subset_res_wide[,2:11],2,as.numeric)
subset_res_wide[,2:11] <- t(apply(subset_res_wide[,2:11],1,scale_minmax))
# 重新设置列名顺序

column_info <- data.frame(id = colnames(subset_res_wide),group = 'Case',width=c(12,rep(1.8,10)),
                          name = c('',show_names),
                          geom = c('text',rep('bar',ncol(subset_res_wide)-1)),palette  = c(NA,original_colnames))

funky_heatmap(subset_res_wide,column_info = column_info, palettes = palettes,scale_column = F)

######
method_sum <- subset_res %>%
  group_by(method) %>%
  summarise(sum = sum(Jaccard_index, na.rm = TRUE))
method_mean <- subset_res %>%
  group_by(method) %>%
  summarise(mean = mean(Jaccard_index, na.rm = TRUE))
rownames(method_sum) <- method_sum$method
method_sum <- method_sum[order,]
ggdotchart(method_sum,
           x = "method",
           y = "sum",
           color = "#A9C0A9",
           sorting = "none",
           add = "segments",
           dot.size = 5,
           add.params = list(color = "lightgray", size = 2.5),
           ggtheme = theme_pubclean()) +
  theme(
    axis.text.y = element_text(face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 0, face = "bold")
  ) +
  labs(y = "The sum of Jaccard index", x = NULL)
########output sourced data
subset_res_output <- subset_res
colnames(subset_res_output) <- c("Cohort","Case","Mean Jaccard index","Mean_subset1_count","Mean_subset2_count","Method")
subset_res_output <- subset_res_output[,c("Method","Cohort","Case","Mean Jaccard index","Mean_subset1_count","Mean_subset2_count")]
subset_res_output$Mean_DA_count <- (subset_res_output$Mean_subset1_count + subset_res_output$Mean_subset2_count)/2
subset_res_output <- subset_res_output[,c("Method","Cohort","Case","Mean Jaccard index","Mean_DA_count")]
colnames(subset_res_output)[5] <- 'Mean DA Count Across replicates'
summary_subset_res_output <- summary_subset_res_j
summary_subset_res_output <- summary_subset_res_output[,c('method',"n","jaccard_mean","jaccard_sd","jaccard_median","count_mean_all","count_median_all")]
colnames(summary_subset_res_output) <- c('Method',"Number of datasets included","Mean Jaccard Index","SD of Jaccard Index","Median Jaccard Index","Mean DA Count Across Datasets","Median DA Count Across Datasets")
summary_subset_res_output$Method <- method_rename(summary_subset_res_output$Method)
subset_res_output$Method <- method_rename(subset_res_output$Method)
write.table(subset_res_output,file = 'Consistency_performance.tsv',sep = '\t',col.names = T,row.names = F,quote = F)
write.table(summary_subset_res_output,file = 'Consistency_Summary_performance.tsv',sep = '\t',col.names = T,row.names = F,quote = F)
