library(dplyr)
library(stringr)

#-----------------------------
# 1. Jaccard Index
#-----------------------------
jaccard_index <- function(set1, set2) {
  set1 <- unique(set1[!is.na(set1) & set1 != ""])
  set2 <- unique(set2[!is.na(set2) & set2 != ""])
  
  union_set <- union(set1, set2)
  if (length(union_set) == 0) return(NA_real_)
  
  inter_set <- intersect(set1, set2)
  length(inter_set) / length(union_set)
}

#-----------------------------
# 2. 从单个结果文件中提取显著特征
#    条件：p.adj.val < 0.05，且杜绝 NA
#-----------------------------
get_sig_features <- function(file) {
  df <- read.delim(file, stringsAsFactors = FALSE, check.names = FALSE)
  
  required_cols <- c("Feature", "p.adj.val")
  missing_cols <- setdiff(required_cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop(
      paste0(
        "File ", basename(file), 
        " is missing required columns: ", 
        paste(missing_cols, collapse = ", ")
      )
    )
  }
  
  sig_features <- df %>%
    filter(
      !is.na(p.adj.val),
      p.adj.val < 0.05,
      !is.na(Feature),
      Feature != ""
    ) %>%
    pull(Feature) %>%
    unique()
  
  sig_features
}

#-----------------------------
# 3. 获取所有符合模式的文件
#    文件格式：
#    单队列策略_res_疾病_repX_splitY_meta策略.tsv
#    例如：
#    edgeR_count_res_T2D_rep9_split1_Fisher.tsv
#    MMUPHin_res_CRC_rep4_split1_rma.uni.tsv
#-----------------------------
file_list <- list.files(
  path = ".",
  pattern = "_res_(CRC|T2D|adenoma|IBD)_rep\\d+_split[12]_.+\\.tsv$",
  full.names = TRUE
)

#-----------------------------
# 4. 解析文件名
#-----------------------------
# 分组说明：
# 1 = Single_strategy
# 2 = Disease
# 3 = Rep编号
# 4 = Split编号
# 5 = Meta_strategy
filename_pattern <- "^(.*)_res_([^_]+)_rep(\\d+)_split([12])_(.+)\\.tsv$"

match_mat <- str_match(basename(file_list), filename_pattern)

file_info <- data.frame(
  filepath = file_list,
  filename = basename(file_list),
  Single_strategy = match_mat[, 2],
  Disease = match_mat[, 3],
  Rep = paste0("rep", match_mat[, 4]),
  Split = paste0("split", match_mat[, 5]),
  Meta_strategy = match_mat[, 6],
  stringsAsFactors = FALSE
) %>%
  filter(!is.na(Single_strategy)) %>%
  mutate(
    Strategy = paste(Single_strategy, Meta_strategy, sep = "_")
  )

# 检查解析结果
print(file_info)

#-----------------------------
# 5. 构造唯一组合
#    按 Disease + Single_strategy + Meta_strategy 来配对 split1/split2
#-----------------------------
strategy_combos <- file_info %>%
  distinct(Disease, Single_strategy, Meta_strategy, Strategy)

#-----------------------------
# 6. 逐组合、逐rep计算 JI
#-----------------------------
results <- list()
idx <- 1

for (i in seq_len(nrow(strategy_combos))) {
  disease_i <- strategy_combos$Disease[i]
  single_i  <- strategy_combos$Single_strategy[i]
  meta_i    <- strategy_combos$Meta_strategy[i]
  strat_i   <- strategy_combos$Strategy[i]
  
  sub_info <- file_info %>%
    filter(
      Disease == disease_i,
      Single_strategy == single_i,
      Meta_strategy == meta_i
    )
  
  rep_ids <- sort(unique(sub_info$Rep))
  
  for (rep_id in rep_ids) {
    split1_file <- sub_info %>%
      filter(Rep == rep_id, Split == "split1") %>%
      pull(filepath)
    
    split2_file <- sub_info %>%
      filter(Rep == rep_id, Split == "split2") %>%
      pull(filepath)
    
    if (length(split1_file) == 1 && length(split2_file) == 1) {
      set1 <- get_sig_features(split1_file)
      set2 <- get_sig_features(split2_file)
      
      ji <- jaccard_index(set1, set2)
      
      results[[idx]] <- data.frame(
        Disease = disease_i,
        Single_strategy = single_i,
        Meta_strategy = meta_i,
        Strategy = strat_i,
        Rep = rep_id,
        JI = ji,
        N_sig_split1 = length(set1),
        N_sig_split2 = length(set2),
        stringsAsFactors = FALSE
      )
      idx <- idx + 1
    } else {
      warning(
        paste(
          "Missing file(s) for:",
          "Disease =", disease_i,
          "| Single_strategy =", single_i,
          "| Meta_strategy =", meta_i,
          "|", rep_id
        )
      )
    }
  }
}

results_df <- bind_rows(results)

#-----------------------------
# 7. 汇总：每个 Disease + Strategy 的平均 JI
#-----------------------------
summary_df <- results_df %>%
  mutate(
    DA_count = (N_sig_split1 + N_sig_split2) / 2
  ) %>%
  group_by(Disease, Single_strategy, Meta_strategy, Strategy) %>%
  summarise(
    Mean_Jaccard_index = if (all(is.na(JI))) NA_real_ else mean(JI, na.rm = TRUE),
    Mean_DA_count = if (all(is.na(DA_count))) NA_real_ else mean(DA_count, na.rm = TRUE),
    N_valid_rep = sum(!is.na(JI)),
    .groups = "drop"
  )