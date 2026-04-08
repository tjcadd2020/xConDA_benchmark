#####构建上下调方向函数
effect_size_dir <- function(result,
                            feature_col = "Feature",
                            padj_col    = "p.adj.val") {
  
  # --- Feature id ---
  feat <- if (feature_col %in% names(result)) as.character(result[[feature_col]]) else rownames(result)
  if (is.null(feat)) stop("Cannot find Feature identifiers (Feature column or rownames).")
  
  # --- padj ---
  padj <- if (padj_col %in% names(result)) result[[padj_col]] else rep(NA_real_, nrow(result))
  
  # --- init dir ---
  dir <- rep(0L, nrow(result))
  names(dir) <- feat
  ok <- !is.na(padj)
  
  # --- 1) effect size column (priority order) ---
  effect_size_candidates <- c("Group.Est","lfc_Group1","effect_size","log2FoldChange",
                              "logFC","log2FC","Effect_size","log2FC","logFC")
  es_col <- intersect(effect_size_candidates, names(result))
  
  if (length(es_col) >= 1) {
    es <- as.numeric(result[[es_col[1]]])
    dir[ok] <- ifelse(is.na(es[ok]) | es[ok] == 0, 0L, ifelse(es[ok] > 0, 1L, -1L))
    return(dir)
  }
  
  # --- 2) Enriched/Enrieched column ---
  enr_col <- names(result)[grepl("^Enri(e)?ched$", names(result), ignore.case = TRUE)]
  if (length(enr_col) >= 1) {
    enr <- as.character(result[[enr_col[1]]])
    dir[ok] <- ifelse(enr[ok] == "Disease", 1L, -1L)
    return(dir)
  }
  
  # --- 3) coef column ---
  if ("coef" %in% names(result)) {
    cf <- as.numeric(result[["coef"]])
    dir[ok] <- ifelse(is.na(cf[ok]) | cf[ok] == 0, 0L, ifelse(cf[ok] > 0, 1L, -1L))
    return(dir)
  }
  
  dir
}

# Initialize an empty data frame to store results for all methods
all_results <- data.frame(cohort = character(0), 
                          case = character(0), 
                          Jaccard_index = numeric(), 
                          mean_subset1_count = numeric(),
                          mean_subset2_count = numeric(),
                          method = character(0))

# Loop over each method
for (key in method_list_keys) {
  if(grepl("None", method_list[[key]][3])){
    method_list[[key]] = method_list[[key]][1:2]
  }
  # Loop over each cohort
  for (single_cohort in single_validation_cohort) {
    
    # Extract the group cases for the current cohort
    group_case <- covariates_information2$Group_val[which(covariates_information2$cohort == single_cohort)]
    
    # Loop over each case in the group
    for (case in group_case) {
      
      Jaccard_index_list <- numeric()  # Initialize the list to store Jaccard indices
      subset1_count_list <- numeric()  # Initialize the list to store subset1 significant count
      subset2_count_list <- numeric()  # Initialize the list to store subset2 significant count
      
      # Loop over the 10 splits
      for (p in 1:10) {
        
        # Construct file names for subset1 and subset2
        subset1_file_name <- paste0(paste(method_list[[key]], collapse = '_'), '_res_', single_cohort, '_', case, 
                                    '_splits', p, '_subset1_.tsv')
        subset2_file_name <- paste0(paste(method_list[[key]], collapse = '_'), '_res_', single_cohort, '_', case, 
                                    '_splits', p, '_subset2_.tsv')
        
        # Read the data from the files
        subset1 <- read.table(file = subset1_file_name, sep = '\t', header = TRUE)
        subset2 <- read.table(file = subset2_file_name, sep = '\t', header = TRUE)
        
        # Identify significant features in each subset
        if(grepl("maaslin2", method_list[[key]][1])){
          subset1 <- subset1[which(subset1$metadata == 'Group'),]
          subset2 <- subset2[which(subset2$metadata == 'Group'),]
        }
        #top_percent <- 0.05 # 取前10%的显著特征
        #subset1_threshold <- quantile(cohort_res$pval, probs = top_percent, na.rm = TRUE)
        dir1 <- effect_size_dir(subset1) 
        dir2 <- effect_size_dir(subset2)
        
        sig1 <- subset1$Feature[subset1$p.adj.val < 0.05 & !is.na(subset1$p.adj.val)]
        sig2 <- subset2$Feature[subset2$p.adj.val < 0.05 & !is.na(subset2$p.adj.val)]
        
        # 只保留方向明确的（dir!=0）
        sig1 <- sig1[!is.na(dir1[sig1]) & dir1[sig1] != 0]
        sig2 <- sig2[!is.na(dir2[sig2]) & dir2[sig2] != 0]
        
        set1 <- paste0(sig1, ":", ifelse(dir1[sig1] > 0, "up", "down"))
        set2 <- paste0(sig2, ":", ifelse(dir2[sig2] > 0, "up", "down"))
        if(length(set1) == 1){
          if(set1 == ':'){
            set1 <- NULL
          }
        }
        if(length(set2) == 1){
          if(set2 == ':'){
            set2 <- NULL
          }
        }
        
        den <- length(union(set1, set2))
        Jaccard_index <- if (den == 0) 0 else length(intersect(set1, set2)) / den
        
        # Store the number of significant features in each subset
        subset1_count <- length(sig1)
        subset2_count <- length(sig2)
        
        # Append to lists
        Jaccard_index_list <- c(Jaccard_index_list, Jaccard_index)
        subset1_count_list <- c(subset1_count_list, subset1_count)
        subset2_count_list <- c(subset2_count_list, subset2_count)
      }
      
      # Compute mean values
      Jaccard_index_mean <- mean(Jaccard_index_list)
      mean_subset1_count <- mean(subset1_count_list)
      mean_subset2_count <- mean(subset2_count_list)
      
      # Append the result to the all_results data frame
      all_results <- rbind(all_results, data.frame(cohort = single_cohort, 
                                                   case = case, 
                                                   Jaccard_index = Jaccard_index_mean, 
                                                   mean_subset1_count = mean_subset1_count,
                                                   mean_subset2_count = mean_subset2_count,
                                                   method = paste0(method_list[[key]], collapse = '_')))
    }
  }
}