evaluation_enrichment2<-function(result_file_path,pval_threshold){
  result <- read.table(file = result_file_path,sep = '\t',header = T)
  effect_size_name <- unique(c('Group.Est','lfc_Group1','effect_size','log2FoldChange','logFC','log2FC','logFC','Effect_size','effect_size'))
  
  if(sum(grepl(colnames(result),pattern = paste0("^",effect_size_name,"$",collapse = "|")))==1){
    
    NA_feature <- result$Feature[is.na(result$p.adj.val)]
    result <- result[!is.na(result$p.adj.val),]
    effect_size <- result[,grep(colnames(result),pattern = paste("^",effect_size_name,"$",sep = "",collapse = "|"))]
    diff_feature = result$Feature[result$p.adj.val<pval_threshold]
    background_feature = result$Feature[result$p.adj.val>=pval_threshold]
    background_feature = c(background_feature,NA_feature)
    diff_up_supra_feature = result$Feature[result$p.adj.val<pval_threshold & effect_size>0]
    diff_down_supra_feature = result$Feature[result$p.adj.val<pval_threshold & effect_size<0]
    anaerobe_feature = Ground_truth$X1[Ground_truth$X2=='Anaerobic']
    aerobe_feature = Ground_truth$X1[Ground_truth$X2=='Aerobic']
    f_anaerobe_feature = Ground_truth$X1[Ground_truth$X2=='F Anaerobic']
    #diff_up_supra_feature_anaerobe <- length(intersect(diff_up_supra_feature,anaerobe_feature))
    diff_up_supra_feature_aerobe <- length(intersect(diff_up_supra_feature,aerobe_feature))
    #上面是a
    ###列联表信息
    #b
    diff_up_supra_feature_other <- length(setdiff(diff_up_supra_feature,aerobe_feature))
    #c
    non_diff_up_supra_feature_aerobe <- length(intersect(union(diff_down_supra_feature,background_feature),aerobe_feature))
    #d
    non_diff_up_supra_feature_other <- length(setdiff(union(diff_down_supra_feature,background_feature),aerobe_feature))
    #diff_up_supra_feature_facultative <- length(intersect(diff_up_supra_feature,f_anaerobe_feature))
    diff_down_supra_feature_anaerobe <- length(intersect(diff_down_supra_feature,anaerobe_feature))
    #上面是a
    ###列联表信息
    #b
    diff_down_supra_feature_other <- length(setdiff(diff_down_supra_feature,anaerobe_feature))
    #c
    non_diff_down_supra_feature_anaerobe <- length(intersect(union(diff_up_supra_feature,background_feature),anaerobe_feature))
    #d
    non_diff_down_supra_feature_other <- length(setdiff(union(diff_up_supra_feature,background_feature),anaerobe_feature))
    #diff_down_supra_feature_aerobe <- length(intersect(diff_down_supra_feature,aerobe_feature))
    #diff_down_supra_feature_facultative <- length(intersect(diff_down_supra_feature,f_anaerobe_feature))
    
  }else{
    if(sum(grepl(colnames(result) , pattern = 'Enrieched'))==1){
      Enrich = result$Enrieched
      NA_feature <- result$Feature[is.na(result$p.adj.val)]
      result <- result[!is.na(result$p.adj.val),]
      diff_feature = result$Feature[result$p.adj.val<pval_threshold]
      background_feature = result$Feature[result$p.adj.val>=pval_threshold]
      background_feature = c(background_feature,NA_feature)
      diff_up_supra_feature = result$Feature[result$p.adj.val<pval_threshold & Enrich=='Disease']
      diff_down_supra_feature = result$Feature[result$p.adj.val<pval_threshold & Enrich=='Ctrl']
      anaerobe_feature = Ground_truth$X1[Ground_truth$X2=='Anaerobic']
      aerobe_feature = Ground_truth$X1[Ground_truth$X2=='Aerobic']
      f_anaerobe_feature = Ground_truth$X1[Ground_truth$X2=='F Anaerobic']
      #diff_up_supra_feature_anaerobe <- length(intersect(diff_up_supra_feature,anaerobe_feature))
      diff_up_supra_feature_aerobe <- length(intersect(diff_up_supra_feature,aerobe_feature))
      #上面是a
      ###列联表信息
      #b
      diff_up_supra_feature_other <- length(setdiff(diff_up_supra_feature,aerobe_feature))
      #c
      non_diff_up_supra_feature_aerobe <- length(intersect(union(diff_down_supra_feature,background_feature),aerobe_feature))
      #d
      non_diff_up_supra_feature_other <- length(setdiff(union(diff_down_supra_feature,background_feature),aerobe_feature))
      #diff_up_supra_feature_facultative <- length(intersect(diff_up_supra_feature,f_anaerobe_feature))
      diff_down_supra_feature_anaerobe <- length(intersect(diff_down_supra_feature,anaerobe_feature))
      #上面是a
      ###列联表信息
      #b
      diff_down_supra_feature_other <- length(setdiff(diff_down_supra_feature,anaerobe_feature))
      #c
      non_diff_down_supra_feature_anaerobe <- length(intersect(union(diff_up_supra_feature,background_feature),anaerobe_feature))
      #d
      non_diff_down_supra_feature_other <- length(setdiff(union(diff_up_supra_feature,background_feature),anaerobe_feature))
      #diff_down_supra_feature_aerobe <- length(intersect(diff_down_supra_feature,aerobe_feature))
      #diff_down_supra_feature_facultative <- length(intersect(diff_down_supra_feature,f_anaerobe_feature))
    }else if(sum(grepl(colnames(result) , pattern = '^coef$'))==1){
      result <- result[which(result$metadata == 'Group'),]
      NA_feature <- result$Feature[is.na(result$p.adj.val)]
      result <- result[!is.na(result$p.adj.val),]
      diff_feature = result$Feature[result$p.adj.val<pval_threshold]
      background_feature = result$Feature[result$p.adj.val>=pval_threshold]
      background_feature = c(background_feature,NA_feature)
      diff_up_supra_feature = result$Feature[result$p.adj.val<pval_threshold & result$coef>0]
      diff_down_supra_feature = result$Feature[result$p.adj.val<pval_threshold & result$coef<0]
      anaerobe_feature = Ground_truth$X1[Ground_truth$X2=='Anaerobic']
      aerobe_feature = Ground_truth$X1[Ground_truth$X2=='Aerobic']
      f_anaerobe_feature = Ground_truth$X1[Ground_truth$X2=='F Anaerobic']
      #diff_up_supra_feature_anaerobe <- length(intersect(diff_up_supra_feature,anaerobe_feature))
      diff_up_supra_feature_aerobe <- length(intersect(diff_up_supra_feature,aerobe_feature))
      #上面是a
      ###列联表信息
      #b
      diff_up_supra_feature_other <- length(setdiff(diff_up_supra_feature,aerobe_feature))
      #c
      non_diff_up_supra_feature_aerobe <- length(intersect(union(diff_down_supra_feature,background_feature),aerobe_feature))
      #d
      non_diff_up_supra_feature_other <- length(setdiff(union(diff_down_supra_feature,background_feature),aerobe_feature))
      #diff_up_supra_feature_facultative <- length(intersect(diff_up_supra_feature,f_anaerobe_feature))
      diff_down_supra_feature_anaerobe <- length(intersect(diff_down_supra_feature,anaerobe_feature))
      #上面是a
      ###列联表信息
      #b
      diff_down_supra_feature_other <- length(setdiff(diff_down_supra_feature,anaerobe_feature))
      #c
      non_diff_down_supra_feature_anaerobe <- length(intersect(union(diff_up_supra_feature,background_feature),anaerobe_feature))
      #d
      non_diff_down_supra_feature_other <- length(setdiff(union(diff_up_supra_feature,background_feature),anaerobe_feature))
      #diff_down_supra_feature_aerobe <- length(intersect(diff_down_supra_feature,aerobe_feature))
      #diff_down_supra_feature_facultative <- length(intersect(diff_down_supra_feature,f_anaerobe_feature))
    }else{
      NA_feature <- result$Feature[is.na(result$p.adj.val)]
      result <- result[!is.na(result$p.adj.val),]
      diff_feature = result$Feature[result$p.adj.val<pval_threshold]
      background_feature = result$Feature[result$p.adj.val>=pval_threshold]
      background_feature = c(background_feature,NA_feature)
      anaerobe_feature = Ground_truth$X1[Ground_truth$X2=='Anaerobic']
      aerobe_feature = Ground_truth$X1[Ground_truth$X2=='Aerobic']
      f_anaerobe_feature = Ground_truth$X1[Ground_truth$X2=='F Anaerobic']
      #diff_up_supra_feature_anaerobe <- length(intersect(diff_up_supra_feature,anaerobe_feature))
      diff_up_supra_feature_aerobe <- NA
      #上面是a
      ###列联表信息
      #b
      diff_up_supra_feature_other <- NA
      #c
      non_diff_up_supra_feature_aerobe <- NA
      #d
      non_diff_up_supra_feature_other <- NA
      #diff_up_supra_feature_facultative <- length(intersect(diff_up_supra_feature,f_anaerobe_feature))
      diff_down_supra_feature_anaerobe <- NA
      #上面是a
      ###列联表信息
      #b
      diff_down_supra_feature_other <- NA
      #c
      non_diff_down_supra_feature_anaerobe <- NA
      #d
      non_diff_down_supra_feature_other <- NA
      #diff_down_supra_feature_aerobe <- length(intersect(diff_down_supra_feature,aerobe_feature))
      #diff_down_supra_feature_facultative <- length(intersect(diff_down_supra_feature,f_anaerobe_feature))
    }
    
  }
  enrichment_count <-
    data.frame(
      supra = c(
        diff_up_supra_feature_aerobe,
        diff_up_supra_feature_other,
        non_diff_up_supra_feature_aerobe,
        non_diff_up_supra_feature_other
      ),
      sub = c(
        diff_down_supra_feature_anaerobe,
        diff_down_supra_feature_other,
        non_diff_down_supra_feature_anaerobe,
        non_diff_down_supra_feature_other
      )
    )
  if(anyNA(enrichment_count)){
    return(c(NA,NA))
  }else{
    supra_feature_enrich_sig<-fisher.test(matrix(enrichment_count$supra,byrow = T,nrow = 2),alternative = 'greater')$p.value
    sub_feature_enrich_sig<-fisher.test(matrix(enrichment_count$sub,byrow = T,nrow = 2),alternative = 'greater')$p.value
    
    return(c(supra_feature_enrich_sig,
             sub_feature_enrich_sig))
  }
  
  
}
evaluation_enrichment_score<-function(result_file_path){
  result <- read.table(file = result_file_path,sep = '\t',header = T)
  effect_size_name <- unique(c('Group.Est','lfc_Group1','effect_size','log2FoldChange','logFC','log2FC','logFC','Effect_size','effect_size'))
  
  if(sum(grepl(colnames(result),pattern = paste0("^",effect_size_name,"$",collapse = "|")))==1){
    
    #NA_feature <- result$Feature[is.na(result$p.adj.val)]
    #result <- result[!is.na(result$p.adj.val),]
    effect_size <- result[,grep(colnames(result),pattern = paste("^",effect_size_name,"$",sep = "",collapse = "|"))]
    #diff_feature = result$Feature[result$p.adj.val<pval_threshold]
    #background_feature = result$Feature[result$p.adj.val>=pval_threshold]
    #background_feature = c(background_feature,NA_feature)
    #diff_up_supra_feature = result$Feature[result$p.adj.val<pval_threshold & effect_size>0]
    #diff_down_supra_feature = result$Feature[result$p.adj.val<pval_threshold & effect_size<0]
    anaerobe_feature = Ground_truth$X1[Ground_truth$X2=='Anaerobic']
    aerobe_feature = Ground_truth$X1[Ground_truth$X2=='Aerobic']
    f_anaerobe_feature = Ground_truth$X1[Ground_truth$X2=='F Anaerobic']
    df <-
      data.frame(
        feature = result$Feature,
        dir = rep(NA, length(result$Feature)),
        pval = rep(NA, length(result$Feature)),
        prior = rep(NA, length(result$Feature))
      )
    df$dir = ifelse(is.na(result$p.adj.val),0,ifelse(effect_size>0,1,-1))
    df$pval = ifelse(is.na(result$p.adj.val),1,result$p.adj.val)
    df$prior = ifelse(df$feature %in% aerobe_feature,'up',ifelse(df$feature %in% anaerobe_feature,'down',NA))
  }else{
    if(sum(grepl(colnames(result) , pattern = 'Enrieched'))==1){
      Enrich = result$Enrieched
      #NA_feature <- result$Feature[is.na(result$p.adj.val)]
      #result <- result[!is.na(result$p.adj.val),]
      #diff_feature = result$Feature[result$p.adj.val<pval_threshold]
      #background_feature = result$Feature[result$p.adj.val>=pval_threshold]
      #background_feature = c(background_feature,NA_feature)
      #diff_up_supra_feature = result$Feature[result$p.adj.val<pval_threshold & Enrich=='Disease']
      #diff_down_supra_feature = result$Feature[result$p.adj.val<pval_threshold & Enrich=='Ctrl']
      anaerobe_feature = Ground_truth$X1[Ground_truth$X2=='Anaerobic']
      aerobe_feature = Ground_truth$X1[Ground_truth$X2=='Aerobic']
      f_anaerobe_feature = Ground_truth$X1[Ground_truth$X2=='F Anaerobic']
      df <-
        data.frame(
          feature = result$Feature,
          dir = rep(NA, length(result$Feature)),
          pval = rep(NA, length(result$Feature)),
          prior = rep(NA, length(result$Feature))
        )
      df$dir = ifelse(is.na(result$p.adj.val),0,ifelse(Enrich == 'Disease',1,-1))
      df$pval = ifelse(is.na(result$p.adj.val),1,result$p.adj.val)
      df$prior = ifelse(df$feature %in% aerobe_feature,'up',ifelse(df$feature %in% anaerobe_feature,'down',NA))
    }else if(sum(grepl(colnames(result) , pattern = '^coef$'))==1){
      result <- result[which(result$metadata == 'Group'),]
      #NA_feature <- result$Feature[is.na(result$p.adj.val)]
      #result <- result[!is.na(result$p.adj.val),]
      #diff_feature = result$Feature[result$p.adj.val<pval_threshold]
      #background_feature = result$Feature[result$p.adj.val>=pval_threshold]
      #background_feature = c(background_feature,NA_feature)
      #diff_up_supra_feature = result$Feature[result$p.adj.val<pval_threshold & result$coef>0]
      #diff_down_supra_feature = result$Feature[result$p.adj.val<pval_threshold & result$coef<0]
      anaerobe_feature = Ground_truth$X1[Ground_truth$X2=='Anaerobic']
      aerobe_feature = Ground_truth$X1[Ground_truth$X2=='Aerobic']
      f_anaerobe_feature = Ground_truth$X1[Ground_truth$X2=='F Anaerobic']
      df <-
        data.frame(
          feature = result$Feature,
          dir = rep(NA, length(result$Feature)),
          pval = rep(NA, length(result$Feature)),
          prior = rep(NA, length(result$Feature))
        )
      df$dir = ifelse(is.na(result$p.adj.val),0,ifelse(result$coef>0,1,-1))
      df$pval = ifelse(is.na(result$p.adj.val),1,result$p.adj.val)
      df$prior = ifelse(df$feature %in% aerobe_feature,'up',ifelse(df$feature %in% anaerobe_feature,'down',NA))
    }else{
      # NA_feature <- result$Feature[is.na(result$p.adj.val)]
      # result <- result[!is.na(result$p.adj.val),]
      # diff_feature = result$Feature[result$p.adj.val<pval_threshold]
      # background_feature = result$Feature[result$p.adj.val>=pval_threshold]
      # background_feature = c(background_feature,NA_feature)
      # anaerobe_feature = Ground_truth$X1[Ground_truth$X2=='Anaerobic']
      # aerobe_feature = Ground_truth$X1[Ground_truth$X2=='Aerobic']
      # f_anaerobe_feature = Ground_truth$X1[Ground_truth$X2=='F Anaerobic']
      df <-
        data.frame(
          feature = result$Feature,
          dir = rep(NA, length(result$Feature)),
          pval = rep(NA, length(result$Feature)),
          prior = rep(NA, length(result$Feature))
        )
      
    }
    
  }
  
  return(df)
  
  
}
##############
current_dir <- getwd()

# 使用 list.files 函数查找文件，pattern 参数指定匹配包含 '_r2_' 的文件
files <- list.files(path = current_dir, pattern = "_r2_", full.names = TRUE)
######这里采用不同阈值对富集结果进行计算
enrichment_result_005 <- data.frame(
  Method = rep(NA,64),
  diff_up_supra_feature_aerobe = rep(NA, 64),
  diff_down_supra_feature_anaerobe = rep(NA, 64)
)
k=1
for(file in files){
  enrichment_res <- evaluation_enrichment2(file,0.05)
  modified_file <- gsub("_res_r2", "", basename(file))
  modified_file <- gsub(".tsv$", "", basename(modified_file))
  enrichment_res <- c(modified_file,enrichment_res)
  enrichment_result_005[k,] = enrichment_res
  k = k+1
} 
enrichment_result_005 <- enrichment_result_005 %>%
  mutate_at(vars(-Method),as.numeric)
##############
enrichment_result_001 <- data.frame(
  Method = rep(NA,64),
  diff_up_supra_feature_aerobe = rep(NA, 64),
  diff_down_supra_feature_anaerobe = rep(NA, 64)
)
k=1
for(file in files){
  enrichment_res <- evaluation_enrichment2(file,0.01)
  modified_file <- gsub("_res_r2", "", basename(file))
  modified_file <- gsub(".tsv$", "", basename(modified_file))
  enrichment_res <- c(modified_file,enrichment_res)
  enrichment_result_001[k,] = enrichment_res
  k = k+1
} 
enrichment_result_001 <- enrichment_result_001 %>%
  mutate_at(vars(-Method),as.numeric)
###############
enrichment_result_01 <- data.frame(
  Method = rep(NA,64),
  diff_up_supra_feature_aerobe = rep(NA, 64),
  diff_down_supra_feature_anaerobe = rep(NA, 64)
)
k=1
for(file in files){
  enrichment_res <- evaluation_enrichment2(file,0.1)
  modified_file <- gsub("_res_r2", "", basename(file))
  modified_file <- gsub(".tsv$", "", basename(modified_file))
  enrichment_res <- c(modified_file,enrichment_res)
  enrichment_result_01[k,] = enrichment_res
  k = k+1
} 
enrichment_result_01 <- enrichment_result_01 %>%
  mutate_at(vars(-Method),as.numeric)

write.table(enrichment_result_001,file='evaluation_result_enrichment_2sig_oral-GT-based_0.01.tsv',sep = '\t',row.names = F,col.names = T)
write.table(enrichment_result_005,file='evaluation_result_enrichment_2sig_oral-GT-based_0.05.tsv',sep = '\t',row.names = F,col.names = T)
write.table(enrichment_result_01,file='evaluation_result_enrichment_2sig_oral-GT-based_0.1.tsv',sep = '\t',row.names = F,col.names = T)
#########################
#################count计算
#########################
enrichment_result_2sig_0.05 <- na.omit(read.table(file = 'evaluation_result_enrichment_2sig_oral-GT-based_0.05.tsv',sep = '\t',header = T))
enrichment_result_2sig_0.01 <- na.omit(read.table(file = 'evaluation_result_enrichment_2sig_oral-GT-based_0.01.tsv',sep = '\t',header = T))
enrichment_result_2sig_0.1 <- na.omit(read.table(file = 'evaluation_result_enrichment_2sig_oral-GT-based_0.1.tsv',sep = '\t',header = T))

.calc_dir_counts <- function(df, dir_group = c("up", "down"), alpha = 0.05){
  dir_group <- match.arg(dir_group)
  
  # pval：转数值；NA 按 1（不显著）
  p <- suppressWarnings(as.numeric(df$pval))
  p[is.na(p)] <- 1
  sig <- p < alpha
  
  # dir：兼容数值/字符
  dir_raw <- df$dir
  if (is.factor(dir_raw)) dir_raw <- as.character(dir_raw)
  
  if (is.numeric(dir_raw)) {
    dir_num <- as.numeric(dir_raw)
  } else {
    dir_chr <- tolower(as.character(dir_raw))
    dir_num <- dplyr::case_when(
      dir_chr %in% c("up", "+1", "1") ~  1,
      dir_chr %in% c("down", "-1")    ~ -1,
      TRUE ~ suppressWarnings(as.numeric(dir_chr))
    )
  }
  
  # 上调组：dir >= 0（包含0）；下调组：dir < 0
  in_dir <- if (dir_group == "up") {
    !is.na(dir_num) & dir_num >= 0
  } else {
    !is.na(dir_num) & dir_num < 0
  }
  
  # prior：把 "NA" 字符也当成 NA
  prior <- df$prior
  if (is.factor(prior)) prior <- as.character(prior)
  prior <- dplyr::na_if(prior, "NA")
  
  # sig 要求匹配；nonsig 只要非NA
  need_label <- if (dir_group == "up") "up" else "down"
  sig_prior_flag    <- !is.na(prior) & prior == need_label
  nonsig_prior_flag <- !is.na(prior)
  
  sig_prior     <- sum(in_dir & sig  & sig_prior_flag, na.rm = TRUE)
  sig_nonprior  <- sum(in_dir & sig,  na.rm = TRUE) - sig_prior
  
  nsig_prior    <- sum(in_dir & !sig & nonsig_prior_flag, na.rm = TRUE)
  nsig_nonprior <- sum(in_dir & !sig, na.rm = TRUE) - nsig_prior
  
  c(sig_prior = sig_prior,
    sig_nonprior = sig_nonprior,
    nonsig_prior = nsig_prior,
    nonsig_nonprior = nsig_nonprior)
}


enrichment_count_df_0.05 <- imap_dfr(enrichment_list, function(df, method_name){
  
  up   <- .calc_dir_counts(df, dir_group = "up",   alpha = 0.05)
  down <- .calc_dir_counts(df, dir_group = "down", alpha = 0.05)
  
  tibble(
    Method = method_name,
    
    up_sig_prior        = up["sig_prior"],
    up_sig_nonprior     = up["sig_nonprior"],
    up_nonsig_prior     = up["nonsig_prior"],
    up_nonsig_nonprior  = up["nonsig_nonprior"],
    
    down_sig_prior        = down["sig_prior"],
    down_sig_nonprior     = down["sig_nonprior"],
    down_nonsig_prior     = down["nonsig_prior"],
    down_nonsig_nonprior  = down["nonsig_nonprior"]
  )
})
enrichment_count_df_0.01 <- imap_dfr(enrichment_list, function(df, method_name){
  
  up   <- .calc_dir_counts(df, dir_group = "up",   alpha = 0.01)
  down <- .calc_dir_counts(df, dir_group = "down", alpha = 0.01)
  
  tibble(
    Method = method_name,
    
    up_sig_prior        = up["sig_prior"],
    up_sig_nonprior     = up["sig_nonprior"],
    up_nonsig_prior     = up["nonsig_prior"],
    up_nonsig_nonprior  = up["nonsig_nonprior"],
    
    down_sig_prior        = down["sig_prior"],
    down_sig_nonprior     = down["sig_nonprior"],
    down_nonsig_prior     = down["nonsig_prior"],
    down_nonsig_nonprior  = down["nonsig_nonprior"]
  )
})
enrichment_count_df_0.1 <- imap_dfr(enrichment_list, function(df, method_name){
  
  up   <- .calc_dir_counts(df, dir_group = "up",   alpha = 0.1)
  down <- .calc_dir_counts(df, dir_group = "down", alpha = 0.1)
  
  tibble(
    Method = method_name,
    
    up_sig_prior        = up["sig_prior"],
    up_sig_nonprior     = up["sig_nonprior"],
    up_nonsig_prior     = up["nonsig_prior"],
    up_nonsig_nonprior  = up["nonsig_nonprior"],
    
    down_sig_prior        = down["sig_prior"],
    down_sig_nonprior     = down["sig_nonprior"],
    down_nonsig_prior     = down["nonsig_prior"],
    down_nonsig_nonprior  = down["nonsig_nonprior"]
  )
})
############################################################
###########################RRES分数计算######################
############################################################
####计算recall
##############################
library(dplyr)
library(purrr)
library(tibble)

calc_metrics_one <- function(df,
                             p_cutoff = 0.05,
                             prior_col = "prior",
                             pval_col  = "pval",
                             dir_col   = "dir",
                             up_label = "up",
                             down_label = "down") {
  pval  <- df[[pval_col]]
  dir   <- df[[dir_col]]
  prior <- df[[prior_col]]
  
  # 显著 & 方向预测
  sig <- !is.na(pval) & (pval < p_cutoff)
  pred_up   <- sig & !is.na(dir) & (dir > 0)
  pred_down <- sig & !is.na(dir) & (dir < 0)
  
  # 先验集合
  prior_up   <- !is.na(prior) & (prior == up_label)
  prior_down <- !is.na(prior) & (prior == down_label)
  
  # TP
  tp_up   <- sum(prior_up & pred_up)
  tp_down <- sum(prior_down & pred_down)
  
  # 分母
  n_prior_up   <- sum(prior_up)
  n_prior_down <- sum(prior_down)
  n_pred_up    <- sum(pred_up)
  n_pred_down  <- sum(pred_down)
  
  safe_div <- function(n, d) if (d == 0) 0 else n / d
  
  up_precision <- safe_div(tp_up, n_pred_up)
  up_recall    <- safe_div(tp_up, n_prior_up)
  up_F1 <- if ((up_precision + up_recall) == 0) 0 else
    2 * up_precision * up_recall / (up_precision + up_recall)
  
  down_precision <- safe_div(tp_down, n_pred_down)
  down_recall    <- safe_div(tp_down, n_prior_down)
  down_F1 <- if ((down_precision + down_recall) == 0) 0 else
    2 * down_precision * down_recall / (down_precision + down_recall)
  
  tibble(
    up_precision = up_precision,
    up_recall    = up_recall,
    up_F1        = up_F1,
    down_precision = down_precision,
    down_recall    = down_recall,
    down_F1        = down_F1
  )
}

# ===== 主输出：每个方法一行 =====
p_cutoff <- 0.05  # <- 阈值在这里改

metrics_df <- imap_dfr(enrichment_list, function(df, method_name) {
  out <- calc_metrics_one(df, p_cutoff = p_cutoff)
  tibble(Method = as.character(method_name)) %>% bind_cols(out)
})

metrics_df

#####排序分数
library(dplyr)

rank_without_precision <- function(perf_df, p_df,
                                   p_up_col   = "diff_up_supra_feature_aerobe",
                                   p_down_col = "diff_down_supra_feature_anaerobe",
                                   weight_mode = c("soft", "hard"),
                                   T = 5,         # soft 权重：-log10(q)/T 封顶到 1（常用 5 或 10）
                                   alpha = 0.05,  # hard 权重：q < alpha 则权重=1，否则=0
                                   combine = c("geom", "mean")) {
  
  weight_mode <- match.arg(weight_mode)
  combine <- match.arg(combine)
  
  df <- perf_df %>%
    inner_join(
      p_df %>% transmute(
        Method,
        p_up   = .data[[p_up_col]],
        p_down = .data[[p_down_col]]
      ),
      by = "Method"
    ) %>%
    mutate(
      # 1) 多重校正（在方法之间校正）
      q_up   = p.adjust(p_up,   method = "BH"),
      q_down = p.adjust(p_down, method = "BH"),
      
      # 2) 显著性权重（避免 p=0 导致 Inf）
      w_up = if (weight_mode == "soft") {
        pmin(1, -log10(pmax(q_up, 1e-300)) / T)
      } else {
        as.numeric(q_up < alpha)
      },
      w_down = if (weight_mode == "soft") {
        pmin(1, -log10(pmax(q_down, 1e-300)) / T)
      } else {
        as.numeric(q_down < alpha)
      },
      
      # 3) 不用 precision：用 recall 作为“命中能力”
      s_up   = up_recall   * w_up,
      s_down = down_recall * w_down,
      
      # 4) 合并 up/down：geom 更严格（两边都要好）；mean 更宽松
      score = if (combine == "geom") sqrt(s_up * s_down) else (s_up + s_down) / 2
    ) %>%
    arrange(dplyr::desc(score)) %>%
    mutate(rank = dplyr::dense_rank(dplyr::desc(score))) %>%
    select(Method, rank, score,
           up_recall, down_recall,
           p_up, p_down, q_up, q_down, w_up, w_down, s_up, s_down)
  
  df
}

rank_recall <- rank_without_precision(metrics_df, enrichment_result_2sig_0.05,
                                      weight_mode = "soft",
                                      T = 5,
                                      combine = "geom")
