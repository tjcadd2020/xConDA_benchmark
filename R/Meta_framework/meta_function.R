#####prevalence筛选函数
filter_by_prevalence <- function(matrix, prevalence) {
  # 计算每列的非零比例
  nonzero_proportion <- colSums(matrix != 0) / nrow(matrix)
  
  # 筛选出非零比例大于 prevalence 的列
  filtered_matrix <- matrix[, which(nonzero_proportion > prevalence), drop = FALSE]
  
  # 返回筛选后的矩阵
  return(filtered_matrix)
}
#####方法运行函数
#CLR标准化函数
library(compositions)
clr_transform<-function(data){
  #要求输入行是样本，列是特征,count
  data_clr <- data.frame(t(compositions::clr(t(data))))
  return(data_clr)
}
#CSS标准化函数
library(metagenomeSeq)
css_transform<-function(data,log_transform=TRUE){
  #要求输入行是样本，列是特征,count
  metaSeqObject = newMRexperiment(t(data))
  metaSeqObject_CSS <- cumNorm( metaSeqObject , p=cumNormStatFast(metaSeqObject) )
  data_CSS = data.frame(t(data.frame(MRcounts(metaSeqObject_CSS, norm=TRUE, log=log_transform))))
  #确定是否返回log转换的数据
  return(data_CSS)
}
#TMM标准化函数
TMM_transform <- function(features) {
  #输入列是特征，行是样本，count
  # Convert to Matrix from Data Frame
  features_norm = as.matrix(features)
  dd <- colnames(features_norm)
  
  # TMM Normalizing the Data
  X <- t(features_norm)
  
  libSize = edgeR::calcNormFactors(X, method = "TMM")
  eff.lib.size = colSums(X) * libSize
  
  ref.lib.size = mean(eff.lib.size)
  #Use the mean of the effective library sizes as a reference library size
  X.output = sweep(X, MARGIN = 2, eff.lib.size, "/") * ref.lib.size
  #Normalized read counts
  
  # Convert back to data frame
  features_TMM <- as.data.frame(t(X.output))
  
  # Rename the True Positive Features - Same Format as Before
  colnames(features_TMM) <- dd
  
  
  # Return as list
  return(features_TMM)
}
##log转换函数
LOG_transform <- function(x) {
  y <- replace(x, x == 0, min(x[x>0]) / 2)
  return(data.frame(log2(y)))
}
##输出指定的标准化数据函数
data_normalization <- function(simulated_data_all,normalization){
  if(normalization == 'CLR'){
    normalized_data<-clr_transform(simulated_data_all$simulated_data_abs)
  }
  if(normalization == 'TMM'){
    normalized_data<-TMM_transform(simulated_data_all$simulated_data_abs)
  }
  if(normalization == 'CSS'){
    normalized_data<-css_transform(simulated_data_all$simulated_data_abs)
  }
  if(normalization == 'LOG'){
    normalized_data<-LOG_transform(simulated_data_all$simulated_data_abs)
  }
  if(normalization == 'count'){
    normalized_data<-simulated_data_all$simulated_data_abs
  }
  if(normalization == 'TSS'){
    normalized_data<-simulated_data_all$simulated_data_rela
  }
  return(normalized_data)
}
#对协变量部分进行更改，因为目前筛选协变量已经包含在metadata中
fastANCOM_analysis<-function(simulated_data_all,covariates=NULL,categorical_variable_name = NULL,normalization){
  ###covariates以向量形式输入
  library(fastANCOM)
  simulated_data <- data_normalization(simulated_data_all,normalization = normalization)
  # batch <- simulated_data_all$metadata$Batch
  # batch_data <- model.matrix(~ batch -1)
  # simulated_data_all$metadata <- cbind(simulated_data_all$metadata,batch_data)
  # covariates <- c(setdiff(covariates,'Batch'),colnames(batch_data))
  if(setequal(categorical_variable_name,NULL)){}else{
    for(i in 1:length(categorical_variable_name)){
      simulated_data_all$metadata[[categorical_variable_name[i]]] <- as.numeric(as.factor(simulated_data_all$metadata[[categorical_variable_name[i]]]))
      if(length(unique(simulated_data_all$metadata[[categorical_variable_name[i]]])) > 2){
        another_name <- paste0(categorical_variable_name[i], '_temp')
        assign(another_name, simulated_data_all$metadata[[categorical_variable_name[i]]])
        
        # 使用 as.formula 来动态构造公式
        formula <- as.formula(paste("~", another_name, "-1"))
        temp_data <- model.matrix(formula, data = simulated_data_all$metadata)
        
        simulated_data_all$metadata <- cbind(simulated_data_all$metadata, temp_data)
        covariates <- c(setdiff(covariates, categorical_variable_name[i]), colnames(temp_data))
      }
      
    }
  }
  
  
  if(is.null(covariates)){
    fa_out<-fastANCOM(Y=as.matrix(simulated_data),x=simulated_data_all$metadata$Group,zero_cut = 1)  
  }else{
    fa_out<-fastANCOM(Y=as.matrix(simulated_data),x=simulated_data_all$metadata$Group,Z=as.matrix(simulated_data_all$metadata[,covariates]),zero_cut = 1)
  }
  fa_res<-fa_out$results$final
  fa_res$Feature <- rownames(fa_res)
  colnames(fa_res)[3] <- 'pval'
  colnames(fa_res)[4] <- 'p.adj.val'
  return(fa_res)
}
ANCOMBC_analysis<-function(simulated_data_all,covariates=NULL,normalization){
  
  library(ANCOMBC)
  library(tidyverse)
  library(TreeSummarizedExperiment)
  rownames <- rownames(simulated_data_all$simulated_data_abs)
  simulated_data <- data_normalization(simulated_data_all,normalization = normalization)
  rownames(simulated_data) <- rownames
  se <- TreeSummarizedExperiment(
    assays = list(counts = t(simulated_data)), 
    colData = simulated_data_all$metadata
  )
  
  # 捕获错误并处理
  result <- tryCatch({
    if (is.null(covariates)) {
      covariates <- append(covariates, 'Group')
      ancombc_out <- ancombc2(
        data = se, assay_name = "counts",
        fix_formula = covariates,
        p_adj_method = "BH", prv_cut = 0, lib_cut = 100,
        group = "Group", struc_zero = TRUE, alpha = 0.05, global = TRUE, verbose = TRUE
      )
    } else {
      covariates <- append(covariates, 'Group')
      formula_str <- paste(covariates, collapse = ' + ')
      ancombc_out <- ancombc2(
        data = se, assay_name = "counts",
        fix_formula = formula_str,
        p_adj_method = "BH", prv_cut = 0, lib_cut = 100,
        group = "Group", struc_zero = TRUE, alpha = 0.05, global = TRUE, verbose = TRUE
      )
    }
    ancombc_out
  }, error = function(e) {
    # 提取错误信息中的特征名
    error_message <- conditionMessage(e)
    print(error_message)
    features_to_remove <- str_extract_all(error_message, "Feature\\d+")[[1]]
    
    if (length(features_to_remove) > 0) {
      cat("Removing features with zero variance:", paste(features_to_remove, collapse = ", "), "\n")
      # 从原始数据中移除这些特征
      simulated_data_all$simulated_data_abs <- simulated_data_all$simulated_data_abs[, !colnames(simulated_data_all$simulated_data_abs) %in% features_to_remove]
      # 递归调用自身重新运行
      #print(simulated_data_all$simulated_data_abs)
      simulated_data <- data_normalization(simulated_data_all,normalization = normalization)
      se <- TreeSummarizedExperiment(
        assays = list(counts = t(simulated_data)), 
        colData = simulated_data_all$metadata
      )
      if (is.null(covariates)) {
        covariates <- append(covariates, 'Group')
        ancombc_out <- ancombc2(
          data = se, assay_name = "counts",
          fix_formula = covariates,
          p_adj_method = "BH", prv_cut = 0, lib_cut = 100,
          group = "Group", struc_zero = TRUE, alpha = 0.05, global = TRUE, verbose = TRUE
        )
      } else {
        covariates <- append(covariates, 'Group')
        formula_str <- paste(covariates, collapse = ' + ')
        ancombc_out <- ancombc2(
          data = se, assay_name = "counts",
          fix_formula = formula_str,
          p_adj_method = "BH", prv_cut = 0, lib_cut = 100,
          group = "Group", struc_zero = TRUE, alpha = 0.05, global = TRUE, verbose = TRUE
        )
      }
      ancombc_out
    } else {
      stop("Error occurred but no features identified to remove.")
    }
  })
  
  # 处理结果
  ancombc_res <- result$res
  ancombc_res_process <- ancombc_res
  colnames(ancombc_res_process)[1] <- 'Feature'
  #colnames(ancombc_res_process)[((length(covariates) + 1) * 4 + 2 + length(covariates))] <- 'p.adj.val'
  #colnames(ancombc_res_process)[((length(covariates) + 1) * 3 + 2 + length(covariates))] <- 'pval'
  ancombc_res_process <- ancombc_res_process %>%
    rename_with(~ gsub('^q_Group.*', 'p.adj.val', .), everything()) %>%
    rename_with(~ gsub('^p_Group.*', 'pval', .), everything())
  
  rownames(ancombc_res_process) <- ancombc_res_process$Feature
  na_feature <- setdiff(colnames(simulated_data_all$simulated_data_rela), ancombc_res_process$Feature)
  if (length(na_feature) > 0) {
    filter_feature <- as.data.frame(matrix(NA, nrow = length(na_feature), ncol = ncol(ancombc_res_process)))
    rownames(filter_feature) <- na_feature
    colnames(filter_feature) <- colnames(ancombc_res_process)
    ancombc_res_process <- rbind(ancombc_res_process, filter_feature)
  }
  ancombc_res_process$Feature <- rownames(ancombc_res_process)
  return(ancombc_res_process)
}
edgeR_analysis<-function(simulated_data_all,covariates = NULL,normalization){
  library(edgeR)
  simulated_data <- data_normalization(simulated_data_all,normalization = normalization)
  col_names <- covariates
  covariates <- c("Group",covariates)
  if(length(covariates)>1){
    covariates <- paste(covariates,collapse = ' + ')
  }
  formula_str <- paste0('~ ',covariates)
  Group  <- factor(simulated_data_all$metadata$Group)
  y <- DGEList(counts=t(simulated_data),group=Group)
  design <- model.matrix(as.formula(formula_str),data = simulated_data_all$metadata)
  y <- estimateDisp(y,design)
  
  fit <- glmQLFit(y,design)
  qlf <- glmQLFTest(fit,coef = 2)
  edgeR_res <- topTags(qlf,n = ncol(simulated_data))
  edgeR_res <- edgeR_res$table
  colnames(edgeR_res)[5] <- 'p.adj.val' 
  colnames(edgeR_res)[4] <- 'pval' 
  edgeR_res$Feature <- rownames(edgeR_res)
  return(edgeR_res)
}
limma_analysis<-function(simulated_data_all,covariates=NULL,normalization){
  simulated_data <- data_normalization(simulated_data_all,normalization = normalization)
  col_names <- covariates
  covariates <- c("Group",covariates)
  if(length(covariates)>1){
    covariates <- paste(covariates,collapse = ' + ')
  }
  formula_str <- paste0('~0+',covariates)
  #print(formula_str)
  library(limma)
  library(edgeR)
  simulated_data_all$metadata$Group <- as.factor(simulated_data_all$metadata$Group)
  mm<-model.matrix(as.formula(formula_str),data = simulated_data_all$metadata)
  colnames(mm)[1:2]<- c('group0','group1')
  fit<-lmFit(t(simulated_data),mm)
  contr <- makeContrasts(group1 - group0, levels = colnames(coef(fit)))
  tmp <- contrasts.fit(fit, contr)
  tmp <- eBayes(tmp)
  top.table <- topTable(tmp, sort.by = "P", n = Inf)
  limma_res <- top.table
  limma_res$Feature <- rownames(limma_res)
  colnames(limma_res)[5] <- 'p.adj.val'
  colnames(limma_res)[4] <- 'pval'
  return(limma_res)
}
LM_random_effect<-function(simulated_data_all,covariates,categorical_variable_name,normalization){
  library(lmerTest)
  library(MuMIn)
  simulated_data <- data_normalization(simulated_data_all,normalization = normalization)
  lmem_res <- data.frame(Feature = rep(NA,ncol(simulated_data_all$simulated_data_rela)),Pvalue=rep(NA,ncol(simulated_data_all$simulated_data_rela)),Effect_size = rep(NA,ncol(simulated_data_all$simulated_data_rela)),stringsAsFactors = F)
  for(k in 1:ncol(simulated_data_all$simulated_data_rela)) {
    # 创建临时数据
    lmem_res$Feature <- as.character(lmem_res$Feature)
    temp_data <- data.frame(cbind(simulated_data[,k], simulated_data_all$metadata))
    colnames(temp_data)[1] <- 'Feature'
    # 使用 paste 构建随机效应的字符串
    if(length(setdiff(covariates,categorical_variable_name))!=0 & !is.null(categorical_variable_name)){
      random_effects_str <- paste0("(1|",categorical_variable_name, ")")
      fixed_effects_str <- paste(setdiff(covariates,categorical_variable_name),collapse = " + ")
      # 构建完整的公式字符串
      formula_str <- paste("Feature ~ Group", paste(random_effects_str,collapse = " + "), fixed_effects_str,sep = " + ")
    }else if(!is.null(categorical_variable_name)){
      random_effects_str <- paste0("(1|",categorical_variable_name, ")")
      # 构建完整的公式字符串
      formula_str <- paste("Feature ~ Group", paste(random_effects_str,collapse = " + "),sep = " + ")
    }else{
      #random_effects_str <- paste0("(1|",categorical_variable_name, ")")
      fixed_effects_str <- paste(setdiff(covariates,categorical_variable_name),collapse = " + ")
      # 构建完整的公式字符串
      formula_str <- paste("Feature ~ Group", fixed_effects_str,sep = " + ")
    }
    #print(formula_str)
    # 使用 tryCatch 捕获错误
    tryCatch({
      # 运行 lmer 模型
      lm_re <- lmer(as.formula(formula_str), data = temp_data)
      #print(r.squaredGLMM(lm_re))
      # 提取 p-value 和估计效应
      # 获取以 "Group" 开头的系数
      group_coeffs <- grep("^Group", rownames(summary(lm_re)$coefficients), value = TRUE)
      
      if (length(group_coeffs) > 0) {
        pval <- summary(lm_re)$coefficients[group_coeffs, 'Pr(>|t|)']
        estimate_effect <- summary(lm_re)$coefficients[group_coeffs, 'Estimate']
        lmem_res[k, ] <- c(colnames(simulated_data)[k], pval, estimate_effect)
      }
      # pval <- summary(lm_re)$coefficients['Group', 'Pr(>|t|)']
      # #print(pval)
      # estimate_effect <- summary(lm_re)$coefficients['Group', 'Estimate']
      # #print(estimate_effect)
      # # 存储结果
      # lmem_res[k,] <- c(colnames(simulated_data)[k], pval, estimate_effect)
    }, error = function(e) {
      # 如果出现错误，则将结果设为 NA
      lmem_res[k,1] <- colnames(simulated_data)[k]
    }, finally = {
      # 确保赋值总是执行
      if (is.na(lmem_res[k,1])) {
        lmem_res[k,1] <- colnames(simulated_data)[k]
      }
    })
  }
  p.adj.val  = p.adjust(lmem_res$Pvalue,method = 'BH')
  lmem_res$p.adj.val = p.adj.val
  colnames(lmem_res)[2] <- 'pval'
  return(lmem_res)
}
LM_fixed_effect<-function(simulated_data_all,covariates,normalization){
  simulated_data <- data_normalization(simulated_data_all,normalization = normalization)
  covariates_str <- paste(covariates, collapse = " + ")
  # 构建完整的公式字符串
  formula_str <- paste("Feature ~ Group +", covariates_str)
  lmf_res <- data.frame(Feature = rep(NA,ncol(simulated_data_all$simulated_data_rela)),Pvalue=rep(NA,ncol(simulated_data_all$simulated_data_rela)),Effect_size = rep(NA,ncol(simulated_data_all$simulated_data_rela)))
  for(k in 1:ncol(simulated_data_all$simulated_data_rela)){
    temp_data <- data.frame(cbind(simulated_data[,k],simulated_data_all$metadata))
    colnames(temp_data)[1]<-'Feature'
    lmf <- lm(as.formula(formula_str),data = temp_data)
    # 获取以 "Group" 开头的系数
    group_coeffs <- grep("^Group", rownames(summary(lmf)$coefficients), value = TRUE)
    
    if (length(group_coeffs) > 0) {
      pval <- summary(lmf)$coefficients[group_coeffs, 'Pr(>|t|)']
      estimate_effect <- summary(lmf)$coefficients[group_coeffs, 'Estimate']
      lmf_res[k, ] <- c(colnames(simulated_data)[k], pval, estimate_effect)
    }
  }
  p.adj.val  = p.adjust(lmf_res$Pvalue,method = 'BH')
  lmf_res$p.adj.val = p.adj.val
  colnames(lmf_res)[2] <- 'pval'
  return(lmf_res)
}
Maaslin2_LM_analysis<-function(simulated_data_all,covariates=NULL,categorical_variable_name = NULL,normalization,transform){
  covariates <- append(covariates,"Group")
  reference = NULL
  if(setequal(categorical_variable_name,NULL)){}else{
    for(i in 1:length(categorical_variable_name)){
      simulated_data_all$metadata[[categorical_variable_name[i]]] <- as.numeric(as.factor(simulated_data_all$metadata[[categorical_variable_name[i]]]))
      if(length(unique(simulated_data_all$metadata[[categorical_variable_name[i]]])) > 2){
        reference <- c(reference,paste(categorical_variable_name[i],unique(simulated_data_all$metadata[[categorical_variable_name[i]]])[1],sep=','))
      }
      
    }
  }
  library(Maaslin2)
  meta <- simulated_data_all$metadata
  data <- simulated_data_all$simulated_data_abs
  filename <- "data_maaslin2_lm"
  rownames(meta)<-meta$Sample
  fit_data = Maaslin2(
    input_data = data, 
    input_metadata = meta, 
    min_prevalence = 0,
    normalization = normalization,analysis_method = "LM",
    output = filename,transform = transform,plot_heatmap = FALSE,plot_scatter = FALSE,
    fixed_effects = covariates,
    random_effects = NULL)
  maa_lm_res<-fit_data$results
  colnames(maa_lm_res)[1] <- 'Feature'
  colnames(maa_lm_res)[8] <- 'p.adj.val'
  return(maa_lm_res)
}
Maaslin2_NEGBIN_analysis<-function(simulated_data_all,covariates=NULL,categorical_variable_name = NULL,normalization,transform){
  covariates <- append(covariates,"Group")
  random_effect_variable <- NULL
  if('subject_id' %in% categorical_variable_name){
    random_effect_variable = 'subject_id'
  }
  reference = NULL
  if(setequal(categorical_variable_name,NULL)){}else{
    for(i in 1:length(categorical_variable_name)){
      simulated_data_all$metadata[[categorical_variable_name[i]]] <- as.numeric(as.factor(simulated_data_all$metadata[[categorical_variable_name[i]]]))
      if(length(unique(simulated_data_all$metadata[[categorical_variable_name[i]]])) > 2){
        reference <- c(reference,paste(categorical_variable_name[i],unique(simulated_data_all$metadata[[categorical_variable_name[i]]])[1],sep=','))
      }
      
    }
  }
  library(Maaslin2)
  meta <- simulated_data_all$metadata
  data <- simulated_data_all$simulated_data_abs
  filename <- "data_maaslin2_negbin"
  rownames(meta)<-meta$Sample
  fit_data = Maaslin2(
    input_data = data, 
    input_metadata = meta, 
    min_prevalence = 0,reference = reference,
    normalization = normalization,analysis_method = "NEGBIN",
    output = filename,transform = transform,plot_heatmap = FALSE,plot_scatter = FALSE,
    fixed_effects = covariates,
    random_effects = random_effect_variable)
  maa_lm_res<-fit_data$results
  colnames(maa_lm_res)[1] <- 'Feature'
  colnames(maa_lm_res)[8] <- 'p.adj.val'
  return(maa_lm_res)
}
Maaslin2_ZINB_analysis<-function(simulated_data_all,covariates=NULL,categorical_variable_name = NULL,normalization,transform){
  covariates <- append(covariates,"Group")
  random_effect_variable <- NULL
  if('subject_id' %in% categorical_variable_name){
    random_effect_variable = 'subject_id'
  }
  reference = NULL
  if(setequal(categorical_variable_name,NULL)){}else{
    for(i in 1:length(categorical_variable_name)){
      simulated_data_all$metadata[[categorical_variable_name[i]]] <- as.numeric(as.factor(simulated_data_all$metadata[[categorical_variable_name[i]]]))
      if(length(unique(simulated_data_all$metadata[[categorical_variable_name[i]]])) > 2){
        reference <- c(reference,paste(categorical_variable_name[i],unique(simulated_data_all$metadata[[categorical_variable_name[i]]])[1],sep=','))
      }
      
    }
  }
  library(Maaslin2)
  meta <- simulated_data_all$metadata
  data <- simulated_data_all$simulated_data_abs
  filename <- "data_maaslin2_zinb"
  rownames(meta)<-meta$Sample
  fit_data = Maaslin2(
    input_data = data, 
    input_metadata = meta, 
    min_prevalence = 0,reference = reference,
    normalization = normalization,analysis_method = "ZINB",
    output = filename,transform = transform,plot_heatmap = FALSE,plot_scatter = FALSE,
    fixed_effects = covariates,
    random_effects = random_effect_variable)
  maa_lm_res<-fit_data$results
  colnames(maa_lm_res)[1] <- 'Feature'
  colnames(maa_lm_res)[8] <- 'p.adj.val'
  return(maa_lm_res)
}
DESeq2_analysis<-function(simulated_data_all,covariates=NULL,normalization){
  
  simulated_data <- data_normalization(simulated_data_all,normalization = normalization)
  col_names <- covariates
  covariates <- c("Group",covariates)
  if(length(covariates)>1){
    covariates <- paste(covariates,collapse = ' + ')
  }
  formula_str <- paste0('~ ',covariates)
  library(DESeq2)
  library(dplyr)
  simulated_data_all$metadata$Group <- as.factor(simulated_data_all$metadata$Group)
  #simulated_data_all$metadata <- simulated_data_all$metadata %>%
  #mutate(across(contains("Ca"), as.factor))
  dds <- DESeqDataSetFromMatrix(countData=t(simulated_data) + 1, 
                                colData=simulated_data_all$metadata, 
                                design=as.formula(formula_str), tidy = FALSE)
  
  dds <- DESeq(dds)
  
  DESeq2_res <- data.frame(results(dds,contrast = c('Group',1,0)))
  DESeq2_res$Feature <- rownames(DESeq2_res)
  colnames(DESeq2_res)[which(colnames(DESeq2_res)=='padj')] <- 'p.adj.val'
  #head(results(dds,contrast = c('Group',1,0)))
  return(DESeq2_res)
}
Maaslin2_CPLM_analysis<-function(simulated_data_all,covariates=NULL,categorical_variable_name = NULL,normalization,transform){
  covariates <- append(covariates,"Group")
  random_effect_variable <- NULL
  if('subject_id' %in% categorical_variable_name){
    random_effect_variable = 'subject_id'
  }
  reference = NULL
  if(setequal(categorical_variable_name,NULL)){}else{
    for(i in 1:length(categorical_variable_name)){
      simulated_data_all$metadata[[categorical_variable_name[i]]] <- as.numeric(as.factor(simulated_data_all$metadata[[categorical_variable_name[i]]]))
      if(length(unique(simulated_data_all$metadata[[categorical_variable_name[i]]])) > 2){
        reference <- c(reference,paste(categorical_variable_name[i],unique(simulated_data_all$metadata[[categorical_variable_name[i]]])[1],sep=','))
      }
      
    }
  }
  library(Maaslin2)
  meta <- simulated_data_all$metadata
  data <- simulated_data_all$simulated_data_abs
  filename <- "data_maaslin2_cplm"
  rownames(meta)<-meta$Sample
  fit_data = Maaslin2(
    input_data = data, 
    input_metadata = meta, 
    min_prevalence = 0,reference = reference,
    normalization = normalization,analysis_method = "CPLM",
    output = filename,transform = transform,plot_heatmap = FALSE,plot_scatter = FALSE,
    fixed_effects = covariates,
    random_effects = random_effect_variable)
  maa_lm_res<-fit_data$results
  colnames(maa_lm_res)[1] <- 'Feature'
  colnames(maa_lm_res)[8] <- 'p.adj.val'
  return(maa_lm_res)
}
###meta函数
library(metap)
meta_combine_p <- function(result_list, sample_size,method='Stouffer') {
  # 获取所有唯一的 Feature
  all_features <- unique(unlist(lapply(result_list, function(df) df$Feature)))
  
  # 按照 all_features 更新每个 dataframe
  result_list <- lapply(result_list, function(df) {
    # 创建一个包含所有 Feature 的新 dataframe
    updated_df <- data.frame(Feature = all_features)
    # 合并原始数据，缺失部分填充 NA
    updated_df <- merge(updated_df, df, by = "Feature", all.x = TRUE)
    # 返回更新后的 dataframe
    return(updated_df)
  })
  #print(result_list)
  # 计算 cohort 数量
  cohort_num <- length(result_list)
  
  # 初始化合并的结果 dataframe
  combined_p_data <- data.frame(
    Feature = all_features,
    pval = rep(NA, length(all_features)),
    p.adj.val = rep(NA, length(all_features))
  )
  
  # 计算每个样本的权重
  sample_size_weights <- sqrt(sample_size)
  
  # 计算合并的 p 值
  for (i in 1:length(all_features)) {
    # 提取当前 Feature 的所有 p 值
    p <- unlist(lapply(result_list, function(df) {
      df$pval[which(df$Feature == all_features[i])]
    }))
    
    # 删除 NA 值和对应的权重
    valid_indices <- !is.na(p)
    p <- p[valid_indices]
    weights <- sample_size_weights[valid_indices]
    # 检查有效 P 值数量
    valid_p_count <- sum(valid_indices)
    
    # 如果有效 P 值少于 2，跳过计算
    if (valid_p_count < 2) {
      combined_p_data$pval[i] <- NA
    } else {
      # 处理 P 值不全为 NA 的情况
      if(method == 'Fisher'){
        combined_p_val <- sumlog(p = as.numeric(p), log.p = FALSE,log.input = FALSE)$p
        combined_p_data$pval[i] <- combined_p_val
      }else if(method == 'Stouffer'){
        combined_p_val <- sumz(p = p, weights = weights)$p
        combined_p_data$pval[i] <- combined_p_val
      }
    }
  }
  
  # 计算调整后的 p 值
  combined_p_data$p.adj.val <- p.adjust(combined_p_data$pval, method = 'BH')
  
  return(combined_p_data)
}

meta_combine<-function(effect_merge_data,se_merge_data,cohort_num,
                       rma_conv = 1e-10,method_type = method_type,
                       rma_maxit = 1000){
  ind_feature <- !is.na(effect_merge_data) & !is.na(se_merge_data)
  effect_merge_data <- effect_merge_data[rowSums(ind_feature)>1,]
  se_merge_data <- se_merge_data[rowSums(ind_feature)>1,]
  ind_feature<-ind_feature[rowSums(ind_feature)>1,]
  #ate_merge_data<-ate_merge_data[rowSums(ind_feature)==cohort_num,]
  #ate_se_merge_data<-ate_se_merge_data[rowSums(ind_feature)==cohort_num,]
  #ind_feature<-ind_feature[rowSums(ind_feature)==cohort_num,]
  rownames(ind_feature)<-rownames(effect_merge_data[rowSums(ind_feature)>1,])
  result <- data.frame(matrix(NA,
                              nrow = nrow(effect_merge_data),
                              ncol = 9 + cohort_num))
  colnames(result) <- c("coef",
                        "stderr",
                        "pval",
                        "k",
                        "tau2",
                        "stderr.tau2",
                        "pval.tau2",
                        "I2",
                        "H2",
                        paste0("weight_",
                               colnames(effect_merge_data)))
  rownames(result)<- rownames(effect_merge_data)
  for(feature in rownames(effect_merge_data)){
    #print(feature)
    yi_now  <- unlist(c(effect_merge_data[feature, ind_feature[feature, ]]))
    sei_now <- unlist(c(se_merge_data[feature, ind_feature[feature, ]]))
    
    all_names <- colnames(effect_merge_data)[ind_feature[feature, ]]
    
    valid_idx <- which(
      is.finite(yi_now) & !is.na(yi_now) &
        is.finite(sei_now) & !is.na(sei_now) &
        sei_now > 0
    )
    
    rma_fit <- NULL
    
    if (length(valid_idx) >= 2) {
      yi_use   <- yi_now[valid_idx]
      sei_use  <- sei_now[valid_idx]
      name_use <- all_names[valid_idx]
      
      rma_fit <- tryCatch(
        metafor::rma.uni(
          yi      = yi_use,
          sei     = sei_use,
          slab    = name_use,
          method  = method_type,
          verbose = FALSE,
          control = list(
            threshold = rma_conv,
            stepadj   = 0.5,
            tau2.max  = 1e6,
            maxiter   = rma_maxit
          )
        ),
        error = function(e) NULL
      )
    }
    
    if (!is.null(rma_fit)) {
      wts <- metafor::weights.rma.uni(rma_fit)
      
      ## 优先使用模型对象自身保存的标签
      if (!is.null(rma_fit$slab) && length(rma_fit$slab) == length(wts)) {
        names(wts) <- rma_fit$slab
      }
      
      result[feature, c("coef",
                        "stderr",
                        "pval",
                        "k",
                        "tau2",
                        "stderr.tau2",
                        "pval.tau2",
                        "I2",
                        "H2")] <- unlist(rma_fit[c("beta",
                                                   "se",
                                                   "pval",
                                                   "k",
                                                   "tau2",
                                                   "se.tau2",
                                                   "QEp",
                                                   "I2",
                                                   "H2")])
      
      result[feature, paste0("weight_", names(wts))] <- wts
    }
  }
  result$Bonferroni_pval<-p.adjust(result$pval, method = "bonf")
  result$BH_pval<-p.adjust(result$pval, method = "BH")
  result<-cbind(rownames(result),result)
  colnames(result)[1]<-"Feature"
  return(result)
}
safe_intersect_or_null <- function(target, candidates) {
  if (length(candidates) == 0) return(NULL)
  res <- intersect(target, candidates)
  if (length(res) == 0) return(NULL)
  return(res)
}
remove_low_feature_samples <- function(mat, verbose = TRUE) {
  # 每一行中非零元素的个数
  non_zero_count <- rowSums(mat != 0)
  
  # 保留非零特征数 ≥ 2 的样本
  keep_idx <- non_zero_count >= 2
  
  if (verbose) {
    removed <- sum(!keep_idx)
    message(removed, " samples removed with fewer than 2 non-zero features.")
  }
  
  mat[keep_idx, , drop = FALSE]
}
filter_features_in_multiple_batches <- function(data, min_batch_count = 2, verbose = TRUE) {
  meta <- data$metadata
  abs_mat <- data$simulated_data_abs
  rela_mat <- data$simulated_data_rela
  
  # 确保 SampleID 对齐
  stopifnot(rownames(abs_mat) == meta$SampleID)
  
  # 获取所有唯一的 Batch
  batches <- unique(meta$Batch)
  
  # 构建一个向量，记录每个特征在多少个 batch 中有非零值
  feature_batch_count <- setNames(rep(0, ncol(abs_mat)), colnames(abs_mat))
  
  for (batch in batches) {
    samples <- meta$SampleID[meta$Batch == batch]
    if (length(samples) == 0) next
    
    batch_data <- abs_mat[samples, , drop = FALSE]
    feature_present <- colSums(batch_data != 0) > 0  # 特征在该 batch 中是否出现
    feature_batch_count[feature_present] <- feature_batch_count[feature_present] + 1
  }
  
  # 筛选在至少 min_batch_count 个 batch 中出现的特征
  keep_features <- names(feature_batch_count[feature_batch_count >= min_batch_count])
  
  if (verbose) {
    removed <- length(feature_batch_count) - length(keep_features)
    message(removed, " features removed; ", length(keep_features), " retained.")
  }
  
  # 保留筛选后的特征列
  abs_mat_filtered <- abs_mat[, keep_features, drop = FALSE]
  rela_mat_filtered <- rela_mat[, keep_features, drop = FALSE]
  
  # 返回新的 list
  return(list(
    metadata = meta,
    simulated_data_abs = abs_mat_filtered,
    simulated_data_rela = rela_mat_filtered
  ))
}




run_meta <- function(method_name, normalization, meta_method, data, transform = NULL, data_index) {
  ##### helper functions #####
  add_unique <- function(x, val) {
    unique(as.character(c(x, val)))
  }
  
  drop_allzero_samples <- function(mat, min_nonzero_features = 2) {
    mat <- as.matrix(mat)
    
    # 基础检查
    if (any(!is.finite(mat))) stop("Input matrix contains NA/Inf/NaN.")
    if (any(mat < 0)) stop("Input matrix contains negative values.")
    
    # 每个样本（行）的“非零特征数”
    nnz_feat <- rowSums(mat > 0)
    
    # 删掉非零特征数 < min_nonzero_features 的样本
    keep_samp <- nnz_feat >= min_nonzero_features
    mat <- mat[keep_samp, , drop = FALSE]
    
    if (nrow(mat) == 0) {
      stop(sprintf(
        "All samples were removed (each had < %d non-zero features after filtering).",
        min_nonzero_features
      ))
    }
    
    mat
  }
  
  prepare_design_info <- function(meta) {
    meta$SampleID <- as.character(meta$SampleID)
    if ("Batch" %in% colnames(meta)) meta$Batch <- as.character(meta$Batch)
    
    # 清理分类变量 unused levels
    for (nm in setdiff(colnames(meta), "SampleID")) {
      x <- meta[[nm]]
      if (is.factor(x) || is.character(x) || is.logical(x)) {
        meta[[nm]] <- droplevels(factor(x))
      }
    }
    
    # Group 必须至少2水平
    if (!("Group" %in% colnames(meta))) return(NULL)
    if (!is.factor(meta$Group)) meta$Group <- droplevels(factor(meta$Group))
    if (nlevels(meta$Group) < 2) return(NULL)
    
    covariates_now <- setdiff(colnames(meta), c("Group", "SampleID", "Batch"))
    
    if (length(covariates_now) > 0) {
      keep_covariates <- covariates_now[
        vapply(covariates_now, function(v) {
          x <- meta[[v]]
          if (is.factor(x)) {
            nlevels(x) >= 2
          } else {
            length(unique(x[!is.na(x)])) >= 2
          }
        }, logical(1))
      ]
    } else {
      keep_covariates <- character(0)
    }
    
    keep_cols <- c("SampleID", "Group", intersect("Batch", colnames(meta)), keep_covariates)
    meta2 <- meta[, keep_cols, drop = FALSE]
    
    categorical_variable_name_now <- NULL
    if (length(keep_covariates) > 0) {
      categorical_variable_name_now <- keep_covariates[
        vapply(keep_covariates, function(v) {
          x <- meta2[[v]]
          is.factor(x) || is.character(x) || is.logical(x)
        }, logical(1))
      ]
      if (length(categorical_variable_name_now) == 0) categorical_variable_name_now <- NULL
    }
    
    if (length(keep_covariates) == 0) keep_covariates <- NULL
    
    list(
      metadata = meta2,
      covariates = keep_covariates,
      categorical_variable_name = categorical_variable_name_now
    )
  }
  
  # ---------- full-rank helpers ----------
  rank_ok <- function(meta, terms) {
    terms <- unique(as.character(terms))
    terms <- terms[terms %in% colnames(meta)]
    
    if (length(terms) == 0) return(TRUE)
    
    mm <- model.matrix(reformulate(terms), data = meta)
    qr(mm)$rank == ncol(mm)
  }
  
  drop_covariates_to_full_rank_minimal <- function(meta,
                                                   covariates,
                                                   protect = character(0),
                                                   always_keep = "Group",
                                                   max_drop = 3,
                                                   allow_drop_protect = FALSE,
                                                   verbose = TRUE,
                                                   label = NULL) {
    covariates <- unique(as.character(covariates))
    covariates <- covariates[covariates %in% colnames(meta)]
    
    protect <- intersect(protect, colnames(meta))
    always_keep <- intersect(always_keep, colnames(meta))
    
    base_terms <- unique(c(always_keep, protect, setdiff(covariates, c(always_keep, protect))))
    base_terms <- base_terms[base_terms %in% colnames(meta)]
    
    if (rank_ok(meta, base_terms)) {
      if (verbose) {
        msg_prefix <- if (!is.null(label)) paste0("[", label, "] ") else ""
        message(msg_prefix, "Design matrix is already full rank with terms: ",
                paste(base_terms, collapse = ", "))
      }
      return(list(
        covariates = setdiff(base_terms, always_keep),
        dropped = character(0),
        full_rank = TRUE,
        terms_used = base_terms
      ))
    }
    
    removable <- setdiff(base_terms, c(always_keep, protect))
    
    max_k <- min(max_drop, length(removable))
    for (k in seq_len(max_k)) {
      combs <- combn(removable, k, simplify = FALSE)
      for (drop_set in combs) {
        kept_terms <- setdiff(base_terms, drop_set)
        if (rank_ok(meta, kept_terms)) {
          if (verbose) {
            msg_prefix <- if (!is.null(label)) paste0("[", label, "] ") else ""
            message(msg_prefix, "Dropped to restore full rank: ",
                    paste(drop_set, collapse = ", "))
            message(msg_prefix, "Terms kept for modeling: ",
                    paste(kept_terms, collapse = ", "))
          }
          return(list(
            covariates = setdiff(kept_terms, always_keep),
            dropped = drop_set,
            full_rank = TRUE,
            terms_used = kept_terms
          ))
        }
      }
    }
    
    if (allow_drop_protect && length(protect) > 0) {
      kept_terms <- setdiff(base_terms, protect)
      if (rank_ok(meta, kept_terms)) {
        if (verbose) {
          msg_prefix <- if (!is.null(label)) paste0("[", label, "] ") else ""
          message(msg_prefix, "Dropped protected terms to restore full rank: ",
                  paste(protect, collapse = ", "))
          message(msg_prefix, "Terms kept for modeling: ",
                  paste(kept_terms, collapse = ", "))
        }
        return(list(
          covariates = setdiff(kept_terms, always_keep),
          dropped = protect,
          full_rank = TRUE,
          terms_used = kept_terms
        ))
      }
    }
    
    stop(
      if (!is.null(label)) paste0("[", label, "] ") else "",
      "Design matrix not full rank. Minimal-drop search up to max_drop=",
      max_drop,
      " failed."
    )
  }
  
  # ---------- batch 内部检查 + 修复 ----------
  check_and_fix_batch_design <- function(batch_name,
                                         batch_meta,
                                         covariates,
                                         categorical_variable_name = NULL,
                                         max_drop = 2,
                                         verbose = TRUE) {
    msg_prefix <- paste0("[Batch ", batch_name, "] ")
    
    if (verbose) {
      message(msg_prefix, "Start design check")
      message(msg_prefix, "n samples = ", nrow(batch_meta))
    }
    
    # Group 检查
    if (!("Group" %in% colnames(batch_meta))) {
      if (verbose) message(msg_prefix, "Skip: no Group column.")
      return(NULL)
    }
    if (!is.factor(batch_meta$Group)) {
      batch_meta$Group <- droplevels(factor(batch_meta$Group))
    }
    if (nlevels(batch_meta$Group) < 2) {
      if (verbose) {
        message(msg_prefix, "Skip: Group has < 2 levels after filtering.")
      }
      return(NULL)
    }
    if (verbose) {
      message(msg_prefix, "Group levels = ",
              paste(levels(batch_meta$Group), collapse = ", "))
    }
    
    covariates <- unique(as.character(covariates))
    covariates <- covariates[covariates %in% colnames(batch_meta)]
    
    # 报告原始协变量
    if (verbose) {
      if (length(covariates) == 0) {
        message(msg_prefix, "No covariates to check.")
      } else {
        message(msg_prefix, "Candidate covariates = ",
                paste(covariates, collapse = ", "))
      }
    }
    
    # 再次检查单水平/常数（虽然 prepare_design_info 已做过，但这里显式报告）
    dropped_singleton <- character(0)
    keep_covs <- character(0)
    
    if (length(covariates) > 0) {
      for (v in covariates) {
        x <- batch_meta[[v]]
        ok <- if (is.factor(x) || is.character(x) || is.logical(x)) {
          nlevels(droplevels(factor(x))) >= 2
        } else {
          length(unique(x[!is.na(x)])) >= 2
        }
        
        if (ok) {
          keep_covs <- c(keep_covs, v)
        } else {
          dropped_singleton <- c(dropped_singleton, v)
        }
      }
    }
    
    if (verbose) {
      if (length(dropped_singleton) > 0) {
        message(msg_prefix, "Dropped singleton/constant covariates: ",
                paste(dropped_singleton, collapse = ", "))
      } else {
        message(msg_prefix, "No singleton/constant covariates were dropped.")
      }
    }
    
    # 满秩检查（batch 内不带 Batch）
    fix <- drop_covariates_to_full_rank_minimal(
      meta = batch_meta,
      covariates = keep_covs,
      protect = character(0),
      always_keep = "Group",
      max_drop = max_drop,
      allow_drop_protect = FALSE,
      verbose = verbose,
      label = paste0("Batch ", batch_name)
    )
    
    covariates_use2 <- fix$covariates
    categorical_use2 <- intersect(as.character(categorical_variable_name), as.character(covariates_use2))
    if (length(categorical_use2) == 0) categorical_use2 <- NULL
    
    if (verbose) {
      if (length(covariates_use2) > 0) {
        message(msg_prefix, "Final covariates used = ",
                paste(covariates_use2, collapse = ", "))
      } else {
        message(msg_prefix, "Final covariates used = NULL")
      }
      
      if (!is.null(categorical_use2)) {
        message(msg_prefix, "Final categorical covariates = ",
                paste(categorical_use2, collapse = ", "))
      } else {
        message(msg_prefix, "Final categorical covariates = NULL")
      }
      
      message(msg_prefix, "Design check finished")
    }
    
    list(
      metadata = batch_meta,
      covariates = covariates_use2,
      categorical_variable_name = categorical_use2,
      dropped_singleton = dropped_singleton,
      dropped_rank = fix$dropped,
      terms_used = fix$terms_used
    )
  }
  
  run_batch_method_simple <- function(batch_data_list, selected_method_func, normalization) {
    lapply(batch_data_list, function(obj) {
      selected_method_func(obj$data, obj$covariates, normalization)
    })
  }
  
  run_batch_method_maaslin <- function(batch_data_list, selected_method_func, normalization, transform) {
    lapply(batch_data_list, function(obj) {
      selected_method_func(
        obj$data,
        obj$covariates,
        obj$categorical_variable_name,
        normalization,
        transform
      )
    })
  }
  
  ##### method map #####
  output_file_name <- paste(method_name, normalization, sep = "_")
  
  method_function_map <- list(
    edgeR = edgeR_analysis,
    lmem = LM_random_effect,
    limma = limma_analysis,
    lfem = LM_fixed_effect,
    ANCOMBC = ANCOMBC_analysis,
    Maaslin2_lm = Maaslin2_LM_analysis,
    fastANCOM = fastANCOM_analysis,
    Maaslin2_zinb = Maaslin2_ZINB_analysis,
    Maaslin2_negbin = Maaslin2_NEGBIN_analysis,
    Maaslin2_cplm = Maaslin2_CPLM_analysis,
    DESeq2 = DESeq2_analysis
  )
  
  selected_method_func <- method_function_map[[method_name]]
  if (is.null(selected_method_func)) stop("Unknown method_name: ", method_name)
  
  ##### Step 1: 先按 batch 过滤 feature #####
  data <- filter_features_in_multiple_batches(data, min_batch_count = 2, verbose = TRUE)
  data$metadata$SampleID <- as.character(data$metadata$SampleID)
  if ("Batch" %in% colnames(data$metadata)) data$metadata$Batch <- as.character(data$metadata$Batch)
  
  ##### Step 2: pooled data 先过滤并对齐 #####
  data$simulated_data_abs <- remove_low_feature_samples(data$simulated_data_abs)
  data$simulated_data_abs <- filter_by_prevalence(data$simulated_data_abs, 0.1)
  data$simulated_data_abs <- drop_allzero_samples(data$simulated_data_abs)
  
  data$simulated_data_rela <- remove_low_feature_samples(data$simulated_data_rela)
  data$simulated_data_rela <- filter_by_prevalence(data$simulated_data_rela, 0.1)
  data$simulated_data_rela <- drop_allzero_samples(data$simulated_data_rela)
  
  common_ids <- Reduce(intersect, list(
    rownames(data$simulated_data_abs),
    rownames(data$simulated_data_rela),
    data$metadata$SampleID
  ))
  
  common_features <- Reduce(intersect, list(
    colnames(data$simulated_data_abs),
    colnames(data$simulated_data_rela)
  ))
  
  if (length(common_ids) == 0) stop("No common samples left after pooled filtering/alignment.")
  if (length(common_features) == 0) stop("No common features left after pooled filtering/alignment.")
  
  data$simulated_data_abs  <- data$simulated_data_abs[common_ids, common_features, drop = FALSE]
  data$simulated_data_rela <- data$simulated_data_rela[common_ids, common_features, drop = FALSE]
  data$metadata <- data$metadata[match(common_ids, data$metadata$SampleID), , drop = FALSE]
  
  stopifnot(nrow(data$simulated_data_abs) == nrow(data$metadata))
  stopifnot(all(rownames(data$simulated_data_abs) == data$metadata$SampleID))
  stopifnot(nrow(data$simulated_data_rela) == nrow(data$metadata))
  stopifnot(all(rownames(data$simulated_data_rela) == data$metadata$SampleID))
  
  ##### Step 3: pooled 设计变量清理 #####
  pooled_design_info <- prepare_design_info(data$metadata)
  if (is.null(pooled_design_info)) stop("Pooled data has fewer than 2 Group levels after filtering.")
  
  data$metadata <- pooled_design_info$metadata
  pooled_covariates <- pooled_design_info$covariates
  pooled_categorical_variable_name <- pooled_design_info$categorical_variable_name
  
  ##### Step 4: 构建 batch-specific data + design #####
  if (!("Batch" %in% colnames(data$metadata))) stop("Batch column is required for meta-analysis workflow.")
  
  batch_levels <- unique(as.character(data$metadata$Batch))
  batch_data_list <- list()
  
  # batch 内 minimal-drop 配置
  batch_minimal_max_drop <- 2
  
  for (batch_name in batch_levels) {
    message("========== Processing Batch: ", batch_name, " ==========")
    
    sample_ids <- as.character(data$metadata$SampleID[as.character(data$metadata$Batch) == batch_name])
    if (length(sample_ids) == 0) next
    
    batch_abs  <- data$simulated_data_abs[sample_ids, , drop = FALSE]
    batch_rela <- data$simulated_data_rela[sample_ids, , drop = FALSE]
    batch_meta <- data$metadata[match(sample_ids, data$metadata$SampleID), , drop = FALSE]
    
    message("[Batch ", batch_name, "] Initial sample size = ", nrow(batch_meta),
            ", feature size = ", ncol(batch_abs))
    
    # batch 内再次过滤
    batch_abs  <- remove_low_feature_samples(batch_abs)
    batch_abs  <- filter_by_prevalence(batch_abs, 0.1)
    batch_abs <- drop_allzero_samples(batch_abs)
    batch_rela <- remove_low_feature_samples(batch_rela)
    batch_rela <- filter_by_prevalence(batch_rela, 0.1)
    batch_rela <- drop_allzero_samples(batch_rela)
    
    keep_ids <- Reduce(intersect, list(rownames(batch_abs), rownames(batch_rela), batch_meta$SampleID))
    keep_features <- Reduce(intersect, list(colnames(batch_abs), colnames(batch_rela)))
    if (length(keep_ids) == 0 || length(keep_features) == 0) {
      message("[Batch ", batch_name, "] Skip: no common samples or features left after batch filtering.")
      next
    }
    
    batch_abs  <- batch_abs[keep_ids, keep_features, drop = FALSE]
    batch_rela <- batch_rela[keep_ids, keep_features, drop = FALSE]
    batch_meta <- batch_meta[match(keep_ids, batch_meta$SampleID), , drop = FALSE]
    
    message("[Batch ", batch_name, "] After batch filtering: sample size = ", nrow(batch_meta),
            ", feature size = ", ncol(batch_abs))
    
    batch_data <- list(
      simulated_data_rela = batch_rela,
      simulated_data_abs = batch_abs,
      metadata = batch_meta
    )
    
    batch_design_info <- prepare_design_info(batch_data$metadata)
    if (is.null(batch_design_info)) {
      message("[Batch ", batch_name, "] Skip: fewer than 2 Group levels after design cleanup.")
      next
    }
    
    # 新增：batch 内部显式检查是否满秩，并最小删除协变量
    batch_check <- check_and_fix_batch_design(
      batch_name = batch_name,
      batch_meta = batch_design_info$metadata,
      covariates = batch_design_info$covariates,
      categorical_variable_name = batch_design_info$categorical_variable_name,
      max_drop = batch_minimal_max_drop,
      verbose = TRUE
    )
    if (is.null(batch_check)) next
    
    batch_data$metadata <- batch_check$metadata
    
    stopifnot(nrow(batch_data$simulated_data_abs) == nrow(batch_data$metadata))
    stopifnot(all(rownames(batch_data$simulated_data_abs) == batch_data$metadata$SampleID))
    stopifnot(nrow(batch_data$simulated_data_rela) == nrow(batch_data$metadata))
    stopifnot(all(rownames(batch_data$simulated_data_rela) == batch_data$metadata$SampleID))
    
    batch_data_list[[batch_name]] <- list(
      data = batch_data,
      covariates = batch_check$covariates,
      categorical_variable_name = batch_check$categorical_variable_name
    )
  }
  
  if (length(batch_data_list) == 0) stop("No valid batch left after batch-specific filtering.")
  sample_sizes <- sapply(batch_data_list, function(obj) nrow(obj$data$metadata))
  
  ##### Step 5: 按方法执行 #####
  # not_used_meta 的 minimal-drop 配置（只删协变量，不删Batch）
  minimal_max_drop <- 2
  minimal_allow_drop_batch <- FALSE
  
  if (method_name == "ANCOMBC") {
    
    if (meta_method == "rma.uni") {
      result_list <- run_batch_method_simple(batch_data_list, selected_method_func, normalization)
      
      all_feature <- Reduce(union, lapply(result_list, function(res) res$Feature))
      effect_merge_data <- data.frame(Feature = all_feature)
      se_merge_data <- data.frame(Feature = all_feature)
      
      batch_names <- names(batch_data_list)
      for (i in seq_along(result_list)) {
        res <- result_list[[i]]
        batch <- batch_names[i]
        effect_merge_data[[paste0("effect", batch)]] <- res$lfc_Group1[match(all_feature, res$Feature)]
        se_merge_data[[paste0("se", batch)]] <- res$se_Group1[match(all_feature, res$Feature)]
      }
      
      rownames(effect_merge_data) <- all_feature
      rownames(se_merge_data) <- all_feature
      effect_merge_data <- effect_merge_data[, -1, drop = FALSE]
      se_merge_data <- se_merge_data[, -1, drop = FALSE]
      
      method_res_REML <- meta_combine(effect_merge_data, se_merge_data, length(batch_data_list),
                                      rma_conv = 1e-10, method_type = "REML", rma_maxit = 1000)
      method_res_EB <- meta_combine(effect_merge_data, se_merge_data, length(batch_data_list),
                                    rma_conv = 1e-10, method_type = "EB", rma_maxit = 1000)
      method_res_PM <- meta_combine(effect_merge_data, se_merge_data, length(batch_data_list),
                                    rma_conv = 1e-10, method_type = "PM", rma_maxit = 1000)
      
      write.table(method_res_REML, file = paste0(output_file_name, "_res_", data_index, "_", meta_method, "_REML.tsv"),
                  sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
      write.table(method_res_EB, file = paste0(output_file_name, "_res_", data_index, "_", meta_method, "_EB.tsv"),
                  sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
      write.table(method_res_PM, file = paste0(output_file_name, "_res_", data_index, "_", meta_method, "_PM.tsv"),
                  sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
      
    } else if (meta_method == "Stouffer") {
      result_list <- run_batch_method_simple(batch_data_list, selected_method_func, normalization)
      method_res <- meta_combine_p(result_list = result_list, sample_size = sample_sizes, method = "Stouffer")
      write.table(method_res, file = paste0(output_file_name, "_res_", data_index, "_", meta_method, ".tsv"),
                  sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
      
    } else if (meta_method == "Fisher") {
      result_list <- run_batch_method_simple(batch_data_list, selected_method_func, normalization)
      method_res <- meta_combine_p(result_list = result_list, sample_size = sample_sizes, method = "Fisher")
      write.table(method_res, file = paste0(output_file_name, "_res_", data_index, "_", meta_method, ".tsv"),
                  sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
      
    } else if (meta_method == "not_used_meta") {
      covariates_use <- add_unique(pooled_covariates, "Batch")
      fix <- drop_covariates_to_full_rank_minimal(
        meta = data$metadata,
        covariates = covariates_use,
        protect = "Batch",
        always_keep = "Group",
        max_drop = minimal_max_drop,
        allow_drop_protect = minimal_allow_drop_batch,
        verbose = TRUE,
        label = "not_used_meta"
      )
      covariates_use2 <- fix$covariates
      
      method_res <- selected_method_func(data, covariates_use2, normalization)
      write.table(method_res, file = paste0(output_file_name, "_res_", data_index, "_no_meta.tsv"),
                  sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
    }
    
  } else if (grepl("Maaslin2", method_name)) {
    
    if (meta_method == "rma.uni") {
      result_list <- run_batch_method_maaslin(batch_data_list, selected_method_func, normalization, transform)
      result_list <- lapply(result_list, function(res) res[which(res$metadata == "Group"), ])
      
      all_feature <- Reduce(union, lapply(result_list, function(res) res$Feature))
      effect_merge_data <- data.frame(Feature = all_feature)
      se_merge_data <- data.frame(Feature = all_feature)
      
      batch_names <- names(batch_data_list)
      for (i in seq_along(result_list)) {
        res <- result_list[[i]]
        batch <- batch_names[i]
        effect_merge_data[[paste0("effect_", batch)]] <- res$coef[match(all_feature, res$Feature)]
        se_merge_data[[paste0("se_", batch)]] <- res$stderr[match(all_feature, res$Feature)]
      }
      
      rownames(effect_merge_data) <- all_feature
      rownames(se_merge_data) <- all_feature
      effect_merge_data <- effect_merge_data[, -1, drop = FALSE]
      se_merge_data <- se_merge_data[, -1, drop = FALSE]
      
      method_res_REML <- meta_combine(effect_merge_data, se_merge_data, length(batch_data_list),
                                      rma_conv = 1e-10, method_type = "REML", rma_maxit = 1000)
      method_res_EB <- meta_combine(effect_merge_data, se_merge_data, length(batch_data_list),
                                    rma_conv = 1e-10, method_type = "EB", rma_maxit = 1000)
      method_res_PM <- meta_combine(effect_merge_data, se_merge_data, length(batch_data_list),
                                    rma_conv = 1e-10, method_type = "PM", rma_maxit = 1000)
      
      write.table(method_res_REML, file = paste0(output_file_name, "_res_", data_index, "_", meta_method, "_REML.tsv"),
                  sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
      write.table(method_res_EB, file = paste0(output_file_name, "_res_", data_index, "_", meta_method, "_EB.tsv"),
                  sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
      write.table(method_res_PM, file = paste0(output_file_name, "_res_", data_index, "_", meta_method, "_PM.tsv"),
                  sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
      
    } else if (meta_method == "Stouffer") {
      result_list <- run_batch_method_maaslin(batch_data_list, selected_method_func, normalization, transform)
      result_list <- lapply(result_list, function(res) res[which(res$metadata == "Group"), ])
      method_res <- meta_combine_p(result_list = result_list, sample_size = sample_sizes, method = "Stouffer")
      write.table(method_res, file = paste0(output_file_name, "_res_", data_index, "_", meta_method, ".tsv"),
                  sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
      
    } else if (meta_method == "Fisher") {
      result_list <- run_batch_method_maaslin(batch_data_list, selected_method_func, normalization, transform)
      result_list <- lapply(result_list, function(res) res[which(res$metadata == "Group"), ])
      method_res <- meta_combine_p(result_list = result_list, sample_size = sample_sizes, method = "Fisher")
      write.table(method_res, file = paste0(output_file_name, "_res_", data_index, "_", meta_method, ".tsv"),
                  sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
      
    } else if (meta_method == "not_used_meta") {
      covariates_use <- add_unique(pooled_covariates, "Batch")
      categorical_use <- add_unique(pooled_categorical_variable_name, "Batch")
      
      fix <- drop_covariates_to_full_rank_minimal(
        meta = data$metadata,
        covariates = covariates_use,
        protect = "Batch",
        always_keep = "Group",
        max_drop = minimal_max_drop,
        allow_drop_protect = minimal_allow_drop_batch,
        verbose = TRUE,
        label = "not_used_meta"
      )
      covariates_use2 <- fix$covariates
      
      categorical_use2 <- intersect(as.character(categorical_use), as.character(covariates_use2))
      if (length(categorical_use2) == 0) categorical_use2 <- NULL
      
      method_res <- selected_method_func(data, covariates_use2, categorical_use2, normalization, transform)
      method_res <- method_res[which(method_res$metadata == "Group"), ]
      
      write.table(method_res, file = paste0(output_file_name, "_res_", data_index, "_no_meta.tsv"),
                  sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
    }
    
  } else if (method_name == "DESeq2") {
    
    if (meta_method == "rma.uni") {
      result_list <- run_batch_method_simple(batch_data_list, selected_method_func, normalization)
      
      all_feature <- Reduce(union, lapply(result_list, function(res) res$Feature))
      effect_merge_data <- data.frame(Feature = all_feature)
      se_merge_data <- data.frame(Feature = all_feature)
      
      batch_names <- names(batch_data_list)
      for (i in seq_along(result_list)) {
        res <- result_list[[i]]
        batch <- batch_names[i]
        effect_merge_data[[paste0("effect", batch)]] <- res$log2FoldChange[match(all_feature, res$Feature)]
        se_merge_data[[paste0("se", batch)]] <- res$lfcSE[match(all_feature, res$Feature)]
      }
      
      rownames(effect_merge_data) <- all_feature
      rownames(se_merge_data) <- all_feature
      effect_merge_data <- effect_merge_data[, -1, drop = FALSE]
      se_merge_data <- se_merge_data[, -1, drop = FALSE]
      
      method_res_REML <- meta_combine(effect_merge_data, se_merge_data, length(batch_data_list),
                                      rma_conv = 1e-10, method_type = "REML", rma_maxit = 1000)
      method_res_EB <- meta_combine(effect_merge_data, se_merge_data, length(batch_data_list),
                                    rma_conv = 1e-10, method_type = "EB", rma_maxit = 1000)
      method_res_PM <- meta_combine(effect_merge_data, se_merge_data, length(batch_data_list),
                                    rma_conv = 1e-10, method_type = "PM", rma_maxit = 1000)
      
      write.table(method_res_REML, file = paste0(output_file_name, "_res_", data_index, "_", meta_method, "_REML.tsv"),
                  sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
      write.table(method_res_EB, file = paste0(output_file_name, "_res_", data_index, "_", meta_method, "_EB.tsv"),
                  sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
      write.table(method_res_PM, file = paste0(output_file_name, "_res_", data_index, "_", meta_method, "_PM.tsv"),
                  sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
      
    } else if (meta_method == "Stouffer") {
      result_list <- run_batch_method_simple(batch_data_list, selected_method_func, normalization)
      method_res <- meta_combine_p(result_list = result_list, sample_size = sample_sizes, method = "Stouffer")
      write.table(method_res, file = paste0(output_file_name, "_res_", data_index, "_", meta_method, ".tsv"),
                  sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
      
    } else if (meta_method == "Fisher") {
      result_list <- run_batch_method_simple(batch_data_list, selected_method_func, normalization)
      method_res <- meta_combine_p(result_list = result_list, sample_size = sample_sizes, method = "Fisher")
      write.table(method_res, file = paste0(output_file_name, "_res_", data_index, "_", meta_method, ".tsv"),
                  sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
      
    } else if (meta_method == "not_used_meta") {
      covariates_use <- add_unique(pooled_covariates, "Batch")
      fix <- drop_covariates_to_full_rank_minimal(
        meta = data$metadata,
        covariates = covariates_use,
        protect = "Batch",
        always_keep = "Group",
        max_drop = minimal_max_drop,
        allow_drop_protect = minimal_allow_drop_batch,
        verbose = TRUE,
        label = "not_used_meta"
      )
      covariates_use2 <- fix$covariates
      
      method_res <- selected_method_func(data, covariates_use2, normalization)
      write.table(method_res, file = paste0(output_file_name, "_res_", data_index, "_no_meta.tsv"),
                  sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
    }
    
  } else {
    
    if (meta_method == "Stouffer") {
      result_list <- run_batch_method_simple(batch_data_list, selected_method_func, normalization)
      method_res <- meta_combine_p(result_list = result_list, sample_size = sample_sizes, method = "Stouffer")
      write.table(method_res, file = paste0(output_file_name, "_res_", data_index, "_", meta_method, ".tsv"),
                  sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
      
    } else if (meta_method == "Fisher") {
      result_list <- run_batch_method_simple(batch_data_list, selected_method_func, normalization)
      method_res <- meta_combine_p(result_list = result_list, sample_size = sample_sizes, method = "Fisher")
      write.table(method_res, file = paste0(output_file_name, "_res_", data_index, "_", meta_method, ".tsv"),
                  sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
      
    } else if (meta_method == "not_used_meta") {
      covariates_use <- add_unique(pooled_covariates, "Batch")
      fix <- drop_covariates_to_full_rank_minimal(
        meta = data$metadata,
        covariates = covariates_use,
        protect = "Batch",
        always_keep = "Group",
        max_drop = minimal_max_drop,
        allow_drop_protect = minimal_allow_drop_batch,
        verbose = TRUE,
        label = "not_used_meta"
      )
      covariates_use2 <- fix$covariates
      
      method_res <- selected_method_func(data, covariates_use2, normalization)
      write.table(method_res, file = paste0(output_file_name, "_res_", data_index, "_no_meta.tsv"),
                  sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
    }
  }
}

