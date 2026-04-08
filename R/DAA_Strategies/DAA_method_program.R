Corncob_analysis <- function(simulated_data_all, covariates=NULL, normalization,
                             timeout_sec = 300,     # 每个特征最多跑60秒（硬切断）
                             min_nz_per_group = 0  # 统一过滤阈值；不想过滤就设0
){
  library(corncob)
  
  # ---- 硬切断依赖：callr（跨平台、真正 kill 子进程）----
  if(!requireNamespace("callr", quietly = TRUE)){
    stop("需要安装 callr 才能硬切断：install.packages('callr')")
  }
  
  simulated_data <- data_normalization(simulated_data_all, normalization = normalization)
  
  covariates_full <- c("Group", covariates)
  covariates_str  <- paste(covariates_full, collapse = " + ")  # 修正：始终 collapse
  phi.formula_str <- paste0("~ ", covariates_str)
  formula_str     <- paste0("cbind(W, M - W) ~ ", covariates_str)
  
  pval_vector <- c()
  Feature <- c()
  effect_size_vector <- c()
  
  M_vec <- colSums(t(simulated_data_all$simulated_data_abs))  # 长度=样本数
  
  # 子进程里执行的函数（必须是“自包含”的）
  fit_one <- function(formula_str, phi_str, one_feature_data){
    library(corncob)
    
    fit <- bbdml(
      formula     = as.formula(formula_str),
      phi.formula = as.formula(phi_str),
      data        = one_feature_data
    )
    
    ct <- summary(fit)$coefficients
    
    # 兼容 Group 既可能是数值(行名 mu.Group)，也可能是因子(行名 mu.Groupxxx)
    rn <- if("mu.Group" %in% rownames(ct)) "mu.Group" else grep("^mu\\.Group", rownames(ct), value = TRUE)[1]
    if(is.na(rn) || length(rn) == 0) stop("Cannot find Group coefficient in model summary.")
    
    pval <- ct[rn, 4]
    eff  <- ct[rn, 1]
    c(pval, eff)
  }
  
  for(p in 1:ncol(simulated_data)){
    W_vec <- as.numeric(t(simulated_data)[p, ])
    one_feature_data <- cbind(simulated_data_all$metadata, W = W_vec, M = as.numeric(M_vec))
    
    # 可选：统一规则过滤（建议至少把 nz==0 这种 separation 直接 NA，避免无意义拟合）
    if(min_nz_per_group > 0){
      nz_by_group <- tapply(one_feature_data$W > 0, one_feature_data$Group, sum)
      if(any(nz_by_group < min_nz_per_group)){
        Feature <- append(Feature, colnames(simulated_data)[p])
        pval_vector <- append(pval_vector, NA)
        effect_size_vector <- append(effect_size_vector, NA)
        next
      }
    }
    
    # ---- 硬切断：超时会杀掉子进程 ----
    result <- tryCatch({
      callr::r(
        func = fit_one,
        args = list(formula_str = formula_str,
                    phi_str = phi.formula_str,
                    one_feature_data = one_feature_data),
        timeout = timeout_sec
      )
    }, error = function(e){
      message(sprintf("Timeout/Error at p = %d (%s): %s", p, colnames(simulated_data)[p], e$message))
      c(NA, NA)
    })
    
    Feature <- append(Feature, colnames(simulated_data)[p])
    pval_vector <- append(pval_vector, result[1])
    effect_size_vector <- append(effect_size_vector, result[2])
  }
  
  p.adj.val <- p.adjust(pval_vector, method = "BH")
  data.frame(Feature = Feature, p.adj.val = p.adj.val, effect_size = effect_size_vector)
}
# Corncob_analysis<-function(simulated_data_all,covariates=NULL,normalization){
#   library(corncob)
#   simulated_data <- data_normalization(simulated_data_all,normalization = normalization)
#   col_names <- covariates
#   covariates <- c("Group",covariates)
#   if(length(covariates)>1){
#     covariates <- paste(covariates,collapse = ' + ')
#   }
#   phi.formula_str <- paste0('~ ',covariates)
#   formula_str <- paste0('cbind(W, M - W) ~ ',covariates)
#   pval_vector <- c()
#   Feature <- c()
#   effect_size_vector <- c()
#   for(p in 1:ncol(simulated_data)){
#     one_feature_data <- cbind(simulated_data_all$metadata, 
#                               W = unlist(t(simulated_data)[p, ]),
#                               M = colSums(t(simulated_data_all$simulated_data_abs)))
#     
#     # 使用 tryCatch 捕获错误
#     result <- tryCatch({
#       corncob <- bbdml(formula = as.formula(formula_str),
#                        phi.formula = as.formula(phi.formula_str),
#                        data = one_feature_data)
#       pval <- data.frame(summary(corncob)$coefficients)['mu.Group', 4]
#       effect_size <- data.frame(summary(corncob)$coefficients)['mu.Group', 1]
#       c(pval,effect_size)  # 返回 pval 作为结果
#     }, error = function(e) {
#       message(sprintf("Error at p = %d: %s", p, e$message))  # 打印错误信息
#       c(NA,NA)  # 如果发生错误，返回 NA
#     })
#     Feature <- append(Feature,colnames(simulated_data)[p])
#     pval_vector <- append(pval_vector,result[1])
#     effect_size_vector <- append(effect_size_vector,result[2])
#   }
#   p.adj.val = p.adjust(pval_vector,method = 'BH')
#   Corncob_res <- data.frame(Feature = Feature, p.adj.val = p.adj.val,effect_size = effect_size_vector)
#   return(Corncob_res)
# }
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
ALDEx2_analysis <- function(simulated_data_all,covariates=NULL,normalization){
  library(ALDEx2)
  simulated_data <- data_normalization(simulated_data_all,normalization = normalization)
  col_names <- covariates
  covariates <- c("Group",covariates)
  if(length(covariates)>1){
    covariates <- paste(covariates,collapse = ' + ')
  }
  formula_str <- paste0('~ ',covariates)
  mm <- model.matrix(as.formula(formula_str), simulated_data_all$metadata)
  x.glm <- aldex.clr(t(simulated_data), mm, denom="all", verbose=F)
  glm.test <- aldex.glm(x.glm, mm, fdr.method='BH')
  aldex2_res <- glm.test
  aldex2_res$Feature  <- rownames(aldex2_res)
  colnames(aldex2_res)[which(colnames(aldex2_res) == 'Group:pval.padj')] <- 'p.adj.val'
  return(aldex2_res)
}
Maaslin2_LM_analysis<-function(simulated_data_all,covariates=NULL,categorical_variable_name = NULL,normalization,transform){
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
  filename <- "data_maaslin2_lm"
  rownames(meta)<-meta$Sample
  fit_data = Maaslin2(
    input_data = data, 
    input_metadata = meta, 
    min_prevalence = 0,reference = reference,
    normalization = normalization,analysis_method = "LM",
    output = filename,transform = transform,plot_heatmap = FALSE,plot_scatter = FALSE,
    fixed_effects = covariates,
    random_effects = random_effect_variable)
  maa_lm_res<-fit_data$results
  colnames(maa_lm_res)[1] <- 'Feature'
  colnames(maa_lm_res)[8] <- 'p.adj.val'
  return(maa_lm_res)
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
ZicoSeq_analysis <- function(simulated_data_all,covariates=NULL,normalization){
  library(GUniFrac)
  simulated_data <- data_normalization(simulated_data_all,normalization = normalization)
  
  #count数据输入
  ZicoSeq.obj <- tryCatch({
    feature_name = c()
    if(normalization == 'count'){
      ZicoSeq.res <- ZicoSeq(meta.dat = simulated_data_all$metadata, feature.dat = t(simulated_data), 
                             grp.name = 'Group', adj.name = covariates, feature.dat.type = "count",
                             # Winsorization to replace outliers
                             is.winsor = TRUE, outlier.pct = 0.03, winsor.end = 'top',
                             # Posterior sampling 
                             is.post.sample = TRUE, post.sample.no = 25, 
                             # Use the square-root transformation
                             link.func = list(function (x) x^0.5), stats.combine.func = max,
                             # Permutation-based multiple testing correction
                             perm.no = 999,  strata = NULL, 
                             # Reference-based multiple stage normalization
                             ref.pct = 0.5, stage.no = 6, excl.pct = 0.2,
                             # Family-wise error rate control
                             is.fwer = FALSE, verbose = FALSE, return.feature.dat = FALSE)
    }
    #相对丰度数据输入
    if(normalization == 'TSS'){
      ZicoSeq.res <- ZicoSeq(meta.dat = simulated_data_all$metadata, feature.dat = t(simulated_data), 
                             grp.name = 'Group', adj.name = covariates, feature.dat.type = "proportion",
                             # Winsorization to replace outliers
                             is.winsor = TRUE, outlier.pct = 0.03, winsor.end = 'top',
                             # Posterior sampling 
                             is.post.sample = FALSE, post.sample.no = 25, 
                             # Use the square-root transformation
                             link.func = list(function (x) x^0.5), stats.combine.func = max,
                             # Permutation-based multiple testing correction
                             perm.no = 999,  strata = NULL, 
                             # Reference-based multiple stage normalization
                             ref.pct = 0.5, stage.no = 6, excl.pct = 0.2,
                             # Family-wise error rate control
                             is.fwer = FALSE, verbose = FALSE, return.feature.dat = FALSE)
    }
    #其他类型数据输入
    if((normalization != 'TSS') && (normalization != 'count')){
      ZicoSeq.res <- ZicoSeq(meta.dat = simulated_data_all$metadata, feature.dat = t(simulated_data), 
                             grp.name = 'Group', adj.name = covariates, feature.dat.type = "other",
                             # Winsorization to replace outliers
                             is.winsor = FALSE,
                             # Posterior sampling 
                             is.post.sample = FALSE, post.sample.no = 25, 
                             # Use the square-root transformation
                             link.func = list(function (x) x^0.5), stats.combine.func = max,
                             # Permutation-based multiple testing correction
                             perm.no = 999,  strata = NULL, 
                             # Reference-based multiple stage normalization
                             ref.pct = 0.5, stage.no = 6, excl.pct = 0.2,
                             # Family-wise error rate control
                             is.fwer = FALSE, verbose = FALSE, return.feature.dat = FALSE)
    }
    list(ZicoSeq.res,feature_name)
  },error = function(e){
    # 提取错误信息中的特征名
    error_message <- conditionMessage(e)
    features_nums <- as.numeric(str_extract_all(error_message, "\\d+")[[1]])
    
    cat("Removing feature", paste(features_nums, collapse = ", "), "\n")
    # 从原始数据中移除这些特征
    feature_name = colnames(simulated_data_all$simulated_data_abs)[features_nums]
    #print(feature_name)
    simulated_data_all$simulated_data_abs <- simulated_data_all$simulated_data_abs[,-features_nums]
    simulated_data_all$simulated_data_rela <- simulated_data_all$simulated_data_rela[,-features_nums]
    simulated_data <- data_normalization(simulated_data_all,normalization = normalization)
    
    
    if(normalization == 'count'){
      ZicoSeq.res <- ZicoSeq(meta.dat = simulated_data_all$metadata, feature.dat = t(simulated_data), 
                             grp.name = 'Group', adj.name = covariates, feature.dat.type = "count",
                             # Winsorization to replace outliers
                             is.winsor = TRUE, outlier.pct = 0.03, winsor.end = 'top',
                             # Posterior sampling 
                             is.post.sample = TRUE, post.sample.no = 25, 
                             # Use the square-root transformation
                             link.func = list(function (x) x^0.5), stats.combine.func = max,
                             # Permutation-based multiple testing correction
                             perm.no = 999,  strata = NULL, 
                             # Reference-based multiple stage normalization
                             ref.pct = 0.5, stage.no = 6, excl.pct = 0.2,
                             # Family-wise error rate control
                             is.fwer = FALSE, verbose = FALSE, return.feature.dat = FALSE)
    }
    #相对丰度数据输入
    if(normalization == 'TSS'){
      ZicoSeq.res <- ZicoSeq(meta.dat = simulated_data_all$metadata, feature.dat = t(simulated_data), 
                             grp.name = 'Group', adj.name = covariates, feature.dat.type = "proportion",
                             # Winsorization to replace outliers
                             is.winsor = TRUE, outlier.pct = 0.03, winsor.end = 'top',
                             # Posterior sampling 
                             is.post.sample = FALSE, post.sample.no = 25, 
                             # Use the square-root transformation
                             link.func = list(function (x) x^0.5), stats.combine.func = max,
                             # Permutation-based multiple testing correction
                             perm.no = 999,  strata = NULL, 
                             # Reference-based multiple stage normalization
                             ref.pct = 0.5, stage.no = 6, excl.pct = 0.2,
                             # Family-wise error rate control
                             is.fwer = FALSE, verbose = FALSE, return.feature.dat = FALSE)
    }
    #其他类型数据输入
    if((normalization != 'TSS') && (normalization != 'count')){
      ZicoSeq.res <- ZicoSeq(meta.dat = simulated_data_all$metadata, feature.dat = t(simulated_data), 
                             grp.name = 'Group', adj.name = covariates, feature.dat.type = "other",
                             # Winsorization to replace outliers
                             is.winsor = FALSE,
                             # Posterior sampling 
                             is.post.sample = FALSE, post.sample.no = 25, 
                             # Use the square-root transformation
                             link.func = list(function (x) x^0.5), stats.combine.func = max,
                             # Permutation-based multiple testing correction
                             perm.no = 999,  strata = NULL, 
                             # Reference-based multiple stage normalization
                             ref.pct = 0.5, stage.no = 6, excl.pct = 0.2,
                             # Family-wise error rate control
                             is.fwer = FALSE, verbose = FALSE, return.feature.dat = FALSE)
    }
    list(ZicoSeq.res,feature_name)
  })
  #print(ZicoSeq.obj[[2]])
  #print(as.numeric(data.frame(ZicoSeq.obj[[1]]$coef.list)['Group',]))
  ZicoSeq_res<-data.frame(Feature = c(names(ZicoSeq.obj[[1]]$p.adj.fdr),ZicoSeq.obj[[2]]),
                          effect_size = c(as.numeric(data.frame(ZicoSeq.obj[[1]]$coef.list)['Group',]),rep(NA,length(ZicoSeq.obj[[2]]))),
                          pval = c(ZicoSeq.obj[[1]]$p.raw,rep(NA,length(ZicoSeq.obj[[2]]))),
                          p.adj.val = c(ZicoSeq.obj[[1]]$p.adj.fdr,rep(NA,length(ZicoSeq.obj[[2]]))))
  return(ZicoSeq_res)
}
fastANCOM_analysis<-function(simulated_data_all,covariates=NULL,categorical_variable_name = NULL,normalization){
  ###covariates以向量形式输入
  library(fastANCOM)
  simulated_data <- data_normalization(simulated_data_all,normalization = normalization)
  if('subject_id' %in% categorical_variable_name){
    categorical_variable_name <- setdiff(categorical_variable_name,'subject_id')
  }
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
    msg2 <- gsub("\\\\n", "\n", error_message)
    
    taxa_block <- sub(
      "(?s).*following taxa:\\s*\\n?(.*?)\\nPlease remove.*",
      "\\1",
      msg2,
      perl = TRUE
    )
    
    taxa_vec <- strsplit(trimws(taxa_block), ",\\s*", perl = TRUE)[[1]]
    taxa_vec <- taxa_vec[nzchar(taxa_vec)]
    features_to_remove <- taxa_vec
    
    if (length(features_to_remove) > 0) {
      cat("Removing features with zero variance:", paste(features_to_remove, collapse = ", "), "\n")
      # 从原始数据中移除这些特征
      colnames <- colnames(data_normalization(simulated_data_all,normalization = normalization))
      colnames(simulated_data_all$simulated_data_abs) <- colnames
      simulated_data_all$simulated_data_abs <- simulated_data_all$simulated_data_abs[, !colnames(simulated_data_all$simulated_data_abs) %in% features_to_remove]
      # 递归调用自身重新运行
      rownames <- rownames(simulated_data_all$simulated_data_abs)
      simulated_data <- data_normalization(simulated_data_all,normalization = normalization)
      rownames(simulated_data) <- rownames
      #print(rownames(simulated_data_all$metadata))
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
  if(normalization == 'TSS'){
    min.row.sum = 1
  }else{
    min.row.sum = 5
  }
  
  y <- estimateDisp(y,design,min.row.sum = min.row.sum)
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
LM_random_effect<-function(simulated_data_all,covariates=NULL,categorical_variable_name=NULL,normalization){
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
LM_fixed_effect<-function(simulated_data_all,covariates=NULL,normalization){
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
LM_analysis <- function(simulated_data_all, normalization) {
  simulated_data <- data_normalization(simulated_data_all, normalization = normalization)
  lm_res <- data.frame(Feature = rep(NA, ncol(simulated_data)),
                       Pvalue = rep(NA, ncol(simulated_data)),
                       Effect_size = rep(NA, ncol(simulated_data)))
  
  for (s in 1:ncol(simulated_data)) {
    temp_data <- data.frame(cbind(simulated_data[, s], simulated_data_all$metadata))
    colnames(temp_data)[1] <- 'Feature'
    lm_v <- lm(Feature ~ Group, data = temp_data)
    
    # 获取以 "Group" 开头的系数
    group_coeffs <- grep("^Group", rownames(summary(lm_v)$coefficients), value = TRUE)
    
    if (length(group_coeffs) > 0) {
      pval <- summary(lm_v)$coefficients[group_coeffs, 'Pr(>|t|)']
      estimate_effect <- summary(lm_v)$coefficients[group_coeffs, 'Estimate']
      lm_res[s, ] <- c(colnames(simulated_data)[s], pval, estimate_effect)
    }
  }
  
  p.adj.val <- p.adjust(lm_res$Pvalue, method = 'BH')
  lm_res$p.adj.val <- p.adj.val
  return(lm_res)
}
LM_interaction_analysis<-function(simulated_data_all,covariates=NULL,normalization){
  simulated_data <- data_normalization(simulated_data_all,normalization = normalization)
  covariates_str <- paste(covariates, collapse = " + ")
  interactions_str <- paste("Group:", covariates, collapse = " + ")
  # 构建完整的公式字符串
  formula_str <- paste("Feature ~ Group +", covariates_str, "+", interactions_str)
  lminter_res <- data.frame(Feature = rep(NA,ncol(simulated_data_all$simulated_data_rela)),Pvalue=rep(NA,ncol(simulated_data_all$simulated_data_rela)),Effect_size = rep(NA,ncol(simulated_data_all$simulated_data_rela)))
  for(k in 1:ncol(simulated_data_all$simulated_data_rela)){
    temp_data <- data.frame(cbind(simulated_data[,k],simulated_data_all$metadata))
    colnames(temp_data)[1]<-'Feature'
    lm_inter <- lm(as.formula(formula_str),data = temp_data)
    # 获取以 "Group" 开头的系数
    group_coeffs <- grep("^Group", rownames(summary(lm_inter)$coefficients), value = TRUE)
    if (length(group_coeffs) > 0) {
      pval <- summary(lm_inter)$coefficients[group_coeffs,'Pr(>|t|)']
      estimate_effect <- summary(lm_inter)$coefficients[group_coeffs,'Estimate']
      lminter_res[k,]<-c(colnames(simulated_data)[k],pval,estimate_effect)
    }
    # pval <- summary(lm_inter)$coefficients['Group','Pr(>|t|)']
    # estimate_effect <- summary(lm_inter)$coefficients['Group','Estimate']
    # lminter_res[k,]<-c(colnames(simulated_data)[k],pval,estimate_effect)
  }
  p.adj.val  = p.adjust(lminter_res$Pvalue,method = 'BH')
  lminter_res$p.adj.val = p.adj.val
  return(lminter_res)
}
median_discretize <- function(data, continuous_vars) {
  for (var in continuous_vars) {
    if (!var %in% names(data)) {
      stop(paste("Variable", var, "not found in the dataset."))
    }
    
    # 计算中位数并划分
    median_value <- median(data[[var]], na.rm = TRUE)
    data <- data %>%
      mutate(
        !!paste0(var, "_cat") := ifelse(
          !!sym(var) <= median_value, "low", "high"
        )
      )
  }
  return(data)
}
merge_categorical_and_continuous <- function(data, covariates, continuous_vars,categorical_variable_name) {
  # 对连续变量离散化
  data <- as.data.frame(data)
  if(is.null(continuous_vars)){
    # 找到以 "Ca" 开头的列
    ca_columns <- setdiff(covariates, continuous_vars)
    # 对这些列逐行拼接
    data$Ca <- apply(data[, ca_columns], 1, paste0, collapse = "")
    return(data)
  }else{
    data <- median_discretize(data, continuous_vars)
    # 动态生成分类变量名
    continuous_cats <- paste0(continuous_vars, "_cat")
    
    ca_columns <- which(colnames(data) %in% categorical_variable_name)
    cat_columns <- which(colnames(data) %in% continuous_cats)
    selected_columns <- c(ca_columns, cat_columns)
    
    # 对这些列逐行拼接
    X <- as.matrix(data[, selected_columns, drop = FALSE])
    data$Ca <- apply(X, 1, paste0, collapse = "")
    return(data)
  }
  
}
Wilcoxon_analysis <- function(simulated_data_all, covariates=NULL, categorical_variable_name=NULL,normalization) {
  library(coin)
  library(dplyr)

  simulated_data <- data_normalization(simulated_data_all, normalization = normalization)
  wilcox_res <- data.frame(
    Feature = rep(NA, ncol(simulated_data_all$simulated_data_rela)), 
    Pvalue = rep(NA, ncol(simulated_data_all$simulated_data_rela)), 
    Effect_size = rep(NA, ncol(simulated_data_all$simulated_data_rela))
  )
  
  if (is.null(covariates)) {
    for (k in 1:ncol(simulated_data_all$simulated_data_rela)) {
      temp_data <- data.frame(cbind(simulated_data[, k], simulated_data_all$metadata))
      colnames(temp_data)[1] <- 'Feature'
      x <- temp_data[which(temp_data$Group == 1), 'Feature']
      y <- temp_data[which(temp_data$Group == 0), 'Feature']
      lfc <- log2(mean(as.numeric(x)) / mean(as.numeric(y)))
      wilcox <- wilcox.test(x = x, y = y)
      pval <- wilcox$p.value
      wilcox_res[k, ] <- c(colnames(simulated_data_all$simulated_data_rela)[k], pval, lfc)
    }
    p.adj.val <- p.adjust(wilcox_res$Pvalue, method = 'BH')
    wilcox_res$p.adj.val <- p.adj.val
    return(wilcox_res)
  } else {
    continuous_vars <- setdiff(covariates,categorical_variable_name)
    simulated_data_all$metadata <- merge_categorical_and_continuous(simulated_data_all$metadata, covariates, continuous_vars,categorical_variable_name)
    #print(simulated_data_all$metadata$Ca)
    formula_str <- "Feature ~ Group | Ca"

    for (k in 1:ncol(simulated_data_all$simulated_data_rela)) {
      tryCatch({
        temp_data <- data.frame(cbind(simulated_data[, k], simulated_data_all$metadata))
        colnames(temp_data)[1] <- 'Feature'
        temp_data$Group <- as.factor(temp_data$Group)
        temp_data$Ca <- as.factor(temp_data$Ca)
        
        
        # 过滤掉 `Ca` 变量水平小于 2 的样本
        ca_counts <- table(temp_data$Ca)
        valid_ca_levels <- names(ca_counts[ca_counts >= 2])
        temp_data <- temp_data[temp_data$Ca %in% valid_ca_levels, ]
        
        # 再次检查 `Ca` 是否还有足够的水平
        if (length(unique(temp_data$Ca)) < 2) {
          formula_str <- "Feature ~ Group | Ca"
        }
        
        wilcox <- wilcox_test(as.formula(formula_str), data = temp_data)
        pval <- pvalue(wilcox)
        x <- temp_data[which(temp_data$Group == 1), 'Feature']
        y <- temp_data[which(temp_data$Group == 0), 'Feature']
        lfc <- log2(mean(as.numeric(x)) / mean(as.numeric(y)))
        wilcox_res[k, ] <- c(colnames(simulated_data_all$simulated_data_rela)[k], pval, lfc)
        
      }, error = function(e) {
        cat(sprintf("block var failed: %s\n", e$message))
        
        if (grepl("less than two observations", e$message)) {
          print("Removing samples with rare `Ca` levels")
          
          # 重新计算 `Ca` 变量的分布，删除 `Ca` 过于稀少的样本
          ca_counts <- table(temp_data$Ca)
          valid_ca_levels <- names(ca_counts[ca_counts >= 2])
          temp_data <- temp_data[temp_data$Ca %in% valid_ca_levels, ]
          
          # 如果仍然不满足 Wilcoxon 统计要求，则跳过
          if (length(unique(temp_data$Ca)) < 2) {
            print("Skipping feature due to insufficient Ca levels")
            next
          }
          
          # 重新执行 Wilcoxon 检验
          wilcox <- wilcox_test(as.formula(formula_str), data = temp_data)
          pval <- pvalue(wilcox)
          x <- temp_data[which(temp_data$Group == 1), 'Feature']
          y <- temp_data[which(temp_data$Group == 0), 'Feature']
          lfc <- log2(mean(as.numeric(x)) / mean(as.numeric(y)))
          wilcox_res[k, ] <- c(colnames(simulated_data_all$simulated_data_rela)[k], pval, lfc)
          
        } else {
          stop(e)  # 其他错误直接抛出
        }
      })
    }
    
    p.adj.val <- p.adjust(wilcox_res$Pvalue, method = 'BH')
    wilcox_res$p.adj.val <- p.adj.val
    return(wilcox_res)
  }
}
VTwins_analysis<-function(simulated_data_all,normalization){
  simulated_data <- data_normalization(simulated_data_all,normalization = normalization)
  library(VTwins)
  library(dplyr)
  pheno_data <- data.frame(id = simulated_data_all$metadata$SampleID,grp = simulated_data_all$metadata$Group)
  rownames(pheno_data) <- pheno_data$id
  
  pheno_data[,2]<-gsub(pheno_data[,2],pattern = 1,replacement = "grp2")
  pheno_data[,2]<-gsub(pheno_data[,2],pattern = 0,replacement = "grp1")
  vtwins_res <- pair_find(data=data.frame(simulated_data),
                          phenodata=pheno_data[rownames(simulated_data_all$simulated_data_rela),],
                          k="euclidean",
                          Cut_pair=25, 
                          method_choose="Permutation",
                          SavePath = "./",
                          ShuffleWstat = "ShuffleWstat", 
                          BoundarySample = "BoundarySample",
                          BoundaryPair = "BoundaryPair",
                          ShuffleTime=1000,
                          DownPercent = 0.2,
                          Uppercent=0.8,PvalueCutoff=0.05)
  if(is.null(vtwins_res)){
    vtwins_res <- data.frame(Feature = colnames(simulated_data_all$simulated_data_rela),
                             p.adj.val = rep(NA,ncol(simulated_data_all$simulated_data_rela)))
  }else{
    for(k in 1:nrow(vtwins_res)){
      vtwins_res$p.adj.val[k] <- ifelse(vtwins_res$Decre.aveRank.P.FDR[k] < vtwins_res$Incre.aveRank.P.FDR[k],vtwins_res$Decre.aveRank.P.FDR[k],vtwins_res$Incre.aveRank.P.FDR[k])
    }
    colnames(vtwins_res)[1] <- 'Feature'
  }
  return(vtwins_res)
}
