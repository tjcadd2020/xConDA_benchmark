###生成具有批次效应的数据
library(SparseDOSSA2)
#不指定ground-truth,不对ground-truth进行生成
sample_num = 10000
#注意：组别和协变量需要存在关联才是confounder
Group = sample(c(rep(0,sample_num/2),rep(1,sample_num/2)))
set.seed(289)
Con1 = ifelse(Group == 0, 
              rnorm(sample_num/2, mean = 2.5, sd = 1), 
              rnorm(sample_num/2, mean = 1, sd = 0.3))
set.seed(623)
Con2 = ifelse(Group == 0, 
              rnorm(sample_num/2, mean = 1, sd = 0.5), 
              rnorm(sample_num/2, mean = 3, sd = 1))

Ca1 = sapply(Group, function(x) {
  if (x == 0) {
    sample(c(0, 1), 1, prob = c(0.3, 0.7))
  } else {
    sample(c(0, 1), 1, prob = c(0.6, 0.4))
  }
})
Batch <- sapply(Group, function(x) {
  if (x == 0) {
    # Group = 0 时的批次分布
    sample(c("A", "B", "C", "D", "E"), 1,
           prob = c(0.08, 0.18, 0.34, 0.25, 0.15))
  } else {
    # Group = 1 时的批次分布
    sample(c("A", "B", "C", "D", "E"), 1,
           prob = c(0.20, 0.10, 0.15, 0.30, 0.25))
  }
})

Batch = as.factor(Batch)
Batch_data <- model.matrix(~ Batch -1)
#Batch = factor(Batch,levels = c(1,2,3))
metadata_matrix <- cbind(as.matrix(data.frame(Group = Group,
                                              Con1 = Con1,Con2 = Con2, Ca1 = Ca1)),Batch_data)
rownames(metadata_matrix) <- paste0('Sample',rownames(metadata_matrix))
IBD_simulation_for_meta <- SparseDOSSA2(template = "IBD", 
                                        n_sample = sample_num, 
                                        new_features = TRUE,n_feature = 1000,median_read_depth = 3000000,metadata_effect_size = 4,perc_feature_spiked_metadata = 0.2,
                                        spike_metadata = "abundance",metadata_matrix = metadata_matrix,
                                        verbose = TRUE)
spike_metadata <- IBD_simulation_for_meta$spike_metadata$feature_metadata_spike_df
###设置协变量的效应值为2
spike_metadata$effect_size[spike_metadata$metadata_datum != 1] <- 2
spike_metadata$effect_size[spike_metadata$metadata_datum == 5] <- 9
spike_metadata$effect_size[spike_metadata$metadata_datum == 6] <- 10
spike_metadata$effect_size[spike_metadata$metadata_datum == 7] <- 7
spike_metadata$effect_size[spike_metadata$metadata_datum == 8] <- 8
spike_metadata$effect_size[spike_metadata$metadata_datum == 9] <- 9

IBD_simulation_for_meta_high_he <- SparseDOSSA2(template = "IBD", 
                                                n_sample = sample_num, 
                                                new_features = TRUE,n_feature = 1000,median_read_depth = 3000000,spike_metadata = spike_metadata,
                                                metadata_matrix = metadata_matrix,
                                                verbose = TRUE)
simulated_data_for_meta_high_he <- IBD_simulation_for_meta_high_he$simulated_data
simulated_data_for_meta_high_he <- simulated_data_for_meta_high_he[which((rowSums(simulated_data_for_meta_high_he!=0)/sample_num)>0.1),]
simulated_data_for_meta_high_he <- simulated_data_for_meta_high_he[,which(colSums(simulated_data_for_meta_high_he!=0)>10)]
ground_truth_for_meta_high_he <- extract_ground_truth(IBD_simulation_for_meta_high_he,simulated_data_for_meta_high_he)
spike_metadata$effect_size[spike_metadata$metadata_datum != 1] <- 2


IBD_simulation_for_meta_high_he$spike_metadata$metadata_matrix <- as.data.frame(IBD_simulation_for_meta_high_he$spike_metadata$metadata_matrix)
IBD_simulation_for_meta_high_he$spike_metadata$metadata_matrix$Batch<-Batch
IBD_simulation_for_meta_high_he$spike_metadata$metadata_matrix <- IBD_simulation_for_meta_high_he$spike_metadata$metadata_matrix[colnames(simulated_data_for_meta_high_he),]


##绝对丰度转化相对丰度函数
abd_to_rel<-function(feat_abd){
  mat <- apply(feat_abd,2, as.numeric)
  rel_abd<-t(apply(mat,1, function(x){
    if(sum(x) == 0){
      return(x)
    }else{
      return(x/sum(x))
    }
  }))
  rownames(rel_abd)<-rownames(feat_abd)
  return(rel_abd)
}
####提取ground_truth函数
extract_ground_truth<-function(simulation,simulated_data){
  spike_metadata<-simulation$spike_metadata$feature_metadata_spike_df
  spike_metadata <- spike_metadata[spike_metadata$feature_spiked %in% rownames(simulated_data),]
  ground_truth <- spike_metadata
  return(ground_truth)
}
####随机抽样函数
sample_selection<-function(metadata,case_sample_num,control_sample_num){
  ###output selected samples
  case <- rownames(metadata)[metadata[,1]==1]
  control <- rownames(metadata)[metadata[,1]==0]
  case_selected <- sample(case,size = case_sample_num,replace = F)
  control_selected <- sample(control,size = control_sample_num,replace = F)
  selected_sample <- append(case_selected,control_selected)
  
  return(selected_sample)
}
####输出模拟数据函数
output_simulated_data<-function(simulation,simulated_data,ground_truth,data_index,case_sample_num,control_sample_num){
  metadata <- simulation$spike_metadata$metadata_matrix
  rownames(metadata) <- colnames(simulated_data)
  selected_sample<-sample_selection(metadata,case_sample_num,control_sample_num)
  metadata <- metadata[selected_sample,]
  simulated_data_selected_abs <- simulated_data[,selected_sample]
  simulated_data_selected <- abd_to_rel(t(simulated_data_selected_abs))
  metadata <- as.data.frame(metadata)
  metadata$SampleID <- rownames(metadata)
  write.table(ground_truth,file = paste0('ground_truth',data_index,'.tsv'),row.names = F,col.names = T,quote = F,sep = '\t')
  write.table(simulated_data_selected,file = paste0('simulated_data_',data_index,'.tsv'),row.names = T,col.names = T,sep = "\t",quote = F)
  write.table(t(simulated_data_selected_abs),file = paste0('simulated_data_abs_',data_index,'.tsv'),row.names = T,col.names = T,sep = "\t",quote = F)
  write.table(metadata,file = paste0('simulated_metadata',data_index,'.tsv'),row.names = F,col.names = T,quote = F,sep = '\t')
  result = list(ground_truth = ground_truth,
                simulated_data_rela = simulated_data_selected,
                simulated_data_abs = t(simulated_data_selected_abs),
                metadata = metadata)
  return(result)
}
