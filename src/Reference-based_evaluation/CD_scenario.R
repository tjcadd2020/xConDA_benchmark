p_threshold = 0.05
eva_correct_count <- eva_significant[,1:5]
eva_correct_count[,2:5]<-eva_correct_count[,2:5]<p_threshold
eva_HMP_correct_count <- eva_HMP_significant[,1:5]
eva_HMP_correct_count[,2:5]<-eva_HMP_correct_count[,2:5]<p_threshold
eva_HallAB_correct_count <- eva_HallAB_significant[,1:5]
eva_HallAB_correct_count[,2:5]<-eva_HallAB_correct_count[,2:5]<p_threshold
eva_NielsenHB_correct_count <- eva_NielsenHB_significant[,1:5]
eva_NielsenHB_correct_count[,2:5]<-eva_NielsenHB_correct_count[,2:5]<p_threshold
rownames(eva_correct_count) <- eva_correct_count$Method
rownames(eva_HMP_correct_count) <- eva_HMP_correct_count$Method
rownames(eva_HallAB_correct_count) <- eva_HallAB_correct_count$Method
rownames(eva_NielsenHB_correct_count) <- eva_NielsenHB_correct_count$Method
method_order = eva_correct_count$Method
eva_HMP_correct_count <- eva_HMP_correct_count[method_order,]
eva_HallAB_correct_count <- eva_HallAB_correct_count[method_order,]
eva_NielsenHB_correct_count <- eva_NielsenHB_correct_count[method_order,]
all_correct_count <- eva_correct_count[,2:5] + eva_HMP_correct_count[,2:5] + eva_HallAB_correct_count[,2:5] + eva_NielsenHB_correct_count[,2:5]
all_correct_count$Method <- method_order
all_correct_count$Sum <- rowSums(all_correct_count[,1:4])
all_correct_count <- all_correct_count[order(all_correct_count$Sum,decreasing = T),]

method_rename <- function(x){
  x <- gsub(x,pattern = 'maaslin2_lm_res',replacement = 'MaAslin2_GLM',ignore.case = FALSE)
  x <- gsub(x,pattern = 'maaslin2_lm',replacement = 'MaAslin2_GLM',ignore.case = FALSE)
  x <- gsub(x,pattern = 'maaslin2',replacement = 'MaAslin2',ignore.case = FALSE)
  x <- gsub(x,pattern = 'cplm',replacement = 'CPLM',ignore.case = FALSE)
  x <- gsub(x,pattern = 'negbin',replacement = 'NEGBIN',ignore.case = FALSE)
  x <- gsub(x,pattern = 'megbin',replacement = 'NEGBIN',ignore.case = FALSE)
  x <- gsub(x,pattern = 'zinb',replacement = 'ZINB',ignore.case = FALSE)
  x <- gsub(x,pattern = 'LM-fixed',replacement = 'LFEM',ignore.case = FALSE)
  x <- gsub(x,pattern = 'lmem',replacement = 'LMEM',ignore.case = FALSE)
  x <- gsub(x,pattern = 'ANCOMBC',replacement = 'ANCOM-BC2',ignore.case = FALSE)
  return(x)
}
all_correct_count$Method <- method_rename(all_correct_count$Method)
rownames(all_correct_count) <- all_correct_count$Method
#######
library(pheatmap)
library(paletteer)


pheatmap(t(all_correct_count[,1:4]),display_numbers = T,cluster_cols = F,cluster_rows = F,color = colorRampPalette(c("#F8F8F8", "#CC5151"))(5),
         fontsize = 13,fontsize_number = 12,width = 4,height = 4,legend = F,border_color = "white",angle_col = 90,
         number_format = "%.0f", number_color = "grey30")

all_correct_count$Method <- factor(all_correct_count$Method,levels = all_correct_count$Method)
ggplot(all_correct_count[,c('Method','Sum')], aes(x = Method, y = Sum)) +
  geom_col(width = 0.7, alpha = 0.85, fill = "#73A4CA") +
  scale_y_continuous(limits = c(0, 12), expand = expansion(mult = c(0, 0.12))) + 
  labs(x = NULL, y = "Count", title = NULL) +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
######
# 把 all_correct_count 的第1到4列：>=3 记为1，否则0
SRS_count <- all_correct_count
SRS_count[, 1:4] <- as.data.frame(lapply(SRS_count[, 1:4], function(x) as.integer(x >= 3)))
SRS_count$Sum <- rowSums(SRS_count[, 1:4])
SRS_count <- SRS_count[order(SRS_count$Sum,decreasing = T),]
###############################
library(ggplot2)
library(pheatmap)
library(patchwork)

## 1) 保持方法顺序（按当前顺序）
SRS_count$Method <- factor(SRS_count$Method, levels = SRS_count$Method)
SRS_count <- SRS_count[order(SRS_count$Method), ]