# ============================================================================
# 基因相关性分析
# ============================================================================
# 本文件包含基因相关性分析相关的代码
# ============================================================================

# 加载function
source("scripts/utils.R")

# 关注的基因并计算基因组间比较
glist<-c("VCAM1","PECAM1","ACTA2","PTPRC","TGFB1","MKI67")
tmp<-dat[glist,]
library(dplyr)
library(tidyr)
groups <- unique(tmp@meta.data$Group_seperate)
results <- data.frame()
for (gene in glist) {
  tmp@meta.data[[paste0(gene, "_expression")]] <- FetchData(tmp, vars = gene)
  for (pair in combn(groups, 2, simplify = FALSE)) {
    group1 <- pair[[1]]
    group2 <- pair[[2]]
    data1 <- tmp@meta.data %>%
      filter(Group_seperate == group1) %>%
      pull(paste0(gene, "_expression")) 
    data2 <- tmp@meta.data %>%
      filter(Group_seperate == group2) %>%
      pull(paste0(gene, "_expression"))
    mean1<-mean(data1[,1])
    mean2<-mean(data2[,1])
    t_test <- t.test(data1, data2)
    wilcox_test <- wilcox.test(data1[,1], data2[,1])
    results <- rbind(results, data.frame(
      Gene = gene,
      Group1 = group1,
      Group2 = group2,
      Group1.mean=mean1,
      Group2.mean=mean2,
      T_Test_P_Value = t_test$p.value,
      Wilcox_Test_P_Value = wilcox_test$p.value
    ))
  }
}
write.csv(results, file = "gene_GroupSep_comparison_results.csv", row.names = FALSE)







