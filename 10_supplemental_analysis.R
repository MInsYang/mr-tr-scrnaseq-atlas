# ============================================================================
# 补充分析
# ============================================================================
# 本文件包含补充分析相关的代码，包括试验条件排除，细胞比例统计分析调整等内容
# ============================================================================

# 加载function
source("scripts/utils.R")

library(dplyr)
library(tidyr)
library(broom)


setwd("/data2/Project/F_Group_Analysis/7_Analysis/fanxiu")
dat<-readRDS("../Figure_1104/dat_celltype.rds")
dat$AF<-ifelse(dat$orig.ident %in% c("CE-V24_Mit","CE-V24_Tri","CE-V30_Tri","CE-V40_Mit","CE-V17_Mit"),"AF","SR")
dat$Type<-ifelse(dat$tissue=="Mit","FMR","FTR")
dat$Sample_class<-paste0(dat$Type,"_",dat$AF)
p<-CellDimPlot(dat,group.by="celltype",split.by="Sample_class",bg_color = "white",raster = F)
pdf("1_umap.pdf",width=12,height=10)
print(p)
dev.off()
dat$Grade<-gsub("MV-|TV-","",dat$Group)
dat$Grade<-gsub("moderate","mild",dat$Grade)
dat$Sample_class_Grade<-paste0(dat$Type,"_",dat$AF,"_",dat$Grade)
dat$Sample_class_Grade<-factor(as.character(dat$Sample_class_Grade),
                               levels=c("FMR_AF_normal","FMR_AF_mild","FMR_AF_severe",
                                        "FMR_SR_normal","FMR_SR_mild","FMR_SR_severe",
                                                        "FTR_AF_mild",
                                        "FTR_SR_normal","FTR_SR_mild","FTR_SR_severe"))
dat$Sample_class_Grade<-droplevels(dat$Sample_class_Grade)
saveRDS(dat,file="1_majorCelltype.rds")
p<-CellDimPlot(dat,group.by="celltype",split.by="Sample_class",bg_color = "white",raster = F)
pdf("1_umap.pdf",width=12,height=10)
print(p)
dev.off()
p<-CellDimPlot(dat,group.by="celltype",split.by="Sample_class_Grade",ncol=5,bg_color = "white",raster = F)
pdf("1_umap_2.pdf",width=20,height=10) #这个宽可能需要改，看着差不都就行了
p
dev.off()

dat$Grade<-gsub("mild","Moderate",dat$Grade)
dat$Grade<-gsub("normal","Mild",dat$Grade)
dat$Grade<-gsub("severe","Severe",dat$Grade)


#--------------------------------------------
# Function: run_cluster_anova
#--------------------------------------------
run_cluster_anova <- function(dat) {
  # 检查输入数据
  required_cols <- c("patient", "majorCluster", "loc")
  if (!all(required_cols %in% colnames(dat))) {
    stop("输入数据必须包含以下列: patient, majorCluster, loc")
  }
  
  # Step 1: 计算每个样本中各 majorCluster 的比例
  cluster_prop <- dat %>%
    group_by(patient, loc, majorCluster) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(patient, loc) %>%
    mutate(prop = n / sum(n)) %>%
    ungroup()
  
  # Step 2: 对每个 majorCluster 做 ANOVA
  anova_results <- cluster_prop %>%
    group_by(majorCluster) %>%
    do({
      model <- lm(prop ~ loc, data = .)
      tidy(aov(model)) %>% 
        filter(term == "loc") %>% 
        dplyr::select(term, df, statistic, p.value)
    }) %>%
    ungroup()
  
  # Step 3: 整理结果表格
  anova_results <- anova_results %>%
    rename(Cluster = majorCluster,
           Term = term,
           DF = df,
           F_value = statistic,
           P_value = p.value) %>%
    arrange(P_value)
  
  # Step 4: 返回结果
  return(anova_results)
}



library(Startrac)

vic<-readRDS("/data2/Project/Fuwai_songGrp/7_Analysis/5.VIC_subtype.rds")
vic$fibrotic<-ifelse(vic$subtype %in% c("VIC5","VIC3"),"Anti-fibrotic VIC","Neutral VIC")
vic$fibrotic<-ifelse(vic$subtype %in% c("VIC7","VIC8","VIC9"),"Pro-fibrotic VIC",as.character(vic$fibrotic))
vic$fibrotic<-factor(vic$fibrotic,levels=c("Neutral VIC","Anti-fibrotic VIC","Pro-fibrotic VIC"))
ratio_plot(vic,sample.by="Sample",anno.by="fibrotic",condition.by = "Group",save.prefix = "fibrotic六组",plot_type="box",
           condition.col=k6,strip.col = NULL,vjust.x = 0.5,hjust.x = 1,ord.condition = c("MV-normal",  "MV-moderate", "MV-severe",  "TV-normal",  "TV-moderate", "TV-severe" ))
ratio_plot(vic,sample.by="Sample",anno.by="fibrotic",condition.by = "Group_seperate",save.prefix = "fibrotic四组",plot_type="box",
           condition.col=k6,strip.col = NULL,vjust.x = 0.5,hjust.x = 1,ord.condition = c("MV-Control", "MV-Case","TV-Control","TV-Case"))
#turkey ratio_plot
dir.create("tukey")
ratio_plot(vic,sample.by="Sample",anno.by="fibrotic",condition.by = "Group_seperate",save.prefix = "tukey/fibrotic四组",plot_type="box",
           condition.col=k6,strip.col = NULL,vjust.x = 0.5,hjust.x = 1,compare.method = "tukey",
           ord.condition = c("MV-Control", "MV-Case","TV-Control","TV-Case"))
ratio_plot(vic,sample.by="Sample",anno.by="fibrotic",condition.by = "Group",save.prefix = "tukey/fibrotic六组",plot_type="box",
           condition.col=k6,strip.col = NULL,vjust.x = 0.5,hjust.x = 1,compare.method="tukey",
           ord.condition = c("MV-normal",  "MV-moderate", "MV-severe",  "TV-normal",  "TV-moderate", "TV-severe" ))
ratio_plot(vic,sample.by="Sample",anno.by="subtype",condition.by = "Group",save.prefix = "tukey/VIC_subtype六组",plot_type="box",
           condition.col=k6,strip.col = NULL,vjust.x = 0.5,hjust.x = 1,compare.method="tukey",
           ord.condition = c("MV-normal",  "MV-moderate", "MV-severe",  "TV-normal",  "TV-moderate", "TV-severe" ))



dat  <- vic@meta.data
dat <- dat[,c("Sample","fibrotic","Group","Group_seperate")]
colnames(dat) <- c("patient","majorCluster","loc","loc2")
dat$Cell_Name <- c("VIC")
dat$patient <- as.character(dat$patient)
dat$majorCluster <- as.character(dat$majorCluster)
dat$loc <- as.character(dat$loc)

Roe <- calTissueDist(dat,
                     byPatient = F,
                     colname.cluster = "majorCluster", # 不同细胞亚群
                     colname.patient = "patient", # 不同样本
                     colname.tissue = "loc", # 不同组织
                     method = "chisq", # "chisq", "fisher", and "freq" 
                     min.rowSum = 0) 
Roe
col_fun <- colorRamp2(c(min(Roe, na.rm = TRUE), 1, max(Roe, na.rm = TRUE)), c("blue", "white", "red"))
Roe<-Roe[,c("MV-normal","MV-moderate","MV-severe","TV-normal","TV-moderate","TV-severe")]
Roe<-Roe[c("Neutral VIC","Anti-fibrotic VIC","Pro-fibrotic VIC"),]
colnames(Roe)<-c("MV-Mild","MV-Moderate","MV-Severe","TV-Normal","TV-Moderate","TV-Severe")
p<-Heatmap(as.matrix(Roe),
        show_heatmap_legend = TRUE, 
        cluster_rows = FALSE, 
        cluster_columns = FALSE,
        row_names_side = 'right', 
        show_column_names = TRUE,
        show_row_names = TRUE,
        col = col_fun,
        row_names_gp = gpar(fontsize = 10),
        column_names_gp = gpar(fontsize = 10),
        heatmap_legend_param = list(
          title = "Ro/e Index",  # 自定义图注名称
          at = seq(0.5, 1.5, by = 0.5), # 例刻度的位置/自己的数据必须修改一下！
          labels = seq(0.5, 1.5, by = 0.5) # 每个刻度的标签/自己的数据必须修改一下！
        ),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.2f", Roe[i, j]), x, y, gp = gpar(fontsize = 10, col = "black"))
        }
)
pdf("VIC_fibrotic_ratio_Roe.pdf",width=5,height=3)
print(p)
dev.off()
result_table <- run_cluster_anova(dat)
print(result_table)
write.xlsx(result_table, "VIC_fibrotic_ANOVA_results.xlsx", row.names = FALSE)
## +++, Ro/e > 1;
## ++, 0.8 < Ro/e ≤ 1;
## +, 0.2 ≤ Ro/e ≤ 0.8;
## +/−, 0 < Ro/e < 0.2;
## −, Ro/e = 0



#2
dat  <- vic@meta.data
dat <- dat[,c("Sample","subtype","Group","Group_seperate")]
colnames(dat) <- c("patient","majorCluster","loc","loc2")
dat$Cell_Name <- c("VIC")
dat$patient <- as.character(dat$patient)
dat$majorCluster <- as.character(dat$majorCluster)
dat$loc <- as.character(dat$loc)

Roe <- calTissueDist(dat,
                     byPatient = F,
                     colname.cluster = "majorCluster", # 不同细胞亚群
                     colname.patient = "patient", # 不同样本
                     colname.tissue = "loc", # 不同组织
                     method = "chisq", # "chisq", "fisher", and "freq" 
                     min.rowSum = 0) 
Roe
col_fun <- colorRamp2(c(min(Roe, na.rm = TRUE), 1, max(Roe, na.rm = TRUE)), c("blue", "white", "red"))
Roe<-Roe[,c("MV-normal","MV-moderate","MV-severe","TV-normal","TV-moderate","TV-severe")]
colnames(Roe)<-c("MV-Mild","MV-Moderate","MV-Severe","TV-Normal","TV-Moderate","TV-Severe")
p<-Heatmap(as.matrix(Roe),
           show_heatmap_legend = TRUE, 
           cluster_rows = FALSE, 
           cluster_columns = FALSE,
           row_names_side = 'right', 
           show_column_names = TRUE,
           show_row_names = TRUE,
           col = col_fun,
           row_names_gp = gpar(fontsize = 10),
           column_names_gp = gpar(fontsize = 10),
           heatmap_legend_param = list(
             title = "Ro/e Index",  # 自定义图注名称
             at = seq(0.5, 1.5, by = 0.5), # 例刻度的位置/自己的数据必须修改一下！
             labels = seq(0.5, 1.5, by = 0.5) # 每个刻度的标签/自己的数据必须修改一下！
           ),
           cell_fun = function(j, i, x, y, width, height, fill) {
             grid.text(sprintf("%.2f", Roe[i, j]), x, y, gp = gpar(fontsize = 10, col = "black"))
           }
)
pdf("VIC_subtype_ratio_Roe.pdf",width=5,height=4)
print(p)
dev.off()
result_table <- run_cluster_anova(dat)
print(result_table)
write.xlsx(result_table, "VIC_subtype_ANOVA_results.xlsx", row.names = FALSE)
## +++, Ro/e > 1;
## ++, 0.8 < Ro/e ≤ 1;
## +, 0.2 ≤ Ro/e ≤ 0.8;
## +/−, 0 < Ro/e < 0.2;
## −, Ro/e = 0

#3.VEC
vec<-readRDS("../5.VEC_subtype.rds")
dat  <- vec@meta.data
dat <- dat[,c("Sample","subtype","Group","Group_seperate")]
colnames(dat) <- c("patient","majorCluster","loc","loc2")
dat$patient <- as.character(dat$patient)
dat$majorCluster <- as.character(dat$majorCluster)
ratio_plot(vec,sample.by="Sample",anno.by="subtype",condition.by = "Group",save.prefix = "tukey/VEC_subtype六组",plot_type="box",
           condition.col=k6,strip.col = NULL,vjust.x = 0.5,hjust.x = 1,compare.method="tukey",
           ord.condition = c("MV-normal",  "MV-moderate", "MV-severe",  "TV-normal",  "TV-moderate", "TV-severe" ))


dat$loc <- as.character(dat$loc)
Roe <- calTissueDist(dat,
                     byPatient = F,
                     colname.cluster = "majorCluster", # 不同细胞亚群
                     colname.patient = "patient", # 不同样本
                     colname.tissue = "loc", # 不同组织
                     method = "chisq", # "chisq", "fisher", and "freq" 
                     min.rowSum = 0) 
Roe
col_fun <- colorRamp2(c(min(Roe, na.rm = TRUE), 1, max(Roe, na.rm = TRUE)), c("blue", "white", "red"))
Roe<-Roe[,c("MV-normal","MV-moderate","MV-severe","TV-normal","TV-moderate","TV-severe")]
colnames(Roe)<-c("MV-Mild","MV-Moderate","MV-Severe","TV-Normal","TV-Moderate","TV-Severe")
p<-Heatmap(as.matrix(Roe),
           show_heatmap_legend = TRUE, 
           cluster_rows = FALSE, 
           cluster_columns = FALSE,
           row_names_side = 'right', 
           show_column_names = TRUE,
           show_row_names = TRUE,
           col = col_fun,
           row_names_gp = gpar(fontsize = 10),
           column_names_gp = gpar(fontsize = 10),
           heatmap_legend_param = list(
             title = "Ro/e Index",  # 自定义图注名称
             at = seq(0.5, 1.5, by = 0.5), # 例刻度的位置/自己的数据必须修改一下！
             labels = seq(0.5, 1.5, by = 0.5) # 每个刻度的标签/自己的数据必须修改一下！
           ),
           cell_fun = function(j, i, x, y, width, height, fill) {
             grid.text(sprintf("%.2f", Roe[i, j]), x, y, gp = gpar(fontsize = 10, col = "black"))
           }
)
pdf("VEC_subtype_ratio_Roe.pdf",width=5,height=4)
print(p)
dev.off()

result_table <- run_cluster_anova(dat)
print(result_table)
write.xlsx(result_table, "VEC_ANOVA_results.xlsx", row.names = FALSE)

#5.transitional 
trans<-readRDS("../VIC_monocle/20250209/viec_trans_label.rds")
#ratio_plot(trans,sample.by="Sample",anno.by="trans",condition.by = "Group",save.prefix = "trans_test",plot_type="box",
#           condition.col=k6,strip.col = NULL,vjust.x = 0.5,hjust.x = 1,ord.condition = c("MV-normal",  "MV-moderate", "MV-severe",  "TV-normal",  "TV-moderate", "TV-severe" ))
dat  <- trans@meta.data
dat <- dat[,c("Sample","trans","Group","Group_seperate")]
colnames(dat) <- c("patient","majorCluster","loc","loc2")
dat$patient <- as.character(dat$patient)
dat$majorCluster <- as.character(dat$majorCluster)
ratio_plot(trans,sample.by="Sample",anno.by="trans",condition.by = "Group",save.prefix = "tukey/trans六组",plot_type="box",
           condition.col=k6,strip.col = NULL,vjust.x = 0.5,hjust.x = 1,compare.method="tukey",
           ord.condition = c("MV-normal",  "MV-moderate", "MV-severe",  "TV-normal",  "TV-moderate", "TV-severe" ))


dat$loc <- as.character(dat$loc)
Roe <- calTissueDist(dat,
                     byPatient = F,
                     colname.cluster = "majorCluster", # 不同细胞亚群
                     colname.patient = "patient", # 不同样本
                     colname.tissue = "loc", # 不同组织
                     method = "chisq", # "chisq", "fisher", and "freq" 
                     min.rowSum = 0) 
Roe
col_fun <- colorRamp2(c(min(Roe, na.rm = TRUE), 1, max(Roe, na.rm = TRUE)), c("blue", "white", "red"))
Roe<-Roe[,c("MV-normal","MV-moderate","MV-severe","TV-normal","TV-moderate","TV-severe")]
colnames(Roe)<-c("MV-Mild","MV-Moderate","MV-Severe","TV-Normal","TV-Moderate","TV-Severe")
p<-Heatmap(as.matrix(Roe),
           show_heatmap_legend = TRUE, 
           cluster_rows = FALSE, 
           cluster_columns = FALSE,
           row_names_side = 'right', 
           show_column_names = TRUE,
           show_row_names = TRUE,
           col = col_fun,
           row_names_gp = gpar(fontsize = 10),
           column_names_gp = gpar(fontsize = 10),
           heatmap_legend_param = list(
             title = "Ro/e Index",  # 自定义图注名称
             at = seq(0.5, 1.5, by = 0.5), # 例刻度的位置/自己的数据必须修改一下！
             labels = seq(0.5, 1.5, by = 0.5) # 每个刻度的标签/自己的数据必须修改一下！
           ),
           cell_fun = function(j, i, x, y, width, height, fill) {
             grid.text(sprintf("%.2f", Roe[i, j]), x, y, gp = gpar(fontsize = 10, col = "black"))
           }
)
pdf("trans_ratio_Roe.pdf",width=5,height=2)
print(p)
dev.off()

result_table <- run_cluster_anova(dat)
print(result_table)
write.xlsx(result_table, "trans_ANOVA_results.xlsx", row.names = FALSE)




#4.mc subtype
obj<-readRDS("../subRD/20250109/MC_subtype.rds")
head(obj@meta.data)
dat  <- obj@meta.data
dat <- dat[,c("Sample","subtype","Group","Group_seperate")]
colnames(dat) <- c("patient","majorCluster","loc","loc2")
dat$patient <- as.character(dat$patient)
dat$majorCluster <- as.character(dat$majorCluster)
ratio_plot(obj,sample.by="Sample",anno.by="subtype",condition.by = "Group",save.prefix = "tukey/MC_subtype六组",plot_type="box",
           condition.col=k6,strip.col = NULL,vjust.x = 0.5,hjust.x = 1,compare.method="tukey",
           ord.condition = c("MV-normal",  "MV-moderate", "MV-severe",  "TV-normal",  "TV-moderate", "TV-severe" ))


dat$loc <- as.character(dat$loc)
Roe <- calTissueDist(dat,
                     byPatient = F,
                     colname.cluster = "majorCluster", # 不同细胞亚群
                     colname.patient = "patient", # 不同样本
                     colname.tissue = "loc", # 不同组织
                     method = "chisq", # "chisq", "fisher", and "freq" 
                     min.rowSum = 0) 
Roe
col_fun <- colorRamp2(c(min(Roe, na.rm = TRUE), 1, max(Roe, na.rm = TRUE)), c("blue", "white", "red"))
Roe<-Roe[,c("MV-normal","MV-moderate","MV-severe","TV-normal","TV-moderate","TV-severe")]
colnames(Roe)<-c("MV-Mild","MV-Moderate","MV-Severe","TV-Normal","TV-Moderate","TV-Severe")
p<-Heatmap(as.matrix(Roe),
           show_heatmap_legend = TRUE, 
           cluster_rows = FALSE, 
           cluster_columns = FALSE,
           row_names_side = 'right', 
           show_column_names = TRUE,
           show_row_names = TRUE,
           col = col_fun,
           row_names_gp = gpar(fontsize = 10),
           column_names_gp = gpar(fontsize = 10),
           heatmap_legend_param = list(
             title = "Ro/e Index",  # 自定义图注名称
             at = seq(0.5, 1.5, by = 0.5), # 例刻度的位置/自己的数据必须修改一下！
             labels = seq(0.5, 1.5, by = 0.5) # 每个刻度的标签/自己的数据必须修改一下！
           ),
           cell_fun = function(j, i, x, y, width, height, fill) {
             grid.text(sprintf("%.2f", Roe[i, j]), x, y, gp = gpar(fontsize = 10, col = "black"))
           }
)
pdf("MC_subtype_ratio_Roe.pdf",width=5,height=4)
print(p)
dev.off()
result_table <- run_cluster_anova(dat)
print(result_table)
write.xlsx(result_table, "MC_ANOVA_results.xlsx", row.names = FALSE)

#6.lc subtype
obj<-readRDS("../subRD/20250109/LC_subtype.rds")
ratio_plot(obj,sample.by="Sample",anno.by="subtype",condition.by = "Group",save.prefix = "tukey/LC_subtype六组",plot_type="box",
           condition.col=k6,strip.col = NULL,vjust.x = 0.5,hjust.x = 1,compare.method="tukey",
           ord.condition = c("MV-normal",  "MV-moderate", "MV-severe",  "TV-normal",  "TV-moderate", "TV-severe" ))


head(obj@meta.data)
dat  <- obj@meta.data
dat <- dat[,c("Sample","subtype","Group","Group_seperate")]
colnames(dat) <- c("patient","majorCluster","loc","loc2")
dat$patient <- as.character(dat$patient)
dat$majorCluster <- as.character(dat$majorCluster)
dat$loc <- as.character(dat$loc)
Roe <- calTissueDist(dat,
                     byPatient = F,
                     colname.cluster = "majorCluster", # 不同细胞亚群
                     colname.patient = "patient", # 不同样本
                     colname.tissue = "loc", # 不同组织
                     method = "chisq", # "chisq", "fisher", and "freq" 
                     min.rowSum = 0) 
Roe
col_fun <- colorRamp2(c(min(Roe, na.rm = TRUE), 1, max(Roe, na.rm = TRUE)), c("blue", "white", "red"))
Roe<-Roe[,c("MV-normal","MV-moderate","MV-severe","TV-normal","TV-moderate","TV-severe")]
colnames(Roe)<-c("MV-Mild","MV-Moderate","MV-Severe","TV-Normal","TV-Moderate","TV-Severe")
p<-Heatmap(as.matrix(Roe),
           show_heatmap_legend = TRUE, 
           cluster_rows = FALSE, 
           cluster_columns = FALSE,
           row_names_side = 'right', 
           show_column_names = TRUE,
           show_row_names = TRUE,
           col = col_fun,
           row_names_gp = gpar(fontsize = 10),
           column_names_gp = gpar(fontsize = 10),
           heatmap_legend_param = list(
             title = "Ro/e Index",  # 自定义图注名称
             at = seq(0.5, 1.5, by = 0.5), # 例刻度的位置/自己的数据必须修改一下！
             labels = seq(0.5, 1.5, by = 0.5) # 每个刻度的标签/自己的数据必须修改一下！
           ),
           cell_fun = function(j, i, x, y, width, height, fill) {
             grid.text(sprintf("%.2f", Roe[i, j]), x, y, gp = gpar(fontsize = 10, col = "black"))
           }
)
pdf("LC_subtype_ratio_Roe.pdf",width=5,height=4)
print(p)
dev.off()
result_table <- run_cluster_anova(dat)
print(result_table)
write.xlsx(result_table, "LC_ANOVA_results.xlsx", row.names = FALSE)

#杨老师，还需要补个图，横坐标是6分组，
#纵坐标是HTR2B和HTR2A，比较不同分组VIC细胞中这俩基因的表达差异，
#需要进行ANOVA和事后Turkey检验

library(Seurat)
library(dplyr)
library(ggplot2)
library(reshape2)
library(ggpubr)

genes <- c("HTR2A", "HTR2B")

# 取 VIC 细胞
vic<-readRDS("/data2/Project/Fuwai_songGrp/7_Analysis/5.VIC_subtype.rds")

# 提取表达
expr <- FetchData(vic, vars = genes)

# 合并 metadata（Group_seperate）
expr$Group <- vic@meta.data$Group
expr$Cell <- rownames(expr)

# 转为长格式
df <- melt(expr, id.vars = c("Cell", "Group"),
           variable.name = "Gene",
           value.name = "Expression")
df$Group<-gsub("normal","mild",df$Group)
df$Group<-factor(as.character(df$Group),levels=c("MV-mild","MV-moderate","MV-severe",
                                                 "TV-mild","TV-moderate","TV-severe"))
p_violin <- ggplot(df, aes(x = Group, y = Expression, fill = Group)) +
  geom_violin(trim = TRUE, scale = "width", alpha = 0.9, width = 1.2) +  # ← 加宽
  facet_wrap(~Gene, ncol = 1, scales = "free_y") +
  theme_pubr() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.position = "none",
    strip.text = element_text(size = 14, face = "bold")
  ) +
  labs(
    x = "",
    y = "Expression",
    title = "HTR2A / HTR2B Expression in VIC (6 Groups)"
  )+scale_fill_manual(values=c("#50a3a4","#fcaf38","#b60e00","#9fb4c6","#ffdd1c","#f95335"))
                      
print(p_violin)

# 保存 PDF
ggsave(
  "HTR2A_HTR2B_6groups_violin.pdf",
  p_violin,
  width = 5,
  height = 6
)

p_box <- ggplot(df, aes(x = Group, y = Expression, fill = Group)) +
  geom_boxplot(outlier.shape = NA, width = 0.7, alpha = 0.9) +  # 不显示点
  facet_wrap(~Gene, ncol = 1, scales = "free_y") +
  theme_pubr() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.position = "none",
    strip.text = element_text(size = 14, face = "bold")
  ) +
  labs(
    x = "",
    y = "Expression",
    title = "HTR2A / HTR2B Expression in VIC (6 Groups)"
  )+scale_fill_manual(values=c("#50a3a4","#fcaf38","#b60e00","#9fb4c6","#ffdd1c","#f95335"))


print(p_box)

# 保存文件
ggsave(
  "HTR2A_HTR2B_6groups_boxplot.pdf",
  p_box,
  width = 5,
  height = 6
)
# 保存 PDF

anova_results <- df %>%
  group_by(Gene) %>%
  do(tidy(aov(Expression ~ Group, data = .)))

anova_results

library(broom)

tukey_results <- df %>%
  group_by(Gene) %>%
  do({
    fit <- aov(Expression ~ Group, data = .)
    tk <- TukeyHSD(fit)
    broom::tidy(tk)
  })

write.xlsx(list(anova=anova_results,tukey=tukey_results), "HTR2A_HTR2B_TukeyHSD.xlsx", row.names = FALSE)



