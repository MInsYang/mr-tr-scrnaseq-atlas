# ============================================================================
# fibrotic / transitional 判断及分析
# ============================================================================
# 本文件包含fibrotic定义及transitional VIC、VEC定义的相关的代码
# ============================================================================

# 加载function
source("scripts/utils.R")

# fibrotic 差异基因计算
for(i in unique(vic$Group)){
  tmp<-subset(vic,Group==i)
  deg<-FindMarkers(tmp,group.by = "fibrotic",ident.1="Anti-fibrotic VIC",ident.2="Pro-fibrotic VIC",only.pos = FALSE,
                   logfc.threshold = -Inf)
  deg<-deg[order(deg$avg_log2FC,decreasing=TRUE),]
  write.xlsx(deg,file=paste0(i,"_fibrotic_AntivsPro_deg.xlsx"),rowNames=TRUE)
}
for(i in unique(vic$Group_seperate)){
  tmp<-subset(vic,Group_seperate==i)
  deg<-FindMarkers(tmp,group.by = "fibrotic",ident.1="Anti-fibrotic VIC",ident.2="Pro-fibrotic VIC",only.pos = FALSE,
                   logfc.threshold = -Inf)
  deg<-deg[order(deg$avg_log2FC,decreasing=TRUE),]
  write.xlsx(deg,file=paste0(i,"_fibrotic_AntivsPro_deg.xlsx"),rowNames=TRUE)
}
for(i in unique(vic$Group_overall)){
  tmp<-subset(vic,Group_overall==i)
  deg<-FindMarkers(tmp,group.by = "fibrotic",ident.1="Anti-fibrotic VIC",ident.2="Pro-fibrotic VIC",only.pos = FALSE,
                   logfc.threshold = -Inf)
  deg<-deg[order(deg$avg_log2FC,decreasing=TRUE),]
  write.xlsx(deg,file=paste0(i,"_fibrotic_AntivsPro_deg.xlsx"),rowNames=TRUE)
}

# transitional 计算
mdat<-readRDS("VIEC.rds")
deg<-mergeXLSX2("celltype_deg.xlsx")
deg<-deg[deg$p_val_adj<0.01,]
top20<-deg %>% group_by(group1) %>% top_n(50,avg_log2FC) 
mdat<-AddModuleScore(mdat,features=list(top20$gene[top20$group1=="VIC"]),name = "VIC_top20")
mdat<-AddModuleScore(mdat,features=list(top20$gene[top20$group1=="VEC"]),name="VEC_top20")
meta<-mdat@meta.data[,c("VIC_top201","VEC_top201","subtype")]
p1<-ggplot(meta,aes(x=VIC_top201,y=VEC_top201,color=subtype))+geom_point()+theme_scp()+labs(x="VIC_top20 score",y="VEC_top20 score")
top50<-deg %>% group_by(group1) %>% top_n(50,avg_log2FC)
mdat<-AddModuleScore(mdat,features=list(top50$gene[top50$group1=="VIC"]),name = "VIC_top50")
mdat<-AddModuleScore(mdat,features=list(top50$gene[top50$group1=="VEC"]),name="VEC_top50")
meta<-mdat@meta.data[,c("VIC_top501","VEC_top501","subtype","celltype")]
p2<-ggplot(meta,aes(x=VIC_top501,y=VEC_top501,color=subtype))+geom_point()+theme_scp()+labs(x="VIC_top50 score",y="VEC_top50 score")
p2
top100<-deg %>% group_by(group1) %>% top_n(100,avg_log2FC)
mdat<-AddModuleScore(mdat,features=list(top100$gene[top100$group1=="VIC"]),name = "VIC_top100")
mdat<-AddModuleScore(mdat,features=list(top100$gene[top100$group1=="VEC"]),name="VEC_top100")
meta<-mdat@meta.data[,c("VIC_top1001","VEC_top1001","subtype")]
p3<-ggplot(meta,aes(x=VIC_top1001,y=VEC_top1001,color=subtype))+geom_point()+theme_scp()+labs(x="VIC_top100 score",y="VEC_top100 score")
# 测试不同top gene数量的效果
p1+p2+p3

setwd("/data2/Project/F_Group_Analysis/7_Analysis/score")
meta<-readRDS("VIEC_top50_score.rds")
p2<-ggplot(meta,aes(x=VIC_top501,y=VEC_top501,color=subtype))+geom_point()+theme_scp()+labs(x="VIC_top50 score",y="VEC_top50 score")+
  scale_color_manual(values=as.character(SCP::palette_scp(n=length(unique(meta$subtype)))))
ggExtra::ggMarginal(
  p2,
  type = "density",        # 边际密度图
  groupColour = TRUE,      # 按组颜色绘制
  groupFill = TRUE         # 按组填充绘制
)
group1<-meta$VIC_top501[meta$celltype=="VIC"]
group2<-meta$VIC_top501[meta$celltype=="VEC"]
density1 <- density(group1)
density2 <- density(group2)
# 找到密度曲线的共同 x 值范围
common_x <- seq(max(min(density1$x), min(density2$x)), min(max(density1$x), max(density2$x)), length.out = 1000)
# 在共同范围内计算两条曲线的 y 值
y1 <- approx(density1$x, density1$y, xout = common_x)$y
y2 <- approx(density2$x, density2$y, xout = common_x)$y
# 找到交点
intersection_idx <- which(diff(sign(y1 - y2)) != 0)
intersection_points <- common_x[intersection_idx]
intersection_points
idx <- which(density1$x < 0.3551256) #vic
idx2 <- which(density2$x > 0.3551256) #vec
vic_df<-meta[meta$celltype=="VIC",]
vec_df<-meta[meta$celltype=="VEC",]
table(vic_df$VIC_top501<0.3551256)
table(vec_df$VIC_top501>0.3551256)
trans_vec<-rownames(vec_df)[vec_df$VIC_top501>0.3551256]


# 对每组数据计算均值和标准误差
mean1 <- mean(group1)
se1 <- sd(group1) / sqrt(length(group1))
mean2 <- mean(group2)
se2 <- sd(group2) / sqrt(length(group2))
# 计算 95% 置信区间
ci1 <- c(mean1 - qt(0.975, df = length(group1) - 1) * se1, mean1 + qt(0.975, df = length(group1) - 1) * se1)
ci2 <- c(mean2 - qt(0.975, df = length(group2) - 1) * se2, mean2 + qt(0.975, df = length(group2) - 1) * se2)
list(ci1 = ci1, ci2 = ci2)
library(ggplot2)
#画图数据
df <- data.frame(
  x = c(density1$x, density2$x),
  y = c(density1$y, density2$y),
  group = rep(c("VIC", "VEC"), each = length(density1$x))
)
# 绘制密度图
p<-ggplot(df, aes(x = x, y = y, color = group)) +
  geom_line(size = 1) +
  geom_vline(xintercept = intersection_points, linetype = "dashed", color = "red") +
  #geom_vline(xintercept = ci1, linetype = "dotted", color = "blue") +
  #geom_vline(xintercept = ci2, linetype = "dotted", color = "green") +
  theme_minimal() +
  labs(title = "Density Curves with Intersection and Confidence Intervals",
       x = "VIC_score",
       y = "Density",
       color = "Group")
p+theme_scp()
pdf("VIC_score_density.pdf",width=6,height=4)
print(p+theme_scp())
dev.off()

mean_function <- function(data, indices) {
  sample_data <- data[indices]  # 根据bootstrap索引抽样
  return(mean(sample_data))
}
library(boot)
# 对每组数据进行 bootstrap 重抽样
boot1 <- boot(group1, mean_function, R = 1000)  # R 是重抽样次数
boot2 <- boot(group2, mean_function, R = 1000)
# 计算置信区间（95%）
ci1 <- boot.ci(boot1, type = "perc")$percent[4:5]  # 从 bootstrap 结果中提取百分位置信区间
ci2 <- boot.ci(boot2, type = "perc")$percent[4:5]
# 打印置信区间
list(ci1 = ci1, ci2 = ci2)


group1<-meta$VEC_top501[meta$celltype=="VIC"]
group2<-meta$VEC_top501[meta$celltype=="VEC"]
density1 <- density(group1)
density2 <- density(group2)
# 找到密度曲线的共同 x 值范围
common_x <- seq(max(min(density1$x), min(density2$x)), min(max(density1$x), max(density2$x)), length.out = 1000)
# 在共同范围内计算两条曲线的 y 值
y1 <- approx(density1$x, density1$y, xout = common_x)$y
y2 <- approx(density2$x, density2$y, xout = common_x)$y
# 找到交点
intersection_idx <- which(diff(sign(y1 - y2)) != 0)
intersection_points <- common_x[intersection_idx]
intersection_points
# 对每组数据计算均值和标准误差
mean1 <- mean(group1)
se1 <- sd(group1) / sqrt(length(group1))
mean2 <- mean(group2)
se2 <- sd(group2) / sqrt(length(group2))
# 计算 95% 置信区间
ci1 <- c(mean1 - qt(0.975, df = length(group1) - 1) * se1, mean1 + qt(0.975, df = length(group1) - 1) * se1)
ci2 <- c(mean2 - qt(0.975, df = length(group2) - 1) * se2, mean2 + qt(0.975, df = length(group2) - 1) * se2)
list(ci1 = ci1, ci2 = ci2)
library(ggplot2)
#画图数据
df <- data.frame(
  x = c(density1$x, density2$x),
  y = c(density1$y, density2$y),
  group = rep(c("VIC", "VEC"), each = length(density1$x))
)
# 绘制密度图
ggplot(df, aes(x = x, y = y, color = group)) +
  geom_line(size = 1) +
  geom_vline(xintercept = intersection_points, linetype = "dashed", color = "red") +
  geom_vline(xintercept = ci1, linetype = "dotted", color = "blue") +
  geom_vline(xintercept = ci2, linetype = "dotted", color = "green") +
  theme_minimal() +
  labs(title = "Density Curves with Intersection and Confidence Intervals",
       x = "VEC score",
       y = "Density",
       color = "Group")
intersection_points
idx <- which(density1$x < 0.1838277) #vic
idx2 <- which(density2$x > 0.1838277) #vec
vic_df<-meta[meta$celltype=="VIC",]
vec_df<-meta[meta$celltype=="VEC",]
table(vic_df$VEC_top501<0.1838277)
table(vec_df$VEC_top501>0.1838277)
trans_vic<-rownames(vic_df)[vic_df$VEC_top501>0.1838277]
saveRDS(trans_vic,file="trans_vic_cid.rds")
saveRDS(trans_vec,file="trans_vec_cid.rds")

meta$trans<-ifelse(rownames(meta) %in% c(trans_vic,trans_vec),"transitional","conventional")
meta$label<-paste0(meta$trans,"_",meta$celltype)
ggplot(meta,aes(x=VIC_top501,y=VEC_top501))+geom_point(aes(color=label))+theme_scp()
df<-table(meta$subtype,meta$label)
fq<-prop.table(df,1)*100
df<-melt(fq)
df<-df[grepl("trans",df$Var2),]
write.xlsx(fq,file="trans_prop.xlsx",rowNames=TRUE)


q4<-c("#8dd2c5","#7fb2d5","#f47f72","#bfbcda")
names(q4)<-c("PECAM1-COL1A1+", "PECAM1-COL1A1-", "PECAM1+COL1A1+" , "PECAM1+COL1A1-")
vic<-readRDS("../5.VIC_subtype.rds")
trans_vic<-readRDS("../20250102/trans_vic_cid.rds")
vic$trans<-"conventional"
vic@meta.data[trans_vic,"trans"]<-"transitional"
vic@meta.data$PECAM1<-vic@assays$RNA@data["PECAM1",rownames(vic@meta.data)]
vic$PECAM1_group<-ifelse(vic$PECAM1>0,"PECAM1+","PECAM1-")
vic@meta.data$COL1A1<-vic@assays$RNA@data["COL1A1",rownames(vic@meta.data)]
vic$COL1A1_group<-ifelse(vic$COL1A1>0,"COL1A1+","COL1A1-")
vic$gene_label<-paste0(vic$PECAM1_group,vic$COL1A1_group)
vic$gene_label<-factor(as.character(vic$gene_label),levels=c("PECAM1-COL1A1+", "PECAM1-COL1A1-", "PECAM1+COL1A1+" , "PECAM1+COL1A1-"))
df<-table(vic$gene_label,vic$trans)
fq<-prop.table(table(vic$gene_label,vic$trans),2)*100
sheets<-list(count=df,ratio=fq)
write.xlsx(sheets,file="vic_transitional_genelabel_ratio.xlsx",rowNames=TRUE)
CellDimPlot(vic,group.by="trans",cells.highlight = rownames(vic@meta.data)[vic$gene_label=="PECAM1+COL1A1+"],raster=FALSE)
tmp<-subset(vic,trans=="transitional")
p<-CellDimPlot(tmp,group.by="gene_label",cells.highlight = rownames(tmp@meta.data)[tmp$gene_label=="PECAM1+COL1A1+"],
            palcolor=q4,raster=FALSE)
pdf("vic_transitional_genelabel.pdf",width=5,height=4)
print(p)
dev.off()
p<-CellDimPlot(vic,group.by="gene_label",cells.highlight = rownames(vic@meta.data)[vic$gene_label=="PECAM1+COL1A1+"],
            palcolor = q4,raster=FALSE)
pdf("vic_all_genelabel.pdf",width=7,height=5)
print(p)
dev.off()

vec<-readRDS("../5.VEC_subtype.rds")
trans_vec<-readRDS("../20250102/trans_vec_cid.rds")
vec$trans<-"conventional"
vec@meta.data[trans_vec,"trans"]<-"transitional"
vec@meta.data$PECAM1<-vec@assays$RNA@data["PECAM1",rownames(vec@meta.data)]
vec$PECAM1_group<-ifelse(vec$PECAM1>0,"PECAM1+","PECAM1-")
vec@meta.data$COL1A1<-vec@assays$RNA@data["COL1A1",rownames(vec@meta.data)]
vec$COL1A1_group<-ifelse(vec$COL1A1>0,"COL1A1+","COL1A1-")
vec$gene_label<-paste0(vec$PECAM1_group,vec$COL1A1_group)
vec$gene_label<-factor(as.character(vec$gene_label),levels=c("PECAM1-COL1A1+", "PECAM1-COL1A1-", "PECAM1+COL1A1+" , "PECAM1+COL1A1-"))
df<-table(vec$gene_label,vec$trans)
fq<-prop.table(table(vec$gene_label,vec$trans),2)*100
sheets<-list(count=df,ratio=fq)
write.xlsx(sheets,file="vec_transitional_genelabel_ratio.xlsx",rowNames=TRUE)
CellDimPlot(vec,group.by="trans",cells.highlight = rownames(vec@meta.data)[vec$gene_label=="PECAM1+COL1A1+"],raster=FALSE)
tmp<-subset(vec,trans=="transitional")
p<-CellDimPlot(tmp,group.by="gene_label",cells.highlight = rownames(tmp@meta.data)[tmp$gene_label=="PECAM1+COL1A1+"],
               palcolor=list(q4),raster=FALSE)
pdf("vec_transitional_genelabel.pdf",width=5,height=4)
print(p)
dev.off()
p<-CellDimPlot(vec,group.by="gene_label",cells.highlight = rownames(vec@meta.data)[vec$gene_label=="PECAM1+COL1A1+"],
               palcolor = q4,raster=FALSE)
pdf("vec_all_genelabel.pdf",width=7,height=5)
print(p)
dev.off()
tmp<-merge(vic,vec)
tmp<-subset(tmp,trans %in% c("transitional_VIC","transitional_VEC"))
saveRDS(tmp,file="transitional_viec.rds")

meta<-readRDS("../score/VIEC_top50_score.rds")
meta$trans<-"conventional"
meta$trans<-paste0(meta$trans,"_",meta$celltype)
meta[rownames(tmp@meta.data),"trans"]<-tmp$trans
meta$gene_label<-tmp@meta.data[rownames(meta),"gene_label"]
saveRDS(meta,file="viec_score_genelable.rds")
meta<-readRDS("viec_score_genelable.rds")
p1<-ggplot(meta,aes(x=VIC_top501,y=VEC_top501))+geom_point(aes(color=trans),size = 1.5, alpha = 0.7)+
  geom_point(
  data = meta[meta$gene_label == "PECAM1+COL1A1+", ],
  aes(x = VIC_top501, y = VEC_top501,fill=gene_label,color=trans),
  size = 1.5, alpha=0.3, shape = 21)+
  geom_point(
    data = meta[meta$gene_label == "PECAM1+COL1A1+" & meta$trans %in% c("transitional_VEC","transitional_VIC"), ],
    aes(x = VIC_top501, y = VEC_top501,fill=gene_label),
    size = 1.5, color="black",alpha=0.3, shape = 21)+labs(x = "VIC score", y = "VEC score", color = "Transition Type") +theme_scp()+
  scale_color_manual(values=q4)+scale_fill_manual(values="#ffdd14")
p1
smeta<-meta[meta$trans %in% c("transitional_VEC","transitional_VIC"),]
p2<-ggplot(smeta,aes(x=VIC_top501,y=VEC_top501))+geom_point(aes(color=trans),size = 3, alpha = 0.7)+geom_point(
  data = smeta[smeta$gene_label == "PECAM1+COL1A1+", ],
  aes(x = VIC_top501, y = VEC_top501,fill=gene_label),
  size = 3.5, color = "black", shape = 21) +labs(x = "VIC score", y = "VEC score", color = "Transition Type") +theme_scp()+
  scale_color_manual(values=q4[c(3,4)])+scale_fill_manual(values="#ffdd14")
p2
pdf("score_dot.pdf",width=7,height=5)
print(p1)
print(p2)
dev.off()

fq1<-read.xlsx("vic_transitional_genelabel_ratio.xlsx",sheet=2)
df<-melt(fq1)
colnames(df)[1]<-"gene_label"
df<-df[df$gene_label=="PECAM1+COL1A1+",]
p3<-ggplot(df,aes(x=variable,y=value))+geom_bar(stat="identity",fill="#ffdd14",color="black")+
  theme_scp()+theme(axis.text.x=element_text(angle=45,vjust=1,hjust=1))+labs(x="",y="(PECAM1+COL1A1+) % in VIC")
p3
fq1<-read.xlsx("vec_transitional_genelabel_ratio.xlsx",sheet=2)
df<-melt(fq1)
colnames(df)[1]<-"gene_label"
df<-df[df$gene_label=="PECAM1+COL1A1+",]
p4<-ggplot(df,aes(x=variable,y=value))+geom_bar(stat="identity",fill="#ffdd14",color="black")+
  theme_scp()+theme(axis.text.x=element_text(angle=45,vjust=1,hjust=1))+labs(x="",y="(PECAM1+COL1A1+) % in VEC")
p4
pdf("gene_label_ratio.pdf",width=5,height=5)
print(p3+p4)
dev.off()