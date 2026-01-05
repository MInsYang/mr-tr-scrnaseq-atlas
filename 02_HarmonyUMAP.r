# ============================================================================
# 02 harmony UMAP celltype annotation
# ============================================================================
# 本文件包含02 harmony+umap+celltype annotation相关的代码
# ============================================================================

# 加载function
source("scripts/utils.R")


# Harmony去批次+降维+聚类+差异分析
sampleInfo<-read.xlsx("~/0_createInfo.xlsx",rowNames=TRUE)
dat<-readRDS("dat_celltype.rds")
dat<-Integration_SCP(dat,batch = "Sample",integration_method = "Harmony",
                       cluster_algorithm = "leiden",cluster_resolution = seq(0.3,1.5,by=0.1))
dat@meta.data[,c("Sample","Group","Group_overall","Group_seperate")]<-sampleInfo[dat$orig.ident,
                                                                                 c("Sample_new","Group","Group_overall","Group_seperate")]
saveRDS(dat,file="4.har.rds")
# 取resulotion=0.3的细胞分群结果
dat<-RunDEtest(dat,group_by="Harmony_SNN_res.0.3",only.pos = TRUE,fc.threshold = 1)
deg<-dat@tools$DEtest_Harmony_SNN_res.0.3$AllMarkers_wilcox
deg<-deg[order(deg$avg_log2FC,decreasing = TRUE),]
deg %>% split(deg$group1) %>% write.xlsx("4.har_clsDEG.xlsx")
CellStatPlot(dat,stat.by="Harmony_SNN_res.0.3",group.by="Sample",plot_type="dot")
# 根据marker基因做细胞类型注释
Idents(dat)<-dat$Harmony_SNN_res.0.3
dat<-RenameIdents(dat,"1"=	"VIC","2"="VIC","3"="VIC","5"="VIC","7"="VIC","10"="VIC","9"="VIC",
                  "8"="VEC",
                  "11"="Mast cell",
                  "4"="Myeloid cell",
                  "6"="Lymphocyte")
dat$celltype<-dat@active.ident
saveRDS(dat,file="dat_celltype.rds")

CellDimPlot(dat,group.by="celltype",label = TRUE,label_insitu = TRUE,palcolor=c1)
CellStatPlot(dat,stat.by="celltype",group.by=c("Sample","Group","Group_overall","Group_seperate"),
             plot_type="dot")
p1<-CellDimPlot(dat,group.by="celltype",palcolor = colall,raster=FALSE) 
p2<-list()
for(i in 1:length(unique(dat$tissue))){
  tp<-CellDimPlot(dat[,rownames(dat@meta.data)[dat$tissue==unique(dat$tissue)[i]]],group.by="celltype",split.by="tissue",palcolor = colall,raster=FALSE)
  p2[[i]]<-tp
}
#p2<-CellDimPlot(dat,group.by="celltype",split.by="tissue",raster=FALSE)
pdf("1_celltypeTissue.pdf",width=14,height=5)
cowplot::plot_grid(p1,p2[[1]],p2[[2]],nrow=1)
dev.off()
p1<-CellDimPlot(dat,group.by="Group_overall",raster=FALSE,palcolor=colgo) 
p2<-list()
for(i in 1:length(unique(dat$tissue))){
  tp<-CellDimPlot(dat[,rownames(dat@meta.data)[dat$tissue==unique(dat$tissue)[i]]],group.by="Group_overall",palcolor=colall,split.by="tissue",raster=FALSE)
  p2[[i]]<-tp
}
pdf("2_celltypeGroupOverall.pdf",width=14,height=5)
cowplot::plot_grid(p1,p2[[1]],p2[[2]],nrow=1)
dev.off()
p1<-CellDimPlot(dat,group.by="VR",raster=FALSE,palcolor=colall) 
p2<-list()
for(i in 1:length(unique(dat$tissue))){
  tp<-CellDimPlot(dat[,rownames(dat@meta.data)[dat$tissue==unique(dat$tissue)[i]]],group.by="VR",palcolor=colall,split.by="tissue",raster=FALSE)
  p2[[i]]<-tp
}
pdf("3_celltypeVR.pdf",width=14,height=5)
cowplot::plot_grid(p1,p2[[1]],p2[[2]],nrow=1)
dev.off()

# marker查看
mk<-data.frame(gene=c("COL1A1","COL3A1","DCN",
                      "PECAM1","CDH5","POSTN",
                      "KIT","TPSAB1","TPSB2",
                      "CD68","CD14","LYZ",
                      "CD3D","CD3E","NKG7"
                      ),
               cluster=rep(c("VIC","VEC","Mast cell","Myeloid cell","Lymphocyte"),c(3,3,3,3,3)))
dat$celltype_rev<-factor(as.character(dat$celltype),levels=rev(levels(dat$celltype)))
Idents(dat)<-dat$celltype_rev
p<-ggscploty(object = dat,features = mk$gene,gene.order=mk$gene,
         featuresAnno = mk$cluster,featuresAnno.order = levels(dat$celltype),
         reduction="HarmonyUMAP2D",pct.exp.var="celltypes_rev",
         mapping = aes(x = gene_name,y = celltype_rev, fill = mean_exp,exp = value),pct.order = levels(dat$celltype_rev))+ geom_scViolin() +
  facet_new(facet_col = "featureAnno",scales = "free_x",
            strip.col = colall,
            x.angle = 90,x.label.hjust = 1)
p

pdf("4_marker_vlnplot.pdf",width=10,height=4)
p+labs(y="celltype")
dev.off()
p1<-GroupHeatmapy(dat,features=as.character(mk$gene),group_palcolor = list(colall),
                  group.by="celltype",nlabel = 0,cluster_columns = FALSE,show_row_names = TRUE,add_bg=FALSE,add_dot=TRUE,show_column_names = TRUE,
                  flip=TRUE,heatmap_palette = "YlGn",dot_size=unit(6,"mm"),exp_legend_title = "zscore",row_title = "",column_title = "",cluster_rows = FALSE)
pdf("5_marker_dotplot.pdf",width=10,height=3)
p1$plot
dev.off()

df<-data.frame(table(dat$celltype))
colnames(df)<-c("celltype","count")
df <- df %>%
  mutate(proportion = count / sum(count) * 100) %>%
  arrange(desc(count))
p<-ggplot(df, aes(x = reorder(celltype, count), y = proportion, fill = celltype)) +
  geom_bar(stat = "identity", width = 0.7) +
  coord_flip() +
  theme_classic() +
  scale_fill_manual(values=colall) +
  labs(x = "celltype", y = "nCellsStat", fill = "Cell Type") +
  geom_text(aes(label = paste0(count, "\n(", round(proportion, 1), "%)")),
            hjust = -0.01, size = 3) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+theme(axis.text=element_text(color="black",size=12))
pdf("6_celltypeStat.pdf",width=5,height=4)
print(p+ylim(c(0,100)))
dev.off()
# 细胞比例检测及画图
ratio_plot(dat,sample.by = "Sample",anno.by = "celltype",condition.by = "Group_overall",save.prefix = "7_GroupOverall",ord.condition = c("Control","Case"),
           plot_type="box",strip.col = colall,facet.ncol=3,width=6,height=6,hjust=1,vjust.x=0.5)
ratio_plot(dat,sample.by = "Sample",anno.by = "celltype",condition.by = "Group_seperate",save.prefix = "8_GroupSeperate",
           ord.condition = c("MV-Control","MV-Case","TV-Control","TV-Case"),
           plot_type="box",strip.col = colall,facet.ncol=3,width=6,height=6,hjust=1,vjust.x=0.5)
meta<-dat@meta.data
meta$celltype_plt<-factor(as.character(meta$celltype),levels=rev(levels(dat$celltype)))
p1<-CellStatPlot(dat,group.by="celltype",stat.by="Group_overall",flip=TRUE,legend.position = "top",palcolor = colall)
p2<-CellStatPlot(dat,group.by="celltype",stat.by="Group_seperate",flip=TRUE,legend.position = "top",palcolor=colall)
p3<-CellStatPlot(dat,stat.by="celltype",group.by="celltype",stat_type="count",flip=TRUE,legend.position = "top",palcolor=colall)
p3
p4<-ggplot(meta,aes(y=celltype_plt,x=nCount_RNA))+geom_boxplot(fill="skyblue")+
  theme_scp()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+labs(y="celltype")
p4
library(patchwork)
pdf("layout1.pdf",width=12,height=5)
(p1 | p2 | p3 | p4) + plot_layout(widths = c(1, 1, 1, 1))
dev.off()

p1<-CellStatPlot(dat,group.by="celltype",stat.by="VR",flip=TRUE,legend.position = "top",palcolor=colall)
p2<-CellStatPlot(dat,group.by="celltype",stat.by="Group",flip=TRUE,legend.position = "top",palcolor=colall)
p3<-CellStatPlot(dat,stat.by="celltype",group.by="celltype",stat_type="count",flip=TRUE,legend.position = "top",palcolor=colall)
p3
p4<-ggplot(meta,aes(y=celltype_plt,x=nCount_RNA))+geom_boxplot(fill="skyblue")+
  theme_scp()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+labs(y="celltype")
p4
library(patchwork)
pdf("layout2.pdf",width=12,height=5)
(p1 | p2 | p3 | p4) + plot_layout(widths = c(1, 1, 1, 1))
dev.off()

ratio_plot(dat,sample.by = "Sample",anno.by = "celltype",condition.by = "VR",save.prefix = "9_VR",ord.condition = c("normal","moderate","severe"),
           plot_type="box",strip.col = colall,facet.ncol=3,width=6,height=6,hjust=1,vjust.x=0.5)
ratio_plot(dat,sample.by = "Sample",anno.by = "celltype",condition.by = "Group",save.prefix = "10_Group",
           ord.condition = c("MV-normal","MV-moderate","MV-severe","TV-normal","TV-moderate","TV-severe"),
           plot_type="box",strip.col = colall,facet.ncol=3,width=8,height=6,condition.col =colall,hjust=1)

col7<-c("#6cb86a","#CE1E5D","#7AC177","#C5ACB5","#EFBC87","#F9E199","#8F5A32","#A1B0AB","#42519F")
col1<-col7[c(1,4,2,5,6)]
colgp<-c("#7AC177","#C5ACB5","#CE1E5D","#EFBC87")
names(colgp)<-c("MV-Case","TV-Control","MV-Control","TV-Case")
col72<-c("#AEC4BF","#DCE5F0","#BD130F","#D4553E","#42519F")
colgo<-c("#42519F","#EFBC87")
colvr<-c("#7AC177","#F9E199","#CE1E5D")
#b13<-c("#9fb4c6","#066781","#50a3a4","#fcaf38","#ffdd1c","#c08520","#f4dfca","#b60e00","#f95335","#674a40","#c59e90","#485176","#7d807a")
b13<-c("#9fb4c6","#ffdd1c","#f95335","#50a3a4","#fcaf38","#b60e00")
names(b13)<-c( "TV-normal" ,  "TV-moderate" ,  "TV-severe" , "MV-normal" ,"MV-moderate",  "MV-severe")
colall<-c(col1,colgo,colvr)
names(colall)<-c("VIC","VEC","Mast cell","Myeloid cell","Lymphocyte","Case","Control","normal","moderate","severe")
colall<-c(colall,b13)
p1<-CellDimPlot(dat,group.by="Group_overall",raster=FALSE,palcolor=col_groupoverall) 
p2<-CellDimPlot(dat,group.by="Group_overall",raster=FALSE,palcolor=c("#42519F","#EFBC87")) 
p3<-CellDimPlot(dat,group.by="VR",raster=FALSE,palcolor=c("#7AC177","#F9E199","#CE1E5D")) 
p1+p2+p3
