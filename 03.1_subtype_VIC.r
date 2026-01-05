# ============================================================================
# 03 细胞亚群分析--- VIC
# ============================================================================
# 本文件包含 VIC 亚群分析相关的代码
# ============================================================================

# 加载function
source("scripts/utils.R")

# VIC亚群分析
vic<-subset(dat,celltype=="VIC")
har<-Integration_SCP(vic,batch = "Sample",integration_method = "Harmony",do_HVF_finding = TRUE,force_linear_reduction = TRUE,force_nonlinear_reduction = TRUE,
                     cluster_algorithm = "leiden",neighbor_k = 30,cluster_resolution = seq(0.3,1.5,by=0.1))
saveRDS(har,file="VIC_har.rds")
dat<-readRDS("VIC_har.rds")
clustree(dat,prefix="Harmony_SNN_res.")
CellDimPlot(dat,group.by="Harmony_SNN_res.0.3",label=TRUE,label_insitu = TRUE,add_mark = TRUE, mark_linetype = 2)
register(MulticoreParam(workers = 5, progressbar = TRUE))
har<-RunDEtest(har,group_by="Harmony_SNN_res.0.3",only.pos=TRUE,fc.threshold = 1.5)
saveRDS(har,file="VIC_har.rds")
deg <- har@tools$DEtest_Harmony_SNN_res.0.3$AllMarkers_wilcox
deg <-deg[order(deg$avg_log2FC,decreasing = TRUE),]
deg %>% split(deg$group1) %>% write.xlsx("VIC_res0.3_fc1.5.xlsx")
deg<-deg[deg$p_val<0.05,]
top2<- deg %>% group_by(group1) %>% top_n(2,avg_log2FC)
top2<-sortMTX(top2,by = "group1",order = seq(1,9))
for( i in unique(top2$group1)){
   p<-FeatureDimPlot(har,features = top2$gene[top2$group1==i],raster=FALSE)
   pdf(paste0("vic2_cls",i,".pdf"),width=10,height=5)
   print(p)
   dev.off()
}
har$Harmony_SNN_res.0.3<-factor(as.character(har$Harmony_SNN_res.0.3),levels=seq(1:9))
p<-GroupHeatmapy(har,features=top2$gene,group.by="Harmony_SNN_res.0.3",column_title = "",row_title="",nlabel = 0,show_row_names = TRUE,
              show_column_names = TRUE,exp_legend_title = "zscore",add_dot = TRUE,add_bg = TRUE)
#p$plot
pdf("vic3.pdf",width=6,height=8)
p$plot
dev.off()
top10<- deg %>% group_by(group1) %>% top_n(10,avg_log2FC)
top10<-sortMTX(top10,by = "group1",order = seq(1,9))
p<-FeatureHeatmap(har,features=top10$gene,max_cells = 10000,group.by="Harmony_SNN_res.0.3",column_title = "",row_title="",
                  nlabel = 30,show_row_names = FALSE,cluster_rows = FALSE,cluster_columns = TRUE,
                  show_column_names = FALSE,exp_legend_title = "zscore",use_raster = FALSE)
#p$plot
Cairo::CairoTIFF("vic4.tiff",width=1000,height=1000)
p$plot
dev.off()
ratio_plot(har,sample.by = "Sample",anno.by = "Harmony_SNN_res.0.3",condition.by = "Group",save.prefix = "vic7_Group",
           ord.condition = c("MV-normal","MV-moderate","MV-severe","TV-normal","TV-moderate","TV-severe"),
           plot_type="box",strip.col = colall,facet.ncol=3,width=8,height=8,hjust=1,vjust.x=0.5,condition.col =as.vector(colxk) )
ratio_plot(har,sample.by = "Sample",anno.by = "Harmony_SNN_res.0.3",condition.by = "Group_seperate",save.prefix = "vic7_GroupSeperate",
           ord.condition = c("MV-Control","MV-Case","TV-Control","TV-Case"),condition.col = as.vector(colxk),
           plot_type="box",strip.col = colall,facet.ncol=3,width=8,height=8,hjust=1,vjust.x=0.5)
ratio_plot(har,sample.by = "Sample",anno.by = "Harmony_SNN_res.0.3",condition.by = "Group",save.prefix = "vic7_Group",
           ord.condition = c("MV-normal","MV-moderate","MV-severe","TV-normal","TV-moderate","TV-severe"),compare.method = "wilcox.test",
           plot_type="box",strip.col = colall,facet.ncol=3,width=8,height=8,hjust=1,vjust.x=0.5,condition.col =as.vector(colxk) )
ratio_plot(har,sample.by = "Sample",anno.by = "Harmony_SNN_res.0.3",condition.by = "Group_seperate",save.prefix = "vic7_GroupSeperate",
           ord.condition = c("MV-Control","MV-Case","TV-Control","TV-Case"),condition.col = as.vector(colxk),compare.method = "wilcox.test",
           plot_type="box",strip.col = colall,facet.ncol=3,width=8,height=8,hjust=1,vjust.x=0.5)
meta<-har@meta.data
meta$celltype_plt<-factor(as.character(meta$Harmony_SNN_res.0.3),levels=rev(levels(har$Harmony_SNN_res.0.3)))
p1<-CellStatPlot(har,group.by="Harmony_SNN_res.0.3",stat.by="Group_overall",flip=TRUE,legend.position = "top",palcolor = colall)
p2<-CellStatPlot(har,group.by="Harmony_SNN_res.0.3",stat.by="Group_seperate",flip=TRUE,legend.position = "top",palcolor=colall)
p4<-ggplot(meta,aes(y=celltype_plt,x=nCount_RNA))+geom_boxplot(fill="skyblue")+
  theme_scp()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+labs(y="celltype")
p4
library(patchwork)
pdf("vic6_1.pdf",width=8,height=5)
(p1 | p2 |  p4) + plot_layout(widths = c(1, 1, 1))
dev.off()
har$Group<-gsub(":","-",har$Group)
har$Group<-gsub("moderete","moderate",har$Group)
har$Group<-gsub("Mit","MV",har$Group)
har$Group<-gsub("Tri","TV",har$Group)
p1<-CellStatPlot(har,group.by="Harmony_SNN_res.0.3",stat.by="VR",flip=TRUE,legend.position = "top",palcolor = colall)
p2<-CellStatPlot(har,group.by="Harmony_SNN_res.0.3",stat.by="Group",flip=TRUE,legend.position = "top",palcolor=colall)
p4<-ggplot(meta,aes(y=celltype_plt,x=nCount_RNA))+geom_boxplot(fill="skyblue")+
  theme_scp()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+labs(y="celltype")
p4
library(patchwork)
pdf("vic6_2.pdf",width=8,height=5)
(p1 | p2 |  p4) + plot_layout(widths = c(1, 1, 1))
dev.off()
har<-RunEnrichment(har,group_by="Harmony_SNN_res.0.3",DE_threshold = "avg_log2FC>0 & p_val<0.05",db = "GO_BP",species = "Homo_sapiens")
mtx<-har@tools$Enrichment_Harmony_SNN_res.0.3_wilcox$enrichment
mtx %>% split(mtx$Groups) %>% write.xlsx("vic5_enrichGOBP.xlsx")
p<-EnrichmentPlot(har,db="GO_BP",group_by="Harmony_SNN_res.0.3",plot_type="bar",ncol=3)
pdf("vic5_bar.pdf",width=20,height=14)
p
dev.off()
# VIC亚群合并及注释
Idents(har)<-har$Harmony_SNN_res.0.5
har<-RenameIdents(har,"4"="1",
                  "11"="3","12"="3")
har$subtype<-paste0("VIC",har@active.ident)
Idents(har)<-har$subtype
har<-RenameIdents(har,"VIC1"="VIC1","VIC2"="VIC2","VIC3"="VIC3","VIC5"="VIC4","VIC6"="VIC5","VIC7"="VIC6","VIC8"="VIC7","VIC9"="VIC8","VIC10"="VIC9")
har$subtype<-har@active.ident
saveRDS(har,file="vic_subtype.rds")

