# ============================================================================
# monocle相关分析
# ============================================================================
# 本文件包含01 data preprocessing相关的代码
# ============================================================================

# 加载function
source("scripts/utils.R")


mcds<-readRDS("/data2/Project/F_Group_Analysis/7_Analysis/VIC_monocle/20250209/p7_viec_trans_reverse.rds")
BEAM_res=BEAM(mcds,branch_point = 2,cores = 1,progenitor_method="duplicate")
saveRDS(BEAM_res, file = "p7_viec_trans_bp2_BEAM_res.rds")
.libPaths(c( "/home/yangmy/R/x86_64-pc-linux-gnu-library/4.3", "/opt/R/4.3.2/lib/R/library"  ,
             "/data2/usr/yangmy_conda/r_env/lib/R/library" ,
             "/home/yangmy/.local/share/r-miniconda/envs/decontX_env/lib/R/library"))
library(monocle)
setwd("/data2/Project/F_Group_Analysis/7_Analysis/VIC_monocle")
mcd<-readRDS("~/monocle_Rds/p7_vic_TV2.Rds")
plot_cell_trajectory(mcd,color_by="State")
BEAM_res=BEAM(mcd,branch_point = 1,cores = 1,progenitor_method="duplicate")
saveRDS(BEAM_res, file = "p7_vic_TV2_bp1_BEAM_res.rds")

dir.create("20250209")
setwd("20250209")
mcds<-readRDS("../vic_MV_AntiFibrotic.Rds")
beam<-BEAM(mcds,branch_point =4,progenitor_method="duplicate")
saveRDS(beam,file="vic_MV_AntiFibrotic_bp4_beam.rds")

beam<-BEAM(mcds,branch_point = 2,progenitor_method="duplicate")
saveRDS(beam,file="vic_MV_AntiFibrotic_bp2_beam.rds")
beam<-BEAM(mcds,branch_point =5,progenitor_method="duplicate")
saveRDS(beam,file="vic_MV_AntiFibrotic_bp5_beam.rds")

pData(mcds)$subtype<-factor(as.character(pData(mcds)$subtype),levels=c("VIC1","VIC2","VIC4","VIC3","VIC5","VIC6"))
p1<-plot_cell_trajectory(mcds,color_by="fibrotic",show_branch_points = FALSE)+facet_wrap(~fibrotic,nrow=1)+
  scale_color_manual(values=c("#3D52A4","#F9B97E",as.character(SCP::palette_scp(n=6))))
mtx1<-pData(mcds)[,c("fibrotic","Pseudotime")]
colnames(mtx1)[1]<-"type"
mtx2<-pData(mcds)[,c("subtype","Pseudotime")]
colnames(mtx2)[1]<-"type"
mtx<-rbind(mtx1,mtx2)
mtx$type<-factor(as.character(mtx$type),levels=c("Anti-fibrotic VIC","Neutral VIC","VIC1","VIC2","VIC4","VIC3","VIC5","VIC6"))
p3<-ggplot(mtx,aes(x=Pseudotime,color=type,fill=type))+geom_density(aes(y=..density..),alpha=0.3)+facet_wrap(~type,nrow=1)+
  scale_color_manual(values=c("#3D52A4","#F9B97E",as.character(SCP::palette_scp(n=6))))+
  scale_fill_manual(values=c("#3D52A4","#F9B97E",as.character(SCP::palette_scp(n=6))))+theme_scp()
p2<-plot_cell_trajectory(mcds,color_by="subtype",show_branch_points = FALSE)+facet_wrap(~subtype,nrow=1)+
  scale_color_manual(values=c(as.character(SCP::palette_scp(n=6))))
p1<-(p1 & theme_void() & theme(legend.position="none",strip.text=element_blank()) )
p2<-(p2 & theme_void() & labs(y="") & theme(legend.position="none",strip.text = element_blank(),
                                           axis.text.y=element_blank(),axis.ticks.y=element_blank()))
library(cowplot)

# 设置 p1 和 p2 的组合
p1 <- p1 + theme(plot.margin = margin(0, 0, 0, 20))  # 右边空隙减小 index2
p2 <- p2 + theme(plot.margin = margin(0, 0, 0, 0)) # 左边空隙减小 index4
p1_p2 <- cowplot::plot_grid(
  p1, 
  p2 , 
  rel_widths = c(2.2, 6),  # 设置 p1 和 p2 的宽度比例
  nrow = 1, 
  align = "h",  # 水平对齐
  axis="tb"
)
p1_p2
# 将 p3 和它的图例分开
p3_legend <- cowplot::get_legend(p3)  # 提取图例
p3 <- p3 + theme(legend.position = "none")  # 移除 p3 的图例
p3<-p3+theme(plot.margin=margin(0,0,0,0))
# 整体拼接：p1_p2 在上，p3 在下，图例单独放右侧
final_plot <- cowplot::plot_grid(
  cowplot::plot_grid(p1_p2, p3, ncol = 1, rel_heights = c(1, 1)),  # 上下两部分
  p3_legend,  # 图例单独放置
  ncol = 2, 
  rel_widths = c(8, 1),  # 设置主图和图例的宽度比例
  align="vh",axis="lr"
)
final_plot
pdf("vic_MV_AntiFibrotic.pdf",width=16,height=5)
print(final_plot)
dev.off()

mcds<-readRDS("~/monocle_Rds/p7_viec_trans.Rds")
mcds<-orderCells(mcds,reverse = TRUE)
saveRDS(mcds,file="p7_viec_trans_reverse.rds")
beam<-readRDS("~/monocle_Rds/p7_viec_trans_bp2_BEAM_res.rds")
write.xlsx(beam,file="viec_trans_bp2_beam.xlsx")
pData(mcds)$Group<-factor(as.character(pData(mcds)$Group),levels=c("MV-normal","MV-moderate","MV-severe","TV-normal","TV-moderate","TV-severe"))
pData(mcds)$trans<-ifelse(pData(mcds)$trans=="transitional_VEC","Transition VEC","Transition VIC")
at2<-c("#3D52A4","#F9B97E")
b13<-c("#50a3a4","#fcaf38","#b60e00","#9fb4c6","#ffdd1c","#f95335")
c13<-c("#9BCFE6","#0087C0","#71C192","#4CBF46","#73AD30","#FFAC3D","#FF7800","#FF9189","#FD112C","#E1718E",
       "#A379BE","#9D769B","#F6E478","#BE5214")
p1<-plot_cell_trajectory(mcds,color_by="trans",show_branch_points = FALSE)+facet_wrap(~trans,nrow=1)+
  scale_color_manual(values=c(at2))
mtx1<-pData(mcds)[,c("trans","Pseudotime")]
colnames(mtx1)[1]<-"type"
mtx1$subtype<-mtx1$type
mtx2<-pData(mcds)[,c("Group","subtype","Pseudotime")]
colnames(mtx2)[1]<-"type"
mtx<-rbind(mtx1,mtx2)
mtx$type<-factor(as.character(mtx$type),levels=c("Transition VEC","Transition VIC","MV-normal","MV-moderate","MV-severe","TV-normal","TV-moderate","TV-severe"))
p3<-ggplot(mtx,aes(x=Pseudotime,color=subtype,fill=subtype))+geom_density_ridges(aes(y=..density..),alpha=0.3)+facet_wrap(~type,nrow=1)+
  scale_color_manual(values=c(at2,b13))+scale_x_reverse() +
  scale_fill_manual(values=c(at2,b13))+theme_scp()
p2<-plot_cell_trajectory(mcds,color_by="subtype",show_branch_points = FALSE)+facet_wrap(~Group,nrow=1)+
  scale_color_manual(values=c(b13))
p1<-(p1 & theme_void() & theme(legend.position="none",strip.text=element_blank()) )
p2<-(p2 & theme_void() & labs(y="") & theme(legend.position="none",strip.text = element_blank(),
                                            axis.text.y=element_blank(),axis.ticks.y=element_blank()))
library(cowplot)

# 设置 p1 和 p2 的组合
p1 <- p1 + theme(plot.margin = margin(0, 0, 0, 20))  # 右边空隙减小 index2
p2 <- p2 + theme(plot.margin = margin(0, 0, 0, 0)) # 左边空隙减小 index4
p1_p2 <- cowplot::plot_grid(
  p1, 
  p2 , 
  rel_widths = c(2.2, 6),  # 设置 p1 和 p2 的宽度比例
  nrow = 1, 
  align = "h",  # 水平对齐
  axis="tb"
)
p1_p2
# 将 p3 和它的图例分开
p3_legend <- cowplot::get_legend(p3)  # 提取图例
p3 <- p3 + theme(legend.position = "none")  # 移除 p3 的图例
p3<-p3+theme(plot.margin=margin(0,0,0,0))
# 整体拼接：p1_p2 在上，p3 在下，图例单独放右侧
final_plot <- cowplot::plot_grid(
  cowplot::plot_grid(p1_p2, p3, ncol = 1, rel_heights = c(1, 1)),  # 上下两部分
  p3_legend,  # 图例单独放置
  ncol = 2, 
  rel_widths = c(8, 1),  # 设置主图和图例的宽度比例
  align="vh",axis="lr"
)
final_plot
pdf("viec_Trans.pdf",width=15,height=5)
print(final_plot)
dev.off()

setwd("/data2/Project/F_Group_Analysis/7_Analysis/VIC_monocle/20250209")
mcds<-readRDS("p7_viec_trans_reverse.rds")
f<-read.xlsx("Trans_vennData.xlsx")
glist<-f$`Beam|DEGs`[!is.na(f$`Beam|DEGs`)]
p<-plot_genes_branched_heatmap(mcds[glist,],branch_point = 2,num_clusters = 3,return_heatmap = T)
res<-c()
for( i in unique(p$annotation_row$Cluster)){
  tmp<-p$annotation_row
  gn<-rownames(tmp[tmp$Cluster==i,,drop=FALSE])
  eg<-bitr(gn,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")
  ego<-enrichGO(eg$ENTREZID,OrgDb="org.Hs.eg.db",ont="BP",pvalueCutoff = 1,qvalueCutoff = 1,readable = TRUE)
  ot<-as.data.frame(ego)
  ot$group=paste0("Cluster",i)
  if(i==1){
    res=ot
  }else{res=rbind(res,ot)}
}
sheet<-list(gene_cluster=p$annotation_row,gene_go=res)
write.xlsx(sheet,file="Trans_vennData_cls.xlsx",rowNames=TRUE)
pdf("Trans_bp2_heatmap.pdf",width=8,height=8)
plot_genes_branched_heatmap(mcds[glist,],branch_point = 2,num_clusters = 3)
dev.off()

mcds<-readRDS("../vic_MV_AntiFibrotic.Rds")
f<-read.xlsx("MV_bp4_vennData.xlsx")
glist<-f[,3][!is.na(f[,3])]
p<-plot_genes_branched_heatmap(mcds[glist,],branch_point = 4,num_clusters = 3,return_heatmap = T)
res<-c()
for( i in unique(p$annotation_row$Cluster)){
  tmp<-p$annotation_row
  gn<-rownames(tmp[tmp$Cluster==i,,drop=FALSE])
  eg<-bitr(gn,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")
  ego<-enrichGO(eg$ENTREZID,OrgDb="org.Hs.eg.db",ont="BP",pvalueCutoff = 1,qvalueCutoff = 1,readable = TRUE)
  ot<-as.data.frame(ego)
  ot$group=paste0("Cluster",i)
  if(i==1){
    res=ot
  }else{res=rbind(res,ot)}
}
sheet<-list(gene_cluster=p$annotation_row,gene_go=res)
write.xlsx(sheet,file="MV_bp4_vennData_cls.xlsx",rowNames=TRUE)
pdf("MV_bp4_heatmap.pdf",width=8,height=8)
plot_genes_branched_heatmap(mcds[glist,],branch_point = 4,num_clusters = 3)
dev.off()

f<-read.xlsx("MV_bp5_vennData.xlsx")
glist<-f[,3][!is.na(f[,3])]
p<-plot_genes_branched_heatmap(mcds[glist,],branch_point = 5,num_clusters = 3,return_heatmap = T)
res<-c()
for( i in unique(p$annotation_row$Cluster)){
  tmp<-p$annotation_row
  gn<-rownames(tmp[tmp$Cluster==i,,drop=FALSE])
  eg<-bitr(gn,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")
  ego<-enrichGO(eg$ENTREZID,OrgDb="org.Hs.eg.db",ont="BP",pvalueCutoff = 1,qvalueCutoff = 1,readable = TRUE)
  ot<-as.data.frame(ego)
  ot$group=paste0("Cluster",i)
  if(i==1){
    res=ot
  }else{res=rbind(res,ot)}
}
sheet<-list(gene_cluster=p$annotation_row,gene_go=res)
write.xlsx(sheet,file="MV_bp5_vennData_cls.xlsx",rowNames=TRUE)
pdf("MV_bp5_heatmap.pdf",width=8,height=8)
plot_genes_branched_heatmap(mcds[glist,],branch_point = 5,num_clusters = 3)
dev.off()

mcds<-readRDS("~/monocle_Rds/p7_vic_TV2.Rds")
f<-read.xlsx("/data2/Project/F_Group_Analysis/7_Analysis/VIC_monocle/20250209/TV-vennData.xlsx")
glist<-f[,3][!is.na(f[,3])]
p<-plot_genes_branched_heatmap(mcds[glist,],branch_point = 1,num_clusters = 3,return_heatmap = T)
res<-c()
for( i in unique(p$annotation_row$Cluster)){
  tmp<-p$annotation_row
  gn<-rownames(tmp[tmp$Cluster==i,,drop=FALSE])
  eg<-bitr(gn,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")
  ego<-enrichGO(eg$ENTREZID,OrgDb="org.Hs.eg.db",ont="BP",pvalueCutoff = 1,qvalueCutoff = 1,readable = TRUE)
  ot<-as.data.frame(ego)
  ot$group=paste0("Cluster",i)
  if(i==1){
    res=ot
  }else{res=rbind(res,ot)}
}
sheet<-list(gene_cluster=p$annotation_row,gene_go=res)
write.xlsx(sheet,file="TV_vennData_cls.xlsx",rowNames=TRUE)
pdf("TV_heatmap.pdf",width=8,height=8)
plot_genes_branched_heatmap(mcds[glist,],branch_point = 1,num_clusters = 3)
dev.off()
hmcols <- colorRampPalette(RColorBrewer::brewer.pal(11, "PiYG")[c(3:9)])(100)
plot_genes_branched_heatmap(mcds[glist,],branch_point = 1,num_clusters = 3,
                            hmcols = hmcols)
"#8E0152" "#C51B7D" "#DE77AE" "#F1B6DA" "#FDE0EF" "#F7F7F7" "#E6F5D0" "#B8E186" "#7FBC41"
"#4D9221" "#276419"
col_fun <- circlize::colorRamp2(
  c(-3,-2, -1,0,1,2 3),
  c("#7FBC41", "#B8E186","#E6F5D0", "#F7F7F7","#FDE0EF","#F1B6DA" , "#C51B7D")  # 或使用 PiYG 色系中的代表色
)
hmcols_smooth <- col_fun(seq(-3, 3, length.out = 256))
pdf("TV_heatmap.pdf",width=6,height=7)
plot_genes_branched_heatmap(
  mcds[glist,],
  branch_point = 1,
  num_clusters = 3,
  hmcols = col_fun(seq(-3, 3, length.out = 65))  # 映射为颜色向量
)
dev.off()
mcds<-readRDS("../VIC_monocle/vic_MV_AntiFibrotic.Rds")
f<-read.xlsx("/data2/Project/F_Group_Analysis/7_Analysis/VIC_monocle/20250209/MV_bp4_vennData.xlsx")
glist<-f[,3][!is.na(f[,3])]
pdf("MV_bp4_heatmap.pdf",width=6,height=7)
plot_genes_branched_heatmap(mcds[glist,],branch_point = 4,num_clusters = 3,
                            hmcols=col_fun(seq(-3, 3, length.out = 100))
                            )
dev.off()

