# ============================================================================
# 03 亚群分析---免疫细胞
# ============================================================================
# 本文件包含免疫细胞亚群分析的相关的代码 
# ============================================================================

# 加载function
source("scripts/utils.R")

# 免疫细胞亚群
dat<-readRDS(".dat_celltype.rds")  
for(ctype in c("Lymphocyte","Myeloid cell")){ 
  clean<-subset(dat,celltype==ctype)
  clean<-Integration_SCP(clean,batch = "Sample",integration_method = "Harmony",neighbor_k=30,
                         cluster_algorithm = "leiden",cluster_resolution = seq(0.3,0.7,by=0.1))
  saveRDS(clean,file=paste0(ctype,".rds"))
  od<-paste0("RD_",ctype)
  dir.create(od)
  for(i in paste0("Harmony_SNN_res.",seq(0.3,0.7,by=0.1))){
    clean<-RunDEtest(clean,group_by=i,only.pos = TRUE,fc.threshold = 1)
    deg<-clean@tools[[paste0("DEtest_",i)]]$AllMarkers_wilcox
    deg<-deg[order(deg$avg_log2FC,decreasing=TRUE),]
    deg %>% split(deg$group1) %>% write.xlsx(paste0(od,"/clsDEG_",i,".xlsx"))
  }
  clean@meta.data<-clean@meta.data[,!colnames(clean@meta.data) %in% paste0("Harmony_SNN_res.",seq(0.8,1.5,by=0.1))]
  p<-clustree(clean,prefix="Harmony_SNN_res.")
  pdf(paste0(od,"/clustree.pdf"),width=14,height=12)
  print(p)
  dev.off()
  for(i in seq(0.3,0.7,by=0.1)){
    p<-CellDimPlot(clean,group.by=paste0("Harmony_SNN_res.",i),label=TRUE,label_insitu = TRUE)
    pdf(paste0(od,"/umap_",i,".pdf"),width=7,height=5)
    print(p)
    dev.off()
  }
}


mtx<-MC@tools$Enrichment_subtype
library(jjAnno,lib.loc="/home/yangmy/R/junjun")
library(ClusterGVis,lib.loc="/home/yangmy/R/junjun")
#vic$celltype<-factor(vic$celltype,levels=lev)
deg<-mergeXLSX("RD_Myeloid cell/clsDEG_Harmony_SNN_res.0.4.xlsx")
top5<-deg %>% group_by(group1) %>% top_n(5,avg_log2FC)
top5$cluster<-paste0("MC",top5$group1)
#分开取20
result <- data.frame()
# 初始化一个变量保存已选的基因
selected_genes <- character()
# 按 group1 的顺序逐步处理每组数据
for (group in unique(deg$group1)) {
  # 获取当前组的记录，按 avg_log2FC 从大到小排序
  current_group <- deg %>%
    filter(group1 == group) %>%
    arrange(desc(avg_log2FC))
  # 过滤掉已被选中的基因
  current_group <- current_group %>%
    filter(!gene %in% selected_genes)
  # 选择 top 20 的记录
  top20 <- current_group %>%
    slice_head(n = 20)
  # 将当前组的 top20 添加到结果中
  result <- bind_rows(result, top20)
  # 更新已选中的基因列表
  selected_genes <- c(selected_genes, top20$gene)
}
# 查看最终结果
result
top20<-sortMTX(result,by = "group1",order =seq(1,8))
top20$cluster<-paste0("MC",top20$group1)
vic<-readRDS("20250109/MC_subtype.rds")
Seurat::Idents(vic)<-vic$subtype
st.data<-prepareDataFromscRNA(object=vic,diffData=top20,showAverage=TRUE)
enrich <- enrichCluster(object = st.data,
                        OrgDb = org.Hs.eg.db,
                        type = "BP",
                        organism = "hsa",
                        pvalueCutoff = 0.5,
                        topn = 5,
                        seed = 5201314)
enricha50 <- enrichCluster(object = st.data,
                           OrgDb = org.Hs.eg.db,
                           type = "BP",
                           organism = "hsa",
                           pvalueCutoff = 0.5,
                           topn = 50,
                           seed = 5201314)
write.xlsx(enricha50,file="MC_top20_GOBP_top50.xlsx")
markGenes<-top5$gene
colvic<-c("#6cb86a","#CE1E5D","#C5ACB5","#EFBC87","#F9E199","#8F5A32","#A1B0AB","#42519F")
colvic<-colvic[1:8]
pdf('MC_heatmap.pdf',height = 12,width = 12)
visCluster(object = st.data,
           plot.type = "both",markGenes = top5$gene,markGenes.side = "left",
           show_row_dend = F,ctAnno.col = colvic,
           column_names_rot = 45,go.col=rep(colvic,each=5),
           go.size=12,sample.col=colvic,cluster.order = c(1:9),
           annoTerm.data=enrich,line.side = "left")
dev.off()

deg<-mergeXLSX("RD_Lymphocyte/clsDEG_Harmony_SNN_res.0.4.xlsx")
top5<-deg %>% group_by(group1) %>% top_n(5,avg_log2FC)
#分开取20
result <- data.frame()
# 初始化一个变量保存已选的基因
selected_genes <- character()
# 按 group1 的顺序逐步处理每组数据
for (group in unique(deg$group1)) {
  # 获取当前组的记录，按 avg_log2FC 从大到小排序
  current_group <- deg %>%
    filter(group1 == group) %>%
    arrange(desc(avg_log2FC))
  # 过滤掉已被选中的基因
  current_group <- current_group %>%
    filter(!gene %in% selected_genes)
  # 选择 top 20 的记录
  top20 <- current_group %>%
    slice_head(n = 20)
  # 将当前组的 top20 添加到结果中
  result <- bind_rows(result, top20)
  # 更新已选中的基因列表
  selected_genes <- c(selected_genes, top20$gene)
}
# 查看最终结果
result
top20<-sortMTX(result,by = "group1",order =seq(1,8))
top20$cluster<-paste0("LC",top20$group1)
vic<-readRDS("20250109/LC_subtype.rds")
Seurat::Idents(vic)<-vic$subtype
st.data<-prepareDataFromscRNA(object=vic,diffData=top20,showAverage=TRUE)
enrich <- enrichCluster(object = st.data,
                        OrgDb = org.Hs.eg.db,
                        type = "BP",
                        organism = "hsa",
                        pvalueCutoff = 0.5,
                        topn = 5,
                        seed = 5201314)
enricha50 <- enrichCluster(object = st.data,
                           OrgDb = org.Hs.eg.db,
                           type = "BP",
                           organism = "hsa",
                           pvalueCutoff = 0.5,
                           topn = 50,
                           seed = 5201314)
write.xlsx(enricha50,file="LC_top20_GOBP_top50.xlsx")
markGenes<-top5$gene
colvic<-c("#6cb86a","#CE1E5D","#C5ACB5","#EFBC87","#F9E199","#8F5A32","#A1B0AB","#42519F")
colvic<-colvic[1:7]
pdf('LC_heatmap.pdf',height = 12,width = 12)
visCluster(object = st.data,
           plot.type = "both",markGenes = top5$gene,markGenes.side = "left",
           show_row_dend = F,ctAnno.col = colvic,
           column_names_rot = 45,go.col=rep(colvic,each=5),
           go.size=12,sample.col=colvic,cluster.order = c(1:9),
           annoTerm.data=enrich,line.side = "left")
dev.off()