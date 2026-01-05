# ============================================================================
# 基因集打分 细胞周期+EMT
# ============================================================================
# 本文件包含基因集打分相关的代码
# ============================================================================

# 加载function
source("scripts/utils.R")

#cellcycle细胞周期评分
g12s<-c("MIR892B","MIR873","CTDSP2","PTENP1-AS","CDK2","CDK3","CDK4","PSME3","CDK6","CTDSPL","CDK7",
        "CDKN1A","CDK2AP2","CDKN1B","CDKN2A","CDKN2B","CDKN2C","CDKN2D","CDKN3","BTN2A2","GPNMB","KHDRBS1",
        "PLK2","DBF4","RCC1","PIM2","FAM107A","CHEK2","TREX1","ECD","PHB2","PRAP1","DCUN1D3","PLK3","PLK5","IQGAP3",
        "TRIM71","LSM11","CACUL1","E2F7","CYP1A1","SASS6","SDE2","DDX3X","DLG1","E2F1","E2F3","EGFR","ARID2","EIF4E",
        "EIF4EBP1","EIF4G1","AIF1","AKT1","EZH2","STOX1","FGF10","FHL1","MYO16","PHF8","PLCB1","SPDYA","STXBP4","MBLAC1",
        "FBXO7","ZNF324","KANK2","HINFP","ANKRD17","GIGYF2","APPL1","LATS2","MTBP","GLI1","GML","RGCC","BRD7","GSPT1","GPR132",
        "ANXA1","APBB1","APC","PRMT2","HSPA8","HYAL1","ID2","ID4","INHBA","ITGB1","KCNA5","GPR15LG","NANOGP8","MIR10A","MIR133A1",
        "MIR137","MIR15A","MIR15B","MIR16-1","MIR193A","MIR208A","MIR214","MIR221","MIR222","MIR26A1","MIR29A","MIR29B1","MIR29C",
        "MIR30C2","MDM2","MLF1","MAP3K11","MN1","MNAT1","MIR133B","MIR372","MUC1","MYC","ATM","NPAT","DDR2","ATP2B4","CRNN","RPS27L",
        "TFDP3","WAC","DACT1","CRLF3","ACTL6B","TRIAP1","GTSE1","METTL13","CPSF3","PKD1","PKD2","PLCG2","PLRG1","PML","POLE","INO80",
        "PAF1","PPP2CA","RFWD3","PBRM1","APPL2","PHF10","FBXW7","PPP3CA","PIDD1","PPP6C","AMBRA1","CENPJ","KMT2E","PRKDC","SUSD2",
        "MEPCE","TRIM39","TCIM","GJC2","PSME1","PSME2","PTEN","MIR362","MIR451A","MIR495","MIR515-1","MIR520A","MIR519D",
        "MIR520H","MIR503","ARID1B","RPTOR","USP29","USP37","PTPN6","CTDSP1","RHOU","BCAT1","RB1","RBBP8","RBL1",
        "RBL2","CCND1","BCL2","RDX","DPF2","ACTB","BCL7A","RPA2","RPL26","RPS6","RPS6KB1","RRM1","RRM2","CCL2","BID","STIL",
        "SKP2","SMARCA2","SMARCA4","SMARCB1","SMARCC1","DDRGK1","SMARCC2","SMARCD1","SMARCD2","SMARCD3","SMARCE1","SOX2",
        "ADAM17","TAF1","TAF10","TBX2","MIR638","TERT","TFDP1","TP53","UBE2E2","WEE1","CACNB4","ZNF655","PAGR1","CDC73","FBXO31",
        "JADE1","CUL5","DPF3","CAMK2A","FAM83D","TMEM14B","DPF1","ARID1A","CDC7","GFI1B","CASP2","USP26","CUL4B","CUL4A","CUL3","CUL2",
        "CUL1","KLF11","LSM10","PIAS1","CDK10","ACTL6A","CRADD","CCNA2","CCND2","BCL7C","BCL7B","KLF4",
        "CCND3","CCNE1","ACVR1","CCNH","TM4SF5","ACVR1B","LATS1","CCNE2","SLFN11","ADAMTS1","CDK1","CDC6","KIF14","CDC25A","CDC34")
g2m<-c("USH1C","RAD50","CDK2","CDK3","CDK4","CDK6","CDKN1A","AKAP8","NDC80","CENPF","KHDRBS1","NES","ARPP19",
       "TOPBP1","CHEK1","FOXN3","CHEK2","PABIR1","PLK3","HUS1B","NEK10","DNM2","ENSA","STOX1","CCNY","FHL1","ATF5",
       "PAXIP1","FOXM1","FBXL7","PLCB1","USP22","BRD4","SIN3A","SYF2","FBXO5","AKAP8L","KCNH5","VPS4A","CCDC57",
       "BABAM1","GPR132","DONSON","NOP53","HSPA2","HUS1","APP","USP50","USP17L2","MIR195","MIR19B1","FOXO4","MRE11",
       "NBN","ATM","NPM1","ORC1","RRM2B","PBX1","ING4","MRNIP","FZR1","TAOK3","MBTPS2","PPME1","LCMT1","DTL","UIMC1","CDK14",
       "ABCB1","PLK1","ETAA1","ATR","USP47","PKIA","CHFR","RCC2","TRIM39","AVEN","MTA3","TAOK1","BARD1","INIP","RAD17","RAD21",
       "RAD51C","RAD51B","RBBP8","CCND1","RFPL1","RINT1","MIIP","RRM1","CLSPN","BLM","NABP1","SKP2","INTS3","SMARCD3","BRCA1",
       "AURKA","TAF2","TP53","TPD52L1","UBE2A","WEE1","WNT10B","NABP2","FBXL15","BRCC3","RNASEH2B","KDM8","ATAD5","CALM1",
       "CTC1","DBF4B","CDK5RAP3","CALM2","CALM3","CDC7","ABRAXAS1","DYRK3","BRSK1","PPM1D","MASTL","ZFYVE19","CDC14B",
       "MBTPS1","RAB11A","IER3","NAE1","CCNA2","CCNB1","PHOX2B","CCNG1","BRSK2","TICRR","PKMYT1","LATS1","CCNB2","ZNF830",
       "CCNQ","AURKB","CHMP4C","TAOK2","VPS4B","MACROH2A1","BABAM2","CDK1","MELK","CDC6","KIF14","CDC25A","CDC25B","CDC25C"
       )
g12s<-c("CDK2","CDKN1A","CDKN1B","PLK2","CHEK2","TREX1","PRAP1","PLK3","E2F7","SDE2","GIGYF2","GML",
        "MDM2","MUC1","ATM","RPS27L","WAC","TRIAP1","GTSE1","PML","RFWD3","PIDD1","PRKDC","CCND1",
        "RPA2","RPL26","TP53","FBXO31","CASP2","CRADD")
g2m<-c("RAD50","CDKN1A","TOPBP1","CHEK1","FOXN3","PABIR1","HUS1B","SYF2","BABAM1",
       "DONSON","NOP53","HUS1","FOXO4","MRE11","NBN","ATM","ORC1","MRNIP","FZR1",
       "TAOK3","MBTPS2","DTL","UIMC1","PLK1","ETAA1","ATR","CHFR","TRIM39","TAOK1",
       "BARD1","INIP","RAD17","RBBP8","RINT1","CLSPN","BLM","NABP1","INTS3","BRCA1",
       "NABP2","BRCC3","CDK5RAP3","ABRAXAS1","BRSK1","CDC14B","MBTPS1","IER3",
       "NAE1","CCNG1","TICRR","ZNF830","TAOK2","BABAM2","CDK1","CDC6")
g2m<-g2m[g2m %in% rownames(vic)]
g12s<-g12s[g12s %in% rownames(vic)]
deg<-vic@tools$DEtest_Group_overall$AllMarkers_wilcox
deg<-deg[deg$p_val<0.05,]
nrow(deg)

g0<-c("CFLAR","CALCOCO1","YPEL3","CST3","SERINC1",
"CLIP4","PCYOX1","TMEM59","RGS2","YPEL5","CD63","CDH13",
"GSN","MR1","CYB5R1","AZGP1","ZFYVE1","DMXL1","EPS8L2","PTTG1IP",
"MIR22HG","PSAP","GOLGA8B","NEAT1","TXNIP")#"KIAA1109","MTRNR2L12" not found
g0<-read.table("~/G0.txt",header=TRUE)
g0<-g0$Veen[!is.na(g0$Veen)]
g0<-g0[!g0 %in% c(cc.genes.updated.2019$s.genes,cc.genes.updated.2019$g2m.genes)]

flist<-list(g0=g0[1:10],g2m=cc.genes.updated.2019$g2m.genes,s=cc.genes.updated.2019$s.genes)
vic<-AddModuleScore(vic,features=flist,name=c("g0.score","g2m.score","s.score"),search=TRUE)
FeatureStatPlot(vic,stat.by=c( "g0.score1",   "g2m.score2",     "s.score3"))
vic$G0.score<-vic$g0.score1
vic$G2M.Score<-vic$g2m.score2
vic$S.score<-vic$s.score3
meta<-vic@meta.data

flist<-list(g2m=cc.genes.updated.2019$g2m.genes,s=cc.genes.updated.2019$s.genes)
dat<-AddModuleScore(dat,features=flist,name=c("g0.score","g2m.score","s.score"),search=TRUE)

#colnames(vic@meta.data)[c(35,36,37)]<-c("G0.score","G2M.score","S.score")
cc.scores<-vic@meta.data[,c("G0.score","G2M.score","S.score")]
assign_phase <- function(scores) {
  # 找到每行的最大值对应的列名
  max_col <- apply(scores, MARGIN = 1, FUN = function(row) {
    colnames(scores)[which.max(row)]
  })
  # 将结果作为一个新列添加到数据框
  scores$Phase <- gsub(".score","",max_col)
  return(scores)
}
cc.scores<-assign_phase(cc.scores)
vic$Phase<-cc.scores[rownames(vic@meta.data),"Phase"]
table(vic$Phase)

CellDimPlot(vic,group.by="Phase",raster=FALSE,palette = "tron")
library(plotly)
plot_ly(
  data = cc.scores,
  x = ~G0.score,          # X轴对应 G0.score
  y = ~G2M.score,         # Y轴对应 G2M.score
  z = ~S.score,           # Z轴对应 S.score
  color = ~Phase,         # 点的颜色对应 Phase
  colors = c("red", "blue", "green"), # 可根据需要调整颜色
  type = "scatter3d",     # 设置为三维散点图
  mode = "markers"        # 使用点作为标记
) %>%
  layout(
    scene = list(
      xaxis = list(title = "G0.score"),
      yaxis = list(title = "G2M.score"),
      zaxis = list(title = "S.score")
    ),
    title = "Plot of Cell Cycle Scores"
  )
htmlwidgets::saveWidget(as_widget(p), "3D_Scatter_Plot.html")



vic<-readRDS("VIC_subtype.rds")
#vic<-CellCycleScoring(vic,s.features=cc.genes.updated.2019$s.genes,g2m.features=cc.genes.updated.2019$g2m.genes)
g2m<-g2m[g2m %in% deg$gene]
g12s<-g12s[g12s %in% deg$gene]
#g2m<-cc.genes.updated.2019$g2m.genes[cc.genes.updated.2019$g2m.genes %in% deg$gene] 
#g12s<-cc.genes.updated.2019$s.genes[cc.genes.updated.2019$s.genes %in% deg$gene]
length(g2m)
length(g12s)
vic<-CellCycleScoring(vic,s.features=g12s,g2m.features=g2m)
p1<-CellDimPlot(vic,group.by="Phase",raster=FALSE)
p1
meta<-vic@meta.data
p2<-ggplot(meta,aes(x=S.Score,y=G2M.Score))+geom_point(aes(color=Phase))+theme_scp()
p2
pdf("vic_cellcycle.pdf",width=12,height=5)
p1+p2
dev.off()
colvic<-c("#A6CEE3","#1F78B4","#B2DF8A","#34A02C","#FDBF6F","#FF7F01","#FB9A99","#E31A1D","#CAB2D6")
colanti<-c("lightgrey","lightgrey","#B2DF8A","lightgrey","#FDBF6F","#FF7F01","lightgrey","lightgrey","lightgrey")
colpro<-c("lightgrey","lightgrey","lightgrey","lightgrey","lightgrey","lightgrey","#FB9A99","#E31A1D","#CAB2D6")
vic$fibrotic<-ifelse(vic$subtype %in% c("VIC5","VIC3","VIC6"),"Anti-fibrotic VIC","Neutral VIC")
#vic$fibrotic<-ifelse(vic$subtype %in% c("VIC5","VIC3"),"Anti-fibrotic VIC","Neutral VIC")
vic$fibrotic<-ifelse(vic$subtype %in% c("VIC7","VIC8","VIC9"),"Pro-fibrotic VIC",as.character(vic$fibrotic))

p1<-CellDimPlot(vic,cells.highlight = rownames(vic@meta.data)[vic$fibrotic=="Anti-fibrotic VIC"],group.by = "subtype",
            raster=FALSE,show_stat = FALSE,palcolor = list(colanti))+labs(title="Anti-fibrotic")f
p2<-CellDimPlot(vic,cells.highlight = rownames(vic@meta.data)[vic$fibrotic=="Pro-fibrotic VIC"],group.by = "subtype",
            raster=FALSE,show_stat=FALSE,palcolor=list(colpro))+labs(title="Pro-fibrotic")
pdf("VIC_fibrotic.pdf",width=12,height=5)
p1+p2
dev.off()
CellDimPlot(vic,group.by = "fibrotic",raster=FALSE,show_stat = FALSE,palcolor = list(c("lightgrey","#1F78B4","#FF7F01")))

mdat<-merge(vic,vec)
mdat<-Integration_SCP(mdat,batch = "Sample",integration_method = "Harmony",cluster_algorithm = "leiden") 
saveRDS(mdat,file="VIEC.rds")

#
calculate_EMT_score_seurat <- function(seurat_obj, mes_genes, epi_genes, assay = "RNA", slot = "data") {
  # Extract expression matrix
  expr_matrix <- GetAssayData(seurat_obj, assay = assay, slot = slot)
  # Filter genes
  mes_present <- intersect(mes_genes, rownames(seurat_obj))
  epi_present <- intersect(epi_genes, rownames(seurat_obj))
  if (length(mes_present) == 0 || length(epi_present) == 0) {
    stop("No matching genes found in the expression matrix")
  }
  # Warnings for missing genes
  if (length(mes_present) != length(mes_genes)) {
    warning("Some mesenchymal genes are missing: ", paste(setdiff(mes_genes, mes_present), collapse = ", "))
  }
  if (length(epi_present) != length(epi_genes)) {
    warning("Some epithelial genes are missing: ", paste(setdiff(epi_genes, epi_present), collapse = ", "))
  }
  #mes_score <- colMeans(seurat_obj@assays$RNA@data[mes_genes, , drop = FALSE])
  #epi_score <- colMeans(seurat_obj@assays$RNA@data[epi_genes, , drop = FALSE])
  mes_score <- colSums(seurat_obj@assays$RNA@data[mes_present, , drop = FALSE])
  epi_score <- colSums(seurat_obj@assays$RNA@data[epi_present, , drop = FALSE])
  # Calculate EMT score
  emt_score <- mes_score - epi_score
  emt_score[is.na(emt_score)] <- 0 # 处理 NaN
  
  # Ensure vector length matches metadata
  if (length(emt_score) != ncol(expr_matrix)) {
    stop("Calculated EMT scores do not match the number of cells in the Seurat object")
  }
  # Add EMT score to metadata
  seurat_obj <- AddMetaData(seurat_obj, metadata = setNames(emt_score, colnames(expr_matrix)), col.name = "EMT_Score")
  return(seurat_obj)
}

# Example usage:
mes_genes <- c("AGER", "FN1", "MMP2", "SNAI2", "VIM", "ZEB2")
epi_genes <- c("CDH1", "CDH3", "CLDN4", "EPCAM", "MAL2", "ST14")
mdat_emt <- calculate_EMT_score_seurat(mdat, mes_genes, epi_genes)
FeatureStatPlot(mdat_emt,stat.by="EMT_Score",group.by="subtype",add_box = TRUE)