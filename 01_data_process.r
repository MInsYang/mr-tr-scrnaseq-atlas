# ============================================================================
# 01 Data Preprocessing
# ============================================================================
# 本文件包含01 data preprocessing相关的代码
# ============================================================================

# 加载function
source("scripts/utils.R")

#构建要分析的样本数据表
crpath="/data2/Project/F_Group_Analysis/Shihongjie_Valves_CR/"
setwd("/data2/Project/F_Group_Analysis/7_Analysis")
fh1<-list.files(crpath,pattern = "*_ot")
sampleInfo<-data.frame(sample=fh1,multiSpecies="no")
sampleInfo$SID<-paste0("CE-",c("V38_Mit","V38_Tri","V40_Mit","V48_Tri","V49_Mit","V49_Tri","V52_Tri","V59_Mit","V60_Tri","V03_Mit","V03_Tri",
                               "V07_Tri","V17_Mit","V20_Tri","V24_Mit","V24_Tri","V30_Tri","V82_Mit","V82_Tri"))
sampleInfo$sample<-paste0(crpath,sampleInfo$sample,"/outs/filtered_feature_bc_matrix")
write.xlsx(sampleInfo,file="0_createInfo.xlsx",colNames=TRUE)

#读取input cellrange outs 文件夹的mtx文件，构建SeuratObject
for(f in 1:nrow(sampleInfo)){
  mtx<-Read10X(sampleInfo$sample[f])
  message("Read10X done: ",sampleInfo$sample[f]," prefix ",sampleInfo$SID[f])
  colnames(mtx)<-paste0(sampleInfo$SID[f],"_",colnames(mtx))
  tmp<-CreateSeuratObject(mtx,min.cells = 10)
  tmp$orig.ident<-sampleInfo$SID[f]
if(f == 1){
  alldata<-tmp
}else{alldata<-merge(alldata,tmp)}
}
alldata<-PercentageFeatureSet(alldata,pattern="^Mt-|^mt-|^MT-",col.name="percent.mito")
alldata<-PercentageFeatureSet(alldata,pattern="^HB[^(P|E|S)]",col.name = "percent.hb")
st1<-qcstat(alldata,with.hb="percent.hb",with.mt="percent.mito")
saveRDS(alldata,file="0_raw.rds")
#FeatureStatPlot(alldata,stat.by=c("percent.mito","percent.hb","nCount_RNA","nFeature_RNA"),group.by="orig.ident")

#使用DoubletFinder按样本进行双细胞鉴定
remotes::install_github("lzmcboy/DoubletFinder_204_fix",force=TRUE)
library(DoubletFinder)
dir.create("doubletFinder")
setwd("doubetFinder")
for( s in unique(alldata$orig.ident)){
  message("processing ",s)
tmp<-subset(alldata,orig.ident == s)
DoubletRate = ncol(tmp)*8*1e-6
seu <- NormalizeData(tmp) %>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData() %>% RunPCA() %>% RunUMAP( dims = 1:10) %>% FindNeighbors(dims=1:10) %>% FindClusters()
sweep.res.list <- paramSweep(seu, PCs = 1:10, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
pK_bcmvn <- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))
## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
homotypic.prop <- modelHomotypic(seu$seurat_clusters)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(DoubletRate*nrow(seu@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
seu <- doubletFinder(seu, PCs = 1:10,  pK = pK_bcmvn, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
saveRDS(seu,file=paste0(s,"_dbl.rds"))
}

#合并doubletFinder结果
for(rds in list.files("doubletFinder")){
  tmp<-readRDS(paste0("doubletFinder/",rds))
  if(rds==list.files("doubletFinder")[1]){
    res<-data.frame(cid=rownames(tmp@meta.data),
                    DF.classification=tmp@meta.data[,colnames(tmp@meta.data)[grepl("DF.classifications",colnames(tmp@meta.data))]])
    
  }else{
    df<-data.frame(cid=rownames(tmp@meta.data),
                   DF.classification=tmp@meta.data[,colnames(tmp@meta.data)[grepl("DF.classification",colnames(tmp@meta.data))]])
    res<-rbind(res,df)
  }
}
saveRDS(res,file="dblFinder_res.rds")

#将doubletFinder结果添加到SeuratObject的meta.data中
rownames(res)<-res$cid
alldata$DF.class<-res[rownames(alldata@meta.data),"DF.classification"]
saveRDS(alldata,file="0_raw.rds")

alldata <- PercentageFeatureSet(alldata, "^RP[SL]", col.name = "percent_ribo")
# Percentage hemoglobin genes - includes all genes starting with HB except HBP.
alldata <- PercentageFeatureSet(alldata, "^HB[^(P|E|S)]", col.name = "percent_hb")

#质控处理
clean<-subset(alldata,DF.class=="Singlet")
clean<-subset(clean,nFeature_RNA>500 & nFeature_RNA<4000 & percent.mito<10 & percent.hb<10)
clean <- clean[rowSums(clean@assays$RNA@counts > 0) >= 10, ]
sampleInfo<-data.frame(sample=c("CE-V03_Mit", "CE-V03_Tri","CE-V07_Tri","CE-V17_Mit","CE-V20_Tri", "CE-V24_Mit",
                                "CE-V24_Tri", "CE-V30_Tri","CE-V38_Mit","CE-V38_Tri", "CE-V40_Mit","CE-V48_Tri",
                                "CE-V49_Mit", "CE-V49_Tri","CE-V52_Tri", "CE-V59_Mit","CE-V60_Tri", "CE-V82_Mit",
                                "CE-V82_Tri"))
sampleInfo$VR<-c("severe","moderate","severe","severe","normal","normal",
                 "normal","moderate","moderate","normal","moderate","normal",
                 "normal","normal","severe","normal","severe","severe",
                 "moderate")
sampleInfo$tissue<-gsub("CE-V[0-9].*_","",sampleInfo$sample)
rownames(sampleInfo)<-sampleInfo$sample
write.xlsx(sampleInfo,file="sampleInfo.xlsx",rowNames=TRUE)
clean@meta.data[,c("VR","tissue")]<-sampleInfo[clean$orig.ident,c("VR","tissue")]
sampleInfo<-read.xlsx("~/sampleInfo客户.xlsx",rowNames=TRUE)
sampleInfo$Group<-ifelse(sampleInfo$Group=="Tri:moderete","Tri:moderate",as.character(sampleInfo$Group))
clean@meta.data[,c("Sample","Group","Group_overall","Group_seperate")]<-sampleInfo[clean$orig.ident,c("Sample_new","Group","Group_overall","Group_seperate")]
saveRDS(clean,file="clean1028.rds")
clean$Group<-paste0(clean$tissue,":",clean$VR)
clean$VR<-factor(as.character(clean$VR),levels=c("normal","moderate","severe"))
clean$Group<-factor(as.character(clean$Group),
                    levels=c("Mit:normal","Mit:moderate","Mit:severe","Tri:normal","Tri:moderate","Tri:severe"))

df_raw<-qcstat(raw,type="Sample",with.mt = "percent.mito",with.hb = "percent.hb")
colnames(df_raw)<-c("Cell count","Gene count(Mean)","Gene count(Median)",
                    "UMI count(Mean)","UMI count(Median)",
                    "Mitochondria gene percent(Mean)","Mitochondria gene percent(Median)",
                    "Hemoglobin gene percent(Mean)","Hemoglobin gene percent(Median)")
df_clean<-qcstat(clean,type="Sample",with.mt = "percent.mito",with.hb = "percent.hb")
colnames(df_clean)<-c("Cell count","Gene count(Mean)","Gene count(Median)",
                    "UMI count(Mean)","UMI count(Median)",
                    "Mitochondria gene percent(Mean)","Mitochondria gene percent(Median)",
                    "Hemoglobin gene percent(Mean)","Hemoglobin gene percent(Median)")
sheets<-list(raw=df_raw,filtered=df_clean)
write.xlsx(sheets,file="Table1.xlsx",rowNames=TRUE,colNames=TRUE)

#按样本计算pca图并绘制95%置信区间椭圆
pcaCol<-c(  "#3477a9" ,  "#f5b375"  , "#4a9d47" ,  "#e45a5f" ,  "#684797"   ,"#FFC0CB","grey35"  )
sampleInfo$Group<-paste0(sampleInfo$tissue,":",sampleInfo$VR)
agg_matrix <- AggregateExpression(clean, group.by = "Sample")$RNA
pca_result <- prcomp(t(agg_matrix), scale. = TRUE)
pca_data <- data.frame(PC1 = pca_result$x[, 1], 
                       PC2 = pca_result$x[, 2], 
                       Sample = rownames(pca_result$x))
pca_data$Group<-sampleInfo[pca_data$Sample,"Group"]
# 绘制PCA图，并按组添加置信区间
# 全部样本
p<-generate_plot(pca_data)
pdf("Figure/2.PCA1.pdf",width=7,height=5)
print(p+scale_color_manual(values=pcaCol))
dev.off()  

#修改细胞meta信息
clean@meta.data[,c("Sample","Group","Group_overall","Group_seperate")]<-sampleInfo[clean$orig.ident,c("Sample_new","Group","Group_overall","Group_seperate")]
clean$Sample<-factor(as.character(clean$Sample),levels=c("MV1","MV2","MV3","MV4","MV5","MV6","MV7","MV8","TV1","TV2","TV3","TV4","TV5","TV6","TV7","TV8","TV9","TV10","TV11"))
clean$Group<-factor(as.character(clean$Group),levels=c("Mit:normal","Mit:moderate","Mit:severe","Tri:normal","Tri:moderate","Tri:severe"))
clean$Group_overall<-factor(as.character(clean$Group_overall),levels=c("Control","Case"))
clean$Group_seperate<-factor(as.character(clean$Group_seperate),levels=c("MV-Control","MV-Case","TV-Control","TV-Case"))
clean$Group<-gsub("Mit:","MV-",clean$Group)
clean$Group<-gsub("Tri:","TV-",clean$Group)
clean$Group<-gsub("moderete","moderate",clean$Group)
clean$Group<-factor(clean$Group,levels=c("MV-normal","MV-moderate","MV-severe","TV-normal","TV-moderate","TV-severe"))
table(clean$Group)
saveRDS(clean,file="dat_celltype.rds")

p<-FeatureStatPlot(clean,stat.by=c("nFeature_RNA","nCount_RNA","percent.mito"),group.by="Sample") & NoLegend()
pdf("1.QCstat_Sample.pdf",width=16,height=4)
print(p)
dev.off()
p<-FeatureStatPlot(clean,stat.by=c("nFeature_RNA","nCount_RNA","percent.mito"),group.by="Group") & NoLegend()
pdf("1.QCstat_Group.pdf",width=7,height=4)
print(p)
dev.off()
p<-FeatureStatPlot(clean,stat.by=c("nFeature_RNA","nCount_RNA","percent.mito"),group.by="Group_overall") & NoLegend()
pdf("1.QCstat_Group_overall.pdf",width=7,height=4)
print(p)
dev.off()
p<-FeatureStatPlot(clean,stat.by=c("nFeature_RNA","nCount_RNA","percent.mito"),group.by="Group_seperate") & NoLegend()
pdf("1.QCstat_Group_seperate.pdf",width=7,height=4)
print(p)
dev.off()

#读取原始数据并绘制质控图
raw<-readRDS("0_raw.rds")
raw@meta.data[,c("VR","tissue")]<-sampleInfo[raw$orig.ident,c("VR","tissue")]
raw$Sample<-raw$orig.ident
raw$Group<-paste0(raw$tissue,":",raw$VR)
raw$VR<-factor(as.character(raw$VR),levels=c("normal","moderate","severe"))
raw$Group<-factor(as.character(raw$Group),
                    levels=c("Mit:normal","Mit:moderate","Mit:severe","Tri:normal","Tri:moderate","Tri:severe"))
p<-FeatureStatPlot(raw,stat.by=c("nFeature_RNA","nCount_RNA","percent.mito"),group.by="Sample") & NoLegend()
pdf("Figure/0.QCstat_Sample.pdf",width=16,height=4)
print(p)
dev.off()
p<-FeatureStatPlot(raw,stat.by=c("nFeature_RNA","nCount_RNA","percent.mito"),group.by="tissue") & NoLegend()
pdf("Figure/0.QCstat_tissue.pdf",width=7,height=4)
print(p)
dev.off()
p<-FeatureStatPlot(raw,stat.by=c("nFeature_RNA","nCount_RNA","percent.mito"),group.by="VR") & NoLegend()
pdf("Figure/0.QCstat_VR.pdf",width=7,height=4)
print(p)
dev.off()
p<-FeatureStatPlot(raw,stat.by=c("nFeature_RNA","nCount_RNA","percent.mito"),group.by="Group") & NoLegend()
pdf("Figure/0.QCstat_Group.pdf",width=7,height=4)
print(p)
dev.off()
