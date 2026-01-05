#cellphonedb 矩阵及meta txt文件准备
vic<-readRDS("vic.rds")
trans<-readRDS("transitional_viec.rds")
trans$subtype<-trans$trans
lc<-readRDS("../subRD/20250109/LC_subtype.rds")
mc<-readRDS("../subRD/20250109/MC_subtype.rds")
dt<-merge(vic,y=c(trans,lc,mc))
norm_count<-GetAssayData(dt,layer="data")
norm_count<-as.data.frame(norm_count)
norm_count<-cbind(toupper(rownames(norm_count)),norm_count)
colnames(norm_count)[1]<-"Gene"
metadata<-data.frame(Cell=rownames(dt@meta.data),cell_type=dt$subtype)
metadata$cell_type<-gsub(' ','_',metadata$cell_type)
dir.name<-"dt_input"
dir.create(dir.name)
write.table(norm_count,file=paste0(dir.name,"/Normalized_counts.txt"),row.names=F,sep="\t",quote=F)
write.table(metadata,file=paste0(dir.name,"/cellphonedb_meta.txt"),row.names=F,sep="\t",quote=F)
saveRDS(all,file=paste0(dir.name,"/scrna.rds"))
dt1<-subset(dt,tissue=="Tri")
norm_count<-GetAssayData(dt1,layer="data")
norm_count<-as.data.frame(norm_count)
norm_count<-cbind(toupper(rownames(norm_count)),norm_count)
colnames(norm_count)[1]<-"Gene"
metadata<-data.frame(Cell=rownames(dt1@meta.data),cell_type=dt1$subtype)
metadata$cell_type<-gsub(' ','_',metadata$cell_type)
dir.name<-"TV_input"
dir.create(dir.name)
write.table(norm_count,file=paste0(dir.name,"/Normalized_counts.txt"),row.names=F,sep="\t",quote=F)
write.table(metadata,file=paste0(dir.name,"/cellphonedb_meta.txt"),row.names=F,sep="\t",quote=F)
saveRDS(all,file=paste0(dir.name,"/scrna.rds"))
dt2<-subset(dt,tissue=="Mit")
norm_count<-GetAssayData(dt2,layer="data")
norm_count<-as.data.frame(norm_count)
norm_count<-cbind(toupper(rownames(norm_count)),norm_count)
colnames(norm_count)[1]<-"Gene"
metadata<-data.frame(Cell=rownames(dt2@meta.data),cell_type=dt2$subtype)
metadata$cell_type<-gsub(' ','_',metadata$cell_type)
dir.name<-"MV_input"
dir.create(dir.name)
write.table(norm_count,file=paste0(dir.name,"/Normalized_counts.txt"),row.names=F,sep="\t",quote=F)
write.table(metadata,file=paste0(dir.name,"/cellphonedb_meta.txt"),row.names=F,sep="\t",quote=F)
saveRDS(all,file=paste0(dir.name,"/scrna.rds"))

library(dplyr)
dir.name<-"MV_sigLR"
input_name<-"../MV_input"
dir.create(dir.name)
setwd(dir.name)
fhlist<-list.files(path=paste0(input_name,"/cpdb_ot"))
file_m=paste0(input_name,"/cpdb_ot/",fhlist[grepl("statistical_analysis_means",fhlist)])
file_p=paste0(input_name,"/cpdb_ot/",fhlist[grepl("statistical_analysis_pvalues",fhlist)])
file_score=paste0(input_name,"/cpdb_ot/",fhlist[grepl("statistical_analysis_interaction_scores",fhlist)])
get_sigLRv5(file_m,file_p,file_score=file_score)
setwd("../")
flist<-list.files(path=dir.name)
res<-c()
for(f in flist){
  mtx<-read.xlsx(paste0(dir.name,"/",f),rowNames = TRUE)
  res<-rbind(res,mtx)
}
sig<-res[res$p_value<0.05,]
sig <- sig %>%
  tidyr::separate(variable, into = c("source", "target"), sep = "\\|",remove=FALSE) #%>% tidyr::separate(interacting_pair,into=c("ligand","receptor"),sep="_",remove=FALSE)
sig<-sig[sig$source != sig$target,]
write.xlsx(sig,file=paste0(dir.name,"/all_sig_res.xlsx"))
# 按比例相关性对细胞类型进行分块
ratio_vic<-prop.table(table(vic$subtype,vic$Sample),2)*100
ratio_trans<-prop.table(table(trans$trans,trans$Sample),2)*100
ratio_lc<-prop.table(table(lc$subtype,lc$Sample),2)*100
ratio_mc<-prop.table(table(mc$subtype,mc$Sample),2)*100
via<-readRDS("7.VIEC.rds")
via@graphs<-list()
viat<-subset(via,subtype %in% c("transitional_VEC","transitional_VIC"),invert=TRUE)
ratio_con<-prop.table(table(viat$subtype,viat$Sample),2)*100
ratio_all<-as.data.frame.matrix(t(rbind(ratio_vic,ratio_trans,ratio_lc,ratio_mc,ratio_con)))
pdf("cor_celltype.pdf",width=10,height=10)
pheatmap(cor(ratio_all))
dev.off()
res<-list("prop"=ratio_all,"cor"=cor(ratio_all))
saveRDS(res,file="prop_cor.rds")

# 按比例相关性，cellphonedb互作sig LR受配体数量统计
c1<-c("VIC4","transitional_VIC","MC5","VIC1","MC3")
c2<-c("MC2","LC1","MC7")
c3<-c("MC8","LC2","LC7","VIC8","LC6")
c4<-c("VIC7","VIC9","LC4")
c5<-c("LC3","VIC6","MC4")
c6<-c("LC5","VIC2","MC6")
c7<-c("VIC3","MC1","VIC5","transitional_VEC")
mtx1<-sig[sig$source %in% c1 & sig$target %in% c1,]
write.xlsx(mtx1,file=paste0(dir.name,"/c1_res.xlsx"))
mtx2<-sig[sig$source %in% c2 & sig$target %in% c2,]
write.xlsx(mtx2,file=paste0(dir.name,"/c2_res.xlsx"))
mtx3<-sig[sig$source %in% c3 & sig$target %in% c3,]
write.xlsx(mtx3,file=paste0(dir.name,"/c3_res.xlsx"))
mtx4<-sig[sig$source %in% c4 & sig$target %in% c4,]
write.xlsx(mtx4,file=paste0(dir.name,"/c4_res.xlsx"))
mtx5<-sig[sig$source %in% c5 & sig$target %in% c5,]
write.xlsx(mtx5,file=paste0(dir.name,"/c5_res.xlsx"))
mtx6<-sig[sig$source %in% c6 & sig$target %in% c6,]
write.xlsx(mtx6,file=paste0(dir.name,"/c6_res.xlsx"))
mtx7<-sig[sig$source %in% c7 & sig$target %in% c7,]
write.xlsx(mtx7,file=paste0(dir.name,"/c7_res.xlsx"))

dir.name<-"dt_sigLR"
input_name<-"../dt_input"
dir.create(dir.name)
setwd(dir.name)
fhlist<-list.files(path=paste0(input_name,"/cpdb_ot"))
file_m=paste0(input_name,"/cpdb_ot/",fhlist[grepl("statistical_analysis_means",fhlist)])
file_p=paste0(input_name,"/cpdb_ot/",fhlist[grepl("statistical_analysis_pvalues",fhlist)])
file_score=paste0(input_name,"/cpdb_ot/",fhlist[grepl("statistical_analysis_interaction_scores",fhlist)])
get_sigLRv5(file_m,file_p,file_score=file_score)
setwd("../")
flist<-list.files(path=dir.name)
res<-c()
for(f in flist){
  mtx<-read.xlsx(paste0(dir.name,"/",f),rowNames = TRUE)
  res<-rbind(res,mtx)
}
sig<-res[res$p_value<0.05,]
sig <- sig %>%
  tidyr::separate(variable, into = c("source", "target"), sep = "\\|",remove=FALSE) #%>% tidyr::separate(interacting_pair,into=c("ligand","receptor"),sep="_",remove=FALSE)
sig<-sig[sig$source != sig$target,]
write.xlsx(sig,file=paste0(dir.name,"/all_sig_res.xlsx"))


setwd("/data2/Project/F_Group_Analysis/7_Analysis/20250512")
cpdb_anno <- read.csv('/data2/usr/yangmy_conda/cpdb_anno.csv', header = T, row.names = 1)
GO_pvals <- read.delim("dt_input/cpdb_ot/statistical_analysis_pvalues_05_13_2025_161925.txt", check.names = FALSE)
GO_means <- read.delim("dt_input/cpdb_ot/statistical_analysis_means_05_13_2025_161925.txt", check.names = FALSE)
ins<-c( #Transitional VEC_specific
        "FCER2_integrin_aMb2_complex",
        "FCER2_integrin_aXb2_complex",
        "SPP1_integrin_a4b1_complex",
        "PLAUR_integrin_a4b1_complex",
        "FCER2_integrin_aVb3_complex",
        "POMC_MC1R",
        "5alphaDihydroprogesterone_byDHRS9_PGR",
        "EFNB1_EPHB2",
        "EFNB1_EPHB3",
        "FGF20_FGFR3",
        "LGALS3_MERTK",
        "Glutamate_byGLS_and_SLC1A3_GRIA3",
        "LTC4_byLTC4S_CYSLTR2",
        "PDGFB_PDGFRA",
        "PDGFB_PDGFR_complex",
        "PDGFB_PDGFRB",
        "ProstaglandinF2a_byAKR1B1_PTGFR",
        "atRetinoicAcid_byALDH1A1_RAreceptor_RARA",
        "atRetinoicAcid_byALDH1A1_RAreceptor_RARA_RXRA",
        "atRetinoicAcid_byALDH1A1_RAreceptor_RARB",
        "atRetinoicAcid_byALDH1A1_RAreceptor_RARB_RXRB",
        "atRetinoicAcid_byALDH1A1_RAreceptor_RARG",
        "atRetinoicAcid_byALDH1A1_RAreceptor_RXRA",
        "atRetinoicAcid_byALDH1A1_RAreceptor_RXRB",
        "HFE_TFRC",
        "WNT5B_FZD5_LRP5",
        "WNT5B_FZD5_LRP6",
        "WNT5B_FZD6_LRP5",
        "WNT5B_FZD6_LRP6",
        "WNT5B_SFRP5",
        "ESAM_ESAM",
        "FTH1_SCARA5",
        "PECAM1_CD38"
)
table(ins %in% GO_means$interacting_pair)
flt_pvals<-GO_pvals[GO_pvals$interacting_pair %in% ins,]
flt_means<-GO_means[GO_means$interacting_pair %in% ins,]
pdf("EndoMT.pdf",width=9,height=9)
ks_cpdb_bdPlot(file_pvals = flt_pvals,
               file_means = flt_means,
               cpdb_anno = cpdb_anno, 
               source_cells = c("MC1", "transitional_VEC"),
               target_cells = c("MC1", "transitional_VEC"),
               comm_cut = 0)
ks_cpdb_bdPlot(file_pvals = flt_pvals,
               file_means = flt_means,
               cpdb_anno = cpdb_anno, 
               source_cells = c("MC7", "LC1","conventional_VEC"),
               target_cells = c("MC7", "LC1","conventional_VEC"),
               comm_cut = 0)
ks_cpdb_bdPlot(file_pvals = flt_pvals,
               file_means = flt_means,
               cpdb_anno = cpdb_anno, 
               source_cells = c("MC2", "transitional_VIC"),
               target_cells = c("MC2", "transitional_VIC"),
               comm_cut = 0)
ks_cpdb_bdPlot(file_pvals = flt_pvals,
               file_means = flt_means,
               cpdb_anno = cpdb_anno, 
               source_cells = c("MC4", "LC3","conventional_VIC"),
               target_cells = c("MC4", "LC3","conventional_VIC"),
               comm_cut = 0)
dev.off()

library(CCPlotR)
flist<-list.files("dt_sigLR/")
mtx<-c()
for( f in flist[-1]){
  tmp<-read.xlsx(paste0("dt_sigLR/",f))
  mtx<-rbind(mtx,tmp)
}
mtx<-mtx[,-1]
#mtx<-read.xlsx("dt_all_sig_res.xlsx")
#st<-c("MC1", "transitional_VEC","MC7", "LC1","conventional_VEC","MC2", "transitional_VIC","MC4", "LC3","conventional_VIC")
sts<-c("MC1|transitional_VEC","MC7|conventional_VEC","LC1|conventional_VEC","MC2|transitional_VIC",
       "MC4|conventional_VIC","LC3|conventional_VIC")
mtx<-mtx[mtx$interacting_pair %in% ins & mtx$variable %in% sts, ]
mtx$interaction_name_2<-mtx$interacting_pair
mtx <- mtx %>%
  tidyr::separate(variable, into = c("source", "target"), sep = "\\|",remove=FALSE) 
mtx$n.target<-gsub("transitional_","Trans-",mtx$target)
mtx$n.source<-gsub("transitional_","Trans-",mtx$source)
mtx$n.target<-gsub("conventional_","Conven-",mtx$n.target)
mtx$n.source<-gsub("conventional_","Conven-",mtx$n.source)
mtx<-sortMultiMTX(mtx,by="n.source",by.1="n.target",order =c("MC1","MC7","LC1","MC2","MC4","LC3"),
                  order.1 =c("Trans-VEC","Conven-VEC","Trans-VIC","Conven-VIC"))
mtx$n.variable<-paste0(mtx$n.source,"_",mtx$n.target)
mtx$n.variable<-factor(as.character(mtx$n.variable),levels=unique(mtx$n.variable))
mtx$interaction_name_2<-mtx$interacting_pair
mtx_pairs <- left_join(mtx, cpdb_anno, by = "interaction_name_2")
mtx_pairs$pval_group<-factor(as.character(mtx_pairs$pval_group),levels=c("not sig","* 0.05","** 0.01"))
p1<-ggplot(mtx_pairs,aes(x=interaction_name,y=n.variable))+
  geom_point(aes(color=iScore,size=pval_group))+theme_scp()+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))+
  scale_color_distiller(palette="Spectral",direction=-1,type="div")+scale_size_manual(values=c(0.5,3,6))
p2<-ggplot(mtx_pairs,aes(x=interaction_name,y=n.variable))+
  geom_point(aes(color=means,size=pval_group))+theme_scp()+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))+
  scale_color_distiller(palette="Spectral",direction=-1,type="div")+scale_size_manual(values=c(0.5,3,6))
pdf("EndoMT_dot.pdf",width=10,height=5)
print(p1)
print(p2)
dev.off()
write.xlsx(mtx_pairs,file="EndoMT_mtx.xlsx")

test<-separate(mtx_pairs,"interaction_name",into=c("gene1","gene2"),sep=":",remove = FALSE)
split_gene1 <- unlist(strsplit(as.character(unique(test$gene1)), split = "\\+"))
split_gene2 <- unlist(strsplit(as.character(unique(test$gene2)), split = "\\+"))
dt<-readRDS("dt_input/scrna.rds")
dsource<-subset(dt,subtype %in% c("LC1","LC3","MC1","MC2","MC4","MC7"))
dtarget<-subset(dt,subtype %in% c("conventional_VEC","conventional_VIC","transitional_VEC","transitional_VIC"))
p<-GroupHeatmapy(dsource,features=unique(split_gene1),group.by="subtype",exp_legend_title = "zscore",
                 cluster_rows = FALSE,cluster_columns = FALSE,row_title = "",column_title = "",
                 flip=TRUE,nlabel = 0,show_row_names = TRUE,show_column_names = TRUE,
                 add_dot=TRUE)
p$plot


ins<-c("CRTAM_CADM1",
       "SPP1_integrin_a9b1_complex",
       "POMC_MC1R",
       "AGRN_PTPRS",
       "LRFN4_PTPRS",
       "Dihydrotestosterone_bySRD5A1_AR",
       "Dihydrotestosterone_bySRD5A3_AR",
       "Dehydroepiandrosterone_bySTS_ESR1",
       "2arachidonoylglycerol_byDAGLB_CNR1",
       "Cholesterol_byLIPA_RORA",
       "Dehydroepiandrosterone_bySTS_PPARA",
       "5alphaDihydroprogesterone_byDHRS9_PGR",
       "EFNB1_EPHB2",
       "EFNB1_EPHB3",
       "EFNB1_EPHB4",
       "EFNB1_EPHB6",
       "AREG_EGFR",
       "TGFA_EGFR",
       "TENM4_ADGRL1",
       "TENM4_ADGRL2",
       "LGALS9_P4HB",
       "Glutamate_byGLS_and_SLC1A3_GRIA3",
       "HLA-E_VSIR",
       "IL1RN_IL1_receptor",
       "IL1RAP_PTPRF",
       "LRFN4_PTPRF",
       "LRFN4_PTPRD",
       "LeukotrieneE4_byDPEP2_CYSLTR2",
       "LTC4_byLTC4S_CYSLTR2",
       "LeukotrieneC4_byLTC4S_CYSLTR2",
       "5SHpETE_byALOX5_OXER1",
       "PLAU_PLAUR",
       "PDGFB_PDGFRA",
       "PDGFC_PDGFRA",
       "PDGFB_PDGFR_complex",
       "PDGFB_PDGFRB",
       "SELPLG_SELL",
       "ProstaglandinF2a_byAKR1B1_PTGFR",
       "ProstaglandinF2a_byPRXL2B_PTGFR",
       "SEMA4D_PLXNB2",
       "SEMA4D_PLXNB1",
       "ThromboxaneA2_byTBXAS1_TBXA2R",
       "TGFB1_TGFbeta_receptor2",
       "TGFB1_TGFbeta_receptor1",
       "TGFB1_TGFBR3",
       "TNF_TNFRSF1A",
       "TNF_TNFRSF1B",
       "TNFSF10_TNFRSF10D",
       "TNFSF8_TNFRSF8",
       "WNT5B_FZD1_LRP5",
       "WNT2B_FZD1_LRP6",
       "WNT5B_FZD1_LRP6",
       "WNT2B_FZD2_LRP5",
       "WNT5B_FZD2_LRP5",
       "WNT2B_FZD2_LRP6",
       "WNT5B_FZD2_LRP6",
       "WNT2B_FZD3_LRP5",
       "WNT5B_FZD3_LRP5",
       "WNT2B_FZD3_LRP6",
       "WNT5B_FZD3_LRP6",
       "WNT2B_FZD4_LRP5",
       "WNT5B_FZD4_LRP5",
       "WNT2B_FZD4_LRP6",
       "WNT5B_FZD4_LRP6",
       "WNT5B_FZD5_LRP5",
       "WNT5B_FZD5_LRP6",
       "WNT5B_FZD7_LRP5",
       "WNT5B_FZD7_LRP6",
       "WNT2B_SFRP1",
       "WNT5B_SFRP1",
       "CD44_TYROBP",
       "CD93_IFNGR1",
       "FTH1_SCARA5",
       "FTL_SCARA5",
       "PECAM1_CD38",
       "PPIA_BSG",
       "IFNG_Type_II_IFNR"
)
flist<-list.files("dt_sigLR/")
mtx<-c()
for( f in flist[-1]){
  tmp<-read.xlsx(paste0("dt_sigLR/",f))
  mtx<-rbind(mtx,tmp)
}
mtx<-mtx[,-1]
sts<-c("MC1|VIC3","MC1|VIC5","MC4|VIC6","LC3|VIC6","MC5|VIC1","MC3|VIC1","MC6|VIC2","LC5|VIC2",
       "MC2|VIC4","LC4|VIC7","LC4|VIC9","LC6|VIC8")
mtx<-mtx[mtx$interacting_pair %in% ins & mtx$variable %in% sts, ]
mtx$interaction_name_2<-mtx$interacting_pair
mtx <- mtx %>%
  tidyr::separate(variable, into = c("source", "target"), sep = "\\|",remove=FALSE) 
mtx$n.target<-gsub("transitional_","Trans-",mtx$target)
mtx$n.source<-gsub("transitional_","Trans-",mtx$source)
mtx$n.target<-gsub("conventional_","Conven-",mtx$n.target)
mtx$n.source<-gsub("conventional_","Conven-",mtx$n.source)
mtx<-sortMultiMTX(mtx,by="n.source",by.1="n.target",order =c("MC1","MC4","LC3","MC5","MC3","MC6","LC5","MC2","LC4","LC6"),
                  order.1 =c("VIC3","VIC5","VIC6","VIC1","VIC2","VIC4","VIC7","VIC9","VIC8"))
mtx$n.variable<-paste0(mtx$n.source,"_",mtx$n.target)
mtx$n.variable<-factor(as.character(mtx$n.variable),levels=unique(mtx$n.variable))
mtx$interaction_name_2<-mtx$interacting_pair
mtx_pairs <- left_join(mtx, cpdb_anno, by = "interaction_name_2")
mtx_pairs$pval_group<-factor(as.character(mtx_pairs$pval_group),levels=c("not sig","* 0.05","** 0.01"))
p1<-ggplot(mtx_pairs,aes(x=interaction_name,y=n.variable))+
  geom_point(aes(color=iScore,size=pval_group))+theme_scp()+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))+
  scale_color_distiller(palette="Spectral",direction=-1,type="div")+scale_size_manual(values=c(0.5,3,6))
p2<-ggplot(mtx_pairs,aes(x=interaction_name,y=n.variable))+
  geom_point(aes(color=means,size=pval_group))+theme_scp()+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))+
  scale_color_distiller(palette="Spectral",direction=-1,type="div")+scale_size_manual(values=c(0.5,3,6))
pdf("antiFib_dot.pdf",width=18,height=8)
print(p1)
print(p2)
dev.off()
write.xlsx(mtx_pairs,file="antiFib_mtx.xlsx")

flt_pvals<-GO_pvals[GO_pvals$interacting_pair %in% ins,]
flt_means<-GO_means[GO_means$interacting_pair %in% ins,]
pdf("antiFib.pdf",width=20,height=20)
ks_cpdb_bdPlot(file_pvals = flt_pvals,
               file_means = flt_means,
               cpdb_anno = cpdb_anno, 
               source_cells = c("MC1", "VIC3","VIC5"),
               target_cells =c("MC1", "VIC3","VIC5"),
               comm_cut = 0,thresh = )
ks_cpdb_bdPlot(file_pvals = flt_pvals,
               file_means = flt_means,
               cpdb_anno = cpdb_anno, 
               source_cells = c("MC4", "LC3","VIC6"),
               target_cells = c("MC4", "LC3","VIC6"),
               comm_cut = 0)
ks_cpdb_bdPlot(file_pvals = flt_pvals,
               file_means = flt_means,
               cpdb_anno = cpdb_anno, 
               source_cells = c("MC5","MC3","VIC1"),
               target_cells = c("MC5","MC3","VIC1"),
               comm_cut = 0)
ks_cpdb_bdPlot(file_pvals = flt_pvals,
               file_means = flt_means,
               cpdb_anno = cpdb_anno, 
               source_cells = c("MC6","LC5","VIC2"),
               target_cells = c("MC6","LC5","VIC2"),
               comm_cut = 0)
ks_cpdb_bdPlot(file_pvals = flt_pvals,
               file_means = flt_means,
               cpdb_anno = cpdb_anno, 
               source_cells = c("MC2","VIC4"),
               target_cells =c("MC2","VIC4"),
               comm_cut = 0)
ks_cpdb_bdPlot(file_pvals = flt_pvals,
               file_means = flt_means,
               cpdb_anno = cpdb_anno, 
               source_cells = c("LC4","VIC7","VIC9"),
               target_cells = c("LC4","VIC7","VIC9"),
               comm_cut = 0)
ks_cpdb_bdPlot(file_pvals = flt_pvals,
               file_means = flt_means,
               cpdb_anno = cpdb_anno, 
               source_cells = c("LC6","VIC8"),
               target_cells = c("LC6","VIC8"),
               comm_cut = 0)
dev.off()


# CellChat
c1<-c("MC8","LC2","LC7")
c2<-c("LC1","MC7","conventional_VEC")
c3<-c("MC5","VIC1","MC3")
c4<-c("MC2","VIC4","transitional_VIC")
c5<-c("VIC7","VIC9","LC4")
c6<-c("LC3","VIC6","MC4","conventional_VIC")
c7<-c("VIC3","MC1","VIC5","transitional_VEC")
c8<-c("VIC8","LC6")
c9<-c("LC5","VIC2","MC6")
cms<-list(c1,c2,c3,c4,c5,c6,c7,c8,c9)
library(CellChat,lib.loc="~/R/cellchat/")
#cellchat + cellphonedb 手工对应后interaction文件，计算cellchat
trim=0.05
type="thresholdedMean"
distance.use="FALSE"
#MV
db.new<-readRDS("../20250526/CellChatDB_cellphonedb.human_C67_new.rds")
dat<-readRDS("../20250512/MV_input/scrna.rds")
head(dat@meta.data)
for( group in unique(dat$Group)){
  data.input = dat@assays$RNA@data[,rownames(dat@meta.data)[dat$Group==group]]
  identity=data.frame(group=dat@meta.data[colnames(data.input),"subtype"],row.names=colnames(data.input))
  colnames(identity)<-"labels"
  cellchat<-createCellChat(data.input,identity,group.by="labels")
  cellchat<-setIdent(cellchat,ident.use="labels")
  levels(cellchat@idents)
  groupSize <- as.numeric(table(cellchat@idents))
  PPI<-PPI.human
  cellchat@DB <- db.new
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat,do.fast=FALSE,do.DE=FALSE,min.cells=3)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(cellchat,type=type,trim=trim,distance.use=distance.use)
  cellchat <- computeCommunProbPathway(cellchat,thresh=1)
  cellchat <- aggregateNet(cellchat,thresh=1)
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP",thresh=1)
  saveRDS(cellchat,file=paste0("MV_db67_",group,"_cellchat_obj.rds"))
}

db.new<-readRDS("../20250526/CellChatDB_cellphonedb.human_C7_new.rds")
for( group in unique(dat$Group)){
  data.input = dat@assays$RNA@data[,rownames(dat@meta.data)[dat$Group==group]]
  identity=data.frame(group=dat@meta.data[colnames(data.input),"subtype"],row.names=colnames(data.input))
  colnames(identity)<-"labels"
  cellchat<-createCellChat(data.input,identity,group.by="labels")
  cellchat<-setIdent(cellchat,ident.use="labels")
  levels(cellchat@idents)
  groupSize <- as.numeric(table(cellchat@idents))
  PPI<-PPI.human
  cellchat@DB <- db.new
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat,do.fast=FALSE,do.DE=FALSE,min.cells=3)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(cellchat,type=type,trim=trim,distance.use=distance.use)
  cellchat <- computeCommunProbPathway(cellchat,thresh=1)
  cellchat <- aggregateNet(cellchat,thresh=1)
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP",thresh=1)
  message(paste0("Length of pathways: ", length(cellchat@netP$pathways)))
  nPatterns = 5
  cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing",
                                            k = nPatterns, heatmap.show = T)
  cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming",
                                            k = nPatterns, heatmap.show = T)
  saveRDS(cellchat,file=paste0("MV_db7_",group,"_cellchat_obj.rds"))
}
#-----------------------------------TV-----------------------------------
db.new<-readRDS("../20250526/CellChatDB_cellphonedb.human_C67_new.rds")
dat<-readRDS("../20250512/TV_input/scrna.rds")
head(dat@meta.data)
for( group in unique(dat$Group)){
  data.input = dat@assays$RNA@data[,rownames(dat@meta.data)[dat$Group==group]]
  identity=data.frame(group=dat@meta.data[colnames(data.input),"subtype"],row.names=colnames(data.input))
  colnames(identity)<-"labels"
  cellchat<-createCellChat(data.input,identity,group.by="labels")
  cellchat<-setIdent(cellchat,ident.use="labels")
  levels(cellchat@idents)
  groupSize <- as.numeric(table(cellchat@idents))
  PPI<-PPI.human
  cellchat@DB <- db.new
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat,do.fast=FALSE,do.DE=FALSE,min.cells=3)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(cellchat,type=type,trim=trim,distance.use=distance.use)
  cellchat <- computeCommunProbPathway(cellchat,thresh=1)
  cellchat <- aggregateNet(cellchat,thresh=1)
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP",thresh=1)
  saveRDS(cellchat,file=paste0("TV_db67_",group,"_cellchat_obj.rds"))
}

db.new<-readRDS("../20250526/CellChatDB_cellphonedb.human_C7.rds")
for( group in unique(dat$Group)){
  data.input = dat@assays$RNA@data[,rownames(dat@meta.data)[dat$Group==group]]
  identity=data.frame(group=dat@meta.data[colnames(data.input),"subtype"],row.names=colnames(data.input))
  colnames(identity)<-"labels"
  cellchat<-createCellChat(data.input,identity,group.by="labels")
  cellchat<-setIdent(cellchat,ident.use="labels")
  levels(cellchat@idents)
  groupSize <- as.numeric(table(cellchat@idents))
  cellchat@DB <- db.new
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat,do.fast=FALSE,do.DE=FALSE,min.cells=3)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(cellchat,type=type,trim=trim,distance.use=distance.use)
  cellchat <- computeCommunProbPathway(cellchat,thresh=1)
  cellchat <- aggregateNet(cellchat,thresh=1)
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP",thresh=1)
  message(paste0("Length of pathways: ", length(cellchat@netP$pathways)))
  nPatterns = 5
  cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing",
                                            k = nPatterns, heatmap.show = T)
  cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming",
                                            k = nPatterns, heatmap.show = T)
  saveRDS(cellchat,file=paste0("TV_db7_",group,"_cellchat_obj.rds"))
}


cellchat.tv1<-readRDS("TV_db67_TV-normal_cellchat_obj.rds")
cellchat.tv2<-readRDS("TV_db67_TV-moderate_cellchat_obj.rds")
cellchat.tv3<-readRDS("TV_db67_TV-severe_cellchat_obj.rds")
object.list <- list(Mild = cellchat.tv1, Moderate = cellchat.tv2, Severe=cellchat.tv3)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
df1<-rankNet(cellchat, mode = "comparison", comparison=c(3,2,1),measure = "weight",slot.name="net",
             sources.use = NULL, targets.use = NULL, stacked = F, do.flip=F,do.stat = TRUE,return.data=TRUE,
             pairLR=c("EFNB1_EPHB4","Cholesterol-Cholesterol-LIPA_RORA","LTC4-LTC4S_CYSLTR2","IFNG_IFNGR1_IFNGR2",
                      "Dihydrotestosterone-DHT-SRD5A3_AR","2arachidonoylglycerol-2AG-DAGLB_CNR1","5alphaDihydroprogesterone-DHRS9_PGR"))
df3<-rankNet(cellchat, mode = "comparison", comparison=c(1,2,3),measure = "weight",slot.name="net",
             sources.use = NULL, targets.use = NULL, stacked = T, do.flip=T,do.stat = TRUE,return.data=TRUE,
             pairLR=c("EFNB1_EPHB4","Cholesterol-Cholesterol-LIPA_RORA","LTC4-LTC4S_CYSLTR2","IFNG_IFNGR1_IFNGR2",
                      "Dihydrotestosterone-DHT-SRD5A3_AR","2arachidonoylglycerol-2AG-DAGLB_CNR1","5alphaDihydroprogesterone-DHRS9_PGR"))
cellchat.mv1<-readRDS("MV_db67_MV-normal_cellchat_obj.rds")
cellchat.mv2<-readRDS("MV_db67_MV-moderate_cellchat_obj.rds")
cellchat.mv3<-readRDS("MV_db67_MV-severe_cellchat_obj.rds")
object.list <- list(Mild = cellchat.mv1, Moderate = cellchat.mv2, Severe=cellchat.mv3)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
df2<-rankNet(cellchat, mode = "comparison", comparison=c(3,2,1),measure = "weight",slot.name="net",
             sources.use = NULL, targets.use = NULL, stacked = F, do.flip=F,do.stat = TRUE,return.data=TRUE,
             pairLR=c("EFNB1_EPHB4","Cholesterol-Cholesterol-LIPA_RORA","LTC4-LTC4S_CYSLTR2","IFNG_IFNGR1_IFNGR2",
                      "Dihydrotestosterone-DHT-SRD5A3_AR","2arachidonoylglycerol-2AG-DAGLB_CNR1","5alphaDihydroprogesterone-DHRS9_PGR"))
df4<-rankNet(cellchat, mode = "comparison", comparison=c(1,2,3),measure = "weight",slot.name="net",
             sources.use = NULL, targets.use = NULL, stacked = T, do.flip=T,do.stat = TRUE,return.data=TRUE,
             pairLR=c("EFNB1_EPHB4","Cholesterol-Cholesterol-LIPA_RORA","LTC4-LTC4S_CYSLTR2","IFNG_IFNGR1_IFNGR2",
                      "Dihydrotestosterone-DHT-SRD5A3_AR","2arachidonoylglycerol-2AG-DAGLB_CNR1","5alphaDihydroprogesterone-DHRS9_PGR"))
pdf("db67.pdf",width=7,height=5)
print(df1$gg.obj+theme_scp()+theme(axis.text.x=element_text(angle=90,size=10))+labs(title="TV"))
print(df2$gg.obj+theme_scp()+theme(axis.text.x=element_text(angle=90,size=10))+labs(title="MV"))
dev.off()
pdf("db67_stack.pdf",width=8,height=7)
print(df3$gg.obj+theme_scp()+labs(title="TV"))
print(df4$gg.obj+theme_scp()+labs(title="MV"))
dev.off()


cellchat.tv1<-readRDS("TV_db7_TV-normal_cellchat_obj.rds")
cellchat.tv2<-readRDS("TV_db7_TV-moderate_cellchat_obj.rds")
cellchat.tv3<-readRDS("TV_db7_TV-severe_cellchat_obj.rds")
object.list <- list(Mild = cellchat.tv1, Moderate = cellchat.tv2, Severe=cellchat.tv3)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
df1<-rankNet(cellchat, mode = "comparison", comparison=c(3,2,1),measure = "weight",slot.name="net",
             sources.use = NULL, targets.use = NULL, stacked = F, do.flip=F,do.stat = TRUE,return.data=TRUE,
             pairLR=c("PDGFB_PDGFRA","PDGFB_PDGFRB","SPP1_ITGAV_ITGB1"))
df3<-rankNet(cellchat, mode = "comparison", comparison=c(3,2,1),measure = "weight",slot.name="net",
             sources.use = NULL, targets.use = NULL, stacked = T, do.flip=T,do.stat = TRUE,return.data=TRUE,
             pairLR=c("PDGFB_PDGFRA","PDGFB_PDGFRB","SPP1_ITGAV_ITGB1"))
cellchat.mv1<-readRDS("MV_db7_MV-normal_cellchat_obj.rds")
cellchat.mv2<-readRDS("MV_db7_MV-moderate_cellchat_obj.rds")
cellchat.mv3<-readRDS("MV_db7_MV-severe_cellchat_obj.rds")
object.list <- list(Mild = cellchat.mv1, Moderate = cellchat.mv2, Severe=cellchat.mv3)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
df2<-rankNet(cellchat, mode = "comparison", comparison=c(3,2,1),measure = "weight",slot.name="net",
             sources.use = NULL, targets.use = NULL, stacked = F, do.flip=F,do.stat = TRUE,return.data=TRUE,
             pairLR=c("PDGFB_PDGFRA","PDGFB_PDGFRB","SPP1_ITGAV_ITGB1"))
df4<-rankNet(cellchat, mode = "comparison", comparison=c(3,2,1),measure = "weight",slot.name="net",
             sources.use = NULL, targets.use = NULL, stacked = T, do.flip=T,do.stat = TRUE,return.data=TRUE,
             pairLR=c("PDGFB_PDGFRA","PDGFB_PDGFRB","SPP1_ITGAV_ITGB1"))
pdf("db7.pdf",width=6,height=4.5)
print(df1$gg.obj+theme_scp()+theme(axis.text.x=element_text(angle=90,size=10))+labs(title="TV"))
print(df2$gg.obj+theme_scp()+theme(axis.text.x=element_text(angle=90,size=10))+labs(title="MV"))
dev.off()
pdf("db7_stacked.pdf",width=7,height=4)
print(df3$gg.obj+theme_scp()+theme(axis.text.x=element_text(angle=90,size=10))+labs(title="TV"))
print(df4$gg.obj+theme_scp()+theme(axis.text.x=element_text(angle=90,size=10))+labs(title="MV"))
dev.off()

tv13<-c("#9fb4c6","#ffdd1c","#f95335")
mv13<-c("#50a3a4","#fcaf38","#b60e00")
names(tv13)<-c( "Mild" ,  "Moderate" ,  "Severe" )
names(mv13)<-c("Mild" ,"Moderate",  "Severe")

plt<-df1$signaling.contribution
for( nm in unique(plt$name)){
  plt.df<-plt[plt$name==nm,]
  p<-ggplot(plt.df,aes(x=group,y=contribution.scaled))+scale_fill_manual(values=tv13)+
    geom_bar(aes(fill=group),stat="identity",position="dodge")+theme_scp()+labs(y="Information flow")
  pdf(paste0("TV_db7_",nm,".pdf"),width=4,height=3)
  print(p)
  dev.off()
}
tv<-plt[plt$name=="IFNG_IFNGR1_IFNGR2",]


plt<-df2$signaling.contribution
for( nm in unique(plt$name)){
  plt.df<-plt[plt$name==nm,]
  p<-ggplot(plt.df,aes(x=group,y=contribution.scaled))+scale_fill_manual(values=mv13)+
    geom_bar(aes(fill=group),stat="identity",position="dodge")+theme_scp()+labs(y="Information flow")
  pdf(paste0("MV_db7_",nm,".pdf"),width=4,height=3)
  print(p)
  dev.off()
}

plt<-df1$signaling.contribution
for( nm in unique(plt$name)){
  plt.df<-plt[plt$name==nm,]
  p<-ggplot(plt.df,aes(x=group,y=contribution.scaled))+scale_fill_manual(values=tv13)+
    geom_bar(aes(fill=group),stat="identity",position="dodge")+theme_scp()+labs(y="Information flow")
  pdf(paste0("TV_db67_",nm,".pdf"),width=4,height=3)
  print(p)
  dev.off()
}
tv<-plt[plt$name=="IFNG_IFNGR1_IFNGR2",]
tv$group<-paste0("TV ",tv$group)

plt<-df2$signaling.contribution
for( nm in unique(plt$name)){
  plt.df<-plt[plt$name==nm,]
  p<-ggplot(plt.df,aes(x=group,y=contribution.scaled))+scale_fill_manual(values=mv13)+
    geom_bar(aes(fill=group),stat="identity",position="dodge")+theme_scp()+labs(y="Information flow")
  pdf(paste0("MV_db67_",nm,".pdf"),width=4,height=3)
  print(p)
  dev.off()
}
mv<-plt[plt$name=="IFNG_IFNGR1_IFNGR2",]
mv$group<-paste0("MV ",mv$group)

tvmv<-rbind(tv,mv)
p<-ggplot(tvmv,aes(x=group,y=contribution.scaled))+scale_fill_manual(values=as.character(c(mv13,tv13)))+
  geom_bar(aes(fill=group),stat="identity",position="dodge")+theme_scp()+
  labs(y="Information flow")+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
pdf(paste0("TVMV_IFNG.pdf"),width=4,height=4)
print(p)
dev.off()