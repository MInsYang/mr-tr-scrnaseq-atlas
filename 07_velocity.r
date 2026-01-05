# ============================================================================
# velocity
# ============================================================================
# 本文件包含velocity相关的代码
# ============================================================================

# 加载function
source("scripts/utils.R")

#velocytoR环境debug
conda activate velocytoR
export LD_PRELOAD="/data2/usr/yangmy_conda/velocytoR/lib/libstdc++.so.6"
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
/usr/local/bin/R
#dyn.load("/data2/usr/yangmy_conda/velocytoR/lib/libstdc++.so.6")
library(velocyto.R)
loomInfo<-data.frame(loomlist=list.files("/data2/Project/F_Group_Analysis/Shihongjie_Valves_CR/",pattern=".loom",recursive = TRUE))

res<-c()
res<-prepareLoom(file.list=fh$loomlist,rename.prefix.list=fh$SID,old_prefix.list = fh$fname)
mtx1<-res$spliced
mtx2<-res$unspliced
colnames(mtx1)<-gsub(":","_",colnames(mtx1))
colnames(mtx2)<-gsub(":","_",colnames(mtx2))
cid<-rownames(vic@meta.data)[rownames(vic@meta.data) %in% colnames(mtx2)]
vic@assays$spliced<-CreateAssayObject(mtx1[,cid])
vic@assays$unspliced<-CreateAssayObject(mtx2[,cid])
vic<-subset(vic,cells=cid)
#vic <- RunSCVELO(srt = vic, assay_X = "RNA", group_by = "subtype", linear_reduction = "Harmony", nonlinear_reduction = "HarmonyUMAP2D")
lays<-VelocityPlot(vic, reduction = "HarmonyUMAP2D", plot_type = "stream",return_layer = TRUE)
VelocityPlot(vic, reduction = "HarmonyUMAP2D", group_by = "subtype")
library(ggplot2)
p<-CellDimPlot(vic,group.by = "subtype", reduction = "HarmonyUMAP2D",pt.size=2,pt.alpha=0.2,
            label = TRUE, label_insitu = TRUE,
            velocity = "stochastic", velocity_plot_type = "stream", 
            velocity_density = 2, velocity_smooth = 1, streamline_n = 20, 
            legend.position = "none", theme_use = "theme_blank"
)
pdf("vic_velocity.pdf",width=8,height=6)
print(p)
dev.off()

CellDimPlot(pancreas_sub,reduction="UMAP",velocity_plot_type="stream",velocity="stochastic",
            group.by="SubCellType",velocity_density = 2,velocity_smooth = 1, streamline_n = 20)

vic <- RunSCVELO(srt = vic, assay_X = "RNA",n_jobs = 1, group_by = "subtype", linear_reduction = "Harmonypca", nonlinear_reduction = "HarmonyUMAP2D")
saveRDS(vic,"vic_scvelo_res.rds")