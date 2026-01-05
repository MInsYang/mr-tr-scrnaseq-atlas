# ============================================================================
# 03 细胞亚群分析---VEC
# ============================================================================
# 本文件包含03 VEC 亚群相关的代码
# ============================================================================

# 加载function
source("scripts/utils.R")


# VEC同样亚群处理 + 亚群注释
vec$subtype<-paste0("VEC",as.character(vec$Harmony_SNN_res.0.3))
#对cluster进行downsample之后计算monocle2
downsample_ratio <- 0.3
meta_data <- vic@meta.data
meta_data$cell <- rownames(meta_data)
set.seed(123)  # 设置随机种子以确保结果可重复
sampled_cells <- meta_data %>%
  group_by(subtype) %>%
  slice_sample(prop = downsample_ratio) %>%
  dplyr::pull(cell)
dsapl <- subset(vic, cells = sampled_cells)
table(dsapl@meta.data$subtype)
saveRDS(dsapl,file="vic_downsample0.3.rds")
run_monocle(dat,prefix="monocle2_vic",annotation="subtype",condition="Group",cores=1,topn=5000,order.method="default")

#CytoTRACE2
library(CytoTRACE2)
res<-cytotrace2(vic,is_seurat=TRUE,ncores=10,species="human")
saveRDS(res,file="vic_cytotrace2.rds")
vic<-readRDS("../../5.VIC_subtype.rds")
vic@reductions$umap<-vic@reductions$HarmonyUMAP2D
annotation<-data.frame(phenotype=vic$subtype,Row.names=rownames(vic@meta.data)) %>% set_rownames(.,colnames(vic))
plots <- plotData(cytotrace2_result = res, 
                  annotation = annotation,
                  expression_data = vic,is_seurat=TRUE
)
#saveRDS(plots,file="vic_cytotrace2_plots.rds")
library(patchwork)
p1<-plots$CytoTRACE2_UMAP
p2<-plots$CytoTRACE2_Potency
p3<-plots$CytoTRACE2_Relative_UMAP
p4<-plots$CytoTRACE2_Boxplot_byPheno
pdf("all_cytotrace2.pdf",width=12,height=10)
(p1+p2+p3+p4) + patchwork::plot_layout(ncol=2)
dev.off()

labels <- c('Differentiated', 'Unipotent', 'Oligopotent', 'Multipotent', 'Pluripotent', 'Totipotent')
colors <- c("#9E0142", "#F46D43", "#FEE08B", "#E6F598", "#66C2A5", "#5E4FA2")
vic<-readRDS("vic_cytotrace2.rds")
vic@meta.data[["CytoTRACE2_Score_clipped"]] <- 5.5 - 6 * vic@meta.data[["CytoTRACE2_Score"]] 
# Then, we apply the minimum and maximum cutoffs
vic@meta.data[["CytoTRACE2_Score_clipped"]]  <- -pmax(pmin(vic@meta.data[["CytoTRACE2_Score_clipped"]], 5), 0)
p2<-FeaturePlot(vic, "CytoTRACE2_Score_clipped",raster=FALSE) +
  scale_colour_gradientn(colours = rev(colors), na.value = "transparent",
                         # breaks=c(0, 0.08333333, 0.25000000, 0.41666667, 0.58333333, 0.75000000, 0.91666667, 1 ),
                         labels = c(labels),
                         limits=c(-5,0), name = "Potency score \n",
                         guide = guide_colorbar(frame.colour = "black", ticks.colour = "black"))
p1<-CellDimPlot(vic,group.by="CytoTRACE2_Potency",show_stat = FALSE,raster=FALSE,palcolor = list(rev(as.character(colors))))

p3<-FeaturePlot(vic, "CytoTRACE2_Relative",raster=FALSE) +
  scale_colour_gradientn(colours = (c( "#000004FF", "#3B0F70FF", "#8C2981FF", "#DE4968FF", "#FE9F6DFF", "#FCFDBFFF")),
                         na.value = "transparent",
                         limits=c(0,1),
                         breaks = seq(0,1, by = 0.2),
                         labels=c("0.0 (More diff.)", "0.2", "0.4", "0.6", "0.8", "1.0 (Less diff.)"),
                         name = "Relative\norder \n" ,
                         guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) 
vic$Phenotype<-vic$subtype
mtd <- vic@meta.data[c("Phenotype", "CytoTRACE2_Score")]
medians <- mtd %>%
  group_by(Phenotype) %>%
  summarise(CytoTRACE2_median_per_pheno = median(CytoTRACE2_Score, na.rm = TRUE)) %>%
  arrange(desc(CytoTRACE2_median_per_pheno))

# Join the median scores back to the original dataframe for coloring
mtd <- mtd %>%
  inner_join(medians, by = "Phenotype")

# Order Phenotype by median CytoTRACE2_Score for plotting
mtd$Phenotype <- factor(mtd$Phenotype, levels = medians$Phenotype)
p4<-ggplot(mtd[!is.na(mtd$Phenotype), ], aes(x = Phenotype, y = CytoTRACE2_Score)) +
  geom_boxplot(aes(fill = CytoTRACE2_median_per_pheno), width = 0.8, alpha = 0.5, outlier.shape = NA) +
  geom_jitter(aes(fill = CytoTRACE2_median_per_pheno), width = 0.05, height = 0, alpha = 0.5, shape = 21, stroke = 0.1, size = 1) +
  theme_classic() +
  scale_y_continuous(breaks=seq(0, 1, by = 0.2),  limits = c(0,1),
                     sec.axis = sec_axis(trans = ~., breaks = seq(0, 1, by = 1/12),
                                         labels = c("", 'Differentiated', "",'Unipotent', "",'Oligopotent', "",'Multipotent',"", 'Pluripotent', "",'Totipotent', "") )) +
  scale_fill_gradientn(colors = rev(colors),
                       breaks=c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
                       limits = c(0,1), 
                       labels = c(labels))+
  scale_color_gradientn(colors = rev(colors),
                        breaks=c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
                        limits = c(0,1), 
                        labels = c(labels))+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10))  +
  labs(x = "Phenotype", y = 'Potency score') +
  ggtitle("Developmental potential by phenotype") +
  theme(legend.position = "None",
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 12),         
        legend.text = element_text(size = 12),
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5, margin = margin(b = 20)),
        axis.ticks.y.right  = element_line(color = c("black", NA, "black",NA, "black",NA, "black",NA,"black", NA, "black",NA, "black")),
        aspect.ratio = 0.8,
        axis.ticks.length.y.right  = unit(0.3, "cm"))
pdf("all_cytotrace2.pdf",width=12,height=10)
(p1+p2+p3+p4) + patchwork::plot_layout(ncol=2)
dev.off()
