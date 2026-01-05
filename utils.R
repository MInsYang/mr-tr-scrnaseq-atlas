# ============================================================================
# 自定义函数库
# ============================================================================
# 本文件包含自定义函数，供其他脚本调用
# ============================================================================
library(SCP)
library(Seurat)
library(openxlsx)
library(DoubletFinder)
library(dplyr)
library(ggplot2)
library(patchwork)
library(tidyverse)
library(tidyr)
library(reshape2)
library(pheatmap)
library(RColorBrewer)
library(ComplexHeatmap)
library(bitr)
library(ggpubr)
library(ggsci)
library(ggrepel)


# Function: qcstat
qcstat<-function(obj,type="orig.ident",with.mt=NULL,with.rb=NULL){
  types<-data.table(table(obj@meta.data[,type]))$V1
  cellcount<-data.table(table(obj@meta.data[,type]))$N
  meangene<-as.data.frame(tapply(obj@meta.data[,"nFeature_RNA"], INDEX=obj@meta.data[,type], FUN=median))[,1]
  meancount<-as.data.frame(tapply(obj@meta.data[,"nCount_RNA"],INDEX=obj@meta.data[,type],FUN=median))[,1]
  format_value <- function(x) {
    ifelse(x == floor(x),  # 如果数值是整数
           formatC(x, format="d", big.mark=","),  # 格式化为整数并添加千分位逗号
           formatC(x, format="f", digits=3, big.mark=",")  # 否则，保留三位小数并添加千分位逗号
    )
  }
  if(!is.null(with.mt)){
    meanmt<-as.data.frame(tapply(obj@meta.data[,with.mt],INDEX=obj@meta.data[,type],FUN=median))[,1]
    df <- data.frame(
      cellcount = format_value(cellcount),
      meangene = format_value(meangene),
      meancount = format_value(meancount),
      meanpctmt = format_value(meanmt)
    )
  }else{
    df <- data.frame(
      cellcount = format_value(cellcount),
      meangene = format_value(meangene),
      meancount = format_value(meancount))
  }
  if(!is.null(with.rb)){
    meanrb<-as.data.frame(tapply(obj@meta.data[,with.rb],INDEX=obj@meta.data[,type],FUN=median))[,1]
    df$"meanpctrb"<-format_value(meanrb)
  }else{
    df <- df
    }
  rownames(df)<-types
  return(df)
}

# Function: ratio_plot
ratio_plot <- function(
    obj,
    sample.by = "Sample",
    anno.by = "plt2",
    condition.by = "Condition",
    strip.col = NULL,
    compa = TRUE,
    condition.col = c("#334d8f","#5bc4ed","#7ED321","#e6542b"),
    save.prefix = "ratio_plot",
    ord.condition = NULL,
    facet.ncol = 5,
    width = 8,
    height = 4,
    angle.x = 90,
    hjust.x = 0,
    vjust.x = 1,
    individual_y = NULL,
    compare.method = c("t.test", "wilcox", "anova", "tukey"),
    plot_type = c("bar","box")
){
  
  compare.method <- match.arg(compare.method)
  plot_type <- match.arg(plot_type)
  
  message(">>> Building frequency table...")
  
  ## ----------- 1. Data preparation -----------------------
  if (inherits(obj, "data.frame")) {
    df0 <- obj
  } else if (inherits(obj, "Seurat")) {
    df0 <- obj@meta.data
  } else {
    stop("obj must be Seurat or data.frame")
  }
  
  if (!all(c(sample.by, anno.by) %in% colnames(df0))) {
    stop("sample.by or anno.by not found in metadata!")
  }
  
  # frequency table
  fq <- as.data.frame.matrix(
    prop.table(table(df0[,sample.by], df0[,anno.by]), 1) * 100
  )
  
  # add condition
  if (!is.null(condition.by)) {
    cond <- unique(df0[, c(sample.by, condition.by)])
    colnames(cond) <- c("sample","Condition")
    rownames(cond) <- cond$sample
    fq$Condition <- cond[rownames(fq),"Condition"]
  } else {
    fq$Condition <- rownames(fq)
  }
  
  df <- reshape2::melt(fq)
  colnames(df)[1:3] <- c("Condition", "Celltype", "freq")
  
  if (!is.null(ord.condition)) {
    df$Condition <- factor(as.character(df$Condition), levels = ord.condition)
  }
  
  # Save raw freq table
  openxlsx::write.xlsx(fq, file = paste0(save.prefix,"_frequency_table.xlsx"), rowNames = TRUE)
  
  ## ----------- 2. Plotting -----------------------
  message(">>> Creating plot...")
  
  p_bar <- ggplot(df, aes(Condition, freq, fill = Condition)) +
    scale_fill_manual(values = condition.col) +
    stat_summary(geom = "bar", fun = mean, width = 0.5) +
    stat_summary(geom = "errorbar", fun.data = mean_se, width = 0.1) +
    geom_point(color = "black", size = 1, alpha = 0.8) +
    facet_wrap(~Celltype, ncol = facet.ncol, scales = "free_y") +
    labs(x="", y="Percentage(%)") +
    theme_pubr() +
    theme(axis.text.x = element_text(angle = angle.x, hjust = hjust.x, vjust = vjust.x)) +
    guides(fill="none")
  
  p_box <- ggplot(df, aes(Condition, freq, fill = Condition)) +
    scale_fill_manual(values = condition.col) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 1.5, alpha = 0.7) +
    facet_wrap(~Celltype, ncol = facet.ncol, scales = "free_y") +
    labs(x="", y="Percentage(%)") +
    theme_pubr() +
    theme(axis.text.x = element_text(angle = angle.x, hjust = hjust.x, vjust = vjust.x)) +
    guides(fill="none")
  
  p_out <- if (plot_type == "box") p_box else p_bar
  
  
  ## ----------- 3. Statistical comparisons -----------------------
  if (compa) {
    message(">>> Running statistical analysis: ", compare.method)
    
    if (compare.method %in% c("t.test","wilcox")) {
      
      res <- rstatix::compare_means(
        freq ~ Condition,
        group.by = "Celltype",
        data = df,
        method = compare.method
      )
      openxlsx::write.xlsx(res, paste0(save.prefix,"_",compare.method,"_results.xlsx"))
      
    } else if (compare.method == "anova") {
      
      anova_res <- df %>% dplyr::group_by(Celltype) %>% 
        rstatix::anova_test(freq ~ Condition)
      openxlsx::write.xlsx(anova_res, paste0(save.prefix,"_anova_results.xlsx"))
      
    } else if (compare.method == "tukey") {
      
      tukey_list <- lapply(split(df, df$Celltype), function(d) {
        
        d <- droplevels(d)
        if (nlevels(d$Condition) < 2) return(NULL)
        
        fit <- aov(freq ~ Condition, data = d)
        tk  <- TukeyHSD(fit, "Condition")
        
        res <- tk$Condition
        res <- as.data.frame(res)
        res$Comparison <- rownames(res)
        rownames(res) <- NULL
        res$Celltype <- unique(d$Celltype)
        
        colnames(res) <- tolower(colnames(res))
        if (!"p.adj" %in% colnames(res)) {
          if ("p adj" %in% colnames(res)) res$p.adj <- res$`p adj`
          if ("p.value" %in% colnames(res)) res$p.adj <- res$`p.value`
        }
        
        res <- res[, c("celltype","comparison","diff","lwr","upr","p.adj")]
        return(res)
      })
      
      tukey_df <- do.call(rbind, tukey_list)
      openxlsx::write.xlsx(tukey_df, paste0(save.prefix,"_TukeyHSD.xlsx"), rowNames = FALSE)
    }
  }
  
  ## ----------- 4. Output plot -----------------------
  message(">>> Saving PDF...")
  
  pdf(paste0(save.prefix, "_plot.pdf"), width=width, height=height)
  print(p_out)
  dev.off()
  
  return(p_out)
}

# Function: generate_ellipse
generate_ellipse <- function(center, cov_matrix, level = 0.95, n = 100) {
  # 计算椭圆的半长轴和半短轴
  chisq_val <- qchisq(level, df = 2)  # 计算卡方值
  eig <- eigen(cov_matrix)  # 特征值分解
  angle <- atan2(eig$vectors[2, 1], eig$vectors[1, 1])  # 椭圆的角度
  scale <- sqrt(chisq_val * eig$values)  # 椭圆的长短半径
  
  t <- seq(0, 2 * pi, length.out = n)  # 生成角度序列
  ellipse_points <- data.frame(
    x = center[1] + scale[1] * cos(t) * cos(angle) - scale[2] * sin(t) * sin(angle),
    y = center[2] + scale[1] * cos(t) * sin(angle) + scale[2] * sin(t) * cos(angle)
  )
  
  return(ellipse_points)
}


# Function: generate_plot
generate_plot<-function(pca_data){
  p<-ggplot(pca_data, aes(x = PC1, y = PC2, color = Group)) +
    geom_point(size = 1.5) +  # 绘制样本点
    geom_text_repel(aes(label = Sample), size = 2, 
                    box.padding = 0.5,  # 调整文本框的边距
                    point.padding = 0.3,  # 调整点和文本之间的距离
                    segment.size = 0.5,max.overlaps = 100) + # 添加样本名称
    theme_scp() +  # 使用简约主题
    labs(title = "PCA Plot with 95% Confidence Ellipses by Group", 
         x = "Principal Component 1", 
         y = "Principal Component 2") +
    scale_color_discrete(name = "Group") +  # 设置颜色图例
    geom_path(data = do.call(rbind, lapply(split(pca_data, pca_data$Group), function(group_data) {
      # 检查组内样本数
      if (nrow(group_data) < 2) return(NULL)  # 如果样本数小于2，返回NULL
      center <- colMeans(group_data[, c("PC1", "PC2")])  # 计算中心点
      cov_matrix <- cov(group_data[, c("PC1", "PC2")])  # 计算协方差矩阵
      ellipse_points <- generate_ellipse(center, cov_matrix)  # 生成椭圆
      ellipse_points$Group <- group_data$Group[1]  # 添加组信息
      return(ellipse_points)
    })), aes(x = x, y = y, group = Group, color = Group), size = 1, alpha = 0.2)
  return(p)
}


# Function: prepareLoom
prepareLoom<-function(file="~/_BAM/G_DBAW5.loom",file.list=NULL,rename.prefix="Data_G_",
                      old_prefix=NULL,old_prefix.list=NULL,rename.prefix.list=NULL,suffix="-1"){
  oneloom<-function(filename=filename,rename.prefix=rename.prefix,old_prefix=old_prefix,suffix=suffix){
    ldat <- read.loom.matrices(filename)
    if(is.null(old_prefix)){
    old_prefix=paste0(tools::file_path_sans_ext(basename(filename)),":")
    }else{old_prefix=old_prefix}
    colnames(ldat$spliced)<-gsub(old_prefix,rename.prefix,colnames(ldat$spliced))
    colnames(ldat$unspliced)<-gsub(old_prefix,rename.prefix,colnames(ldat$unspliced))
    colnames(ldat$ambiguous)<-gsub(old_prefix,rename.prefix,colnames(ldat$ambiguous))
    if(!is.null(suffix)){
      colnames(ldat$spliced)=gsub("x$",suffix,colnames(ldat$spliced))
      colnames(ldat$unspliced)=gsub("x$",suffix,colnames(ldat$unspliced))
      colnames(ldat$ambiguous)=gsub("x$",suffix,colnames(ldat$ambiguous))
    }
    return(ldat)
  }
  
  result<-c()
  if(is.null(file.list)){
    if(!is.null(file)){stop("set at least one loom file")}
    result<-oneloom(filename=file,rename.prefix=rename.prefix,suffix=suffix)
  }
  if(!is.null(file.list)){
    if(!is.null(file)){warning("when setting 'file.list','file' will not be used")}
    if(!is.null(rename.prefix)){warning("when setting 'rename.prefix.list', 'rename.prefix' will not be used")}
    if(length(file.list) != length(rename.prefix.list)){stop(" 'rename.prefix.list' length must be the same length as 'filelist'")}
    for(i in 1:length(file.list)){
      message("geting ",file.list[i])
      tmp<-oneloom(filename=file.list[i],rename.prefix=rename.prefix.list[i],suffix=suffix,old_prefix=old_prefix.list[i])
      if(i==1){result=tmp}
      if(i>1){
        result$spliced<-cbind(result$spliced,tmp$spliced)
        result$unspliced<-cbind(result$unspliced,tmp$unspliced)
        result$ambiguous<-cbind(result$ambiguous,tmp$ambiguous)
      }
    }
  }
  return(result)
}


# Function: assign_phase
assign_phase <- function(scores) {
  # 找到每行的最大值对应的列名
  max_col <- apply(scores, MARGIN = 1, FUN = function(row) {
    colnames(scores)[which.max(row)]
  })
  # 将结果作为一个新列添加到数据框
  scores$Phase <- gsub(".score","",max_col)
  return(scores)
}


# Function: calculate_EMT_score_seurat
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


# Function: find_first_intersection
find_first_intersection <- function(dens1, dens2) {
  x_common <- dens1$x
  y_diff <- dens1$y - dens2$y
  crossing <- which(diff(sign(y_diff)) != 0) # 找到符号变化的点
  x_intersection <- x_common[crossing[1]]   # 第一个交点
  return(x_intersection)
}


# Function: mean_function
mean_function <- function(data, indices) {
  sample_data <- data[indices]  # 根据bootstrap索引抽样
  return(mean(sample_data))
}


# Function: run_cluster_anova
run_cluster_anova <- function(dat) {
  # 检查输入数据
  required_cols <- c("patient", "majorCluster", "loc")
  if (!all(required_cols %in% colnames(dat))) {
    stop("输入数据必须包含以下列: patient, majorCluster, loc")
  }
  
  # Step 1: 计算每个样本中各 majorCluster 的比例
  cluster_prop <- dat %>%
    group_by(patient, loc, majorCluster) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(patient, loc) %>%
    mutate(prop = n / sum(n)) %>%
    ungroup()
  
  # Step 2: 对每个 majorCluster 做 ANOVA
  anova_results <- cluster_prop %>%
    group_by(majorCluster) %>%
    do({
      model <- lm(prop ~ loc, data = .)
      tidy(aov(model)) %>% 
        filter(term == "loc") %>% 
        dplyr::select(term, df, statistic, p.value)
    }) %>%
    ungroup()
  
  # Step 3: 整理结果表格
  anova_results <- anova_results %>%
    rename(Cluster = majorCluster,
           Term = term,
           DF = df,
           F_value = statistic,
           P_value = p.value) %>%
    arrange(P_value)
  
  # Step 4: 返回结果
  return(anova_results)
}



# modify from SCP
GroupHeatmapy<-function (srt, features = NULL, group.by = NULL, split.by = NULL,legend_group.by=FALSE, 
                         within_groups = FALSE, grouping.var = NULL, numerator = NULL, legend.horizontal=TRUE,
                         cells = NULL, aggregate_fun = base::mean, exp_cutoff = 0, show_annotation_name=FALSE,
                         border = TRUE, flip = FALSE, slot = "counts", assay = NULL, 
                         exp_method = c("zscore", "raw", "fc", "log2fc", "log1p"), 
                         exp_legend_title = NULL, limits = NULL, lib_normalize = identical(slot, 
                                                                                           "counts"), libsize = NULL, feature_split = NULL, feature_split_by = NULL, 
                         n_split = NULL, split_order = NULL, split_method = c("kmeans", 
                                                                              "hclust", "mfuzz"), decreasing = FALSE, fuzzification = NULL, 
                         cluster_features_by = NULL, cluster_rows = FALSE, cluster_columns = FALSE, 
                         cluster_row_slices = FALSE, cluster_column_slices = FALSE, 
                         show_row_names = FALSE, show_column_names = FALSE, row_names_side = ifelse(flip, 
                                                                                                    "left", "right"), column_names_side = ifelse(flip, "bottom", 
                                                                                                                                                 "top"), row_names_rot = 0, column_names_rot = 90, row_title = NULL, 
                         column_title = NULL, row_title_side = "left", column_title_side = "top", 
                         row_title_rot = 0, column_title_rot = ifelse(flip, 90, 0), 
                         anno_terms = FALSE, anno_keys = FALSE, anno_features = FALSE, 
                         terms_width = unit(4, "in"), terms_fontsize = 8, keys_width = unit(2, 
                                                                                            "in"), keys_fontsize = c(6, 10), features_width = unit(2, 
                                                                                                                                                   "in"), features_fontsize = c(6, 10), IDtype = "symbol", 
                         species = "Homo_sapiens", db_update = FALSE, db_version = "latest", 
                         db_combine = FALSE, convert_species = FALSE, Ensembl_version = 103, 
                         mirror = NULL, db = "GO_BP", TERM2GENE = NULL, TERM2NAME = NULL, 
                         minGSSize = 10, maxGSSize = 500, GO_simplify = FALSE, GO_simplify_cutoff = "p.adjust < 0.05", 
                         simplify_method = "Wang", simplify_similarityCutoff = 0.7, 
                         pvalueCutoff = NULL, padjustCutoff = 0.05, topTerm = 5, show_termid = FALSE, 
                         topWord = 20, words_excluded = NULL, nlabel = 20, features_label = NULL, 
                         label_size = 10, label_color = "black", add_bg = FALSE, bg_alpha = 0.5, 
                         add_dot = FALSE, dot_size = unit(8, "mm"), add_reticle = FALSE, 
                         reticle_color = "grey", add_violin = FALSE, fill.by = "feature", 
                         fill_palette = "Dark2", fill_palcolor = NULL, heatmap_palette = "RdBu", 
                         heatmap_palcolor = NULL, group_palette = "Paired", group_palcolor = NULL, 
                         cell_split_palette = "simspec", cell_split_palcolor = NULL, 
                         feature_split_palette = "simspec", feature_split_palcolor = NULL, 
                         cell_annotation = NULL, cell_annotation_palette = "Paired", 
                         cell_annotation_palcolor = NULL, cell_annotation_params = if (flip) list(width = unit(10, 
                                                                                                               "mm")) else list(height = unit(10, "mm")), feature_annotation = NULL, 
                         feature_annotation_palette = "Dark2", feature_annotation_palcolor = NULL, 
                         feature_annotation_params = if (flip) list(height = unit(5, 
                                                                                  "mm")) else list(width = unit(5, "mm")), use_raster = NULL, 
                         raster_device = "png", raster_by_magick = FALSE, height = NULL, 
                         width = NULL, units = "inch", seed = 11, ht_params = list()) 
{
  set.seed(seed)
  if (isTRUE(raster_by_magick)) {
    check_R("magick")
  }
  if (is.null(features)) {
    stop("No feature provided.")
  }
  split_method <- match.arg(split_method)
  data_nm <- c(ifelse(isTRUE(lib_normalize), "normalized", 
                      ""), slot)
  data_nm <- paste(data_nm[data_nm != ""], collapse = " ")
  if (length(exp_method) == 1 && is.function(exp_method)) {
    exp_name <- paste0(as.character(x = formals()$exp_method), 
                       "(", data_nm, ")")
  }
  else {
    exp_method <- match.arg(exp_method)
    exp_name <- paste0(exp_method, "(", data_nm, ")")
  }
  if (!is.null(grouping.var) && exp_method != "log2fc") {
    warning("When 'grouping.var' is specified, 'exp_method' can only be 'log2fc'", 
            immediate. = TRUE)
    exp_method <- "log2fc"
  }
  exp_name <- exp_legend_title %||% exp_name
  if (!is.null(grouping.var)) {
    if (identical(split.by, grouping.var)) {
      stop("'grouping.var' must be different from 'split.by'")
    }
    if (!is.factor(srt@meta.data[[grouping.var]])) {
      srt@meta.data[[grouping.var]] <- factor(srt@meta.data[[grouping.var]], 
                                              levels = unique(srt@meta.data[[grouping.var]]))
    }
    if (is.null(numerator)) {
      numerator <- levels(srt@meta.data[[grouping.var]])[1]
      warning("'numerator' is not specified. Use the first level in 'grouping.var': ", 
              numerator, immediate. = TRUE)
    }
    else {
      if (!numerator %in% levels(srt@meta.data[, grouping.var])) {
        stop("'", numerator, "' is not an element of the '", 
             grouping.var, "'")
      }
    }
    srt@meta.data[["grouping.var.use"]] <- srt@meta.data[[grouping.var]] == 
      numerator
    add_dot <- FALSE
    exp_name <- paste0(numerator, "/", "other\n", exp_method, 
                       "(", data_nm, ")")
  }
  assay <- assay %||% DefaultAssay(srt)
  if (is.null(group.by)) {
    srt@meta.data[["All.groups"]] <- factor("")
    group.by <- "All.groups"
  }
  if (any(!group.by %in% colnames(srt@meta.data))) {
    stop(group.by[!group.by %in% colnames(srt@meta.data)], 
         " is not in the meta data of the Seurat object.")
  }
  if (!is.null(group.by)) {
    for (g in group.by) {
      if (!is.factor(srt@meta.data[[g]])) {
        srt@meta.data[[g]] <- factor(srt@meta.data[[g]], 
                                     levels = unique(srt@meta.data[[g]]))
      }
    }
  }
  if (length(split.by) > 1) {
    stop("'split.by' only support one variable.")
  }
  if (any(!split.by %in% colnames(srt@meta.data))) {
    stop(split.by[!split.by %in% colnames(srt@meta.data)], 
         " is not in the meta data of the Seurat object.")
  }
  if (!is.null(split.by)) {
    if (!is.factor(srt@meta.data[[split.by]])) {
      srt@meta.data[[split.by]] <- factor(srt@meta.data[[split.by]], 
                                          levels = unique(srt@meta.data[[split.by]]))
    }
  }
  group_elements <- unlist(lapply(srt@meta.data[, group.by, 
                                                drop = FALSE], function(x) length(unique(x))))
  if (any(group_elements == 1) && exp_method == "zscore") {
    stop("'zscore' cannot be applied to the group(s) consisting of one element: ", 
         paste0(names(group_elements)[group_elements == 1], 
                collapse = ","))
  }
  if (length(group_palette) == 1) {
    group_palette <- rep(group_palette, length(group.by))
  }
  if (length(group_palette) != length(group.by)) {
    stop("'group_palette' must be the same length as 'group.by'")
  }
  group_palette <- setNames(group_palette, nm = group.by)
  raw.group.by <- group.by
  raw.group_palette <- group_palette
  if (isTRUE(within_groups)) {
    new.group.by <- c()
    new.group_palette <- group_palette
    for (g in group.by) {
      groups <- split(colnames(srt), srt[[g, drop = TRUE]])
      new.group_palette[g] <- list(rep(new.group_palette[g], 
                                       length(groups)))
      for (nm in names(groups)) {
        srt[[make.names(nm)]] <- factor(NA, levels = levels(srt[[g, 
                                                                 drop = TRUE]]))
        srt[[make.names(nm)]][colnames(srt) %in% groups[[nm]], 
        ] <- nm
        new.group.by <- c(new.group.by, make.names(nm))
      }
    }
    group.by <- new.group.by
    group_palette <- unlist(new.group_palette)
  }
  if (!is.null(feature_split) && !is.factor(feature_split)) {
    feature_split <- factor(feature_split, levels = unique(feature_split))
  }
  if (length(feature_split) != 0 && length(feature_split) != 
      length(features)) {
    stop("feature_split must be the same length as features")
  }
  if (is.null(feature_split_by)) {
    feature_split_by <- group.by
  }
  if (any(!feature_split_by %in% group.by)) {
    stop("feature_split_by must be a subset of group.by")
  }
  if (!is.null(cell_annotation)) {
    if (length(cell_annotation_palette) == 1) {
      cell_annotation_palette <- rep(cell_annotation_palette, 
                                     length(cell_annotation))
    }
    if (length(cell_annotation_palcolor) == 1) {
      cell_annotation_palcolor <- rep(cell_annotation_palcolor, 
                                      length(cell_annotation))
    }
    npal <- unique(c(length(cell_annotation_palette), length(cell_annotation_palcolor), 
                     length(cell_annotation)))
    if (length(npal[npal != 0]) > 1) {
      stop("cell_annotation_palette and cell_annotation_palcolor must be the same length as cell_annotation")
    }
    if (any(!cell_annotation %in% c(colnames(srt@meta.data), 
                                    rownames(srt@assays[[assay]])))) {
      stop("cell_annotation: ", paste0(cell_annotation[!cell_annotation %in% 
                                                         c(colnames(srt@meta.data), rownames(srt@assays[[assay]]))], 
                                       collapse = ","), " is not in the Seurat object.")
    }
  }
  if (!is.null(feature_annotation)) {
    if (length(feature_annotation_palette) == 1) {
      feature_annotation_palette <- rep(feature_annotation_palette, 
                                        length(feature_annotation))
    }
    if (length(feature_annotation_palcolor) == 1) {
      feature_annotation_palcolor <- rep(feature_annotation_palcolor, 
                                         length(feature_annotation))
    }
    npal <- unique(c(length(feature_annotation_palette), 
                     length(feature_annotation_palcolor), length(feature_annotation)))
    if (length(npal[npal != 0]) > 1) {
      stop("feature_annotation_palette and feature_annotation_palcolor must be the same length as feature_annotation")
    }
    if (any(!feature_annotation %in% colnames(srt@assays[[assay]]@meta.features))) {
      stop("feature_annotation: ", paste0(feature_annotation[!feature_annotation %in% 
                                                               colnames(srt@assays[[assay]]@meta.features)], 
                                          collapse = ","), " is not in the meta data of the ", 
           assay, " assay in the Seurat object.")
    }
  }
  if (length(width) == 1) {
    width <- rep(width, length(group.by))
  }
  if (length(height) == 1) {
    height <- rep(height, length(group.by))
  }
  if (length(width) >= 1) {
    names(width) <- group.by
  }
  if (length(height) >= 1) {
    names(height) <- group.by
  }
  if (isTRUE(flip)) {
    cluster_rows_raw <- cluster_rows
    cluster_columns_raw <- cluster_columns
    cluster_row_slices_raw <- cluster_row_slices
    cluster_column_slices_raw <- cluster_column_slices
    cluster_rows <- cluster_columns_raw
    cluster_columns <- cluster_rows_raw
    cluster_row_slices <- cluster_column_slices_raw
    cluster_column_slices <- cluster_row_slices_raw
  }
  if (is.null(cells)) {
    cells <- colnames(srt@assays[[1]])
  }
  if (all(!cells %in% colnames(srt@assays[[1]]))) {
    stop("No cells found.")
  }
  if (!all(cells %in% colnames(srt@assays[[1]]))) {
    warning("Some cells not found.", immediate. = TRUE)
  }
  cells <- intersect(cells, colnames(srt@assays[[1]]))
  if (is.null(features)) {
    features <- VariableFeatures(srt, assay = assay)
  }
  index <- features %in% c(rownames(srt@assays[[assay]]), colnames(srt@meta.data))
  features <- features[index]
  features_unique <- make.unique(features)
  if (!is.null(feature_split)) {
    feature_split <- feature_split[index]
    names(feature_split) <- features_unique
  }
  cell_groups <- list()
  for (cell_group in group.by) {
    if (!is.factor(srt@meta.data[[cell_group]])) {
      srt@meta.data[[cell_group]] <- factor(srt@meta.data[[cell_group]], 
                                            levels = unique(srt@meta.data[[cell_group]]))
    }
    cell_groups[[cell_group]] <- setNames(srt@meta.data[cells, 
                                                        cell_group], cells)
    cell_groups[[cell_group]] <- na.omit(cell_groups[[cell_group]])
    cell_groups[[cell_group]] <- factor(cell_groups[[cell_group]], 
                                        levels = levels(cell_groups[[cell_group]])[levels(cell_groups[[cell_group]]) %in% 
                                                                                     cell_groups[[cell_group]]])
    if (!is.null(split.by)) {
      if (!is.factor(srt@meta.data[[split.by]])) {
        srt@meta.data[[split.by]] <- factor(srt@meta.data[[split.by]], 
                                            levels = unique(srt@meta.data[[split.by]]))
      }
      levels <- apply(expand.grid(levels(srt@meta.data[[split.by]]), 
                                  levels(cell_groups[[cell_group]])), 1, function(x) paste0(x[2:1], 
                                                                                            collapse = " : "))
      cell_groups[[cell_group]] <- setNames(paste0(cell_groups[[cell_group]][cells], 
                                                   " : ", srt@meta.data[cells, split.by]), cells)
      cell_groups[[cell_group]] <- factor(cell_groups[[cell_group]], 
                                          levels = levels[levels %in% cell_groups[[cell_group]]])
    }
    if (!is.null(grouping.var)) {
      levels <- apply(expand.grid(c("TRUE", "FALSE"), levels(cell_groups[[cell_group]])), 
                      1, function(x) paste0(x[2:1], collapse = " ; "))
      cell_groups[[cell_group]] <- setNames(paste0(cell_groups[[cell_group]][cells], 
                                                   " ; ", srt@meta.data[cells, "grouping.var.use"]), 
                                            cells)
      cell_groups[[cell_group]] <- factor(cell_groups[[cell_group]], 
                                          levels = levels[levels %in% cell_groups[[cell_group]]])
    }
  }
  gene <- features[features %in% rownames(srt@assays[[assay]])]
  gene_unique <- features_unique[features %in% rownames(srt@assays[[assay]])]
  meta <- features[features %in% colnames(srt@meta.data)]
  mat_raw <- as_matrix(rbind(slot(srt@assays[[assay]], slot)[gene, 
                                                             cells, drop = FALSE], t(srt@meta.data[cells, meta, drop = FALSE])))[features, 
                                                                                                                                 , drop = FALSE]
  rownames(mat_raw) <- features_unique
  if (isTRUE(lib_normalize) && min(mat_raw, na.rm = TRUE) >= 
      0) {
    if (!is.null(libsize)) {
      libsize_use <- libsize
    }
    else {
      libsize_use <- colSums(slot(srt@assays[[assay]], 
                                  "counts")[, colnames(mat_raw), drop = FALSE])
      isfloat <- any(libsize_use%%1 != 0, na.rm = TRUE)
      if (isTRUE(isfloat)) {
        libsize_use <- rep(1, length(libsize_use))
        warning("The values in the 'counts' slot are non-integer. Set the library size to 1.", 
                immediate. = TRUE)
        if (!is.null(grouping.var)) {
          exp_name <- paste0(numerator, "/", "other\n", 
                             exp_method, "(", slot, ")")
        }
        else {
          exp_name <- paste0(exp_method, "(", slot, ")")
        }
      }
    }
    mat_raw[gene_unique, ] <- t(t(mat_raw[gene_unique, , 
                                          drop = FALSE])/libsize_use * median(libsize_use))
  }
  mat_raw_list <- list()
  mat_perc_list <- list()
  for (cell_group in names(cell_groups)) {
    mat_tmp <- t(aggregate(t(mat_raw[features_unique, , drop = FALSE]), 
                           by = list(cell_groups[[cell_group]][colnames(mat_raw)]), 
                           FUN = aggregate_fun))
    colnames(mat_tmp) <- mat_tmp[1, , drop = FALSE]
    mat_tmp <- mat_tmp[-1, , drop = FALSE]
    class(mat_tmp) <- "numeric"
    mat_raw_list[[cell_group]] <- mat_tmp
    mat_perc <- t(aggregate(t(mat_raw[features_unique, , 
                                      drop = FALSE]), by = list(cell_groups[[cell_group]][colnames(mat_raw)]), 
                            FUN = function(x) {
                              sum(x > exp_cutoff)/length(x)
                            }))
    colnames(mat_perc) <- mat_perc[1, , drop = FALSE]
    mat_perc <- mat_perc[-1, , drop = FALSE]
    class(mat_perc) <- "numeric"
    if (isTRUE(flip)) {
      mat_perc <- t(mat_perc)
    }
    mat_perc_list[[cell_group]] <- mat_perc
  }
  mat_list <- list()
  for (cell_group in group.by) {
    mat_tmp <- mat_raw_list[[cell_group]]
    if (is.null(grouping.var)) {
      mat_tmp <- matrix_process(mat_tmp, method = exp_method)
      mat_tmp[is.infinite(mat_tmp)] <- max(abs(mat_tmp[!is.infinite(mat_tmp)]), 
                                           na.rm = TRUE) * ifelse(mat_tmp[is.infinite(mat_tmp)] > 
                                                                    0, 1, -1)
      mat_tmp[is.na(mat_tmp)] <- mean(mat_tmp, na.rm = TRUE)
      mat_list[[cell_group]] <- mat_tmp
    }
    else {
      compare_groups <- strsplit(colnames(mat_tmp), " ; ")
      names_keep <- names(which(table(sapply(compare_groups, 
                                             function(x) x[[1]])) == 2))
      group_keep <- which(sapply(compare_groups, function(x) x[[1]] %in% 
                                   names_keep))
      group_TRUE <- intersect(group_keep, which(sapply(compare_groups, 
                                                       function(x) x[[2]]) == "TRUE"))
      group_FALSE <- intersect(group_keep, which(sapply(compare_groups, 
                                                        function(x) x[[2]]) == "FALSE"))
      mat_tmp <- log2(mat_tmp[, group_TRUE]/mat_tmp[, group_FALSE])
      colnames(mat_tmp) <- gsub(" ; .*", "", colnames(mat_tmp))
      mat_tmp[is.infinite(mat_tmp)] <- max(abs(mat_tmp[!is.infinite(mat_tmp)]), 
                                           na.rm = TRUE) * ifelse(mat_tmp[is.infinite(mat_tmp)] > 
                                                                    0, 1, -1)
      mat_tmp[is.na(mat_tmp)] <- 0
      mat_list[[cell_group]] <- mat_tmp
      cell_groups[[cell_group]] <- factor(gsub(" ; .*", 
                                               "", cell_groups[[cell_group]]), levels = unique(gsub(" ; .*", 
                                                                                                    "", levels(cell_groups[[cell_group]]))))
    }
  }
  mat_split <- do.call(cbind, mat_list[feature_split_by])
  if (is.null(limits)) {
    if (!is.function(exp_method) && exp_method %in% c("zscore", 
                                                      "log2fc")) {
      b <- ceiling(min(abs(quantile(do.call(cbind, mat_list), 
                                    c(0.01, 0.99), na.rm = TRUE)), na.rm = TRUE) * 
                     2)/2
      colors <- colorRamp2(seq(-b, b, length = 100), palette_scp(palette = heatmap_palette, 
                                                                 palcolor = heatmap_palcolor))
    }
    else {
      b <- quantile(do.call(cbind, mat_list), c(0.01, 0.99), 
                    na.rm = TRUE)
      colors <- colorRamp2(seq(b[1], b[2], length = 100), 
                           palette_scp(palette = heatmap_palette, palcolor = heatmap_palcolor))
    }
  }
  else {
    colors <- colorRamp2(seq(limits[1], limits[2], length = 100), 
                         palette_scp(palette = heatmap_palette, palcolor = heatmap_palcolor))
  }
  cell_metadata <- cbind.data.frame(data.frame(row.names = cells, 
                                               cells = cells), cbind.data.frame(srt@meta.data[cells, 
                                                                                              c(group.by, intersect(cell_annotation, colnames(srt@meta.data))), 
                                                                                              drop = FALSE], t(srt@assays[[assay]]@data[intersect(cell_annotation, 
                                                                                                                                                  rownames(srt@assays[[assay]])) %||% integer(), cells, 
                                                                                                                                        drop = FALSE])))
  feature_metadata <- cbind.data.frame(data.frame(row.names = features_unique, 
                                                  features = features, features_uique = features_unique), 
                                       srt@assays[[assay]]@meta.features[features, intersect(feature_annotation, 
                                                                                             colnames(srt@assays[[assay]]@meta.features)), drop = FALSE])
  feature_metadata[, "duplicated"] <- feature_metadata[["features"]] %in% 
    features[duplicated(features)]
  lgd <- list()
  lgd[["ht"]] <- Legend(title = exp_name, col_fun = colors, 
                        border = TRUE)
  if (isTRUE(add_dot)) {
    lgd[["point"]] <- Legend(labels = paste0(seq(20, 100, 
                                                 length.out = 5), "%"), title = "Percent", type = "points", 
                             pch = 21, size = dot_size * seq(0.2, 1, length.out = 5), 
                             grid_height = dot_size * seq(0.2, 1, length.out = 5) * 
                               0.8, grid_width = dot_size, legend_gp = gpar(fill = "grey30"), 
                             border = FALSE, background = "transparent", direction = "vertical")
  }
  ha_top_list <- NULL
  cluster_columns_list <- list()
  column_split_list <- list()
  for (i in seq_along(group.by)) {
    cell_group <- group.by[i]
    cluster_columns_list[[cell_group]] <- cluster_columns
    if (is.null(split.by)) {
      column_split_list[[cell_group]] <- NULL
    }
    else {
      column_split_list[[cell_group]] <- factor(gsub(" : .*", 
                                                     "", levels(cell_groups[[cell_group]])), levels = levels(srt@meta.data[[cell_group]]))
    }
    if (isTRUE(cluster_column_slices) && !is.null(split.by)) {
      if (!isTRUE(cluster_columns)) {
        if (nlevels(column_split_list[[cell_group]]) == 
            1) {
          stop("cluster_column_slices=TRUE can not be used when there is only one group.")
        }
        dend <- cluster_within_group(mat_list[[cell_group]], 
                                     column_split_list[[cell_group]])
        cluster_columns_list[[cell_group]] <- dend
        column_split_list[[cell_group]] <- length(unique(column_split_list[[cell_group]]))
      }
    }
    if (cell_group != "All.groups") {
      funbody <- paste0("\n        grid.rect(gp = gpar(fill = palette_scp(", 
                        paste0("c('", paste0(levels(srt@meta.data[[cell_group]]), 
                                             collapse = "','"), "')"), ",palette = '", group_palette[i], 
                        "',palcolor=c(", paste0("'", paste0(group_palcolor[[i]], 
                                                            collapse = "','"), "'"), "))[nm]))\n      ")
      funbody <- gsub(pattern = "\n", replacement = "", 
                      x = funbody)
      eval(parse(text = paste("panel_fun <- function(index, nm) {", 
                              funbody, "}", sep = "")), envir = environment())
      anno <- list()
      anno[[cell_group]] <- anno_block(align_to = split(seq_along(levels(cell_groups[[cell_group]])), 
                                                        gsub(pattern = " : .*", replacement = "", x = levels(cell_groups[[cell_group]]))), 
                                       panel_fun = getFunction("panel_fun", where = environment()), 
                                       which = ifelse(flip, "row", "column"), show_name = FALSE)
      ha_cell_group <- do.call("HeatmapAnnotation", args = c(anno, 
                                                             which = ifelse(flip, "row", "column"), show_annotation_name = show_annotation_name, 
                                                             annotation_name_side = ifelse(flip, "top", "left"), 
                                                             border = TRUE))
      ha_top_list[[cell_group]] <- ha_cell_group
    }
    if (!is.null(split.by)) {
      funbody <- paste0("\n      grid.rect(gp = gpar(fill = palette_scp(", 
                        paste0("c('", paste0(levels(srt@meta.data[[split.by]]), 
                                             collapse = "','"), "')"), ",palette = '", cell_split_palette, 
                        "',palcolor=c(", paste0("'", paste0(unlist(cell_split_palcolor), 
                                                            collapse = "','"), "'"), "))[nm]))\n    ")
      funbody <- gsub(pattern = "\n", replacement = "", 
                      x = funbody)
      eval(parse(text = paste("panel_fun <- function(index, nm) {", 
                              funbody, "}", sep = "")), envir = environment())
      anno <- list()
      anno[[split.by]] <- anno_block(align_to = split(seq_along(levels(cell_groups[[cell_group]])), 
                                                      gsub(pattern = ".* : ", replacement = "", x = levels(cell_groups[[cell_group]]))), 
                                     panel_fun = getFunction("panel_fun", where = environment()), 
                                     which = ifelse(flip, "row", "column"), show_name = i == 
                                       1)
      ha_split_by <- do.call("HeatmapAnnotation", args = c(anno, 
                                                           which = ifelse(flip, "row", "column"), show_annotation_name = show_annotation_name, 
                                                           annotation_name_side = ifelse(flip, "top", "left"), 
                                                           border = TRUE))
      if (is.null(ha_top_list[[cell_group]])) {
        ha_top_list[[cell_group]] <- ha_split_by
      }
      else {
        ha_top_list[[cell_group]] <- c(ha_top_list[[cell_group]], 
                                       ha_split_by)
      }
    }
  }
  for (i in seq_along(raw.group.by)) {
    cell_group <- raw.group.by[i]
    if (cell_group != "All.groups" & legend_group.by) {
      lgd[[cell_group]] <- Legend(title = cell_group, labels = levels(srt@meta.data[[cell_group]]), 
                                  legend_gp = gpar(fill = palette_scp(levels(srt@meta.data[[cell_group]]), 
                                                                      palette = raw.group_palette[i], palcolor = group_palcolor[[i]])), 
                                  border = TRUE)
    }
  }
  if (!is.null(split.by)) {
    lgd[[split.by]] <- Legend(title = split.by, labels = levels(srt@meta.data[[split.by]]), 
                              legend_gp = gpar(fill = palette_scp(levels(srt@meta.data[[split.by]]), 
                                                                  palette = cell_split_palette, palcolor = cell_split_palcolor)), 
                              border = TRUE)
  }
  if (!is.null(cell_annotation)) {
    subplots_list <- list()
    for (i in seq_along(cell_annotation)) {
      cellan <- cell_annotation[i]
      palette <- cell_annotation_palette[i]
      palcolor <- cell_annotation_palcolor[[i]]
      cell_anno <- cell_metadata[, cellan]
      names(cell_anno) <- rownames(cell_metadata)
      if (!is.numeric(cell_anno)) {
        if (is.logical(cell_anno)) {
          cell_anno <- factor(cell_anno, levels = c(TRUE, 
                                                    FALSE))
        }
        else if (!is.factor(cell_anno)) {
          cell_anno <- factor(cell_anno, levels = unique(cell_anno))
        }
        for (cell_group in group.by) {
          subplots <- CellStatPlot(srt, flip = flip, 
                                   cells = names(cell_groups[[cell_group]]), 
                                   plot_type = "pie", stat.by = cellan, group.by = cell_group, 
                                   split.by = split.by, palette = palette, palcolor = palcolor, 
                                   title = NULL, individual = TRUE, combine = FALSE)
          subplots_list[[paste0(cellan, ":", cell_group)]] <- subplots
          graphics <- list()
          for (nm in names(subplots)) {
            funbody <- paste0("\n              g <- as_grob(subplots_list[['", 
                              cellan, ":", cell_group, "']]", "[['", 
                              nm, "']]  + facet_null() + theme_void() + theme(plot.title = element_blank(),strip.text=element_blank(), plot.subtitle = element_blank(), legend.position = 'none'));\n              g$name <- '", 
                              paste0(cellan, ":", cell_group, "-", nm), 
                              "';\n              grid.draw(g)\n              ")
            funbody <- gsub(pattern = "\n", replacement = "", 
                            x = funbody)
            eval(parse(text = paste("graphics[[nm]] <- function(x, y, w, h) {", 
                                    funbody, "}", sep = "")), envir = environment())
          }
          x_nm <- sapply(strsplit(levels(cell_groups[[cell_group]]), 
                                  " : "), function(x) {
                                    if (length(x) == 2) {
                                      paste0(c(cell_group, x[1], x[2]), collapse = ":")
                                    }
                                    else {
                                      paste0(c(cell_group, x[1], ""), collapse = ":")
                                    }
                                  })
          ha_cell <- list()
          ha_cell[[cellan]] <- anno_customize(x = x_nm, 
                                              graphics = graphics, which = ifelse(flip, 
                                                                                  "row", "column"), border = TRUE, verbose = FALSE)
          if(show_annotation_name){
            anno_args <- c(ha_cell, which = ifelse(flip, 
                                                   "row", "column"), show_annotation_name = cell_group == 
                             group.by[1], annotation_name_side = ifelse(flip, 
                                                                        "top", "left"))
          }else{anno_args <- c(ha_cell, which = ifelse(flip, 
                                                       "row", "column"), show_annotation_name = FALSE)}    
          anno_args <- c(anno_args, cell_annotation_params[setdiff(names(cell_annotation_params), 
                                                                   names(anno_args))])
          ha_top <- do.call(HeatmapAnnotation, args = anno_args)
          if (is.null(ha_top_list[[cell_group]])) {
            ha_top_list[[cell_group]] <- ha_top
          }
          else {
            ha_top_list[[cell_group]] <- c(ha_top_list[[cell_group]], 
                                           ha_top)
          }
        }
        lgd[[cellan]] <- Legend(title = cellan, labels = levels(cell_anno), 
                                legend_gp = gpar(fill = palette_scp(cell_anno, 
                                                                    palette = palette, palcolor = palcolor)), 
                                border = TRUE)
      }
      else {
        for (cell_group in group.by) {
          subplots <- FeatureStatPlot(srt, assay = assay, 
                                      slot = "data", flip = flip, stat.by = cellan, 
                                      cells = names(cell_groups[[cell_group]]), 
                                      group.by = cell_group, split.by = split.by, 
                                      palette = palette, palcolor = palcolor, fill.by = "group", 
                                      same.y.lims = TRUE, individual = TRUE, combine = FALSE)
          subplots_list[[paste0(cellan, ":", cell_group)]] <- subplots
          graphics <- list()
          for (nm in names(subplots)) {
            funbody <- paste0("\n              g <- as_grob(subplots_list[['", 
                              cellan, ":", cell_group, "']]", "[['", 
                              nm, "']]  + facet_null() + theme_void() + theme(plot.title = element_blank(), strip.text=element_blank(),plot.subtitle = element_blank(), legend.position = 'none'));\n              g$name <- '", 
                              paste0(cellan, ":", cell_group, "-", nm), 
                              "';\n              grid.draw(g)\n              ")
            funbody <- gsub(pattern = "\n", replacement = "", 
                            x = funbody)
            eval(parse(text = paste("graphics[[nm]] <- function(x, y, w, h) {", 
                                    funbody, "}", sep = "")), envir = environment())
          }
          x_nm <- sapply(strsplit(levels(cell_groups[[cell_group]]), 
                                  " : "), function(x) {
                                    if (length(x) == 2) {
                                      paste0(c(cellan, cell_group, x[1], x[2]), 
                                             collapse = ":")
                                    }
                                    else {
                                      paste0(c(cellan, cell_group, x[1], ""), 
                                             collapse = ":")
                                    }
                                  })
          ha_cell <- list()
          ha_cell[[cellan]] <- anno_customize(x = x_nm, 
                                              graphics = graphics, which = ifelse(flip, 
                                                                                  "row", "column"), border = TRUE, verbose = FALSE)
          if(show_annotation_name){
            anno_args <- c(ha_cell, which = ifelse(flip, 
                                                   "row", "column"), show_annotation_name = cell_group == 
                             group.by[1], annotation_name_side = ifelse(flip, 
                                                                        "top", "left"))
          }else{
            anno_args <- c(ha_cell, which = ifelse(flip, 
                                                   "row", "column"), show_annotation_name = FALSE)
          }
          anno_args <- c(anno_args, cell_annotation_params[setdiff(names(cell_annotation_params), 
                                                                   names(anno_args))])
          ha_top <- do.call(HeatmapAnnotation, args = anno_args)
          if (is.null(ha_top_list[[cell_group]])) {
            ha_top_list[[cell_group]] <- ha_top
          }
          else {
            ha_top_list[[cell_group]] <- c(ha_top_list[[cell_group]], 
                                           ha_top)
          }
        }
      }
    }
  }
  if (is.null(feature_split)) {
    if (is.null(n_split) || isTRUE(nrow(mat_split) <= n_split)) {
      row_split_raw <- row_split <- feature_split <- NULL
    }
    else {
      if (n_split == 1) {
        row_split_raw <- row_split <- feature_split <- setNames(rep(1, 
                                                                    nrow(mat_split)), rownames(mat_split))
      }
      else {
        if (split_method == "mfuzz") {
          status <- tryCatch(check_R("e1071"), error = identity)
          if (inherits(status, "error")) {
            warning("The e1071 package was not found. Switch split_method to 'kmeans'", 
                    immediate. = TRUE)
            split_method <- "kmeans"
          }
          else {
            mat_split_tmp <- mat_split
            colnames(mat_split_tmp) <- make.unique(colnames(mat_split_tmp))
            mat_split_tmp <- standardise(mat_split_tmp)
            min_fuzzification <- mestimate(mat_split_tmp)
            if (is.null(fuzzification)) {
              fuzzification <- min_fuzzification + 0.1
            }
            else {
              if (fuzzification <= min_fuzzification) {
                warning("fuzzification value is samller than estimated:", 
                        round(min_fuzzification, 2), immediate. = TRUE)
              }
            }
            cl <- e1071::cmeans(mat_split_tmp, centers = n_split, 
                                method = "cmeans", m = fuzzification)
            if (length(cl$cluster) == 0) {
              stop("Clustering with mfuzz failed (fuzzification=", 
                   round(fuzzification, 2), "). Please set a larger fuzzification parameter manually.")
            }
            row_split <- feature_split <- cl$cluster
          }
        }
        if (split_method == "kmeans") {
          km <- kmeans(mat_split, centers = n_split, 
                       iter.max = 10000, nstart = 20)
          row_split <- feature_split <- km$cluster
        }
        if (split_method == "hclust") {
          hc <- hclust(as.dist(dist(mat_split)))
          row_split <- feature_split <- cutree(hc, k = n_split)
        }
      }
      groupmean <- aggregate(t(mat_split), by = list(unlist(lapply(cell_groups[feature_split_by], 
                                                                   levels))), mean)
      maxgroup <- groupmean[, 1][apply(groupmean[, names(row_split)], 
                                       2, which.max)]
      maxgroup <- factor(maxgroup, levels = levels(unlist(cell_groups[feature_split_by])))
      df <- data.frame(row_split = row_split, order_by = maxgroup)
      df_order <- aggregate(df[["order_by"]], by = list(df[, 
                                                           "row_split"]), FUN = function(x) names(sort(table(x), 
                                                                                                       decreasing = TRUE))[1])
      df_order[, "row_split"] <- df_order[, "Group.1"]
      df_order[["order_by"]] <- as.numeric(factor(df_order[["x"]], 
                                                  levels = levels(maxgroup)))
      df_order <- df_order[order(df_order[["order_by"]], 
                                 decreasing = decreasing), , drop = FALSE]
      if (!is.null(split_order)) {
        df_order <- df_order[split_order, , drop = FALSE]
      }
      split_levels <- c()
      for (i in seq_len(nrow(df_order))) {
        raw_nm <- df_order[i, "row_split"]
        feature_split[feature_split == raw_nm] <- paste0("C", 
                                                         i)
        level <- paste0("C", i, "(", sum(row_split == 
                                           raw_nm), ")")
        row_split[row_split == raw_nm] <- level
        split_levels <- c(split_levels, level)
      }
      row_split_raw <- row_split <- factor(row_split, levels = split_levels)
      feature_split <- factor(feature_split, levels = paste0("C", 
                                                             seq_len(nrow(df_order))))
    }
  }
  else {
    row_split_raw <- row_split <- feature_split <- feature_split[row.names(mat_split)]
  }
  if (!is.null(feature_split)) {
    feature_metadata[["feature_split"]] <- feature_split
  }
  else {
    feature_metadata[["feature_split"]] <- NA
  }
  ha_left <- NULL
  if (!is.null(row_split)) {
    if (isTRUE(cluster_row_slices)) {
      if (!isTRUE(cluster_rows)) {
        dend <- cluster_within_group(t(mat_split), row_split_raw)
        cluster_rows <- dend
        row_split <- length(unique(row_split_raw))
      }
    }
    funbody <- paste0("\n      grid.rect(gp = gpar(fill = palette_scp(", 
                      paste0("c('", paste0(levels(row_split_raw), collapse = "','"), 
                             "')"), ",palette = '", feature_split_palette, 
                      "',palcolor=c(", paste0("'", paste0(unlist(feature_split_palcolor), 
                                                          collapse = "','"), "'"), "))[nm]))\n    ")
    funbody <- gsub(pattern = "\n", replacement = "", x = funbody)
    eval(parse(text = paste("panel_fun <- function(index, nm) {", 
                            funbody, "}", sep = "")), envir = environment())
    ha_clusters <- HeatmapAnnotation(features_split = anno_block(align_to = split(seq_along(row_split_raw), 
                                                                                  row_split_raw), panel_fun = getFunction("panel_fun", 
                                                                                                                          where = environment()), width = unit(0.1, "in"), 
                                                                 height = unit(0.1, "in"), show_name = FALSE, which = ifelse(flip, 
                                                                                                                             "column", "row")), which = ifelse(flip, "column", 
                                                                                                                                                               "row"), border = TRUE)
    if (is.null(ha_left)) {
      ha_left <- ha_clusters
    }
    else {
      ha_left <- c(ha_left, ha_clusters)
    }
    lgd[["Cluster"]] <- Legend(title = "Cluster", labels = intersect(levels(row_split_raw), 
                                                                     row_split_raw), legend_gp = gpar(fill = palette_scp(intersect(levels(row_split_raw), 
                                                                                                                                   row_split_raw), type = "discrete", palette = feature_split_palette, 
                                                                                                                         palcolor = feature_split_palcolor, matched = TRUE)), 
                               border = TRUE)
  }
  if (isTRUE(cluster_rows) && !is.null(cluster_features_by)) {
    mat_cluster <- do.call(cbind, mat_list[cluster_features_by])
    if (is.null(row_split)) {
      dend <- as.dendrogram(hclust(as.dist(dist(mat_cluster))))
      dend_ordered <- reorder(dend, wts = colMeans(mat_cluster), 
                              agglo.FUN = mean)
      cluster_rows <- dend_ordered
    }
    else {
      row_split <- length(unique(row_split_raw))
      dend <- cluster_within_group2(t(mat_cluster), row_split_raw)
      cluster_rows <- dend
    }
  }
  cell_group <- group.by[1]
  ht_args <- list(matrix = mat_list[[cell_group]], col = colors, 
                  row_split = row_split, column_split = column_split_list[[cell_group]], 
                  cluster_rows = cluster_rows, cluster_columns = cluster_columns_list[[cell_group]], 
                  cluster_row_slices = cluster_row_slices, cluster_column_slices = cluster_column_slices, 
                  use_raster = TRUE)
  ht_args <- c(ht_args, ht_params[setdiff(names(ht_params), 
                                          names(ht_args))])
  ht_list <- do.call(Heatmap, args = ht_args)
  features_ordered <- rownames(mat_list[[1]])[unlist(suppressWarnings(row_order(ht_list)))]
  feature_metadata[["index"]] <- setNames(object = seq_along(features_ordered), 
                                          nm = features_ordered)[rownames(feature_metadata)]
  if (is.null(features_label)) {
    if (nlabel > 0) {
      if (length(features) > nlabel) {
        index_from <- ceiling((length(features_ordered)/nlabel)/2)
        index_to <- length(features_ordered)
        index <- unique(round(seq(from = index_from, 
                                  to = index_to, length.out = nlabel)))
      }
      else {
        index <- seq_along(features_ordered)
      }
    }
    else {
      index <- NULL
    }
  }
  else {
    index <- which(features_ordered %in% features_label)
    drop <- setdiff(features_label, features_ordered)
    if (length(drop) > 0) {
      warning(paste0(paste0(drop, collapse = ","), "was not found in the features"), 
              immediate. = TRUE)
    }
  }
  if (length(index) > 0) {
    ha_mark <- HeatmapAnnotation(gene = anno_mark(at = which(rownames(feature_metadata) %in% 
                                                               features_ordered[index]), labels = feature_metadata[which(rownames(feature_metadata) %in% 
                                                                                                                           features_ordered[index]), "features"], side = ifelse(flip, 
                                                                                                                                                                                "top", "left"), labels_gp = gpar(fontsize = label_size, 
                                                                                                                                                                                                                 col = label_color), link_gp = gpar(fontsize = label_size, 
                                                                                                                                                                                                                                                    col = label_color), which = ifelse(flip, "column", 
                                                                                                                                                                                                                                                                                       "row")), which = ifelse(flip, "column", "row"), show_annotation_name = show_annotation_name)
    if (is.null(ha_left)) {
      ha_left <- ha_mark
    }
    else {
      ha_left <- c(ha_mark, ha_left)
    }
  }
  ha_right <- NULL
  if (!is.null(feature_annotation)) {
    for (i in seq_along(feature_annotation)) {
      featan <- feature_annotation[i]
      palette <- feature_annotation_palette[i]
      palcolor <- feature_annotation_palcolor[[i]]
      featan_values <- feature_metadata[, featan]
      if (!is.numeric(featan_values)) {
        if (is.logical(featan_values)) {
          featan_values <- factor(featan_values, levels = c(TRUE, 
                                                            FALSE))
        }
        else if (!is.factor(featan_values)) {
          featan_values <- factor(featan_values, levels = unique(featan_values))
        }
        ha_feature <- list()
        ha_feature[[featan]] <- anno_simple(x = as.character(featan_values), 
                                            col = palette_scp(featan_values, palette = palette, 
                                                              palcolor = palcolor), na_col = "transparent", 
                                            which = ifelse(flip, "column", "row"), border = TRUE)
        anno_args <- c(ha_feature, which = ifelse(flip, 
                                                  "column", "row"), show_annotation_name = show_annotation_name, 
                       annotation_name_side = ifelse(flip, "left", 
                                                     "top"), border = TRUE)
        anno_args <- c(anno_args, feature_annotation_params[setdiff(names(feature_annotation_params), 
                                                                    names(anno_args))])
        ha_feature <- do.call(HeatmapAnnotation, args = anno_args)
        if (is.null(ha_right)) {
          ha_right <- ha_feature
        }
        else {
          ha_right <- c(ha_right, ha_feature)
        }
        lgd[[featan]] <- Legend(title = featan, labels = levels(featan_values), 
                                legend_gp = gpar(fill = palette_scp(featan_values, 
                                                                    palette = palette, palcolor = palcolor)), 
                                border = TRUE)
      }
      else {
        col_fun <- colorRamp2(breaks = seq(min(featan_values, 
                                               na.rm = TRUE), max(featan_values, na.rm = TRUE), 
                                           length = 100), colors = palette_scp(palette = palette, 
                                                                               palcolor = palcolor))
        ha_feature <- list()
        ha_feature[[featan]] <- anno_simple(x = featan_values, 
                                            col = col_fun, na_col = "transparent", which = ifelse(flip, 
                                                                                                  "column", "row"), border = TRUE)
        anno_args <- c(ha_feature, which = ifelse(flip, 
                                                  "column", "row"), show_annotation_name = show_annotation_name, 
                       annotation_name_side = ifelse(flip, "left", 
                                                     "top"), border = TRUE)
        anno_args <- c(anno_args, feature_annotation_params[setdiff(names(feature_annotation_params), 
                                                                    names(anno_args))])
        ha_feature <- do.call(HeatmapAnnotation, args = anno_args)
        if (is.null(ha_right)) {
          ha_right <- ha_feature
        }
        else {
          ha_right <- c(ha_right, ha_feature)
        }
        lgd[[featan]] <- Legend(title = featan, col_fun = col_fun, 
                                border = TRUE)
      }
    }
  }
  enrichment <- heatmap_enrichment(geneID = feature_metadata[["features"]], 
                                   geneID_groups = feature_metadata[["feature_split"]], 
                                   feature_split_palette = feature_split_palette, feature_split_palcolor = feature_split_palcolor, 
                                   ha_right = ha_right, flip = flip, anno_terms = anno_terms, 
                                   anno_keys = anno_keys, anno_features = anno_features, 
                                   terms_width = terms_width, terms_fontsize = terms_fontsize, 
                                   keys_width = keys_width, keys_fontsize = keys_fontsize, 
                                   features_width = features_width, features_fontsize = features_fontsize, 
                                   IDtype = IDtype, species = species, db_update = db_update, 
                                   db_version = db_version, db_combine = db_combine, convert_species = convert_species, 
                                   Ensembl_version = Ensembl_version, mirror = mirror, db = db, 
                                   TERM2GENE = TERM2GENE, TERM2NAME = TERM2NAME, minGSSize = minGSSize, 
                                   maxGSSize = maxGSSize, GO_simplify = GO_simplify, GO_simplify_cutoff = GO_simplify_cutoff, 
                                   simplify_method = simplify_method, simplify_similarityCutoff = simplify_similarityCutoff, 
                                   pvalueCutoff = pvalueCutoff, padjustCutoff = padjustCutoff, 
                                   topTerm = topTerm, show_termid = show_termid, topWord = topWord, 
                                   words_excluded = words_excluded)
  res <- enrichment$res
  ha_right <- enrichment$ha_right
  ht_list <- NULL
  vlnplots_list <- NULL
  x_nm_list <- NULL
  y_nm_list <- NULL
  if (fill.by == "group") {
    palette <- group_palette
    palcolor <- group_palcolor
  }
  else {
    palette <- feature_annotation_palette
    palcolor <- feature_annotation_palcolor
  }
  for (cell_group in group.by) {
    if (cell_group == group.by[1]) {
      left_annotation <- ha_left
    }
    else {
      left_annotation <- NULL
    }
    if (cell_group == group.by[length(group.by)]) {
      right_annotation <- ha_right
    }
    else {
      right_annotation <- NULL
    }
    if (isTRUE(add_violin)) {
      vlnplots <- FeatureStatPlot(srt, assay = assay, slot = "data", 
                                  flip = flip, stat.by = rownames(mat_list[[cell_group]]), 
                                  cells = names(cell_groups[[cell_group]]), group.by = cell_group, 
                                  split.by = split.by, palette = fill_palette, 
                                  palcolor = fill_palcolor, fill.by = fill.by, 
                                  same.y.lims = TRUE, individual = TRUE, combine = FALSE)
      lgd[["ht"]] <- NULL
      for (nm in names(vlnplots)) {
        gtable <- as_grob(vlnplots[[nm]] + facet_null() + 
                            theme_void() + theme(legend.position = "none",strip.text=element_blank()))
        gtable$name <- paste0(cell_group, "-", nm)
        vlnplots[[nm]] <- gtable
      }
      vlnplots_list[[paste0("heatmap_group:", cell_group)]] <- vlnplots
      x_nm <- rownames(mat_list[[cell_group]])
      x_nm_list[[paste0("heatmap_group:", cell_group)]] <- x_nm
      y_nm <- sapply(strsplit(levels(cell_groups[[cell_group]]), 
                              " : "), function(x) {
                                if (length(x) == 2) {
                                  paste0(c(cell_group, x[1], x[2]), collapse = ":")
                                }
                                else {
                                  paste0(c(cell_group, x[1], ""), collapse = ":")
                                }
                              })
      y_nm_list[[paste0("heatmap_group:", cell_group)]] <- y_nm
    }
    funbody <- paste0(if (isTRUE(add_dot) || isTRUE(add_violin)) {
      "grid.rect(x, y,\n          width = width, height = height,\n          gp = gpar(col = 'white', lwd = 1, fill = 'white')\n        );"
    }, if (isTRUE(add_bg)) {
      paste0("\n        grid.rect(x, y,\n          width = width, height = height,\n          gp = gpar(col = fill, lwd = 1, fill = adjcolors(fill, ", 
             bg_alpha, "))\n        );\n        ")
    }, if (isTRUE(add_reticle)) {
      paste0("\n        ind_mat = restore_matrix(j, i, x, y);\n        ind_top = ind_mat[1,];\n        ind_left = ind_mat[,1];\n        for(col in seq_len(ncol(ind_mat))){\n          grid.lines(x = unit(rep(x[ind_top[col]],each=2),'npc'),y = unit(c(0,1),'npc'),gp = gpar(col = '", 
             reticle_color, "', lwd = 1.5));\n        };\n        for(row in seq_len(nrow(ind_mat))){\n          grid.lines(x = unit(c(0,1),'npc'),y = unit(rep(y[ind_left[row]],each=2),'npc'),gp = gpar(col = '", 
             reticle_color, "', lwd = 1.5));\n        };\n        ")
    }, if (isTRUE(add_dot)) {
      paste0("perc <- pindex(mat_perc_list[['", cell_group, 
             "']]", ", i, j);\n        grid.points(x, y,\n          pch = 21,\n          size = dot_size*perc,\n          gp = gpar(col = 'black', lwd = 1, fill = fill)\n        );\n        ")
    }, if (isTRUE(add_violin)) {
      if (isTRUE(flip)) {
        paste0("\n        groblist <- extractgrobs(vlnplots = vlnplots_list[[paste0('heatmap_group:', '", 
               cell_group, "')]],\n               x_nm =  x_nm_list[[paste0('heatmap_group:', '", 
               cell_group, "')]],\n               y_nm= y_nm_list[[paste0('heatmap_group:', '", 
               cell_group, "')]],\n               x = j,y = i);\n        grid_draw(groblist, x = x, y = y, width = width, height = height);\n        ")
      }
      else {
        paste0("\n        groblist <- extractgrobs(vlnplots = vlnplots_list[[paste0('heatmap_group:', '", 
               cell_group, "')]],\n               x_nm =  x_nm_list[[paste0('heatmap_group:', '", 
               cell_group, "')]],\n               y_nm= y_nm_list[[paste0('heatmap_group:', '", 
               cell_group, "')]],\n               x = i,y = j);\n        grid_draw(groblist, x = x, y = y, width = width, height = height);\n        ")
      }
    })
    funbody <- gsub(pattern = "\n", replacement = "", x = funbody)
    eval(parse(text = paste("layer_fun <- function(j, i, x, y, width, height, fill) {", 
                            funbody, "}", sep = "")), envir = environment())
    ht_args <- list(name = cell_group, matrix = if (flip) t(mat_list[[cell_group]]) else mat_list[[cell_group]], 
                    col = colors, layer_fun = getFunction("layer_fun", 
                                                          where = environment()), row_title = row_title %||% 
                      if (flip) ifelse(cell_group != "All.groups", 
                                       cell_group, "") else character(0), row_title_side = row_title_side, 
                    column_title = column_title %||% if (flip) character(0) else ifelse(cell_group != 
                                                                                          "All.groups", cell_group, ""), column_title_side = column_title_side, 
                    row_title_rot = row_title_rot, column_title_rot = column_title_rot, 
                    row_split = if (flip) column_split_list[[cell_group]] else row_split, 
                    column_split = if (flip) row_split else column_split_list[[cell_group]], 
                    cluster_rows = if (flip) cluster_columns_list[[cell_group]] else cluster_rows, 
                    cluster_columns = if (flip) cluster_rows else cluster_columns_list[[cell_group]], 
                    cluster_row_slices = if (flip) cluster_column_slices else cluster_row_slices, 
                    cluster_column_slices = if (flip) cluster_row_slices else cluster_column_slices, 
                    show_row_names = show_row_names, show_column_names = show_column_names, 
                    row_names_side = row_names_side, column_names_side = column_names_side, 
                    row_names_rot = row_names_rot, column_names_rot = column_names_rot, 
                    top_annotation = if (flip) left_annotation else ha_top_list[[cell_group]], 
                    left_annotation = if (flip) ha_top_list[[cell_group]] else left_annotation, 
                    bottom_annotation = if (flip) right_annotation else NULL, 
                    right_annotation = if (flip) NULL else right_annotation, 
                    show_heatmap_legend = FALSE, border = border, use_raster = use_raster, 
                    raster_device = raster_device, raster_by_magick = raster_by_magick, 
                    width = if (is.numeric(width[cell_group])) unit(width[cell_group], 
                                                                    units = units) else NULL, height = if (is.numeric(height[cell_group])) unit(height[cell_group], 
                                                                                                                                                units = units) else NULL)
    if (any(names(ht_params) %in% names(ht_args))) {
      warning("ht_params: ", paste0(intersect(names(ht_params), 
                                              names(ht_args)), collapse = ","), " were duplicated and will not be used.", 
              immediate. = TRUE)
    }
    ht_args <- c(ht_args, ht_params[setdiff(names(ht_params), 
                                            names(ht_args))])
    if (isTRUE(flip)) {
      if (is.null(ht_list)) {
        ht_list <- do.call(Heatmap, args = ht_args)
      }
      else {
        ht_list <- ht_list %v% do.call(Heatmap, args = ht_args)
      }
    }
    else {
      ht_list <- ht_list + do.call(Heatmap, args = ht_args)
    }
  }
  if ((!is.null(row_split) && length(index) > 0) || any(c(anno_terms, 
                                                          anno_keys, anno_features)) || !is.null(width) || !is.null(height)) {
    fix <- TRUE
    if (is.null(width) || is.null(height)) {
      message("\nThe size of the heatmap is fixed because certain elements are not scalable.\nThe width and height of the heatmap are determined by the size of the current viewport.\nIf you want to have more control over the size, you can manually set the parameters 'width' and 'height'.\n")
    }
  }
  else {
    fix <- FALSE
  }
  rendersize <- heatmap_rendersize(width = width, height = height, 
                                   units = units, ha_top_list = ha_top_list, ha_left = ha_left, 
                                   ha_right = ha_right, ht_list = ht_list, legend_list = lgd, 
                                   flip = flip)
  width_sum <- rendersize[["width_sum"]]
  height_sum <- rendersize[["height_sum"]]
  if (isTRUE(fix)) {
    fixsize <- heatmap_fixsize(width = width, width_sum = width_sum, 
                               height = height, height_sum = height_sum, units = units, 
                               ht_list = ht_list, legend_list = lgd)
    ht_width <- fixsize[["ht_width"]]
    ht_height <- fixsize[["ht_height"]]
    gTree <- grid.grabExpr({
      draw(ht_list, annotation_legend_list = lgd)
      for (enrich in db) {
        enrich_anno <- names(ha_right)[grep(paste0("_split_", 
                                                   enrich), names(ha_right))]
        if (length(enrich_anno) > 0) {
          for (enrich_anno_element in enrich_anno) {
            enrich_obj <- strsplit(enrich_anno_element, 
                                   "_split_")[[1]][1]
            decorate_annotation(enrich_anno_element, 
                                slice = 1, {
                                  grid.text(paste0(enrich, " (", enrich_obj, 
                                                   ")"), x = unit(1, "npc"), y = unit(1, 
                                                                                      "npc") + unit(2.5, "mm"), just = c("left", 
                                                                                                                         "bottom"))
                                })
          }
        }
      }
    }, width = ht_width, height = ht_height, wrap = TRUE, 
    wrap.grobs = TRUE)
  }
  else {
    ht_width <- unit(width_sum, units = units)
    ht_height <- unit(height_sum, units = units)
    if(legend.horizontal){
      lgd<-packLegend(list=lgd,direction="horizontal")
    }
    gTree <- grid.grabExpr({
      draw(ht_list, annotation_legend_list = lgd)
      for (enrich in db) {
        enrich_anno <- names(ha_right)[grep(paste0("_split_", 
                                                   enrich), names(ha_right))]
        if (length(enrich_anno) > 0) {
          for (enrich_anno_element in enrich_anno) {
            enrich_obj <- strsplit(enrich_anno_element, 
                                   "_split_")[[1]][1]
            decorate_annotation(enrich_anno_element, 
                                slice = 1, {
                                  grid.text(paste0(enrich, " (", enrich_obj, 
                                                   ")"), x = unit(1, "npc"), y = unit(1, 
                                                                                      "npc") + unit(2.5, "mm"), just = c("left", 
                                                                                                                         "bottom"))
                                })
          }
        }
      }
    }, width = ht_width, height = ht_height, wrap = TRUE, 
    wrap.grobs = TRUE)
  }
  if (isTRUE(fix)) {
    p <- panel_fix_overall(gTree, width = as.numeric(ht_width), 
                           height = as.numeric(ht_height), units = units)
  }
  else {
    p <- wrap_plots(gTree)
  }
  return(list(plot = p, matrix_list = mat_list, feature_split = feature_split, 
              cell_metadata = cell_metadata, feature_metadata = feature_metadata, 
              enrichment = res))
}
environment(GroupHeatmapy)<-environment(GroupHeatmap)
