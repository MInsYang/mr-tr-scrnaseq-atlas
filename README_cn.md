# 单细胞分析代码库

本代码库包含单细胞RNA测序数据分析

## 文件结构

```
scripts/
├── README.md                        # 本文件
├── README_cn.md                     # 中文说明
├── utils.R                          # 自定义函数库
├── 01_data_process.r                # 数据读取和预处理
├── 02_HarmonyUMAP.r                 # Harmony整合、UMAP降维和细胞类型注释
├── 03.1_subtype_VIC.r              # VIC细胞亚群分析
├── 03.2_subtype_VEC.r              # VEC细胞亚群分析（含Monocle2、CytoTRACE2）
├── 03.3_subtype_immuneCells.r      # 免疫细胞亚群分析（淋巴细胞、髓系细胞）
├── 04_fibrotic_transitional_analysis.r  # 纤维化和过渡态分析
├── 05_geneset_score.r              # 基因集打分（细胞周期、EMT等）
├── 06_monocle2_analysis.r          # Monocle2轨迹分析
├── 07_velocity.r                   # RNA velocity分析
├── 08_gene_corr_analysisr.r        # 基因相关性分析
├── 09_CCI.r                        # 细胞通讯分析（CellPhoneDB）
└── 10_supplemental_analysis.R      # 补充分析
```

## 使用说明

### 1. 环境准备

确保已安装以下R包（根据实际需要）：

- Seurat
- SCP
- DoubletFinder
- Harmony
- monocle/monocle2
- velocyto.R/SCVELO
- CytoTRACE2
- CellPhoneDB
- CellChat
- openxlsx
- dplyr
- data.table
- ComplexHeatmap
- ClusterGVis
- ggplot2

## 文件说明

### 01_data_process.r

**功能：** 数据读取和预处理

- 构建样本信息表
- 读取10X CellRanger输出数据
- 创建Seurat对象并合并多个样本
- 计算线粒体和血红蛋白基因比例
- 使用DoubletFinder进行双细胞鉴定和过滤
- QC统计和保存


### 02_HarmonyUMAP.r

**功能：** Harmony整合、降维、聚类和细胞类型注释

- Harmony去批次效应整合
- PCA降维和UMAP可视化
- Leiden聚类算法（多种分辨率）
- 差异表达分析
- 基于标记基因的细胞类型注释
- 细胞类型统计和可视化


### 03.1_subtype_VIC.r

**功能：** VIC（瓣膜间质细胞）亚群分析

- VIC细胞提取和重新整合
- 亚群聚类和分辨率选择（clustree分析）
- 亚群差异表达分析
- 标记基因识别和可视化
- 亚群在不同条件下的比例分析
- 热图和特征图展示


### 03.2_subtype_VEC.r

**功能：** VEC（瓣膜内皮细胞）亚群分析

- VEC亚群处理
- 下采样处理用于轨迹分析
- Monocle2伪时间轨迹分析
- CytoTRACE2细胞潜能评估
- 细胞分化状态可视化


### 03.3_subtype_immuneCells.r

**功能：** 免疫细胞亚群分析

- 淋巴细胞（Lymphocyte）和髓系细胞（Myeloid cell）分别分析
- Harmony整合和聚类
- 多种分辨率的聚类结果比较
- 亚群差异表达分析
- 标记基因识别
- 富集分析（GO/KEGG等）


### 04_fibrotic_transitional_analysis.r

**功能：** 纤维化和过渡态分析

- Pro-fibrotic vs Anti-fibrotic VIC差异分析
- 不同分组条件下的纤维化相关基因分析
- VIC和VEC过渡态细胞识别（基于ModuleScore）
- 过渡态细胞在UMAP空间中的分布
- 过渡态细胞比例统计


### 05_geneset_score.r

**功能：** 基因集打分

- 细胞周期评分（G1/S期和G2/M期）
- EMT（上皮-间质转化）评分
- 其他基因集模块评分（AddModuleScore）
- 评分在UMAP空间中的可视化
- 评分在不同条件下的比较分析


### 06_monocle2_analysis.r

**功能：** Monocle2轨迹分析

- Monocle2伪时间轨迹构建
- 轨迹分支点识别
- BEAM（Branched Expression Analysis Modeling）分析
- 不同状态在伪时间上的分布
- 轨迹可视化
- 伪时间密度分布分析


### 07_velocity.r

**功能：** RNA velocity分析

- velocyto.R环境配置
- Loom文件读取和准备
- 剪切和未剪切转录本矩阵处理
- RunSCVELO分析
- 速度向量在UMAP空间中的可视化
- 流场图（stream plot）展示


### 08_gene_corr_analysisr.r

**功能：** 基因相关性分析

- 关注基因的提取（如VCAM1、PECAM1、ACTA2等）
- 基因在不同组间的表达量比较
- 组间统计检验（t检验、Wilcoxon检验）
- 结果汇总和导出


### 09_CCI.r

**功能：** 细胞通讯分析（CellPhoneDB）

- CellPhoneDB输入文件准备（表达矩阵和元数据）
- 分别处理不同组织（Tri、Mit）的数据
- CellPhoneDB结果读取和处理
- 配体-受体相互作用分析
- 信号通路富集分析
- 结果可视化


### 10_supplemental_analysis.R

**功能：** 补充分析

- 样本分类（AF vs SR，FMR vs FTR）
- 试验条件排除分析
- 细胞比例统计分析
- ANOVA分析
- 数据重新分组和可视化
- 其他补充性分析


##  注意事项

1. **路径修改**：所有文件路径都是示例路径，请根据实际情况修改
2. **依赖包**：某些分析需要特定的R包和环境，请确保已正确安装
   - velocyto.R需要特定的conda环境和系统库配置
   - CellPhoneDB需要Python环境和命令行工具
3. **数据文件**：确保所需的数据文件（.rds文件、loom文件等）存在于指定路径
4. **内存需求**：单细胞分析通常需要较大内存，请确保系统资源充足
5. **运行顺序**：建议按照编号顺序运行脚本，因为后续脚本依赖前面脚本的输出
6. **文件命名**：注意`08_gene_corr_analysisr.r`文件名中有拼写错误（analysisr应为analysis）

## 自定义函数

所有自定义函数都在 `utils.R` 中定义，包括：

- `qcstat()` - QC统计函数
- `ratio_plot()` - 比例图绘制
- `run_monocle()` - Monocle轨迹分析封装函数
- `prepareLoom()` - Loom文件准备函数
- `run_cluster_anova()` - 聚类ANOVA分析
- 其他辅助函数

##  版本历史

- 初始版本

如有问题或建议，请联系代码维护者。