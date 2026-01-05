# Single-Cell Analysis Code Repository

This repository contains analysis scripts for single-cell RNA-seq (scRNA-seq) data.

## Repository Structure

```
scripts/
├── README.md                         # This file
├── README_cn.md                      # Chinese version
├── utils.R                           # Custom utility functions
├── 01_data_process.r                 # Data loading and preprocessing
├── 02_HarmonyUMAP.r                  # Harmony integration, UMAP, and cell-type annotation
├── 03.1_subtype_VIC.r                # VIC subtype analysis
├── 03.2_subtype_VEC.r                # VEC subtype analysis (including Monocle2, CytoTRACE2)
├── 03.3_subtype_immuneCells.r        # Immune cell subtype analysis (lymphoid & myeloid)
├── 04_fibrotic_transitional_analysis.r  # Fibrosis and transitional-state analysis
├── 05_geneset_score.r                # Gene set scoring (cell cycle, EMT, etc.)
├── 06_monocle2_analysis.r            # Monocle2 trajectory analysis
├── 07_velocity.r                     # RNA velocity analysis
├── 08_gene_corr_analysisr.r          # Gene correlation analysis
├── 09_CCI.r                          # Cell-cell communication (CellPhoneDB)
└── 10_supplemental_analysis.R        # Supplemental analysis
```

## Usage

### 1. Environment Setup

Make sure the following R packages are installed (as needed):

- Seurat
- SCP
- DoubletFinder
- Harmony
- monocle/monocle2
- velocyto.R / SCVELO
- CytoTRACE2
- CellPhoneDB
- CellChat
- openxlsx
- dplyr
- data.table
- ComplexHeatmap
- ClusterGVis
- ggplot2

## Script Overview

### 01_data_process.r

**Purpose:** Data loading and preprocessing

- Build sample metadata table
- Read 10X Genomics CellRanger outputs
- Create Seurat objects and merge multiple samples
- Compute mitochondrial and hemoglobin gene fractions
- Detect and filter doublets using DoubletFinder
- QC summary statistics and save intermediate objects


### 02_HarmonyUMAP.r

**Purpose:** Harmony integration, dimensionality reduction, clustering, and cell-type annotation

- Batch correction and integration with Harmony
- PCA and UMAP visualization
- Leiden clustering (multiple resolutions)
- Differential expression analysis
- Marker-based cell-type annotation
- Cell-type composition summaries and visualization


### 03.1_subtype_VIC.r

**Purpose:** VIC (valvular interstitial cell) subtype analysis

- Subset VIC cells and re-integrate
- Sub-clustering and resolution selection (clustree)
- Differential expression between subclusters
- Marker identification and visualization
- Subcluster proportion comparison across conditions
- Heatmaps and feature plots


### 03.2_subtype_VEC.r

**Purpose:** VEC (valvular endothelial cell) subtype analysis

- VEC subsetting and processing
- Downsampling for trajectory inference
- Monocle2 pseudotime trajectory analysis
- CytoTRACE2 developmental potential estimation
- Visualization of differentiation states


### 03.3_subtype_immuneCells.r

**Purpose:** Immune cell subtype analysis

- Separate analysis for lymphoid and myeloid cells
- Harmony integration and clustering
- Comparison across multiple clustering resolutions
- Differential expression between subclusters
- Marker identification
- Enrichment analyses (GO/KEGG, etc.)


### 04_fibrotic_transitional_analysis.r

**Purpose:** Fibrosis and transitional-state analysis

- Differential analysis of pro-fibrotic vs anti-fibrotic VIC programs
- Fibrosis-related gene analysis across experimental groups
- Identification of transitional-state cells in VIC/VEC (ModuleScore-based)
- Distribution of transitional cells in UMAP space
- Proportion statistics of transitional cells across conditions


### 05_geneset_score.r

**Purpose:** Gene set scoring

- Cell cycle scoring (G1/S and G2/M)
- EMT (epithelial–mesenchymal transition) scoring
- Other module scores using AddModuleScore
- Visualization of scores in UMAP space
- Group-wise comparisons of scores


### 06_monocle2_analysis.r

**Purpose:** Monocle2 trajectory analysis

- Build Monocle2 pseudotime trajectories
- Identify branching points
- BEAM (Branched Expression Analysis Modeling)
- Distribution of states along pseudotime
- Trajectory visualization
- Pseudotime density analysis


### 07_velocity.r

**Purpose:** RNA velocity analysis

- Configure velocyto.R environment
- Read and prepare loom files
- Process spliced/unspliced matrices
- Run SCVELO analysis
- Visualize velocity vectors on UMAP
- Stream/flow plots


### 08_gene_corr_analysisr.r

**Purpose:** Gene correlation analysis

- Extract genes of interest (e.g., VCAM1, PECAM1, ACTA2)
- Compare gene expression across groups
- Statistical tests (t-test, Wilcoxon test)
- Summarize and export results


### 09_CCI.r

**Purpose:** Cell–cell communication analysis (CellPhoneDB)

- Prepare CellPhoneDB input files (expression matrix and metadata)
- Process data by tissue (Tri, Mit)
- Read and post-process CellPhoneDB outputs
- Ligand–receptor interaction analysis
- Pathway enrichment analysis
- Visualization of results


### 10_supplemental_analysis.R

**Purpose:** Supplemental analyses

- Sample grouping (AF vs SR, FMR vs FTR)
- Exclusion/sensitivity analyses under different criteria
- Cell proportion statistics
- ANOVA
- Regrouping and visualization
- Additional exploratory analyses


## Notes

1. **Update file paths:** All paths are placeholders. Please update them according to your local environment.
2. **Dependencies:** Some analyses require specific environments and system libraries.
   - `velocyto.R` may require a dedicated conda environment and system dependencies.
   - `CellPhoneDB` requires a Python environment and CLI tools.
3. **Input files:** Ensure required inputs (e.g., `.rds`, `.loom`) exist in the specified locations.
4. **Memory requirements:** Single-cell analyses can be memory-intensive; make sure sufficient resources are available.
5. **Execution order:** It is recommended to run scripts in numeric order, as later scripts depend on outputs from earlier steps.
6. **Filename typo:** The file `08_gene_corr_analysisr.r` contains a typo (`analysisr` should be `analysis`).

## Custom Utility Functions

All custom helper functions are defined in `utils.R`, including:

- `qcstat()` - QC summary statistics
- `ratio_plot()` - Proportion plot visualization
- `run_monocle()` - Wrapper for Monocle trajectory analysis
- `prepareLoom()` - Loom preparation helper
- `run_cluster_anova()` - Cluster-level ANOVA
- Other utility helpers

## Version History

- Initial release

For questions or suggestions, please contact the repository maintainer.
