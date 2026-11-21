###############################################################################
# Merge the four annotated lineages into one object (adapted from the public pipeline of Shirazi et al. (https://github.com/SucreLab/HumanNeonatal))
#
# Source      : 1_merge_all_cells.R
# Input       : annotated_epi/endo/meso/imm .rds
# Output      : alldata_neonatal_injury.rds
#
# Code below is unchanged from the original analysis script.
###############################################################################
# -----------------------------------------------------------------------------
# PATHS: this file keeps the absolute paths of the analysis server on which the
# study was run (/data/nas1/xiangxiaoyu_OD/1.BPD+MSCs/1.scRNA/... and
# /home/xxyloveme/study/BPD+MSCs/...). Set PROJ_DIR in config/paths.R and
# replace those prefixes with your own project directory before running.
# -----------------------------------------------------------------------------


# ---
#   title: "scRNA-seq Neonatal Lung Injury 3.1"
# author: Shawyon Shirazi & Nick Negretti
# date: 3/15/24
# output: rmarkdown::github_document
# ---
  
  ## Load libraries and helper functions
  # ```{r, results="hide", message = FALSE}
setwd("/data/nas1/xiangxiaoyu_OD/1.BPD+MSCs/1.scRNA/03_All_cells/")

source("./helper_functions/globals.R")
source("./helper_functions/libraries.R")
source("./helper_functions/cluster.R")
source("./helper_functions/colors.R")

N_WORKERS <- 12
plan("multicore", workers = N_WORKERS)
RhpcBLASctl::blas_set_num_threads(4)
RhpcBLASctl::omp_set_num_threads(4)


# ```{r}
epi <- readRDS("/data/nas1/xiangxiaoyu_OD/1.BPD+MSCs/1.scRNA/02_Cell_Clustering_and_Annotation/1_epithelium/annotated_epi.rds")
DefaultAssay(epi) <- "RNA"
epi <- DietSeurat(epi,
                  assays = c("RNA"))
epi$bulk_celltype <- "Epithelial"

endo <- readRDS("/data/nas1/xiangxiaoyu_OD/1.BPD+MSCs/1.scRNA/02_Cell_Clustering_and_Annotation/2_endothelium/annotated_endo.rds")
DefaultAssay(endo) <- "RNA"
endo <- DietSeurat(endo,
                   assays = c("RNA"))
endo$bulk_celltype <- "Endothelial"

meso <- readRDS("/data/nas1/xiangxiaoyu_OD/1.BPD+MSCs/1.scRNA/02_Cell_Clustering_and_Annotation/3_mesenchyme/annotated_meso.rds")
DefaultAssay(meso) <- "RNA"
meso <- DietSeurat(meso,
                   assays = c("RNA"))
meso$bulk_celltype <- "Mesenchymal"

imm <- readRDS("/data/nas1/xiangxiaoyu_OD/1.BPD+MSCs/1.scRNA/02_Cell_Clustering_and_Annotation/4_immune/annotated_imm.rds")
DefaultAssay(imm) <- "RNA"
imm <- DietSeurat(imm,
                  assays = c("RNA"))
imm$bulk_celltype <- "Immune"


# ```{r}
merged_seurat <- Reduce(function(x, y) {merge(x, y)},
                        list(epi, endo, meso, imm))

merged_seurat$cell_subtype <- Idents(merged_seurat)
merged_seurat$bulk_celltype <- as.factor(merged_seurat$bulk_celltype)

# ```{r}
merged_seurat <- SCTransform(merged_seurat,
                             variable.features.n = 2000,
                             conserve.memory = FALSE,
                             vars.to.regress = c("percent.mt",
                                                 "percent.rb"),
                             method = "glmGamPoi",
                             ncells = NULL)

merged_seurat <- cluster_pca_umap(merged_seurat,
                                  dims_umap = 1:40,
                                  dims_neighbors = 1:40,
                                  cluster_res = 1,
                                  harmony = TRUE,
                                  harmonysigma = 0.4)


# ```{r}
filename <- "./alldata_neonatal_injury.rds"
saveRDS(merged_seurat,
        filename,
        compress = "gzip")
