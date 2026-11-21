###############################################################################
# Re-embed the merged object for the cell browser (adapted from the public pipeline of Shirazi et al. (https://github.com/SucreLab/HumanNeonatal))
#
# Source      : 2_cellxgene_files.R
# Input       : alldata_neonatal_injury.rds
# Output      : re-embedded object
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
#   title: "scRNA-seq Neonatal Lung Injury 3.2"
# author: Shawyon Shirazi & Nick Negretti
# date: 3/15/24
# output: rmarkdown::github_document
# ---
  
  ## Load libraries and helper functions
  # ```{r, results="hide", message = FALSE}
source("./helper_functions/globals.R")
source("./helper_functions/libraries.R")
source("./helper_functions/cluster.R")
source("./helper_functions/colors.R")

N_WORKERS <- 12
plan("multicore", workers = N_WORKERS)

# ```{r}
merged_seurat <- readRDS("./alldata_neonatal_injury.rds")


# ```{r}
Idents(merged_seurat) <- merged_seurat$celltype

# ```{r}
merged_seurat <- Seurat::RunPCA(merged_seurat, verbose = FALSE)

merged_seurat <- RunHarmony(merged_seurat,
                            "dataset",
                            assay.use = "SCT",
                            lambda = 1,
                            sigma = 0.4)

merged_seurat <- Seurat::RunUMAP(merged_seurat,
                                 reduction = "harmony",
                                 dims = 1:40,
                                 verbose = FALSE,
                                 n.neighbors = 30,
                                 metric = "euclidean")

merged_seurat <- Seurat::FindNeighbors(merged_seurat,
                                       reduction = "harmony",
                                       dims = 1:40,
                                       k.param = 10)

merged_seurat <- Seurat::FindClusters(merged_seurat,
                                      resolution = 0.3,
                                      algorithm = 1)
