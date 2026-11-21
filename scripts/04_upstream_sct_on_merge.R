###############################################################################
# SCTransform and Harmony integration of the merged object (adapted from the public pipeline of Shirazi et al. (https://github.com/SucreLab/HumanNeonatal))
#
# Source      : r.00_rawdata.r (lines 303-386)
# Input       : label_merge.rds
# Output      : sct_on_merge.rds
#
# Code below is unchanged from the original analysis script.
###############################################################################
# -----------------------------------------------------------------------------
# PATHS: this file keeps the absolute paths of the analysis server on which the
# study was run (/data/nas1/xiangxiaoyu_OD/1.BPD+MSCs/1.scRNA/... and
# /home/xxyloveme/study/BPD+MSCs/...). Set PROJ_DIR in config/paths.R and
# replace those prefixes with your own project directory before running.
# -----------------------------------------------------------------------------


############3_sct_on_merge#########################
# ---
#   title: "scRNAseq Neonatal Lung Injury 1.3"
# author: Shawyon Shirazi & Nick Negretti
# date: 3/14/24
# output: rmarkdown::github_document
# ---
  
  ## Load libraries and helper functions
  # ```{r, results="hide", message = FALSE, load-libraries}
source("./helper_functions/globals.R")
source("./helper_functions/libraries.R")
source("./helper_functions/cluster.R")
source("./helper_functions/colors.R")
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))

N_WORKERS <- 12
plan("multicore", workers = N_WORKERS)

# Load merged object
# ```{r}
alldata <- readRDS("./label_merge.rds")

# Calculate percent ribosomal genes
# ```{r percent ribo}
alldata <- PercentageFeatureSet(alldata,
                                pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA",
                                col.name = "percent.rb")

# SCTranform
# ```{r}

#报错
# 确保dataset是因子类型
alldata$dataset <- as.factor(alldata$dataset)

alldata <- SCTransform(alldata,
                       variable.features.n = 3000,
                       conserve.memory = FALSE,
                       vars.to.regress = c("percent.mt",
                                           "percent.rb"),
                       method = "glmGamPoi",
                       ncells = NULL)


# Harmonize and cluster the data
# ```{r}
plan("sequential")
alldata <- cluster_pca_umap(alldata,
                            dims_umap = 1:35,
                            dims_neighbors = 1:35,
                            cluster_res = 0.6)
saveRDS(alldata,
        "./sct_on_merge.rds",
        compress = "bzip2")

# Write data as h5ad for explorer
# ```{r}
writeAnnData <- function(seurat_obj, filename, assay = "SCT") {
  anndata::AnnData(X = t(GetAssayData(seurat_obj,
                                      slot = "data",
                                      assay = assay)),
                   var = data.frame(gene = rownames(seurat_obj),
                                    row.names = rownames(seurat_obj)),
                   obs = data.frame(louvain = as.character(Idents(seurat_obj)),
                                    dataset = as.character(seurat_obj$dataset),
                                    seq_run = as.character(seurat_obj$orig.ident),
                                    age = as.character(seurat_obj$age),
                                    condition = as.character(seurat_obj$condition),
                                    broad_condition = as.character(seurat_obj$broad_condition),
                                    row.names = colnames(seurat_obj)),
                   obsm = list(X_umap = matrix(Embeddings(seurat_obj,
                                                          reduction = "umap"),
                                               ncol = 2))
  )$write_h5ad(filename, compression = "gzip")
}

filename <- "./preprocessed_preview.h5ad"
if (!file.exists(filename)) { # Create if doesn't exist
  writeAnnData(alldata, filename)
}
