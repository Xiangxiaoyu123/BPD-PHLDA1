###############################################################################
# Split the integrated object into the four lineages (adapted from the public pipeline of Shirazi et al. (https://github.com/SucreLab/HumanNeonatal))
#
# Source      : r.00_rawdata.r (lines 387-578)
# Input       : sct_on_merge.rds
# Output      : split_epi/meso/endo/imm .rds
#
# Code below is unchanged from the original analysis script.
###############################################################################
# -----------------------------------------------------------------------------
# PATHS: this file keeps the absolute paths of the analysis server on which the
# study was run (/data/nas1/xiangxiaoyu_OD/1.BPD+MSCs/1.scRNA/... and
# /home/xxyloveme/study/BPD+MSCs/...). Set PROJ_DIR in config/paths.R and
# replace those prefixes with your own project directory before running.
# -----------------------------------------------------------------------------


############4_split_celltypes#########################
setwd("/data/nas1/xiangxiaoyu_OD/1.BPD+MSCs/1.scRNA/01_Data_Download/")
# ---
#   title: "scRNAseq Neonatal Lung Injury 1.4"
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

N_WORKERS <- 12
plan("multicore", workers = N_WORKERS)


# Load merged object
# ```{r}
alldata <- readRDS("./sct_on_merge.rds")

# Function to plot both marker genes and umap clusters
# ```{r}
dotdim_preprocess <- function(alldata, filename) {
  marker_genes <- c("EPCAM", #Epithelial
                    "KRT19",
                    "DNAH12",
                    "PTPRC", #Immune
                    "AIF1",
                    "SPI1",
                    "CD247",
                    "VPREB3",
                    "TPSAB1",
                    "PECAM1", #Endothelial
                    "CDH5",
                    "COL1A1", #Mesenchymal
                    "COL1A2")
  dotplot.preprocess <- DotPlot(alldata,
                                features = marker_genes) +
    theme(axis.text.x = element_text(angle = 45,
                                     hjust = 1))
  ggsave(paste(filename,
               "dotplot.pdf",
               sep = "_"),
         plot = dotplot.preprocess)
  dimplot.preprocess <- DimPlot(alldata,
                                label = TRUE)
  ggsave(paste(filename,
               "dimplot.png",
               sep = "_"),
         plot = dimplot.preprocess)
}


# ```{r}
setwd("/data/nas1/xiangxiaoyu_OD/1.BPD+MSCs/1.scRNA/01_Data_Download/")
dotdim_preprocess(alldata, "./4_split_celltypes/annot_lineages")


# Extra step to identify unknown clusters
# ```{r}
# clust_35_markers <- FindMarkers(alldata, 35, NULL)
# clust_35_markers$gene <- rownames(clust_35_markers)  # Workaround for print issue
# clust_35_markers

# Separate out the clusters to be split by broad celltype
# ```{r soupx_norm}
imm_idents <- c(2,3,4,8,11,13,21,23,24,25,29)
epi_idents <- c(12,14,18,27)
meso_idents <- c(1,5,6,9,15,17,28)
endo_idents <- c(0,7,10,16,19,20,22,26)

setdiff(0:29,
        c(epi_idents,
          meso_idents,
          endo_idents,
          immune_idents)) #Any cluster assignment missing?

c(epi_idents,
  meso_idents,
  endo_idents,
  immune_idents)[duplicated(c(epi_idents,
                              meso_idents,
                              endo_idents,
                              immune_idents))] # Any duplicated assignments?


# Subset the data and save a pre-SCTransformed version
# ```{r}
split_epi <- subset(alldata,
                    idents = epi_idents)
split_meso <- subset(alldata,
                     idents = meso_idents)
split_endo <- subset(alldata,
                     idents = endo_idents)
split_imm <- subset(alldata,
                    idents = imm_idents)

saveRDS(split_epi,
        "./4_split_celltypes/pre_sct_split_epi.rds",
        compress = "bzip2")
saveRDS(split_meso,
        "./4_split_celltypes/pre_sct_split_meso.rds",
        compress = "bzip2")
saveRDS(split_endo,
        "./4_split_celltypes/pre_sct_split_endo.rds",
        compress = "bzip2")
saveRDS(split_imm,
        "./4_split_celltypes/pre_sct_split_imm.rds",
        compress = "bzip2")

# Run SCTransform and clustering on each data subset and save the data
# ```{r}
N_WORKERS <- 30
plan("multicore",
     workers = N_WORKERS)
split_epi <- SCTransform(split_epi,
                         variable.features.n = 2000,
                         conserve.memory = TRUE,
                         vars.to.regress = c("percent.mt",
                                             "percent.rb"),
                         method = "glmGamPoi",
                         ncells = NULL)
split_epi <- cluster_pca_umap(split_epi,
                              dims_umap = 1:35,
                              dims_neighbors = 1:35,
                              cluster_res = 0.6)
saveRDS(split_epi,
        "./4_split_celltypes/post_sct_split_epi.rds",
        compress = "bzip2")

# ```{r}
N_WORKERS <- 30
plan("multicore",
     workers = N_WORKERS)
split_meso <- SCTransform(split_meso,
                          variable.features.n = 2000,
                          conserve.memory = TRUE,
                          vars.to.regress = c("percent.mt",
                                              "percent.rb"),
                          method = "glmGamPoi",
                          ncells = NULL)
split_meso <- cluster_pca_umap(split_meso,
                               dims_umap = 1:35,
                               dims_neighbors = 1:35,
                               cluster_res = 0.6)
saveRDS(split_meso,
        "./4_split_celltypes/post_sct_split_meso.rds",
        compress = "bzip2")

# ```{r}
N_WORKERS <- 30
plan("multicore",
     workers = N_WORKERS)
split_endo <- SCTransform(split_endo,
                          variable.features.n = 2000,
                          conserve.memory = TRUE,
                          vars.to.regress = c("percent.mt",
                                              "percent.rb"),
                          method = "glmGamPoi",
                          ncells = NULL)
split_endo <- cluster_pca_umap(split_endo,
                               dims_umap = 1:35,
                               dims_neighbors = 1:35,
                               cluster_res = 0.6)
saveRDS(split_endo,
        "./4_split_celltypes/post_sct_split_endo.rds",
        compress = "bzip2")

# ```{r}

N_WORKERS <- 30
plan("multicore",
     workers = N_WORKERS)
split_imm <- SCTransform(split_imm,
                         variable.features.n = 2000,
                         conserve.memory = TRUE,
                         vars.to.regress = c("percent.mt",
                                             "percent.rb"),
                         method = "glmGamPoi",
                         ncells = NULL)
split_imm <- cluster_pca_umap(split_imm,
                              dims_umap = 1:35,
                              dims_neighbors = 1:35,
                              cluster_res = 0.6)
saveRDS(split_imm,
        "./4_split_celltypes/post_sct_split_imm.rds",
        compress = "gzip")
