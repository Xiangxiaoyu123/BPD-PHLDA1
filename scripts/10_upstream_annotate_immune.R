###############################################################################
# Annotate the immune compartment (adapted from the public pipeline of Shirazi et al. (https://github.com/SucreLab/HumanNeonatal))
#
# Source      : r.02_cellCluster.r (lines 1547-1817)
# Input       : split_imm.rds
# Output      : annotated_imm.rds
#
# Code below is unchanged from the original analysis script.
###############################################################################
# -----------------------------------------------------------------------------
# PATHS: this file keeps the absolute paths of the analysis server on which the
# study was run (/data/nas1/xiangxiaoyu_OD/1.BPD+MSCs/1.scRNA/... and
# /home/xxyloveme/study/BPD+MSCs/...). Set PROJ_DIR in config/paths.R and
# replace those prefixes with your own project directory before running.
# -----------------------------------------------------------------------------


#########################4_immune##################
setwd("/data/nas1/xiangxiaoyu_OD/1.BPD+MSCs/1.scRNA/02_Cell_Clustering_and_Annotation/")
# ---
#   title: "scRNA-seq Neonatal Lung Injury 2.4"
# author: Shawyon Shirazi & Nick Negretti
# date: 3/15/24
# output: rmarkdown::github_document
# ---
  
  # Analysis of the lung immune
  
  ## Load libraries and helper functions
  # ```{r, results="hide", message = FALSE, load-libraries}
source("./helper_functions/globals.R")
source("./helper_functions/libraries.R")
source("./helper_functions/cluster.R")
source("./helper_functions/colors.R")

N_WORKERS <- 12
plan("multicore", workers = N_WORKERS)


# ```{r}
imm <- readRDS("/data/nas1/xiangxiaoyu_OD/1.BPD+MSCs/1.scRNA/01_Data_Download/4_split_celltypes/post_sct_split_imm.rds")


# ```{r}
N_WORKERS <- 20
plan("multicore", workers = N_WORKERS)
imm_sct <- SCTransform(imm,
                       variable.features.n = 2000,
                       conserve.memory = TRUE,
                       vars.to.regress = c("percent.mt",
                                           "percent.rb"),
                       method = "glmGamPoi",
                       ncells = NULL)

# ```{r}
imm_metric <- VlnPlot(imm_sct,
                      features = c("nFeature_RNA",
                                   "nCount_RNA",
                                   "percent.mt"),
                      group.by = "condition")
ggsave("./data/checkpoint/imm/metric_testing.png",
       plot = imm_metric)

## UMAP and clustering pt1
# ```{r, fig.width = 12}
imm_sct <- cluster_pca_umap(imm_sct,
                            k_param = 10,
                            dims_umap = 1:15,
                            dims_neighbors = 1:15,
                            cluster_res = 1.2,
                            harmony = TRUE,
                            harmonysigma = 0.4)

p_cluster <- DimPlot(imm_sct,
                     label = TRUE) +
  umap_theme() +
  theme(aspect.ratio = 1)

imm_cluster <- p_cluster +
  plot_annotation("Immune")
ggsave("./4_immune/cluster1.png",
       plot = imm_cluster)

# ```{r}
marker_genes <- c("EPCAM", "PECAM1", "PTPRC", "COL1A1", # Broad subtype
                  "CD3E", "IL7R", "CD8A", "CCL5", "THEMIS", "CD6",  # T Cell
                  "NCR1", "KLRB1", "NKG7", "FGFBP2", "SH2D1B", "MCTP2", # NK Cell
                  "LYZ", "MARCO", "FCGR1A", "C1QA", "APOC1", # Macrophage
                  "S100A12", "FCN1", "CD14", "VCAN", "ITGAX", "ITGAM", # Monocyte
                  "CLEC4C","FCER1A", "CD1C", "SIPA1L3", "HDAC9", "FCGR2B", "CLEC10A", "FOXN3", # DC
                  "MS4A1", "CD19", "CD79A", "VPREB3", "FCRL1", "EBF1", "BCL11A", "IGHM", "IGHG1", # B Cell
                  "CPA3", "KIT", "TPSD1", "GATA2", "MS4A2", # Mast Cell
                  "GPD1", "PPIC", "AMIGO2", "RBP4", "CES1", "TREM1", "CD52", # Alveolar macrophage
                  "S100A9", "S100A8") # Neutrophil
imm_markers <- DotPlot(imm_sct,
                       features = marker_genes) +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1))
ggsave("./4_immune/markers1.pdf",
       plot = imm_markers,
       width = 17,height = 9)
ggsave("./4_immune/markers1.png",
       plot = imm_markers,
       width = 17,height = 9)


# ```{r metrics}
vlnplot_metrics <- VlnPlot(imm_sct,
                           features = c("nCount_RNA",
                                        "nCount_SCT",
                                        "percent.mt",
                                        "percent.rb"),
                           group.by = "seurat_clusters")
ggsave("./4_immune/metrics1.png",
       plot = vlnplot_metrics,
       width = 20)
DimPlot(imm_sct,
        group.by = "broad_condition")

# ```{r}
imm_sub_sct <- subset(imm_sct,
                      idents = c(11,21,26,23,24,17),
                      invert = TRUE)

# ```{r}
N_WORKERS <- 20
plan("multicore", workers = N_WORKERS)
imm_sub_sct2 <- SCTransform(imm_sub_sct,
                            variable.features.n = 2000,
                            conserve.memory = TRUE,
                            vars.to.regress = c("percent.mt",
                                                "percent.rb"),
                            method = "glmGamPoi",
                            ncells = NULL)

## UMAP and clustering pt2
# ```{r}
var_feat_backup_imm_sub <- VariableFeatures(imm_sub_sct2)

imm_sub_sct2 <- cluster_pca_umap(imm_sub_sct2,
                                 k_param = 10,
                                 dims_umap = 1:15,
                                 dims_neighbors = 1:15,
                                 cluster_res = 1.5,
                                 harmonysigma = 0.4)

# ```{r}
p_cluster <- DimPlot(imm_sub_sct2,
                     label = TRUE) +
  umap_theme() +
  theme(aspect.ratio = 1)

imm_cluster3 <- p_cluster +
  plot_annotation("Immune")
ggsave("./4_immune/cluster2.png",width = 9,height = 7,plot = imm_cluster3)
ggsave("./4_immune/cluster2.pdf",width = 9,height = 7,plot = imm_cluster3)

# ```{r}
marker_genes <- c("EPCAM", "PECAM1", "PTPRC", "COL1A1",
                  "CD3E", "IL7R", "CD8A", "CCL5", "THEMIS", "CD6",
                  "NCR1", "KLRB1", "NKG7", "FGFBP2", "SH2D1B", "MCTP2",
                  "LYZ", "MARCO", "FCGR1A", "C1QA", "APOC1",
                  "S100A12", "FCN1", "CD14", "VCAN", "ITGAX", "ITGAM",
                  "CLEC4C","FCER1A", "CD1C", "SIPA1L3", "HDAC9", "FCGR2B", "CLEC10A", "FOXN3",
                  "MS4A1", "CD19", "CD79A", "VPREB3", "FCRL1", "EBF1", "BCL11A", "IGHM", "IGHG1",
                  "CPA3", "KIT", "TPSD1", "GATA2", "MS4A2",
                  "GPD1", "PPIC", "AMIGO2", "RBP4", "CES1", "TREM1", "CD52",
                  "S100A9", "S100A8")
imm_dotplot3 <- DotPlot(imm_sub_sct2,
                        features = marker_genes) +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1))
ggsave("./4_immune/markers2.pdf",
       width = 15,height = 9,plot = imm_dotplot3)
ggsave("./4_immune/markers2.png",
       width = 15,height = 9,plot = imm_dotplot3)

# ```{r}
vlnplot_metrics <- VlnPlot(imm_sub_sct2,
                           features = c("nCount_RNA",
                                        "nCount_SCT",
                                        "percent.mt",
                                        "percent.rb"),
                           group.by = "seurat_clusters")
ggsave("./4_immune/metrics2.png",
       plot = vlnplot_metrics,
       width = 15)

FeaturePlot(imm_sub_sct2, c("CD2"))
head(FindMarkers(imm_sub_sct2,31,NULL),50)

DimPlot(imm_sub_sct2,
        group.by = "broad_condition")

## Save the DEgenes per cluster
# ```{r}
N_WORKERS <- 20
plan("multicore", workers = N_WORKERS)
filename <- "./4_immune/pre_annotated_imm_for_deg.rds"

imm_sub_sct2 <- PrepSCTFindMarkers(imm_sub_sct2)

if (!file.exists(filename)) {
  markers <- parallel_find_all_markers(imm_sub_sct2)
  
  names(markers) <- levels(Idents(imm_sub_sct2))
  saveRDS(markers,
          filename)
} else {
  markers <- readRDS(filename)
}

wb_markers_names <- levels(Idents(imm_sub_sct2))
wb_markers <- createWorkbook()
for (idx in seq_along(wb_markers_names)) {
  addWorksheet(wb_markers,
               wb_markers_names[idx])
  writeData(wb_markers,
            wb_markers_names[idx],
            markers[[idx]],
            rowNames = TRUE)
}
saveWorkbook(wb_markers,
             file = ".4_immune/imm_cluster_deg.xlsx",
             overwrite = TRUE)

## Annotation
# ```{r}
imm_relabel <- imm_sub_sct2
imm_relabel$louvain <- as.factor(imm_sub_sct2$seurat_clusters)
imm_relabel$celltype <- as.factor(imm_sub_sct2$seurat_clusters)

idxes <- levels(imm_relabel$celltype)
labels <- rep("T Cell", length(levels(imm_relabel$celltype)))
labels[idxes %in% c(11,17,26,31)] <- "NK Cell"
labels[idxes %in% c(7,13,14,25)] <- "NKT Cell"
labels[idxes %in% c(33,32,8)] <- "Alveolar macrophage"
labels[idxes %in% c(12,19,22,29)] <- "Monocyte"
labels[idxes %in% c(15,34)] <- "pDC"
labels[idxes %in% c(23)] <- "cDC"
labels[idxes %in% c(5,21)] <- "B Cell"
labels[idxes %in% c(35)] <- "Mast cell"
labels[idxes %in% c(27)] <- "Basophil"
labels[idxes %in% c(30)] <- "Plasma cell"
labels[idxes %in% c(18)] <- "Neutrophil"

levels(imm_relabel$celltype) <- labels
Idents(imm_relabel) <- imm_relabel$celltype

# ```{r}
imm_relabel$celltype <- ordered(as.factor(imm_relabel$celltype),
                                c("T Cell",
                                  "NK Cell",
                                  "NKT Cell",
                                  "Alveolar macrophage",
                                  "Monocyte",
                                  "pDC",
                                  "cDC",
                                  "B Cell",
                                  "Mast cell",
                                  "Plasma cell",
                                  "Basophil",
                                  "Neutrophil"
                                ))


## Save data for next step merge
# ```{r}
plan("multicore", workers = N_WORKERS)
filename <- "./4_immune/annotated_imm.rds"
if (!file.exists(filename)) {
  saveRDS(imm_relabel,
          filename,
          compress = "gzip")
}

p_cluster <- DimPlot(imm_relabel,
                     label = TRUE) +
  umap_theme() +
  theme(aspect.ratio = 1)

imm_cluster4 <- p_cluster +
  plot_annotation("Immune")
ggsave("./4_immune/celltype_finally.png",width = 9,height = 7,plot = imm_cluster4)
ggsave("./4_immune/celltype_finally.pdf",width = 9,height = 7,plot = imm_cluster4)
