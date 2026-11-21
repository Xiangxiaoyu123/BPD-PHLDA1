###############################################################################
# Annotate the mesenchymal compartment (adapted from the public pipeline of Shirazi et al. (https://github.com/SucreLab/HumanNeonatal))
#
# Source      : r.02_cellCluster.r (lines 1322-1546)
# Input       : split_meso.rds
# Output      : annotated_meso.rds (seven mesenchymal subpopulations)
#
# Code below is unchanged from the original analysis script.
###############################################################################
# -----------------------------------------------------------------------------
# PATHS: this file keeps the absolute paths of the analysis server on which the
# study was run (/data/nas1/xiangxiaoyu_OD/1.BPD+MSCs/1.scRNA/... and
# /home/xxyloveme/study/BPD+MSCs/...). Set PROJ_DIR in config/paths.R and
# replace those prefixes with your own project directory before running.
# -----------------------------------------------------------------------------


#########################3_mesenchyme##################
# ---
#   title: "scRNA-seq Neonatal Lung Injury 2.3"
# author: Shawyon Shirazi & Nick Negretti
# date: 3/15/24
# output: rmarkdown::github_document
# ---
  
  # Analysis of the lung mesenchyme
  
  ## Load libraries and helper functions
  # ```{r, results="hide", message = FALSE}
source("./helper_functions/globals.R")
source("./helper_functions/libraries.R")
source("./helper_functions/cluster.R")
source("./helper_functions/colors.R")

N_WORKERS <- 12
plan("multicore", workers = N_WORKERS)

## Load data from merge and clean pipeline
# ```{r}
meso <- readRDS("/data/nas1/xiangxiaoyu_OD/1.BPD+MSCs/1.scRNA/01_Data_Download/4_split_celltypes/post_sct_split_meso.rds")

## UMAP and clustering pt1
# ```{r}
meso <- cluster_pca_umap(meso,
                         dims_umap = 1:30,
                         dims_neighbors = 1:30,
                         cluster_res = 1,
                         harmonysigma = 0.25,
                         harmonylambda = 2)

# ```{r, fig.width = 12}
p_cluster <- DimPlot(meso,
                     label = TRUE) +
  theme(aspect.ratio = 1)

meso_cluster <- p_cluster +
  plot_annotation("Mesenchyme")
ggsave("./3_mesenchyme/cluster1.png",
       plot = meso_cluster,width = 8,height = 6)
ggsave("./3_mesenchyme/cluster1.pdf",
       plot = meso_cluster,width = 8,height = 6)

# ```{r}
marker_genes <- c("EPCAM", "PECAM1", "PTPRC", # Broad subtype
                  "COL1A1", "COL1A2", # Mesenchymal specific
                  "ACTA2", "PDGFRB", "ACTG2", "MYH11", # SMC
                  "CNN1", "ITGA7", # Vascular SMC
                  "DES", "LGR6", # Airway SMC
                  "TCF21", "WNT2", "NPNT", "INMT", # Alveolar FB
                  "CTHRC1", "POSTN", # Activated FB
                  "SERPINF1", "PI16", "MFAP5", "SCARA5", #Adventitial FB
                  "WNT5A", "FGF18", "PDGFRA", "HHIP", # MyoFB
                  "COX4I2", "CSPG4", "LAMC3", "TRPC6", # Pericyte
                  "MKI67", "TOP2A", "CDC20", "TBX4") # Prolif and extra
meso_markers <- DotPlot(meso,
                        features = marker_genes) +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1))
ggsave("./3_mesenchyme/markers1.pdf",plot = meso_markers,
       width = 13,height = 8)
ggsave("./3_mesenchyme/markers1.png",plot = meso_markers,
       width = 13,height = 8)

# ```{r metrics}
vlnplot_metrics <- VlnPlot(meso,
                           features = c("nCount_RNA",
                                        "nCount_SCT",
                                        "percent.mt",
                                        "percent.rb"),
                           group.by = "seurat_clusters")
ggsave("./3_mesenchyme/metrics1.png")

# ```{r}
meso_sub <- subset(meso,
                   idents = c(12,18,19),
                   invert = TRUE)

# ```{r}
DimPlot(meso_sub)

# ```{r}
N_WORKERS <- 20
plan("multicore", workers = N_WORKERS)
meso_sub_sct <- SCTransform(meso_sub,
                            variable.features.n = 2000,
                            conserve.memory = FALSE,
                            vars.to.regress = c("percent.mt",
                                                "percent.rb"),
                            method = "glmGamPoi",
                            ncells = NULL)
var_feat_backup2 <- VariableFeatures(meso_sub_sct)

## UMAP and clustering pt2
# ```{r, fig.width = 12}
meso_sub_sct <- cluster_pca_umap(meso_sub_sct,
                                 k_param = 30,
                                 dims_umap = 1:50,
                                 dims_neighbors = 1:50,
                                 cluster_res = 1,
                                 harmonysigma = 1)

p_cluster <- DimPlot(meso_sub_sct,label = T) +
  theme(aspect.ratio = 1)

meso_cluster2 <- p_cluster +
  plot_annotation("Mesenchyme")
ggsave("./3_mesenchyme/cluster2.png",plot = meso_cluster2, width = 8,height = 6)
ggsave("./3_mesenchyme/cluster2.pdf",plot = meso_cluster2, width = 8,height = 6)


# ```{r}
marker_genes <- c("EPCAM", "PECAM1", "PTPRC",
                  "COL1A1", "COL1A2",
                  "ACTA2", "PDGFRB", "ACTG2", "MYH11",
                  "CNN1", "ITGA7",
                  "DES", "LGR6",
                  "TCF21", "WNT2", "NPNT", "INMT",
                  "CTHRC1", "POSTN",
                  "SERPINF1", "PI16", "MFAP5", "SCARA5",
                  "WNT5A", "FGF18", "PDGFRA", "HHIP",
                  "COX4I2", "CSPG4", "LAMC3", "TRPC6",
                  "MKI67", "TOP2A", "CDC20", "TBX4")
meso_markers2 <- DotPlot(meso_sub_sct,
                         features = marker_genes) +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1))
ggsave("./3_mesenchyme/markers2.pdf",
       width = 13,height = 8)
ggsave("./3_mesenchyme/markers2.png",
       width = 13,height = 8)

# ```{r metrics}
vlnplot_metrics <- VlnPlot(meso_sub_sct,
                           features = c("nCount_RNA",
                                        "nCount_SCT",
                                        "percent.mt",
                                        "percent.rb"),
                           group.by = "seurat_clusters")
ggsave("./3_mesenchyme/metrics2.png", width = 15)


FeaturePlot(meso_sub_sct,"CTHRC1")

## Save the DEgenes per cluster
# ```{r}
N_WORKERS <- 20
plan("multicore", workers = N_WORKERS)
filename <- "./3_mesenchyme/pre_annotated_meso_for_deg.rds"
meso_sub_sct <- PrepSCTFindMarkers(meso_sub_sct)

if (!file.exists(filename)) {
  markers <- parallel_find_all_markers(meso_sub_sct)
  
  names(markers) <- levels(Idents(meso_sub_sct))
  saveRDS(markers,
          filename)
} else {
  markers <- readRDS(filename)
}


wb_markers_names <- levels(Idents(meso_sub_sct))
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
             file = "./3_mesenchyme/meso_cluster_deg.xlsx",
             overwrite = TRUE)

## Annotation
# ```{r}
meso_relabel <- meso_sub_sct
meso_relabel$louvain <- as.factor(meso_sub_sct$seurat_clusters)
meso_relabel$celltype <- as.factor(meso_sub_sct$seurat_clusters)

idxes <- levels(meso_relabel$celltype)
labels <- rep("Alveolar FB", length(levels(meso_relabel$celltype)))
labels[idxes %in% c(2,14)] <- "Ductal MyoFB"
labels[idxes %in% c(8)] <- "Alveolar MyoFB"
labels[idxes %in% c(5,10,15,17)] <- "VSMC"
labels[idxes %in% c(3,12)] <- "Adventitial FB"
labels[idxes %in% c(9,11)] <- "Pericyte"
labels[idxes %in% c(0,4)] <- "Activated FB"

levels(meso_relabel$celltype) <- labels
Idents(meso_relabel) <- meso_relabel$celltype

# Reorder the clusters
# ```{r}
meso_relabel$celltype <- ordered(as.factor(meso_relabel$celltype),
                                 c("Alveolar FB",
                                   "Activated FB",
                                   "Adventitial FB",
                                   "Alveolar MyoFB",
                                   "Ductal MyoFB",
                                   "Pericyte",
                                   "VSMC"))

# ```{r}
DimPlot(meso_relabel,
        group.by = "celltype")+
  umap_theme() +
  scale_colour_manual(name = "Cluster",
                      values = acd_color_meso) +
  theme(aspect.ratio = 1)
  

## Save data for next step merge
# ```{r}
plan("multicore", workers = N_WORKERS)
filename <- "./3_mesenchyme/annotated_meso.rds"
saveRDS(meso_relabel,
        filename,
        compress = "gzip")
