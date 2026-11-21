###############################################################################
# Annotate the endothelial compartment (adapted from the public pipeline of Shirazi et al. (https://github.com/SucreLab/HumanNeonatal))
#
# Source      : r.02_cellCluster.r (lines 1084-1321)
# Input       : split_endo.rds
# Output      : annotated_endo.rds
#
# Code below is unchanged from the original analysis script.
###############################################################################
# -----------------------------------------------------------------------------
# PATHS: this file keeps the absolute paths of the analysis server on which the
# study was run (/data/nas1/xiangxiaoyu_OD/1.BPD+MSCs/1.scRNA/... and
# /home/xxyloveme/study/BPD+MSCs/...). Set PROJ_DIR in config/paths.R and
# replace those prefixes with your own project directory before running.
# -----------------------------------------------------------------------------


#######################2_endothelium#################
# ---
#   title: "scRNA-seq Neonatal Lung Injury 2.2"
# author: Shawyon Shirazi & Nick Negretti
# date: 3/15/24
# output: rmarkdown::github_document
# ---
  
  # Analysis of the lung endothelium
  
  ## Load libraries and helper functions
  # ```{r, results="hide", message = FALSE, load-libraries}
setwd("/data/nas1/xiangxiaoyu_OD/1.BPD+MSCs/1.scRNA/02_Cell_Clustering_and_Annotation/")
source("./helper_functions/globals.R")
source("./helper_functions/libraries.R")
source("./helper_functions/cluster.R")
source("./helper_functions/colors.R")

N_WORKERS <- 12
plan("multicore", workers = N_WORKERS)


## Load data from merge and clean pipeline
#{r}
endo <- readRDS("/data/nas1/xiangxiaoyu_OD/1.BPD+MSCs/1.scRNA/01_Data_Download/4_split_celltypes/post_sct_split_endo.rds")


## UMAP and clustering pt1
#{r, fig.width = 12}

endo <- cluster_pca_umap(endo,
                         k_param = 10,
                         dims_umap = 1:20,
                         dims_neighbors = 1:20,
                         cluster_res = 1.5,
                         algorithm = 4)

p_cluster <- DimPlot(endo,
                     label = TRUE) +
  umap_theme() +
  theme(aspect.ratio = 1)

endo_cluster <- p_cluster +
  plot_annotation("Endothelium")
ggsave("./2_endothelium/cluster1.png",width = 8,height = 6)
ggsave("./2_endothelium/cluster1.pdf",width = 8,height = 6)

# {r}
marker_genes <- c("EPCAM", "PTPRC", "COL1A1", # Broad subtype
                  "PECAM1", "VWF", # Endo and vessel specific
                  "CCL21", "PROX1", "MMRN1", # Lymphatic
                  "CA4", "EDNRB", "APLN", # aCap
                  "FCN3", "IL7R", "CD36", "APLNR", # gCap
                  "HEY1", "GJA5", "DKK2", # Arterial
                  "COL15A1", "EBF1", "PLVAP", "SPRY1", # Systemic Venous
                  "ACKR1", # Venous
                  "CPE", "C7", "HDAC9", # Pulmonary Venous
                  "MKI67", "TOP2A", "CDC20") # Prolif and extra
endo_markers <- DotPlot(endo,
                        features = marker_genes) +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1))
ggsave("./2_endothelium/markers1.pdf",width = 13,height = 8)
ggsave("./2_endothelium/markers1.png",width = 13,height = 8)

# {r metrics}
# 提取元数据
endo_meta_data <- endo@meta.data
# 转换为长格式用于ggplot2
plot_data <- endo_meta_data %>%
  select(condition, nFeature_RNA, nCount_RNA, percent.mt) %>%
  pivot_longer(cols = c(nFeature_RNA, nCount_RNA, percent.mt), 
               names_to = "metric", values_to = "value")
# 创建小提琴图
endo_metric <- ggplot(plot_data, aes(x = condition, y = value, fill = condition)) +
  geom_violin(scale = "width", trim = TRUE) +
  facet_wrap(~ metric, scales = "free_y", ncol = 1) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("./2_endothelium/metrics1.png",
       plot = endo_metric,width = 7,height = 9)
ggsave("./2_endothelium/metrics1.pdf",
       plot = endo_metric,width = 7,height = 9)


# {r}
N_WORKERS <- 20
plan("multicore", workers = N_WORKERS)
endo_sub_sct <- SCTransform(subset(endo,
                                   idents = c(25),
                                   invert = TRUE),
                            variable.features.n = 2000,
                            conserve.memory = TRUE,
                            vars.to.regress = c("percent.mt",
                                                "percent.rb"),
                            method = "glmGamPoi",
                            ncells = NULL)
var_feat_backup2 <- VariableFeatures(endo_sub_sct)


## UMAP and clustering pt2
# {r, fig.width = 12}
endo_sub_sct <- cluster_pca_umap(endo_sub_sct,
                                 k_param = 5,
                                 dims_umap = 1:15,
                                 dims_neighbors = 1:15,
                                 cluster_res = 1,
                                 algorithm = 4)

p_cluster <- DimPlot(endo_sub_sct,
                     label = TRUE) +
  umap_theme() +
  theme(aspect.ratio = 1)

endo_cluster2 <- p_cluster +
  plot_annotation("Endothelium")
ggsave("./2_endothelium/cluster2.png",
       plot = endo_cluster2,width = 7,height = 6)
ggsave("./2_endothelium/cluster2.pdf",
       plot = endo_cluster2,width = 7,height = 6)

# {r}
marker_genes <- c("EPCAM", "PTPRC", "COL1A1",
                  "PECAM1", "VWF",
                  "CCL21", "PROX1", "MMRN1",
                  "CA4", "EDNRB", "APLN",
                  "FCN3", "IL7R", "CD36", "APLNR",
                  "HEY1", "GJA5", "DKK2",
                  "COL15A1", "EBF1", "PLVAP", "SPRY1",
                  "ACKR1",
                  "CPE", "C7", "HDAC9",
                  "MKI67", "TOP2A", "CDC20", "ANKRD1")
endo_markers2 <- DotPlot(endo_sub_sct,
                         features = marker_genes) +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1))
ggsave("./2_endothelium/markers2.pdf",width = 13,height = 8)
ggsave("./2_endothelium/markers2.png",width = 13,height = 8)


# {r metrics}
endo_meta_data <- endo_sub_sct@meta.data
# 转换为长格式用于ggplot2
plot_data <- endo_meta_data %>%
  select(condition, nFeature_RNA, nCount_RNA, percent.mt) %>%
  pivot_longer(cols = c(nFeature_RNA, nCount_RNA, percent.mt), 
               names_to = "metric", values_to = "value")
# 创建小提琴图
endo_metric <- ggplot(plot_data, aes(x = condition, y = value, fill = condition)) +
  geom_violin(scale = "width", trim = TRUE) +
  facet_wrap(~ metric, scales = "free_y", ncol = 1) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("./2_endothelium/metrics2.png",
       plot = endo_metric,width = 7,height = 9)
ggsave("./2_endothelium/metrics2.pdf",
       plot = endo_metric,width = 7,height = 9)

## Save the DEgenes per cluster
# {r}
N_WORKERS <- 20
plan("multicore", workers = N_WORKERS)
endo_sub_sct <- PrepSCTFindMarkers(endo_sub_sct)

filename <- "./2_endothelium/pre_annotated_endo_for_deg.rds"
if (!file.exists(filename)) {
  markers <- parallel_find_all_markers(endo_sub_sct)
  
  names(markers) <- levels(Idents(endo_sub_sct))
  saveRDS(markers, filename)
} else {
  markers <- readRDS(filename)
}


wb_markers_names <- levels(Idents(endo_sub_sct))
wb_markers <- createWorkbook()
for (idx in seq_along(wb_markers_names)){
  addWorksheet(wb_markers,
               wb_markers_names[idx])
  writeData(wb_markers,
            wb_markers_names[idx],
            markers[[idx]],
            rowNames = TRUE)
}
saveWorkbook(wb_markers,
             file = "./2_endothelium/endo_cluster_deg.xlsx",
             overwrite = TRUE)


## Annotation
# {r}
endo_relabel <- endo_sub_sct
endo_relabel$louvain <- as.factor(endo_sub_sct$seurat_clusters)
endo_relabel$celltype <- as.factor(endo_sub_sct$seurat_clusters)

idxes <- levels(endo_relabel$celltype)
labels <- rep("gCap", length(levels(endo_relabel$celltype)))
labels[idxes %in% c(3,4)] <- "aCap"
labels[idxes %in% c(10)] <- "abCap"
labels[idxes %in% c(6,16)] <- "Arterial EC"
labels[idxes %in% c(2)] <- "Pulmonary venous EC"
labels[idxes %in% c(7)] <- "Lymphatic"
labels[idxes %in% c(13,17,19)] <- "Systemic venous EC"

levels(endo_relabel$celltype) <- labels
endo_relabel$celltype <- ordered(as.factor(endo_relabel$celltype),
                                 c("gCap",
                                   "aCap",
                                   "abCap",
                                   "Arterial EC",
                                   "Pulmonary venous EC",
                                   "Systemic venous EC",
                                   "Lymphatic"
                                 ))

p_cluster <- DimPlot(endo_relabel,
                     group.by = "celltype",
                     label = TRUE) +
  umap_theme() +
  scale_colour_manual(name = "Cluster",
                      values = color_category_20) +
  theme(aspect.ratio = 1)

p_cluster + plot_annotation("Endothelium")


## Save data for next step merge
# {r}
plan("multicore", workers = N_WORKERS)
filename <- "./2_endothelium/annotated_endo.rds"
saveRDS(endo_relabel,
        filename,
        compress = "gzip")
