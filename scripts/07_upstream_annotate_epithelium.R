###############################################################################
# Annotate the epithelial compartment (adapted from the public pipeline of Shirazi et al. (https://github.com/SucreLab/HumanNeonatal))
#
# Source      : r.02_cellCluster.r (lines 641-1083)
# Input       : split_epi.rds
# Output      : annotated_epi.rds
#
# Code below is unchanged from the original analysis script.
###############################################################################
# -----------------------------------------------------------------------------
# PATHS: this file keeps the absolute paths of the analysis server on which the
# study was run (/data/nas1/xiangxiaoyu_OD/1.BPD+MSCs/1.scRNA/... and
# /home/xxyloveme/study/BPD+MSCs/...). Set PROJ_DIR in config/paths.R and
# replace those prefixes with your own project directory before running.
# -----------------------------------------------------------------------------


########################文献代码######################
#######################1_epithelium#################
# ---
#   title: "scRNAseq Neonatal Lung Injury 2.1"
# author: Shawyon Shirazi & Nick Negretti
# date: 3/13/24
# output: rmarkdown::github_document
# ---
  
  # Analysis of the lung epithelium
  
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
# ```{r}
epi <- readRDS("/data/nas1/xiangxiaoyu_OD/1.BPD+MSCs/1.scRNA/01_Data_Download/4_split_celltypes/post_sct_split_epi.rds")

# Plot epithelial data feature counts, RNA counts, and % mitochondrial genes
# ```{r qc metrics violin}
# epi_metric <- VlnPlot(epi,
#                       c("nFeature_RNA", "nCount_RNA", "percent.mt"),
#                       group.by = "condition")
# 方法B：手动ggplot2
library(ggplot2)
library(dplyr)
library(tidyr)
# 提取元数据
epi_meta_data <- epi@meta.data
# 转换为长格式用于ggplot2
plot_data <- epi_meta_data %>%
  select(condition, nFeature_RNA, nCount_RNA, percent.mt) %>%
  pivot_longer(cols = c(nFeature_RNA, nCount_RNA, percent.mt), 
               names_to = "metric", values_to = "value")
# 创建小提琴图
epi_metric <- ggplot(plot_data, aes(x = condition, y = value, fill = condition)) +
  geom_violin(scale = "width", trim = TRUE) +
  facet_wrap(~ metric, scales = "free_y", ncol = 1) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("./1_epithelium/epi_metric_testing1.png",
       plot = epi_metric,width = 7,height = 9)
ggsave("./1_epithelium/epi_metric_testing1.pdf",
       plot = epi_metric,width = 7,height = 9)


## UMAP and clustering part 1
# ```{r, fig.width = 12}
epi <- cluster_pca_umap(epi,
                        k_param = 10,
                        dims_umap = 1:15,
                        dims_neighbors = 1:15,
                        cluster_res = 1.2,
                        harmony = TRUE)

p_cluster <- DimPlot(epi, label = TRUE) +
  umap_theme() +
  theme(aspect.ratio = 1)

epi_cluster <- p_cluster +
  plot_annotation("Epithelium")

ggsave("./1_epithelium/cluster1.png",
       plot = epi_cluster,width = 6,height = 5)
ggsave("./1_epithelium/cluster1.pdf",
       plot = epi_cluster,width = 6,height = 5)

# ```{r}
marker_genes <- c("EPCAM", "PECAM1", "PTPRC", "COL1A1", # Broad subtype
                  "AGER", "RTKN2", #AT1
                  "SFTPC", "ABCA3", "SFTPD", "LAMP3", #AT2
                  "HOPX", # AT1 and AT2
                  "KRT5", "KRT17", "TP63", # Basal
                  "FOXJ1", "TMEM190", "RSPH1", # Multiciliated
                  "SCGB3A2", "SOX4", "MGP", # RASC
                  "SCGB3A1", "MUC5B", # Goblet
                  "MKI67", "TOP2A", "CDC20", "CDK1") # Prolif and extra
epi_markers <- DotPlot(epi,
                       features = marker_genes) +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1))
ggsave("./1_epithelium/markers1.png",
       plot = epi_markers,width = 8,height = 9)
ggsave("./1_epithelium/markers1.pdf",
       plot = epi_markers,width = 8,height = 9)

# ```{r metrics}
# vlnplot_metrics <- VlnPlot(epi,
#                            features = c("nCount_RNA",
#                                         "nCount_SCT",
#                                         "percent.mt",
#                                         "percent.rb"),
#                            group.by = "seurat_clusters")
# 创建 violin plot
library(ggplot2)
library(tidyr)
library(dplyr)
library(patchwork)

# 提取需要的元数据
plot_data_2 <- epi@meta.data %>%
  select(seurat_clusters, nCount_RNA, nCount_SCT, percent.mt, percent.rb)

# 转换为长格式
plot_data_long_2 <- plot_data_2 %>%
  pivot_longer(cols = c(nCount_RNA, nCount_SCT, percent.mt, percent.rb),
               names_to = "metric",
               values_to = "value")

# 确保 seurat_clusters 是因子（为了正确的绘图顺序）
plot_data_long_2$seurat_clusters <- as.factor(plot_data_long_2$seurat_clusters)

vlnplot_metrics <- ggplot(plot_data_long_2, aes(x = seurat_clusters, y = value, fill = seurat_clusters)) +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.7) +
  geom_boxplot(width = 0.2, fill = "white", alpha = 0.5, outlier.size = 0.5) +
  facet_wrap(~ metric, scales = "free_y", ncol = 2) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"  # 隐藏图例，因为x轴已经显示聚类信息
  ) +
  labs(
    title = "Quality Metrics by Cluster",
    x = "Seurat Clusters",
    y = "Value"
  )

ggsave("./1_epithelium/metrics1.png", vlnplot_metrics, width = 12, height = 8)
ggsave("./1_epithelium/metrics1.pdf", vlnplot_metrics, width = 12, height = 8, dpi = 300)

# Remove mesenchymal, Endothelial, epithelial, and immune doublets. 
# Remove non-biological clusters.
# ```{r}
epi_sub <- subset(epi,
                  idents = c(15),
                  invert = TRUE)

## UMAP and clustering part 2
# ```{r, fig.width = 12}
epi_sub <- cluster_pca_umap(epi_sub,
                            k_param = 10,
                            dims_umap = 1:15,
                            dims_neighbors = 1:15,
                            cluster_res = 1.2)

p_cluster <- DimPlot(epi_sub, label = TRUE) +
  umap_theme() +
  theme(aspect.ratio = 1)

epi_cluster2 <- p_cluster +
  plot_annotation("Epithelium")

ggsave("./1_epithelium/cluster2.png",width = 7,height = 9,plot = epi_cluster2)
ggsave("./1_epithelium/cluster2.pdf",width = 7,height = 9,plot = epi_cluster2)

# ```{r}
marker_genes <- c("EPCAM", "PECAM1", "PTPRC", "COL1A1",
                  "AGER", "RTKN2",
                  "SFTPC", "ABCA3", "SFTPD", "LAMP3",
                  "HOPX",
                  "KRT5", "KRT17", "TP63",
                  "FOXJ1", "TMEM190", "RSPH1",
                  "SCGB3A2", "SOX4", "MGP",
                  "SCGB3A1", "MUC5B",
                  "NKX2-1", "CLDN4", "CDKN2A", "ANKRD1", "AREG", "NKD1", "LEF1",
                  "SOX9", "SLC34A2", "LCN2", "CLDN18", "SLPI", "SFTPA1",
                  "HMGA2", "AXIN2", "KRT8", "ID2", "NFKBIA", "CD44")
epi_dotplot2 <- DotPlot(epi_sub,
                        features = marker_genes) +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1))
ggsave("./1_epithelium/markers2.png",width = 12,height = 9,plot = epi_dotplot2)
ggsave("./1_epithelium/markers2.pdf",width = 12,height = 9,plot = epi_dotplot2)

# ```{r metrics2}
# 提取需要的元数据
plot_data_2_sub <- epi_sub@meta.data %>%
  select(seurat_clusters, nCount_RNA, nCount_SCT, percent.mt, percent.rb)

# 转换为长格式
plot_data_long_2_sub <- plot_data_2_sub %>%
  pivot_longer(cols = c(nCount_RNA, nCount_SCT, percent.mt, percent.rb),
               names_to = "metric",
               values_to = "value")

# 确保 seurat_clusters 是因子（为了正确的绘图顺序）
plot_data_long_2_sub$seurat_clusters <- as.factor(plot_data_long_2_sub$seurat_clusters)

vlnplot_metrics_sub <- ggplot(plot_data_long_2_sub, aes(x = seurat_clusters, y = value, fill = seurat_clusters)) +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.7) +
  geom_boxplot(width = 0.2, fill = "white", alpha = 0.5, outlier.size = 0.5) +
  facet_wrap(~ metric, scales = "free_y", ncol = 2) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"  # 隐藏图例，因为x轴已经显示聚类信息
  ) +
  labs(
    title = "Quality Metrics by Cluster",
    x = "Seurat Clusters",
    y = "Value"
  )

ggsave("./1_epithelium/metrics2.png", vlnplot_metrics_sub, width = 12, height = 8)
ggsave("./1_epithelium/metrics2.pdf", vlnplot_metrics_sub, width = 12, height = 8, dpi = 300)


# # ```{r}
# epi_sub2 <- subset(epi_sub,
#                    idents = c(14),
#                    invert = TRUE)
# 
# # ```{r}
# N_WORKERS <- 20
# plan("multicore", workers = N_WORKERS)
# epi_sub2_sct <- SCTransform(epi_sub2,
#                             variable.features.n = 2000,
#                             batch_var = "dataset",
#                             conserve.memory = TRUE,
#                             vars.to.regress = c("percent.mt",
#                                                 "percent.rb"),
#                             method = "glmGamPoi",
#                             ncells = NULL)
# 
# ## UMAP and clustering part 3
# # ```{r soupx norm}
# var_feat_backup <- VariableFeatures(epi_sub2_sct)
# epi_sub2_sct <- cluster_pca_umap(epi_sub2_sct,
#                                  k_param = 10,
#                                  dims_umap = 1:15,
#                                  dims_neighbors = 1:15,
#                                  cluster_res = 1.5)
# 
# # ```{r}
# p_cluster <- DimPlot(epi_sub2_sct, label = TRUE) +
#   umap_theme() +
#   theme(aspect.ratio = 1)
# 
# epi_cluster3 <- p_cluster +
#   plot_annotation("Epithelium")
# ggsave("./1_epithelium/cluster3.png")
# 
# # ```{r}
# marker_genes <- c("EPCAM", "PECAM1", "PTPRC", "COL1A1",
#                   "AGER", "RTKN2",
#                   "SFTPC", "ABCA3", "SFTPD", "LAMP3",
#                   "HOPX",
#                   "KRT5", "KRT17", "TP63",
#                   "FOXJ1", "TMEM190", "RSPH1",
#                   "SCGB3A2", "SOX4", "MGP",
#                   "SCGB3A1", "MUC5B", "MKI67",
#                   "NKX2-1", "CLDN4", "CDKN2A", "ANKRD1", "AREG", "NKD1", "LEF1",
#                   "SOX9", "SLC34A2", "LCN2", "CLDN18", "SLPI", "SFTPA1",
#                   "HMGA2", "AXIN2", "KRT8", "ID2", "NFKBIA", "CD44")
# epi_dotplot3 <- DotPlot(epi_sub2_sct,
#                         features = marker_genes) +
#   theme(axis.text.x = element_text(angle = 45,
#                                    hjust = 1))
# ggsave("./1_epithelium/markers3.pdf",
#        width = 15)
# 
# # ```{r metrics3}
# vlnplot_metrics <- VlnPlot(epi_sub2_sct,
#                            features = c("nCount_RNA",
#                                         "nCount_SCT",
#                                         "percent.mt",
#                                         "percent.rb"),
#                            group.by = "seurat_clusters")
# ggsave("./1_epithelium/metrics3.png",
#        plot = vlnplot_metrics)
# FeaturePlot(epi_sub2_sct,
#             c("SCGB3A2"))
# 
# # ```{r}
# epi_sub3_sct <- subset(epi_sub2_sct,
#                        idents = c(13, 15, 21),
#                        invert = TRUE)
# 
# 
# # ```{r}
N_WORKERS <- 20
plan("multicore", workers = N_WORKERS)
epi_sub3_sct2 <- epi_sub
epi_sub3_sct2 <- SCTransform(epi_sub3_sct2,
                             variable.features.n = 2000,
                             conserve.memory = TRUE,
                             vars.to.regress = c("percent.mt",
                                                 "percent.rb"),
                             method = "glmGamPoi",
                             ncells = NULL)


## UMAP and clustering part 4
# ```{r}
var_feat_backup <- VariableFeatures(epi_sub3_sct2)
epi_sub3_sct2 <- cluster_pca_umap(epi_sub3_sct2,
                                  k_param = 10,
                                  dims_umap = 1:15,
                                  dims_neighbors = 1:15,
                                  cluster_res = 1.2)
# ```{r}
plot_PCA <- DimHeatmap(epi_sub3_sct2,
                       dims = 1:20,
                       nfeatures = 40,
                       fast = FALSE)
ggsave("./1_epithelium/PCA2.png",
       plot_PCA,
       height = 40)

# ```{r}
p_cluster <- DimPlot(epi_sub3_sct2, label = TRUE) +
  umap_theme() +
  theme(aspect.ratio = 1)

epi_cluster3 <- p_cluster +
  plot_annotation("Epithelium")
ggsave("./1_epithelium/cluster4.png",width = 7,height = 8)
ggsave("./1_epithelium/cluster4.pdf",width = 7,height = 8)

# ```{r}
marker_genes <- c("EPCAM", "PECAM1", "PTPRC", "COL1A1",
                  "AGER", "RTKN2",
                  "SFTPC", "ABCA3", "SFTPD", "LAMP3",
                  "HOPX",
                  "KRT5", "KRT17", "TP63",
                  "FOXJ1", "TMEM190", "RSPH1",
                  "SCGB3A2", "SOX4", "MGP",
                  "SCGB3A1", "MUC5B", "MKI67",
                  "NKX2-1", "CLDN4", "CDKN2A", "ANKRD1", "AREG", "NKD1", "LEF1",
                  "SOX9", "SLC34A2", "LCN2", "CLDN18", "SLPI", "SFTPA1",
                  "HMGA2", "AXIN2", "KRT8", "ID2", "NFKBIA", "CD44")
epi_dotplot3 <- DotPlot(epi_sub3_sct2,
                        features = marker_genes) +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1))
ggsave("./1_epithelium/markers4.pdf",
       width = 15,height = 9)
ggsave("./1_epithelium/markers4.png",
       width = 15,height = 9)

## Save the DEgenes per cluster
# ```{r}
N_WORKERS <- 20
plan("multicore", workers = N_WORKERS)
filename <- "./1_epithelium/pre_annotated_epi_for_deg.rds"


# 再次尝试标记基因查找
epi_sub3_sct2 <- PrepSCTFindMarkers(epi_sub3_sct2)

if (!file.exists(filename)) {
  markers <- parallel_find_all_markers(epi_sub3_sct2)
  
  names(markers) <- levels(Idents(epi_sub3_sct2))
  saveRDS(markers, filename)
} else {
  markers <- readRDS(filename)
}

wb_markers_names <- levels(Idents(epi_sub3_sct2))
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
             file = "./1_epithelium/epi_cluster_deg.xlsx",
             overwrite = TRUE)

## Annotation
# ```{r}
epi_relabel <- epi_sub3_sct2
epi_relabel$louvain <- as.factor(epi_sub3_sct2$seurat_clusters)
epi_relabel$celltype <- as.factor(epi_sub3_sct2$seurat_clusters)

idxes <- levels(epi_relabel$celltype)
labels <- rep("AT2", length(levels(epi_relabel$celltype)))

labels[idxes %in% c(0,6,7,15,22)] <- "AT1"
labels[idxes %in% c(5,9)] <- "AT2"
labels[idxes %in% c(10)] <- "Basal"
labels[idxes %in% c(1,3,11,13,14,18,19,21)] <- "Multiciliated"
labels[idxes %in% c(2,4,8,17,20)] <- "RASC"
labels[idxes %in% c(12)] <- "Secretory MUC5B"
labels[idxes %in% c(16)] <- "Secretory -3A1, -3A2"

levels(epi_relabel$celltype) <- labels
Idents(epi_relabel) <- epi_relabel$celltype

# ```{r}
DimPlot(epi_relabel,
        group.by = "celltype") +
  scale_color_manual(values = color_scanpy_default)

# ```{r}
p_cluster <- DimPlot(epi_relabel,
                     group.by = "celltype",
                     label = TRUE) +
  umap_theme() +
  scale_colour_manual(name = "Cluster",
                      values = color_category_20) +
  theme(aspect.ratio = 1)

p_cluster + plot_annotation("Epithelium")


# Reorder the clusters
# ```{r}
epi_relabel$celltype <- ordered(as.factor(epi_relabel$celltype),
                                c("AT1",
                                  "AT2",
                                  "RASC",
                                  "Secretory -3A1, -3A2",
                                  "Secretory MUC5B",
                                  "Basal",
                                  "Multiciliated"))

## Save data for next step merge
# ```{r}
plan("multicore", workers = N_WORKERS)
filename <- "./1_epithelium/annotated_epi.rds"
if (!file.exists(filename)) {
  saveRDS(epi_relabel,
          filename,
          compress = "bzip2")
}
