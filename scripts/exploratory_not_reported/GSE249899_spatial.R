###############################################################################
# Spatial transcriptomic dataset GSE249899
#
# Source      : r.03_Cell_prop.r (lines 8532-8721)
# Input       : GSE249899
# Output      : exploratory figures
# Note        : Not reported in the manuscript.
#
# Code below is unchanged from the original analysis script.
###############################################################################
# -----------------------------------------------------------------------------
# PATHS: this file keeps the absolute paths of the analysis server on which the
# study was run (/data/nas1/xiangxiaoyu_OD/1.BPD+MSCs/1.scRNA/... and
# /home/xxyloveme/study/BPD+MSCs/...). Set PROJ_DIR in config/paths.R and
# replace those prefixes with your own project directory before running.
# -----------------------------------------------------------------------------


#############空转测序验证 GSE249899##########
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(hdf5r)
library(png)
library(jsonlite)
library(jpeg)
library(ggplot2)
library(dplyr)

.libPaths(c("/data/home/xiangxiaoyu_OD/.conda/envs/sc_analysis/lib/R/library/"))

.libPaths(c("/data/home/xiangxiaoyu_OD/.conda/envs/sc_analysis/lib/R/library/",
            "/data/home/xiangxiaoyu_OD/R/x86_64-conda-linux-gnu-library/4.3",
            "/data/nas2/software/miniconda3/envs/R4.3.3/lib/R/library",
            "/data/nas2/software/miniconda/lib/R/library"))

# 1. 读取.h5文件
cat("正在读取空间转录组数据...\n")
# 读取A1样本
setwd("/data/nas1/xiangxiaoyu_OD/1.BPD+MSCs/1.scRNA/09_ST/datasets/A1/spatial/")
a1_data <- Read10X_h5("/data/nas1/xiangxiaoyu_OD/1.BPD+MSCs/1.scRNA/09_ST/datasets/A1/filtered_feature_bc_matrix.h5")
# 创建标准的单细胞Seurat对象，不使用"Spatial" assay
a1_seurat <- CreateSeuratObject(
  counts = a1_data,
  project = "A1",
  assay = "RNA",  # 使用标准的RNA assay
  min.cells = 3,
  min.features = 200)

# 读取D1样本  
d1_data <- Read10X_h5("/data/nas1/xiangxiaoyu_OD/1.BPD+MSCs/1.scRNA/09_ST/datasets/D1/filtered_feature_bc_matrix.h5")
d1_seurat <- CreateSeuratObject(counts = a1_data,
                                project = "D1",
                                assay = "RNA",  # 使用标准的RNA assay
                                min.cells = 3,
                                min.features = 200)
scRNA.combined.No1 <- merge(a1_seurat, y = d1_seurat, add.cell.ids = c("A1", "D1"))
scRNA.combined.No1[["percent.mt"]] <- PercentageFeatureSet(scRNA.combined.No1, pattern = "^MT-")
# 简单质量控制可视化
# 使用基本的ggplot2创建质量控制图，避免使用VlnPlot
qc_data <- data.frame(Group = scRNA.combined.No1$orig.ident,
  nFeature_RNA = scRNA.combined.No1$nFeature_RNA,
  nCount_RNA = scRNA.combined.No1$nCount_RNA,
  percent_mt = scRNA.combined.No1$percent.mt)

# 创建质量控制图
p1 <- ggplot(qc_data, aes(x = Group, y = nFeature_RNA)) +
  geom_violin(aes(fill = Group)) +
  geom_boxplot(width = 0.1) +
  labs(title = "Number of Features per Cell", y = "nFeature_RNA") +
  theme_minimal()

p2 <- ggplot(qc_data, aes(x = Group, y = nCount_RNA)) +
  geom_violin(aes(fill = Group)) +
  geom_boxplot(width = 0.1) +
  labs(title = "UMI Counts per Cell", y = "nCount_RNA") +
  theme_minimal()

p3 <- ggplot(qc_data, aes(x = Group, y = percent_mt)) +
  geom_violin(aes(fill = Group)) +
  geom_boxplot(width = 0.1) +
  labs(title = "Mitochondrial Gene Percentage", y = "percent.mt") +
  theme_minimal()

# 数据过滤（根据实际情况调整阈值）
# scRNA.combined.No1 <- subset(scRNA.combined.No1, 
#                              subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 20)
# 数据标准化
scRNA.combined.No1 <- NormalizeData(scRNA.combined.No1)
# 寻找高变基因
scRNA.combined.No1 <- FindVariableFeatures(scRNA.combined.No1, selection.method = "vst", 
                                           nfeatures = 2000)
# 数据缩放
scRNA.combined.No1 <- ScaleData(scRNA.combined.No1)
# PCA分析
scRNA.combined.No1 <- RunPCA(scRNA.combined.No1, 
                             features = VariableFeatures(object = scRNA.combined.No1))
# 细胞聚类
scRNA.combined.No1 <- FindNeighbors(scRNA.combined.No1, dims = 1:10)
scRNA.combined.No1 <- FindClusters(scRNA.combined.No1, resolution = 0.4)
scRNA.combined.No1 <- FindClusters(object = scRNA.combined.No1, 
                                   resolution = c(seq(0, 1.0, 0.1)))
library(clustree)
clustree(scRNA.combined.No1@meta.data, prefix = "RNA_snn_res.")

#resolution = 0.4
# UMAP降维
scRNA.combined.No1 <- RunUMAP(scRNA.combined.No1, dims = 1:10)
scRNA.combined.No1 <- RunTSNE(scRNA.combined.No1, dims = 1:10)
Idents(scRNA.combined.No1) <- "RNA_snn_res.0.4"
# 可视化聚类结果
DimPlot(scRNA.combined.No1, reduction = "umap", label = TRUE)

##细胞类型注释
setwd("/data/nas1/xiangxiaoyu_OD/1.BPD+MSCs/1.scRNA/09_ST/")
scRNA.combined.No1 <- JoinLayers(scRNA.combined.No1)
cluster_markers <- FindAllMarkers(scRNA.combined.No1, only.pos = TRUE, min.pct = 0.25, 
                                  logfc.threshold = 0.25)
top_markers <- cluster_markers %>% 
  filter(cluster == !!cluster) %>% 
  arrange(desc(avg_log2FC)) %>% 
  head(20) %>% 
  pull(gene)
top_markers <- cluster_markers %>%
  group_by(cluster) %>% #按照clusetr分类
  top_n(n = 10, wt = avg_log2FC) %>% #按照log2fc取前十
  dplyr::arrange(desc(avg_log2FC))%>%
  dplyr::arrange(cluster)


lineage_markers <- list(
  Epithelial = c("KRT19", "SFTPC", "SCGB1A1", "SCGB3A2", "MUC5B", "KRT5"),
  Endothelial = c("PECAM1", "VWF", "CDH5", "EMCN", "ACKR1", "CCL21"),
  Mesenchymal = c("PDGFRB", "ACTA2", "COL1A1", "TCF21", "SERPINF1", "WNT5A", "COX4I2", "TAGLN"),
  Immune = c("PTPRC", "CD3D", "CD79A", "NKG7", "C1QA", "CD14", "CSF2RA", "MS4A1", "CPA3", "MZB1", "S100A8")
)
lineage_markers_genes <- unique(c(lineage_markers$Epithelial,lineage_markers$Endothelial,
                                  lineage_markers$Mesenchymal,lineage_markers$Immune))
DotPlot(scRNA.combined.No1, features = lineage_markers_genes)+
  RotatedAxis()

# 定义您感兴趣的基因
genes_of_interest <- c("TLE5", "ACTG1", "TSC22D2", "HSPA1A", "PHLDA1", "SLC26A7", "MT1A", "NFAT5")
# 检查哪些基因在数据中
available_genes <- genes_of_interest[genes_of_interest %in% rownames(scRNA.combined.No1)]

# 可视化可用的基因
if(length(available_genes) > 0) {
  # 1. UMAP上的基因表达图
  for(gene in available_genes) {
    p <- FeaturePlot(scRNA.combined.No1, features = gene, order = TRUE) + 
      ggtitle(paste(gene, "Expression"))
    print(p)
  }
  
  # 2. 使用ggplot2创建小提琴图
  for(gene in available_genes) {
    # 提取基因表达数据
    expr_data <- FetchData(scRNA.combined.No1, vars = c(gene, "seurat_clusters"))
    colnames(expr_data) <- c("expression", "cluster")
    
    p <- ggplot(expr_data, aes(x = cluster, y = expression, fill = cluster)) +
      geom_violin(scale = "width") +
      geom_boxplot(width = 0.1, fill = "white") +
      labs(title = paste(gene, "Expression by Cluster"),
           x = "Cluster", y = "Expression Level") +
      theme_minimal() +
      theme(legend.position = "none")
    print(p)
  }
  
  # 3. 计算基因表达统计
  expression_stats <- data.frame(
    Gene = available_genes,
    Expressed_Cells = sapply(available_genes, function(gene) {
      sum(GetAssayData(scRNA.combined.No1, slot = "counts")[gene, ] > 0)
    }),
    Total_Cells = ncol(scRNA.combined.No1),
    Expression_Rate = sapply(available_genes, function(gene) {
      round(sum(GetAssayData(scRNA.combined.No1, slot = "counts")[gene, ] > 0) / ncol(scRNA.combined.No1) * 100, 2)
    })
  )
  
  print(expression_stats)
  
  # 可视化表达率
  p <- ggplot(expression_stats, aes(x = reorder(Gene, -Expression_Rate), y = Expression_Rate)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    theme_minimal() +
    labs(title = "基因表达率",
         x = "基因", y = "表达率 (%)") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  print(p)
}

# 保存处理后的数据
saveRDS(scRNA.combined.No1, file = "scRNA.combined.No1_processed.rds")
cat("数据已保存为 scRNA.combined.No1_processed.rds\n")

# 导出基因表达统计数据
if(length(available_genes) > 0) {
  write.csv(expression_stats, "gene_expression_stats.csv", row.names = FALSE)
  cat("基因表达统计数据已保存为 gene_expression_stats.csv\n")
}

getwd()
