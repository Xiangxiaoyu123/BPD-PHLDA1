###############################################################################
# Quality control, SCTransform, Harmony, Louvain clustering and lineage annotation
#
# Source      : r.02_cellCluster.r (lines 1-640)
# Input       : seurat_list.rds
# Output      : QC figures, integrated object, 30 clusters, lineage labels (Fig. 2A-C; Supplementary Fig. 1A-D)
#
# Code below is unchanged from the original analysis script.
###############################################################################
# -----------------------------------------------------------------------------
# PATHS: this file keeps the absolute paths of the analysis server on which the
# study was run (/data/nas1/xiangxiaoyu_OD/1.BPD+MSCs/1.scRNA/... and
# /home/xxyloveme/study/BPD+MSCs/...). Set PROJ_DIR in config/paths.R and
# replace those prefixes with your own project directory before running.
# -----------------------------------------------------------------------------


#############2.单细胞数据处理##############
setwd("/data/nas1/xiangxiaoyu_OD/1.BPD+MSCs/1.scRNA/02_Cell_Clustering_and_Annotation/")
# 首先，我们需要为每个细胞计算线粒体和核糖体基因的百分比。
# 然后，我们将多个样本的Seurat对象合并成一个对象，并添加样本信息（已经存在于每个对象的meta.data中，但我们之前已经添加了"sample"列）。
# 运行SCTransform，回归掉线粒体百分比、核糖体百分比和样本批次（即样本来源）。
# 运行PCA。
# 使用Harmony整合，指定整合的变量为样本来源（因为我们已经回归了批次，所以这一步可能不需要？但文献说用了，所以我们还是做）。
# 然后进行UMAP和聚类。

# 质量控制：过滤掉线粒体基因百分比 > 10%、基因数 < 500 或 > 7000 的细胞
# SCTransform归一化：使用glmGamPoi方法，回归线粒体和核糖体基因的影响，校正样本批次效应
# PCA降维：计算主成分
# Harmony整合：校正样本间的批次效应
# UMAP可视化：基于Harmony降维结果进行非线性降维
# Louvain聚类：识别细胞簇
# 标记基因鉴定：找出每个簇的特异性标记基因
# 细胞类型注释：根据已知标记基因初步鉴定细胞类型

# 加载必要的包
library(Seurat)
library(ggplot2)
library(patchwork)
library(harmony)
library(sctransform)
library(glmGamPoi)

# 1. 质量控制过滤
# 首先合并所有样本为一个Seurat对象
merged_seurat <- merge(seurat_list[[1]], y = seurat_list[-1], 
                       add.cell.ids = names(seurat_list), 
                       project = "BPD_MSC_Study")

# 计算线粒体基因百分比
merged_seurat[["percent.mt"]] <- PercentageFeatureSet(merged_seurat, pattern = "^MT-")

# 计算核糖体基因百分比
merged_seurat[["percent.rb"]] <- PercentageFeatureSet(merged_seurat, pattern = "^RP[SL]")

# 可视化QC指标
qc_plot <- VlnPlot(merged_seurat, 
                   features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), 
                   ncol = 4, pt.size = 0.1) & 
  theme(plot.title = element_text(size = 10))

ggsave("qc_metrics_before_filtering.pdf", qc_plot, width = 12, height = 6)

# 应用QC过滤（根据文献描述）
merged_seurat <- subset(merged_seurat, 
                        subset = nFeature_RNA > 500 & 
                          nFeature_RNA < 7000 & 
                          percent.mt < 10)

# 可视化过滤后的QC指标
qc_plot_filtered <- VlnPlot(merged_seurat, 
                            features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), 
                            ncol = 4, pt.size = 0.1) & 
  theme(plot.title = element_text(size = 10))

ggsave("qc_metrics_after_filtering.pdf", qc_plot_filtered, width = 12, height = 6)

# 2. 使用SCTransform进行归一化和方差稳定化
# 根据文献，使用glmGamPoi方法，回归线粒体和核糖体基因的影响
# 首先检查您的 Seurat 对象结构
print(merged_seurat)
cat("\nAssays available:", paste(Assays(merged_seurat), collapse = ", "), "\n")

# 如果 SCT assay 不存在，我们需要重新运行 SCTransform
# 首先确保默认 assay 是 RNA
DefaultAssay(merged_seurat) <- "RNA"

# 检查sample列的类型和水平
class(merged_seurat$sample)
levels(merged_seurat$sample)
unique(merged_seurat$sample)

merged_seurat$sample <- as.factor(merged_seurat$sample)

merged_seurat <- SCTransform(merged_seurat, 
                             method = "glmGamPoi",
                             vars.to.regress = c("percent.mt", "percent.rb"),
                             # 移除 batch_var = "sample"
                             ncells = NULL,
                             verbose = TRUE)

# 3. 降维和特征选择
# 在使用SCTransform后，我们需要明确指定使用SCT assay的数据
# 首先确保默认assay是SCT
DefaultAssay(merged_seurat) <- "SCT"

# 运行PCA - 明确指定使用SCT assay
merged_seurat <- RunPCA(merged_seurat, 
                        assay = "SCT", 
                        npcs = 50, 
                        verbose = FALSE)
# 可视化PCA结果
pca_plot <- ElbowPlot(merged_seurat, ndims = 50)
ggsave("pca_elbow_plot.pdf", pca_plot, width = 8, height = 6)

# 4. 使用Harmony进行批次效应校正
library(harmony)
merged_seurat <- RunHarmony(merged_seurat,group.by.vars = "sample")
reduction.use <- "harmony"


# 5. 非线性降维（UMAP）
merged_seurat <- RunUMAP(merged_seurat, 
                         reduction = "harmony", 
                         dims = 1:30, 
                         reduction.name = "umap.harmony")

# 6. 聚类分析
merged_seurat <- FindNeighbors(merged_seurat, 
                               reduction = "harmony", 
                               dims = 1:30)

merged_seurat <- FindClusters(merged_seurat, 
                              resolution = seq(from = 0.1, to = 1.0,by = 0.1), 
                              algorithm = 1) # Louvain算法

library(clustree)
my80colors <- unique(c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
                       '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
                       '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
                       '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
                       '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
                       '#968175'))
pdf(file = "05.clustree.pdf",height = 6,width = 15)
clustree(merged_seurat) +
  theme(legend.position = "bottom") + 
  scale_color_manual(values = my80colors)+
  scale_edge_color_continuous(low = "grey80", high = "red")
dev.off()

#确定resoluation=0.6
Idents(merged_seurat) <- "SCT_snn_res.0.6"

# 7. 可视化结果
# UMAP按聚类着色
umap_cluster <- DimPlot(merged_seurat, 
                        reduction = "umap.harmony", 
                        group.by = "SCT_snn_res.0.6",
                        label = TRUE) +
  ggtitle("UMAP by Cluster")


# UMAP按样本着色
umap_sample <- DimPlot(merged_seurat, 
                       reduction = "umap.harmony", 
                       group.by = "sample") +
  ggtitle("UMAP by Sample")


# 合并可视化
combined_umap <- umap_cluster + umap_sample
ggsave("umap_visualization.pdf", combined_umap, width = 14, height = 6)

saveRDS(merged_seurat,file = "merged_seurat.rds")


# 8. 寻找标记基因
# 识别每个簇的标记基因
# 1. 首先识别四大主要细胞类型
# 创建四大细胞类型的标记基因列表
# 1. 首先识别四大主要细胞类型（已完成）
# 创建四大细胞类型的标记基因列表
lineage_markers <- list(
  Epithelial = c("EPCAM", "KRT19", "SFTPC", "SCGB1A1", "SCGB3A2", "MUC5B", "KRT5", "FOXJ1"),
  Endothelial = c("PECAM1", "VWF", "CDH5", "EMCN", "ACKR1", "CCL21"),
  Mesenchymal = c("PDGFRB", "ACTA2", "COL1A1", "TCF21", "SERPINF1", "WNT5A", "COX4I2", "TAGLN"),
  Immune = c("PTPRC", "CD3D", "CD79A", "NKG7", "C1QA", "CD14", "CSF2RA", "MS4A1", "CPA3", "MZB1", "S100A8")
)

# 计算每个细胞的四大细胞类型得分
for (lineage in names(lineage_markers)) {
  merged_seurat <- AddModuleScore(
    merged_seurat,
    features = list(lineage_markers[[lineage]]),
    name = paste0(lineage, "_Score")
  )
}

# 为每个细胞分配主要细胞类型
lineage_scores <- FetchData(merged_seurat, vars = c("Epithelial_Score1", "Endothelial_Score1", 
                                                    "Mesenchymal_Score1", "Immune_Score1"))
merged_seurat$lineage <- apply(lineage_scores, 1, function(x) {
  names(which.max(x))
})

# 2. 对每个主要细胞类型进行迭代细分
# 创建所有细胞亚型的标记基因列表
celltype_markers <- list(
  # Epithelial subtypes
  "AT1" = c("AGER", "RTKN2", "LAMA3", "COL4A3", "PDPN"),
  "AT2" = c("SFTPC", "ABCA3", "NAPSA", "SFTPD", "LAMP3"),
  "RASC" = c("SCGB3A2", "SGMS2", "ERBB4", "CDK14"),
  "Secretory_3A1_3A2" = c("SCGB3A2", "SCGB3A1", "SCGB1A1", "KLK11"),
  "Secretory_MUC5B" = c("MUC5B", "BPIFB1", "MSMB", "CYP2F1"),
  "Basal" = c("KRT17", "KRT15", "KRT5", "MIR205HG", "TP63"),
  "Multiciliated" = c("FOXJ1", "CAPS", "DNAH12", "RSPH1", "CFAP299", "LRRIQ1", "SNTN", "TMEM190", "ERICH3", "HYDIN"),
  
  # Endothelial subtypes
  "aCap" = c("EDNRB", "SOSTDC1", "HPGD", "EMCN", "TSPAN12", "S100A3", "APLN"),
  "gCap" = c("GPIHBP1", "FCN3", "BTNL9", "NTRK2", "BMP6", "IL7R"),
  "abCap" = c("ANKRD1", "CCN1", "CCN2", "SERPINE1", "KITLG", "EDN1", "TNFRSF12A", "GDF15", "APLN"),
  "Arterial_EC" = c("DKK2", "GJA5", "CXCL12", "THSD7A", "SERPINE2", "BMX", "TSPAN2", "EGFL8"),
  "Pulmonary_venous_EC" = c("ACKR1", "HDAC9", "CPE", "PTGIS", "C7", "SELP", "FAM155A"),
  "Systemic_venous_EC" = c("ACKR1", "MCTP1", "CCL14", "PLVAP", "ADGRL4", "SPRY1", "COL15A1"),
  "Lymphatic" = c("CCL21", "MMRN1", "TFF3", "PKHD1L1", "PROX1", "RELN", "TBX1"),
  
  # Mesenchymal subtypes
  "Alveolar_FB" = c("TCF21", "NPNT", "INMT", "WNT2", "BMP5", "FAT3", "MEOX2", "NTNG1", "RSPO2", "FGFR4"),
  "Activated_FB" = c("CTHRC1", "TIMP1", "POSTN"),
  "Adventitial_FB" = c("SERPINF1", "PI16", "MFAP5", "SCARA5", "CCDC80", "FBLN1", "C1R", "C3", "DCLK1", "PRRX1", "IGFBP6", "GSN", "PDGFRL"),
  "Alveolar_MyoFB" = c("WNT5A", "FGF18", "PDGFRA", "CXCL14", "ASPN", "DACH2", "LTBP2", "CCDC68", "COL24A1", "LRP4", "ANO4"),
  "Ductal_MyoFB" = c("WNT5A", "FGF18", "HHIP", "ZNF385D", "HPSE2", "MYH11", "ACTA2", "DOCK3", "DES", "PRUNE2", "SCUBE1", "LGR6"),
  "Pericyte" = c("COX4I2", "CSPG4", "LAMC3", "TRPC6", "HIGD1B", "FAM162B", "KCNK3", "AFF2", "ADARB2", "LRRTM4"),
  "VSMC" = c("ACTA2", "TAGLN", "TPM2", "NTRK3", "ADIRF", "MUSTN1", "RCAN2", "RGS6", "CSMD1", "PLN", "ECRG4"),
  
  # Immune subtypes
  "T_Cell" = c("CD8A", "THEMIS", "CD6", "CAMK4", "BCL11B", "IL7R", "CD2", "LEF1", "ICOS", "CD3D"),
  "NK_Cell" = c("NKG7", "CCL4", "PRF1", "GNLY", "KLRD1", "KLRB1", "KLRF1", "SPON2", "FGFBP2", "MYOM2", "LINGO2"),
  "NKT_Cell" = c("KLRB1", "CD7", "CTSW", "TNFRSF18", "CCL5", "IL7R", "GZMK"),
  "Alveolar_Macrophage" = c("C1QA", "C1QB", "C1QC", "APOC1", "ACP5", "CAPG", "MARCO", "LYZ", "APOE"),
  "Monocyte" = c("S100A12", "FCN1", "VCAN", "KYNU", "CD14"),
  "pDC" = c("RHEX", "IRF4", "PPP1R14B", "PTPRS", "CCDC50", "RUBCN", "IRF4", "IRF7", "PLD4", "LINC01478"),
  "cDC" = c("CSF2RA", "CPVL", "CD86", "FCER1A", "CD1C", "FCGR2B", "CLEC10A", "HLA-DRA", "FPR3"),
  "B_Cell" = c("BANK1", "MS4A1", "IGHM", "CD79A", "CD79B", "EBF1", "VPREB3", "FCRL1", "PAX5", "IGHD", "COL19A1", "TNFRSF13C", "CD22"),
  "Mast_cell" = c("CPA3", "KIT", "TPSD1", "MS4A2", "TPSB2", "TPSAB1", "HPGDS", "CTSG", "CADPS"),
  "Plasma_cell" = c("MZB1", "JCHAIN", "SSR4", "TXNDC5", "SEC11C", "FKBP11", "PRDX4", "IGHG1", "IGHG3", "POU2AF1"),
  "Basophil" = c("ATP10D", "LINC02458", "HDC", "HSH2D", "NLK", "INPP5A", "ENPP3", "TAFA2", "TRPM6", "IL4"),
  "Neutrophil" = c("S100A8", "S100A9", "CSF3R", "NAMPT", "AQP9", "LIMK2", "FBXL5", "LUCAT1", "FCGR3B", "KCNJ15")
)

umap_plot_01 <- DimPlot(merged_seurat, 
        group.by = "SCT_snn_res.0.6", 
        label = TRUE,cols = my80colors)
umap_plot_02 <- DimPlot(merged_seurat, 
        group.by = "lineage", 
        label = TRUE,cols = my80colors)
ggsave("01.umap_res_0.6.pdf", umap_plot_01, width = 10, height = 8)
ggsave("02.umap_celltype_main4.pdf", umap_plot_02, width = 10, height = 8)

##marker基因表达
FeaturePlot(merged_seurat, features = lineage_markers$Epithelial,
            reduction = "umap.harmony")
FeaturePlot(merged_seurat, features = lineage_markers$Endothelial,
            reduction = "umap.harmony")
FeaturePlot(merged_seurat, features = lineage_markers$Mesenchymal,
            reduction = "umap.harmony")
FeaturePlot(merged_seurat, features = lineage_markers$Immune,
            reduction = "umap.harmony")


# 设置默认assay为RNA，并进行归一化
DefaultAssay(merged_seurat) <- "RNA"
merged_seurat <- NormalizeData(merged_seurat, verbose = FALSE)

# 提取所有标记基因
all_markers <- unique(unlist(celltype_markers))
# 只保留在数据中存在的基因
existing_markers <- all_markers[all_markers %in% rownames(merged_seurat)]

# 按细胞类型组织标记基因
marker_groups <- list()
for (celltype in names(celltype_markers)) {
  markers <- celltype_markers[[celltype]]
  existing <- markers[markers %in% existing_markers]
  if (length(existing) > 0) {
    marker_groups[[celltype]] <- existing
  }
}

# 创建气泡图函数
create_dotplot <- function(seurat_obj, markers, group_name, lineage_type) {
  # 提取当前lineage的细胞子集
  lineage_cells <- subset(seurat_obj, subset = lineage == lineage_type)
  
  # 提取当前细胞类型的标记基因
  celltype_markers_sub <- markers[[group_name]]
  
  # 创建气泡图
  p <- DotPlot(lineage_cells, 
               features = celltype_markers_sub,
               group.by = "SCT_snn_res.0.6",
               cols = c("lightgrey", "blue"),
               dot.scale = 3) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title = element_blank()) +
    ggtitle(paste(group_name, "Markers in", lineage_type, "Cells")) +
    ylab("Cluster") +
    xlab("Genes")
  
  return(p)
}

# 为每个主要细胞类型创建气泡图
# 上皮细胞
epithelial_plots <- list()
for (epithelial_type in c("AT1", "AT2", "RASC", "Secretory_3A1_3A2", "Secretory_MUC5B", "Basal", "Multiciliated")) {
  if (epithelial_type %in% names(marker_groups)) {
    epithelial_plots[[epithelial_type]] <- create_dotplot(merged_seurat, marker_groups, epithelial_type, "Epithelial_Score1")
  }
}

# 内皮细胞
endothelial_plots <- list()
for (endothelial_type in c("aCap", "gCap", "abCap", "Arterial_EC", "Pulmonary_venous_EC", "Systemic_venous_EC", "Lymphatic")) {
  if (endothelial_type %in% names(marker_groups)) {
    endothelial_plots[[endothelial_type]] <- create_dotplot(merged_seurat, marker_groups, endothelial_type, "Endothelial_Score1")
  }
}

# 间充质细胞
mesenchymal_plots <- list()
for (mesenchymal_type in c("Alveolar_FB", "Activated_FB", "Adventitial_FB", "Alveolar_MyoFB", "Ductal_MyoFB", "Pericyte", "VSMC")) {
  if (mesenchymal_type %in% names(marker_groups)) {
    mesenchymal_plots[[mesenchymal_type]] <- create_dotplot(merged_seurat, marker_groups, mesenchymal_type, "Mesenchymal_Score1")
  }
}

# 免疫细胞
immune_plots <- list()
for (immune_type in c("T_Cell", "NK_Cell", "NKT_Cell", "Alveolar_Macrophage", "Monocyte", "pDC", "cDC", "B_Cell", "Mast_cell", "Plasma_cell", "Basophil", "Neutrophil")) {
  if (immune_type %in% names(marker_groups)) {
    immune_plots[[immune_type]] <- create_dotplot(merged_seurat, marker_groups, immune_type, "Immune_Score1")
  }
}

# 保存气泡图
# 上皮细胞
for (plot_name in names(epithelial_plots)) {
  ggsave(paste0("epithelial_", plot_name, "_dotplot.pdf"), 
         epithelial_plots[[plot_name]], width = 10, height = 8)
}

# 内皮细胞
for (plot_name in names(endothelial_plots)) {
  ggsave(paste0("endothelial_", plot_name, "_dotplot.pdf"), 
         endothelial_plots[[plot_name]], width = 10, height = 8)
}

# 间充质细胞
for (plot_name in names(mesenchymal_plots)) {
  ggsave(paste0("mesenchymal_", plot_name, "_dotplot.pdf"), 
         mesenchymal_plots[[plot_name]], width = 10, height = 8)
}

# 免疫细胞
for (plot_name in names(immune_plots)) {
  ggsave(paste0("immune_", plot_name, "_dotplot.pdf"), 
         immune_plots[[plot_name]], width = 10, height = 8)
}

# 创建每个lineage的综合气泡图（选择每个亚型的前3个标记基因）
create_summary_dotplot <- function(seurat_obj, markers, lineage_type, celltypes) {
  # 提取当前lineage的细胞子集
  lineage_cells <- subset(seurat_obj, subset = lineage == lineage_type)
  
  # 选择每个亚型的前3个标记基因
  top_markers <- c()
  for (celltype in celltypes) {
    if (celltype %in% names(markers)) {
      top_markers <- c(top_markers, head(markers[[celltype]], 3))
    }
  }
  top_markers <- unique(top_markers)
  
  # 创建气泡图
  p <- DotPlot(lineage_cells, 
               features = top_markers,
               group.by = "SCT_snn_res.0.6",
               cols = c("lightgrey", "blue"),
               dot.scale = 3) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title = element_blank()) +
    ggtitle(paste("Top Markers in", lineage_type, "Cells")) +
    ylab("Cluster") +
    xlab("Genes")
  
  return(p)
}

# 创建每个lineage的综合气泡图
epithelial_summary <- create_summary_dotplot(
  merged_seurat, marker_groups, "Epithelial_Score1",
  c("AT1", "AT2", "RASC", "Secretory_3A1_3A2", "Secretory_MUC5B", "Basal", "Multiciliated")
)

endothelial_summary <- create_summary_dotplot(
  merged_seurat, marker_groups, "Endothelial_Score1",
  c("aCap", "gCap", "abCap", "Arterial_EC", "Pulmonary_venous_EC", "Systemic_venous_EC", "Lymphatic")
)

mesenchymal_summary <- create_summary_dotplot(
  merged_seurat, marker_groups, "Mesenchymal_Score1",
  c("Alveolar_FB", "Activated_FB", "Adventitial_FB", "Alveolar_MyoFB", "Ductal_MyoFB", "Pericyte", "VSMC")
)

immune_summary <- create_summary_dotplot(
  merged_seurat, marker_groups, "Immune_Score1",
  c("T_Cell", "NK_Cell", "NKT_Cell", "Alveolar_Macrophage", "Monocyte", "pDC", "cDC", "B_Cell", "Mast_cell", "Plasma_cell", "Basophil", "Neutrophil")
)

# 保存综合气泡图
ggsave("epithelial_summary_dotplot.pdf", epithelial_summary, width = 14, height = 8)
ggsave("endothelial_summary_dotplot.pdf", endothelial_summary, width = 14, height = 8)
ggsave("mesenchymal_summary_dotplot.pdf", mesenchymal_summary, width = 14, height = 8)
ggsave("immune_summary_dotplot.pdf", immune_summary, width = 16, height = 10)







# 输出每个lineage中SCT_snn_res.0.6分群的细胞数量统计
for (lineage_type in c("Epithelial_Score1", "Endothelial_Score1", "Mesenchymal_Score1", "Immune_Score1")) {
  lineage_cells <- subset(merged_seurat, subset = lineage == lineage_type)
  cluster_stats <- table(lineage_cells$SCT_snn_res.0.6)
  write.csv(cluster_stats, paste0(lineage_type, "_cluster_statistics.csv"))
}


##
# 设置默认assay为RNA，并进行归一化
DefaultAssay(merged_seurat) <- "RNA"
merged_seurat <- NormalizeData(merged_seurat, verbose = FALSE)

# 提取所有标记基因
all_markers <- unique(unlist(celltype_markers))
# 只保留在数据中存在的基因
existing_markers <- all_markers[all_markers %in% rownames(merged_seurat)]

# 按细胞类型组织标记基因
marker_groups <- list()
for (celltype in names(celltype_markers)) {
  markers <- celltype_markers[[celltype]]
  existing <- markers[markers %in% existing_markers]
  if (length(existing) > 0) {
    marker_groups[[celltype]] <- existing
  }
}

# 查看每个lineage在SCT_snn_res.0.6中的分布
for (lineage_type in c("Epithelial_Score1", "Endothelial_Score1", "Mesenchymal_Score1", "Immune_Score1")) {
  lineage_cells <- subset(merged_seurat, subset = lineage == lineage_type)
  cat(paste0("\n", lineage_type, " distribution in SCT_snn_res.0.6:\n"))
  print(table(lineage_cells$SCT_snn_res.0.6))
}


# 创建气泡图函数 - 针对特定细胞亚群
create_dotplot_for_subclusters <- function(seurat_obj, markers, group_name, lineage_type, subclusters) {
  # 提取当前lineage的细胞子集
  lineage_cells <- subset(seurat_obj, subset = lineage == lineage_type)
  
  # 进一步提取特定亚群的细胞
  subcluster_cells <- subset(lineage_cells, subset = SCT_snn_res.0.6 %in% subclusters)
  
  # 提取当前细胞类型的标记基因
  celltype_markers_sub <- markers[[group_name]]
  
  # 创建气泡图
  p <- DotPlot(subcluster_cells, 
               features = celltype_markers_sub,
               group.by = "SCT_snn_res.0.6",
               cols = c("lightgrey", "blue"),
               dot.scale = 3) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title = element_blank()) +
    ggtitle(paste(group_name, "Markers in", lineage_type, "Subclusters:", paste(subclusters, collapse = ", "))) +
    ylab("Cluster") +
    xlab("Genes")
  
  return(p)
}




# 为每个主要细胞类型创建气泡图
# 上皮细胞 - 假设亚群为2,5,7,10,12
epithelial_subclusters <- c(2, 5, 7, 10, 12)
epithelial_plots <- list()
for (epithelial_type in c("AT1", "AT2", "RASC", "Secretory_3A1_3A2", "Secretory_MUC5B", "Basal", "Multiciliated")) {
  if (epithelial_type %in% names(marker_groups)) {
    epithelial_plots[[epithelial_type]] <- create_dotplot_for_subclusters(
      merged_seurat, marker_groups, epithelial_type, "Epithelial_Score1", epithelial_subclusters
    )
  }
}

# 内皮细胞 - 假设亚群为1,3,6,8,11
endothelial_subclusters <- c(1, 3, 6, 8, 11)
endothelial_plots <- list()
for (endothelial_type in c("aCap", "gCap", "abCap", "Arterial_EC", "Pulmonary_venous_EC", "Systemic_venous_EC", "Lymphatic")) {
  if (endothelial_type %in% names(marker_groups)) {
    endothelial_plots[[endothelial_type]] <- create_dotplot_for_subclusters(
      merged_seurat, marker_groups, endothelial_type, "Endothelial_Score1", endothelial_subclusters
    )
  }
}

# 间充质细胞 - 假设亚群为4,9,13,15,17
mesenchymal_subclusters <- c(4, 9, 13, 15, 17)
mesenchymal_plots <- list()
for (mesenchymal_type in c("Alveolar_FB", "Activated_FB", "Adventitial_FB", "Alveolar_MyoFB", "Ductal_MyoFB", "Pericyte", "VSMC")) {
  if (mesenchymal_type %in% names(marker_groups)) {
    mesenchymal_plots[[mesenchymal_type]] <- create_dotplot_for_subclusters(
      merged_seurat, marker_groups, mesenchymal_type, "Mesenchymal_Score1", mesenchymal_subclusters
    )
  }
}

# 免疫细胞 - 假设亚群为0,14,16,18,19,20
immune_subclusters <- c(0, 14, 16, 18, 19, 20)
immune_plots <- list()
for (immune_type in c("T_Cell", "NK_Cell", "NKT_Cell", "Alveolar_Macrophage", "Monocyte", "pDC", "cDC", "B_Cell", "Mast_cell", "Plasma_cell", "Basophil", "Neutrophil")) {
  if (immune_type %in% names(marker_groups)) {
    immune_plots[[immune_type]] <- create_dotplot_for_subclusters(
      merged_seurat, marker_groups, immune_type, "Immune_Score1", immune_subclusters
    )
  }
}

# 保存气泡图
# 上皮细胞
for (plot_name in names(epithelial_plots)) {
  ggsave(paste0("epithelial_", plot_name, "_subclusters_dotplot.pdf"), 
         epithelial_plots[[plot_name]], width = 10, height = 8)
}

# 内皮细胞
for (plot_name in names(endothelial_plots)) {
  ggsave(paste0("endothelial_", plot_name, "_subclusters_dotplot.pdf"), 
         endothelial_plots[[plot_name]], width = 10, height = 8)
}

# 间充质细胞
for (plot_name in names(mesenchymal_plots)) {
  ggsave(paste0("mesenchymal_", plot_name, "_subclusters_dotplot.pdf"), 
         mesenchymal_plots[[plot_name]], width = 10, height = 8)
}

# 免疫细胞
for (plot_name in names(immune_plots)) {
  ggsave(paste0("immune_", plot_name, "_subclusters_dotplot.pdf"), 
         immune_plots[[plot_name]], width = 10, height = 8)
}

# 创建每个lineage的综合气泡图（选择每个亚型的前3个标记基因）
create_summary_dotplot_for_subclusters <- function(seurat_obj, markers, lineage_type, celltypes, subclusters) {
  # 提取当前lineage的细胞子集
  lineage_cells <- subset(seurat_obj, subset = lineage == lineage_type)
  
  # 进一步提取特定亚群的细胞
  subcluster_cells <- subset(lineage_cells, subset = SCT_snn_res.0.6 %in% subclusters)
  
  # 选择每个亚型的前3个标记基因
  top_markers <- c()
  for (celltype in celltypes) {
    if (celltype %in% names(markers)) {
      top_markers <- c(top_markers, head(markers[[celltype]], 3))
    }
  }
  top_markers <- unique(top_markers)
  
  # 创建气泡图
  p <- DotPlot(subcluster_cells, 
               features = top_markers,
               group.by = "SCT_snn_res.0.6",
               cols = c("lightgrey", "blue"),
               dot.scale = 3) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title = element_blank()) +
    ggtitle(paste("Top Markers in", lineage_type, "Subclusters:", paste(subclusters, collapse = ", "))) +
    ylab("Cluster") +
    xlab("Genes")
  
  return(p)
}

# 创建每个lineage的综合气泡图
epithelial_summary <- create_summary_dotplot_for_subclusters(
  merged_seurat, marker_groups, "Epithelial_Score1",
  c("AT1", "AT2", "RASC", "Secretory_3A1_3A2", "Secretory_MUC5B", "Basal", "Multiciliated"),
  epithelial_subclusters
)

endothelial_summary <- create_summary_dotplot_for_subclusters(
  merged_seurat, marker_groups, "Endothelial_Score1",
  c("aCap", "gCap", "abCap", "Arterial_EC", "Pulmonary_venous_EC", "Systemic_venous_EC", "Lymphatic"),
  endothelial_subclusters
)

mesenchymal_summary <- create_summary_dotplot_for_subclusters(
  merged_seurat, marker_groups, "Mesenchymal_Score1",
  c("Alveolar_FB", "Activated_FB", "Adventitial_FB", "Alveolar_MyoFB", "Ductal_MyoFB", "Pericyte", "VSMC"),
  mesenchymal_subclusters
)

immune_summary <- create_summary_dotplot_for_subclusters(
  merged_seurat, marker_groups, "Immune_Score1",
  c("T_Cell", "NK_Cell", "NKT_Cell", "Alveolar_Macrophage", "Monocyte", "pDC", "cDC", "B_Cell", "Mast_cell", "Plasma_cell", "Basophil", "Neutrophil"),
  immune_subclusters
)

# 保存综合气泡图
ggsave("epithelial_summary_subclusters_dotplot.pdf", epithelial_summary, width = 14, height = 8)
ggsave("endothelial_summary_subclusters_dotplot.pdf", endothelial_summary, width = 14, height = 8)
ggsave("mesenchymal_summary_subclusters_dotplot.pdf", mesenchymal_summary, width = 14, height = 8)
ggsave("immune_summary_subclusters_dotplot.pdf", immune_summary, width = 16, height = 10)

# 输出每个lineage中特定亚群的细胞数量统计
for (lineage_type in c("Epithelial_Score1", "Endothelial_Score1", "Mesenchymal_Score1", "Immune_Score1")) {
  lineage_cells <- subset(merged_seurat, subset = lineage == lineage_type)
  
  # 根据lineage类型选择相应的亚群
  if (lineage_type == "Epithelial_Score1") {
    subclusters <- epithelial_subclusters
  } else if (lineage_type == "Endothelial_Score1") {
    subclusters <- endothelial_subclusters
  } else if (lineage_type == "Mesenchymal_Score1") {
    subclusters <- mesenchymal_subclusters
  } else {
    subclusters <- immune_subclusters
  }
  
  subcluster_cells <- subset(lineage_cells, subset = SCT_snn_res.0.6 %in% subclusters)
  cluster_stats <- table(subcluster_cells$SCT_snn_res.0.6)
  write.csv(cluster_stats, paste0(lineage_type, "_subclusters_statistics.csv"))
}

# 输出会话信息
sink("cell_annotation_session_info.txt")
sessionInfo()
sink()
