###############################################################################
# Single-cell dataset GSE137950
#
# Source      : r.03_Cell_prop.r (lines 9379-9802)
# Input       : GSE137950
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


#############单细胞数据 GSE137950_RAW##########
setwd("/data/nas1/xiangxiaoyu_OD/1.BPD+MSCs/1.scRNA/09_ST/scRNA/GSE137950_RAW/")
library(Seurat)
library(data.table) # 用于快速读取大文件，也可以用readr或base R的函数

# 读取数据文件
# 方法1: 使用fread（推荐，速度快）
neo1_data <- fread("GSM4094900_Neo1_out_gene_exon_tagged.dge.txt.gz")
neo3_data <- fread("GSM4094901_Neo3_out_gene_exon_tagged.dge.txt.gz")

# 将数据转换为矩阵格式
# 首先将第一列设置为行名（基因名），然后移除第一列
convert_to_matrix <- function(dt) {
  # 提取基因名
  gene_names <- dt[[1]]
  # 移除第一列（基因名列）
  dt <- dt[, -1]
  # 转换为矩阵
  mat <- as.matrix(dt)
  # 设置行名为基因名
  rownames(mat) <- gene_names
  return(mat)
}

# 转换两个数据集
neo1_matrix <- convert_to_matrix(neo1_data)
neo3_matrix <- convert_to_matrix(neo3_data)

# 查看转换后的矩阵
cat("Neo1矩阵维度:", dim(neo1_matrix), "\n")
cat("Neo3矩阵维度:", dim(neo3_matrix), "\n")
print(head(neo1_matrix[, 1:3]))

# 创建Seurat对象
neo1_seurat <- CreateSeuratObject(
  counts = neo1_matrix,
  project = "Neo1",
  min.cells = 3,      # 基因至少在3个细胞中表达
  min.features = 200  # 细胞至少检测到200个基因
)

neo3_seurat <- CreateSeuratObject(
  counts = neo3_matrix,
  project = "Neo3", 
  min.cells = 3,
  min.features = 200
)

# 查看创建的Seurat对象
print(neo1_seurat)
print(neo3_seurat)

# 合并两个样本（可选）
merged_seurat_No2 <- merge(
  x = neo1_seurat,
  y = neo3_seurat,
  add.cell.ids = c("Neo1", "Neo3"),
  project = "Merged_Neo")

print(merged_seurat_No2t)

merged_seurat_No2[["percent.mt"]] <- PercentageFeatureSet(merged_seurat_No2, pattern = "^MT-")
# 简单质量控制可视化
# 使用基本的ggplot2创建质量控制图，避免使用VlnPlot
qc_data <- data.frame(Group = merged_seurat_No2$orig.ident,
                      nFeature_RNA = merged_seurat_No2$nFeature_RNA,
                      nCount_RNA = merged_seurat_No2$nCount_RNA,
                      percent_mt = merged_seurat_No2$percent.mt)

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
merged_seurat_No2 <- subset(merged_seurat_No2,
                        subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10)
# 数据标准化
merged_seurat_No2 <- NormalizeData(merged_seurat_No2)
# 寻找高变基因
merged_seurat_No2 <- FindVariableFeatures(merged_seurat_No2, selection.method = "vst", 
                                      nfeatures = 2000)
# 数据缩放
merged_seurat_No2 <- ScaleData(merged_seurat_No2)
# PCA分析
merged_seurat_No2 <- RunPCA(merged_seurat_No2, 
                        features = VariableFeatures(object = merged_seurat_No2))
# 细胞聚类
merged_seurat_No2 <- FindNeighbors(merged_seurat_No2, dims = 1:10)
merged_seurat_No2 <- FindClusters(merged_seurat_No2, resolution = 0.5)

#resolution = 0.4
# UMAP降维
merged_seurat_No2 <- RunUMAP(merged_seurat_No2, dims = 1:10)
# merged_seurat_No2 <- RunTSNE(merged_seurat_No2, dims = 1:10)
# Idents(merged_seurat_No2) <- "RNA_snn_res.0.4"
# 可视化聚类结果
DimPlot(merged_seurat_No2, reduction = "umap", label = TRUE)


##细胞类型注释
merged_seurat_No2 <- JoinLayers(merged_seurat_No2)
cluster_markers_No2 <- FindAllMarkers(merged_seurat_No2, only.pos = TRUE, min.pct = 0.25, 
                                  logfc.threshold = 0.25)
top_markers_No2 <- cluster_markers_No2 %>%
  group_by(cluster) %>% #按照clusetr分类
  top_n(n = 10, wt = avg_log2FC) %>% #按照log2fc取前十
  dplyr::arrange(desc(avg_log2FC))%>%
  dplyr::arrange(cluster)

##细胞类型鉴定
lineage_markers <- list(
  Epithelial = c("KRT19", "SFTPC", "SCGB1A1", "SCGB3A2", "MUC5B", "KRT5"),
  Endothelial = c("PECAM1", "VWF", "CDH5", "EMCN", "ACKR1", "CCL21"),
  Mesenchymal = c("PDGFRB", "ACTA2", "COL1A1", "TCF21", "SERPINF1", "WNT5A", "COX4I2", "TAGLN"),
  Immune = c("PTPRC", "CD3D", "CD79A", "NKG7", "C1QA", "CD14", "CSF2RA", "MS4A1", "CPA3", "MZB1", "S100A8")
)
lineage_markers_genes <- unique(c(lineage_markers$Epithelial,lineage_markers$Endothelial,
                                  lineage_markers$Mesenchymal,lineage_markers$Immune))
DotPlot(merged_seurat_No2, features = lineage_markers_genes)+
  RotatedAxis()

#中性粒细胞(Neutrophils)
cat(top_markers_No2[top_markers_No2$cluster == "0",]$gene,sep = "\n")
#Epithelial Cells
cat(top_markers_No2[top_markers_No2$cluster == "1",]$gene,sep = "\n")
#Macrophage
cat(top_markers_No2[top_markers_No2$cluster == "2",]$gene,sep = "\n")
#Blood derived tam
cat(top_markers_No2[top_markers_No2$cluster == "3",]$gene,sep = "\n")
#Mesenchymal Cells
cat(top_markers_No2[top_markers_No2$cluster == "4",]$gene,sep = "\n")

new.cluster.ids <- c("0"="Neutrophils", 
                     "1"="Epithelial Cells", 
                     "2"="Macrophage", 
                     "3"="Blood derived tam", 
                     "4"="Mesenchymal Cells")

merged_seurat_No2 <- RenameIdents(merged_seurat_No2, new.cluster.ids) 
merged_seurat_No2$celltype <- merged_seurat_No2@active.ident

DimPlot(merged_seurat_No2, reduction = "umap", label = TRUE,cols = my36colors,group.by = "celltype")

##注释细胞气泡图
dot_plot_No2 <- DotPlot(merged_seurat_No2, 
                      features = top_markers_No2$gene,
                      cols = c("#fbebdf", "#a30e14"),
                      dot.scale = 4) + 
  RotatedAxis() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 6)) +
  ggtitle("Cell Type Marker Expression - Dot Plot")

# 保存图形
ggsave("celltype_markers_dotplot.pdf", dot_plot_No2, width = 9.9, height = 3.8)


# 定义您感兴趣的基因
genes_of_interest <- c("AES", "ACTG1", "TSC22D2", "HSPA1A", "PHLDA1", "SLC26A7", "NFAT5")
DotPlot(merged_seurat_No2, features = genes_of_interest)+
  RotatedAxis()
FeaturePlot(merged_seurat_No2,features = "PHLDA1",split.by = "orig.ident",pt.size = 0.1)


# 第一部分：基因表达图
# 生成一个存储多个ggplot对象的列表
p_merge_No2 <- list()

# 计算每个基因在所有细胞中的表达范围，用于统一颜色标度
expression_limits_No2 <- list()
for(gene in genes_of_interest) {
  expr_data_No2 <- FetchData(merged_seurat_No2, vars = gene)[,1]
  expression_limits_No2[[gene]] <- c(0, quantile(expr_data_No2[expr_data_No2 > 0], 0.95)) # 使用95%分位数作为上限
}

# 循环绘制每个基因在两个样本中的图
plot_index_No2 <- 1
for(gene in genes_of_interest) {
  # 打印当前基因
  print(gene)
  
  # 获取该基因的表达范围
  gene_limits_No2 <- expression_limits_No2[[gene]]
  
  # 绘制HOX样本的图
  hox_cells_No2 <- WhichCells(merged_seurat_No2, expression = orig.ident == "Neo1")
  p_merge_No2[[plot_index_No2]] <- FeaturePlot(merged_seurat_No2, 
                                               cells = hox_cells_No2,
                                               features = gene, 
                                               order = TRUE, 
                                               pt.size = 0.4) +
    scale_color_gradientn(
      colors = c("#fef9ef", "#f9dbb4", "#f1b782", "#e67c54", "#d13f2d", "#a2292d", "#7c2124"),
      limits = gene_limits_No2,
      name = "Expression"
    ) +
    theme_void() +
    labs(title = paste(gene, "(Neo1)")) +  # 添加基因名和样本标识作为标题
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      plot.background = element_blank(),
      legend.position = "none"  # 每个单独图表不显示图例
    )
  plot_index_No2 <- plot_index_No2 + 1
  
  # 绘制NOX样本的图
  nox_cells_No2 <- WhichCells(merged_seurat_No2, expression = orig.ident == "Neo3")
  p_merge_No2[[plot_index_No2]] <- FeaturePlot(merged_seurat_No2, 
                                               cells = nox_cells_No2,
                                               features = gene, 
                                               order = TRUE, 
                                               pt.size = 0.4) +
    scale_color_gradientn(
      colors = c("#fef9ef", "#f9dbb4", "#f1b782", "#e67c54", "#d13f2d", "#a2292d", "#7c2124"),
      limits = gene_limits_No2,
      name = "Expression"
    ) +
    theme_void() +
    labs(title = paste(gene, "(Neo3)")) +  # 添加基因名和样本标识作为标题
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      plot.background = element_blank(),
      legend.position = "none"  # 每个单独图表不显示图例
    )
  plot_index_No2 <- plot_index_No2 + 1
}

# 创建一个单独的图例（使用第一个基因的第一个样本）
legend_plot_No2 <- FeaturePlot(merged_seurat_No2, 
                               cells = WhichCells(merged_seurat_No2, expression = orig.ident == "Neo1"),
                               features = genes_of_interest[1], 
                               order = TRUE) +
  scale_color_gradientn(
    colors = c("#fef9ef", "#f9dbb4", "#f1b782", "#e67c54", "#d13f2d", "#a2292d", "#7c2124"),
    limits = expression_limits_No2[[genes_of_interest[1]]],
    name = "Expression"
  ) +
  theme(legend.position = "right")

# 提取图例
legend_No2 <- cowplot::get_legend(legend_plot_No2)

# 使用 patchwork 的 wrap_plots() 函数组合图表
# 每个基因占一行，每行有两个子图（HOX和NOX）
plot_grid_No2 <- wrap_plots(p_merge_No2, ncol = 2)  # 2列布局，每个基因一行

# 添加图例到网格
combined_plot_No2 <- cowplot::plot_grid(
  plot_grid_No2, 
  legend_No2,
  rel_widths = c(1, 0.15)  # 调整相对宽度
)

# 显示合并后的图
print(combined_plot_No2)
# 保存基因表达图
ggsave("target_genes_expression_by_sample_No2.pdf", combined_plot_No2, width = 7.9, height = 16)

# 第三部分：更详细的小提琴图，按细胞类型分面
# 提取表达数据
exp_data_No2 <- FetchData(merged_seurat_No2, vars = c(genes_of_interest, "celltype", "orig.ident"))

# 转换数据为长格式
library(tidyr)
exp_data_long_No2 <- pivot_longer(exp_data_No2, 
                                  cols = all_of(genes_of_interest),
                                  names_to = "Gene",
                                  values_to = "Expression")

# 创建分面小提琴图
faceted_violin_No2 <- ggplot(exp_data_long_No2, aes(x = orig.ident, y = Expression, fill = orig.ident)) +
  geom_violin(scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.1, fill = "white", alpha = 0.7, outlier.shape = NA) +
  facet_grid(Gene ~ celltype, scales = "free_y") +
  scale_fill_manual(values = c("Neo1" = "#E69F00", "Neo3" = "#56B4E9")) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 10),
    legend.position = "bottom"
  ) +
  labs(
    x = "Sample",
    y = "Expression Level",
    title = "Gene Expression by Cell Type and Sample - No2",
    fill = "Sample")
ggsave("target_genes_faceted_violin_No2.pdf", faceted_violin_No2, width = 6.8, height = 5.9)

# 第四部分：点图显示平均表达和表达百分比
# 计算每个细胞类型和样本组合的平均表达和表达百分比
avg_exp_No2 <- AverageExpression(merged_seurat_No2, 
                                 features = genes_of_interest,
                                 group.by = c("celltype", "orig.ident"),
                                 assays = "RNA")$RNA

# 重组数据以便绘图
plot_data_No2 <- data.frame()
for(gene in genes_of_interest) {
  for(group in colnames(avg_exp_No2)) {
    # 解析组名 (celltype_sample)
    parts_No2 <- strsplit(group, "_")[[1]]
    celltype <- paste(parts_No2[1:(length(parts_No2)-1)], collapse = "_")
    sample_No2 <- parts_No2[length(parts_No2)]
    
    # 计算表达百分比
    pct_exp_No2 <- mean(FetchData(merged_seurat_No2, 
                                  vars = gene,
                                  cells = WhichCells(merged_seurat_No2, 
                                                     expression = celltype == celltype & orig.ident == sample_No2))[,1] > 0)
    
    plot_data_No2 <- rbind(plot_data_No2, data.frame(
      Gene = gene,
      celltype = celltype,
      Sample = sample_No2,
      AvgExpression = avg_exp_No2[gene, group],
      PercentExpressed = pct_exp_No2 * 100
    ))
  }
}

# 创建点图
dot_plot_No2 <- ggplot(plot_data_No2, aes(x = Sample, y = Gene)) +
  geom_point(aes(size = PercentExpressed, color = AvgExpression)) +
  scale_color_gradientn(values = seq(0,1,0.1),colours = c("#4575b4","#abd9e9","#e0f3f8","#ffffbf","#fdae61","#d73027","#800026"))+
  scale_size(range = c(1, 8), name = "Percent\nExpressed") +
  facet_grid(. ~ celltype) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 8)) +
  labs(
    x = "Sample",
    y = "Gene",
    title = "Gene Expression by Cell Type and Sample - No2")

print(dot_plot_No2)
ggsave("target_genes_dotplot_by_sample_No2.pdf", dot_plot_No2, width = 6.7, height = 3.5)

# 第五部分：统计分析和结果输出
# 执行统计检验（Wilcoxon检验）
stat_results_No2 <- data.frame()

for(gene in genes_of_interest) {
  for(celltype in unique(merged_seurat_No2$celltype)) {
    # 提取该细胞类型中两个样本的表达数据
    Neo1_cells_No2 <- WhichCells(merged_seurat_No2, 
                                expression = celltype == celltype & orig.ident == "Neo1")
    Neo3_cells_No2 <- WhichCells(merged_seurat_No2, 
                                expression = celltype == celltype & orig.ident == "Neo3")
    
    if(length(Neo1_cells_No2) > 3 && length(Neo3_cells_No2) > 3) {  # 确保有足够的细胞进行检验
      Neo1_expr_No2 <- FetchData(merged_seurat_No2, vars = gene, cells = Neo1_cells_No2)[,1]
      Neo3_expr_No2 <- FetchData(merged_seurat_No2, vars = gene, cells = Neo3_cells_No2)[,1]
      
      # 执行Wilcoxon检验
      test_result_No2 <- wilcox.test(Neo1_expr_No2, Neo3_expr_No2)
      
      stat_results_No2 <- rbind(stat_results_No2, data.frame(
        Gene = gene,
        celltype = celltype,
        Neo1_Mean = mean(Neo1_expr_No2),
        Neo3_Mean = mean(Neo3_expr_No2),
        Neo1_Median = median(Neo1_expr_No2),
        Neo3_Median = median(Neo3_expr_No2),
        P_Value = test_result_No2$p.value,
        Significant = test_result_No2$p.value < 0.05
      ))
    }
  }
}

# 添加FDR校正
stat_results_No2$FDR <- p.adjust(stat_results_No2$P_Value, method = "fdr")
stat_results_No2$Significant_FDR <- stat_results_No2$FDR < 0.05

# 输出统计结果
print(stat_results_No2)
write.csv(stat_results_No2, "gene_expression_statistics_No2.csv", row.names = FALSE)

# 创建显著性标记的热图
library(ggplot2)
sig_heatmap_No2 <- ggplot(stat_results_No2, aes(x = celltype, y = Gene, fill = -log10(P_Value))) +
  geom_tile() +
  scale_fill_gradientn(colors = c("white", "yellow", "red"),
                       name = "-log10(P-value)") +
  geom_text(aes(label = ifelse(Significant, "*", "")), 
            size = 6, color = "black") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)) +
  labs(
    x = "Cell Type",
    y = "Gene",
    title = "Significance of Expression Differences (HOX vs NOX) - No2")
print(sig_heatmap_No2)
ggsave("significance_heatmap_No2.pdf", sig_heatmap_No2, width = 4, height = 2.5)

saveRDS(merged_seurat_No2,file = "GSE137950_merged_seurat.rds")

save.image(file = "/data/nas1/xiangxiaoyu_OD/1.BPD+MSCs/1.scRNA/09_ST/scRNA/scRNA.workplace.RData")
