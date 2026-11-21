###############################################################################
# Single-cell analysis of the hyperoxic rat lung dataset
#
# Source      : r.03_Cell_prop.r (lines 8722-9378)
# Input       : GSE242310
# Output      : Fig. 6E-H; Supplementary Fig. 9C-E
#
# Code below is unchanged from the original analysis script.
###############################################################################
# -----------------------------------------------------------------------------
# PATHS: this file keeps the absolute paths of the analysis server on which the
# study was run (/data/nas1/xiangxiaoyu_OD/1.BPD+MSCs/1.scRNA/... and
# /home/xxyloveme/study/BPD+MSCs/...). Set PROJ_DIR in config/paths.R and
# replace those prefixes with your own project directory before running.
# -----------------------------------------------------------------------------


#############单细胞数据 GSE242310_RAW##########
setwd("/data/nas1/xiangxiaoyu_OD/1.BPD+MSCs/1.scRNA/09_ST/scRNA/GSE242310_RAW/")

sapply(c("NOX", "HOX"), dir.create, showWarnings = FALSE)
# 移动文件到对应目录 (请根据实际解压后的文件名调整)
file.rename("GSM7757968_NOX_barcodes.tsv", "NOX/barcodes.tsv")
file.rename("GSM7757968_NOX_features.tsv", "NOX/features.tsv")
file.rename("GSM7757968_NOX_matrix.mtx", "NOX/matrix.mtx")

file.rename("GSM7757969_HOX_barcodes.tsv", "HOX/barcodes.tsv")
file.rename("GSM7757969_HOX_features.tsv", "HOX/features.tsv")
file.rename("GSM7757969_HOX_matrix.mtx", "HOX/matrix.mtx")

library(Seurat)
# 设置样本名称
samples <- c("NOX", "HOX")

# 创建一个列表来存储每个样本的数据
seurat_list <- lapply(samples, function(sample) {
  # 读取10X数据
  data_dir <- file.path(".", sample) # 指向样本子目录
  counts <- Read10X(data.dir = data_dir)
  # 创建Seurat对象，并添加样本标识
  seurat_obj <- CreateSeuratObject(counts = counts, project = sample)
  return(seurat_obj)
})

# 合并所有样本的Seurat对象
merged_seurat_No2 <- merge(seurat_list[[1]], y = seurat_list[-1], 
                       add.cell.ids = samples, project = "Combined_Project")

merged_seurat_No2[["percent.mt"]] <- PercentageFeatureSet(merged_seurat_No2, pattern = "^mt-")
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
merged_seurat <- subset(merged_seurat,
                        subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 20)
# 数据标准化
merged_seurat <- NormalizeData(merged_seurat)
# 寻找高变基因
merged_seurat <- FindVariableFeatures(merged_seurat, selection.method = "vst", 
                                           nfeatures = 2000)
# 数据缩放
merged_seurat <- ScaleData(merged_seurat)
# PCA分析
merged_seurat <- RunPCA(merged_seurat, 
                        features = VariableFeatures(object = merged_seurat))
# 细胞聚类
merged_seurat <- FindNeighbors(merged_seurat, dims = 1:10)
merged_seurat <- FindClusters(merged_seurat, resolution = 0.4)

#resolution = 0.4
# UMAP降维
merged_seurat <- RunUMAP(merged_seurat, dims = 1:10)
# merged_seurat <- RunTSNE(merged_seurat, dims = 1:10)
# Idents(merged_seurat) <- "RNA_snn_res.0.4"
# 可视化聚类结果
DimPlot(merged_seurat, reduction = "umap", label = TRUE)


##细胞类型注释
merged_seurat <- JoinLayers(merged_seurat)
cluster_markers <- FindAllMarkers(merged_seurat, only.pos = TRUE, min.pct = 0.25, 
                                  logfc.threshold = 0.25)
top_markers <- cluster_markers %>%
  group_by(cluster) %>% #按照clusetr分类
  top_n(n = 5, wt = avg_log2FC) %>% #按照log2fc取前十
  dplyr::arrange(desc(avg_log2FC))%>%
  dplyr::arrange(cluster)

#Neutrophil
cat(top_markers[top_markers$cluster == "0",]$gene,sep = "\n")
#Mature alpha-beta T cell
cat(top_markers[top_markers$cluster == "1",]$gene,sep = "\n")
#ARG1+, TREM2+ regulatory myeloid (mreg) cell
cat(top_markers[top_markers$cluster == "2",]$gene,sep = "\n")
#Eosinophil progenitor cell
cat(top_markers[top_markers$cluster == "3",]$gene,sep = "\n")
#T1 B cell
cat(top_markers[top_markers$cluster == "4",]$gene,sep = "\n")
#Macrophage
cat(top_markers[top_markers$cluster == "5",]$gene,sep = "\n")
#Epithelial cell
cat(top_markers[top_markers$cluster == "6",]$gene,sep = "\n")
#Natural helper lymphocyte
cat(top_markers[top_markers$cluster == "7",]$gene,sep = "\n")
#Erythrocyte
cat(top_markers[top_markers$cluster == "8",]$gene,sep = "\n")
#Cycling cell
cat(top_markers[top_markers$cluster == "9",]$gene,sep = "\n")
#Erythroid progenitor cell
cat(top_markers[top_markers$cluster == "10",]$gene,sep = "\n")
#Pro-NK cell
cat(top_markers[top_markers$cluster == "11",]$gene,sep = "\n")
#Reticular cell/Marrow fibroblast
cat(top_markers[top_markers$cluster == "12",]$gene,sep = "\n")
#Lipofibroblast
cat(top_markers[top_markers$cluster == "13",]$gene,sep = "\n")
#Immature gamma-delta T cell
cat(top_markers[top_markers$cluster == "14",]$gene,sep = "\n")
#CD11b-positive dendritic cell
cat(top_markers[top_markers$cluster == "15",]$gene,sep = "\n")

# 小鼠谱系标记基因
lineage_markers <- list(
  Epithelial = c("Krt19", "Sftpc", "Scgb1a1", "Scgb3a2", "Muc5b", "Krt5", "Ager", "Rtkn2", "Lama3", "Col4a3", "Pdpn"),
  Endothelial = c("Pecam1", "Vwf", "Cdh5", "Emcn", "Ackr1", "Ccl21", "Ednrb", "Sostdc1", "Hpgd", "Tspan12", "S100a3", "Apln"),
  Mesenchymal = c("Pdgfrb", "Acta2", "Col1a1", "Tcf21", "Serpinf1", "Wnt5a", "Cox4i2", "Tagln", "Pdgfra", "Cxcl14", "Aspn", "Ltbp2"),
  Immune = c("Ptprc", "Cd3d", "Cd79a", "Nkg7", "C1qa", "Cd14", "Csf2ra", "Ms4a1", "Cpa3", "Mzb1", "S100a8", "Cd8a", "Themis", "Cd6", "Il7r")
)

# 提取所有唯一基因
lineage_markers_genes <- unique(c(lineage_markers$Epithelial, lineage_markers$Endothelial,
                                  lineage_markers$Mesenchymal, lineage_markers$Immune))


# 定义每个细胞类型的标记基因
celltype_markers <- list(
  # 0: Neutrophil
  "Neutrophil" = c("S100a8", "S100a9", "Csf3r", "Fcgr3", "Cxcr2", "Mmp8", "Mmp9"),
  
  # 1: Mature alpha-beta T cell
  "Mature_ab_T_cell" = c("Cd3d", "Cd3e", "Cd4", "Cd8a", "Il7r", "Trac", "Trbc1", "Lef1"),
  
  # 2: ARG1+, TREM2+ regulatory myeloid (mreg) cell
  "Regulatory_myeloid" = c("Arg1", "Trem2", "Mrc1", "Cd163", "Il10", "Tgfb1", "Vegfa"),
  
  # 3: Eosinophil progenitor cell
  "Eosinophil_progenitor" = c("Il5ra", "Siglecf", "Ccr3", "Epx", "Prg2", "Prg3", "Cebpe"),
  
  # 4: T1 B cell
  "T1_B_cell" = c("Cd19", "Cd79a", "Ms4a1", "Ighm", "Igd", "Cd24a", "Cd93"),
  
  # 5: Macrophage
  "Macrophage" = c("Adgre1", "Cd68", "C1qa", "C1qb", "C1qc", "Mertk", "Fcgr1", "Itgam"),
  
  # 6: Epithelial cell
  "Epithelial" = c("Epcam", "Krt8", "Krt18", "Cdh1", "Sftpc", "Scgb1a1", "Scgb3a2"),
  
  # 7: Natural helper lymphocyte
  "Natural_helper_lymphocyte" = c("Il2rb", "Il18r1", "Klrb1c", "Thy1", "Il17rb", "Rora"),
  
  # 8: Erythrocyte
  "Erythrocyte" = c("Hba-a1", "Hba-a2", "Hbb-bs", "Hbb-bt", "Alas2", "Gypa", "Slc4a1"),
  
  # 9: Cycling cell
  "Cycling_cell" = c("Mki67", "Pcna", "Top2a", "Birc5", "Cenpa", "Cenpf", "Ube2c"),
  
  # 10: Erythroid progenitor cell
  "Erythroid_progenitor" = c("Gata1", "Klf1", "Epor", "Tal1", "Lmo2", "Gfi1b", "Cd34"),
  
  # 11: Pro-NK cell
  "Pro_NK_cell" = c("Ncr1", "Klrb1c", "Klrk1", "Klrd1", "Eomes", "Tbx21", "Il2rb", "Cd7"),
  
  # 12: Reticular cell/Marrow fibroblast
  "Reticular_fibroblast" = c("Pdgfrb", "Cspg4", "Lepr", "Nt5e", "Cxcl12", "Kitl", "Vcam1"),
  
  # 13: Lipofibroblast
  "Lipofibroblast" = c("Plin2", "Adipoq", "Lpl", "Fabp4", "Pparg", "Cidec", "Dgat2"),
  
  # 14: Immature gamma-delta T cell
  "Immature_gd_T_cell" = c("Trdc", "Trgc1", "Trgc2", "Cd3d", "Cd3e", "Il2rg", "Rag1", "Rag2"),
  
  # 15: CD11b-positive dendritic cell
  "CD11b_DC" = c("Itgam", "Itgax", "Cd11c", "H2-Ab1", "H2-Aa", "Cd74", "Ccr7", "Flt3")
)

# 提取所有唯一基因
all_markers <- unique(unlist(celltype_markers))

# 检查哪些基因在数据中可用
available_markers <- all_markers[all_markers %in% rownames(merged_seurat)]

cat("可用的标记基因:", length(available_markers), "/", length(all_markers), "\n")

# 创建气泡图
library(ggplot2)

# 方法1: 使用DotPlot函数
dot_plot <- DotPlot(merged_seurat, 
                    features = unique(top_markers$gene),
                    cols = c("blue", "red"),
                    dot.scale = 4) + 
  RotatedAxis() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 6)) +
  ggtitle("Cell Type Marker Expression - Dot Plot")

print(dot_plot)

dot_plot_2 <- DotPlot(merged_seurat, 
                    features = available_markers,
                    cols = c("#fbebdf", "#a30e14"),
                    dot.scale = 4) + 
  RotatedAxis() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 6)) +
  ggtitle("Cell Type Marker Expression - Dot Plot")

# 保存图形
ggsave("celltype_markers_dotplot.pdf", dot_plot, width = 14.9, height = 4.8)
ggsave("celltype_markers_dotplot_2.pdf", dot_plot_2, width = 14.9, height = 4.8)

#细胞类型注释
new.cluster.ids <- c("0"="Neutrophil", 
                     "1"="Mature_ab_T_cell", 
                     "2"="Regulatory_myeloid", 
                     "3"="Eosinophil_progenitor", 
                     "4"="T1_B_cell", 
                     "5"="Macrophage", 
                     "6"="Epithelial", 
                     "7"="Natural_helper_lymphocyte", 
                     "8"="Erythrocyte", 
                     "9"="Cycling_cell", 
                     "10"="Erythroid_progenitor", 
                     "11"="Pro-NK cell", 
                     "12"="Mesenchyme", 
                     "13"="Lipofibroblast", 
                     "14"="Immature_gd_T_cell", 
                     "15"="CD11b_DC")

new.cluster.ids <- c("0"="Imm", 
                     "1"="Imm", 
                     "2"="Imm", 
                     "3"="Imm", 
                     "4"="Imm", 
                     "5"="Imm", 
                     "6"="Epithelial", 
                     "7"="Imm", 
                     "8"="Imm", 
                     "9"="Cycling_cell", 
                     "10"="Imm", 
                     "11"="Imm", 
                     "12"="Mesenchyme", 
                     "13"="Lipofibroblast", 
                     "14"="Imm", 
                     "15"="Imm")

merged_seurat <- RenameIdents(merged_seurat, new.cluster.ids) 
merged_seurat$celltype <- merged_seurat@active.ident
merged_seurat$celltype_big <- merged_seurat@active.ident
  
my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
               '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
               '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
               '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
               '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
               '#968175')
DimPlot(merged_seurat, reduction = "umap", label = TRUE,cols = my36colors,group.by = "celltype_big")
DimPlot(merged_seurat, reduction = "umap", label = TRUE,cols = my36colors)


# 定义您感兴趣的基因
genes_of_interest <- c("AES", "ACTG1", "TSC22D2", "HSPA1A", "PHLDA1", "SLC26A7", "MT1A", "NFAT5")
library(homologene)
#使用homologene函进行转换
#@genelist是要转换的基因列表
#@inTax是输入的基因列表所属的物种号，10090是小鼠
#@outTax是要转换成的物种号，9606是人
genes_of_interest_mouse <-homologene(genes_of_interest, inTax = 9606, outTax = 10090)

# 检查哪些基因在数据中
available_genes <- genes_of_interest_mouse$`10090`[genes_of_interest_mouse$`10090` %in% rownames(merged_seurat)]

DotPlot(merged_seurat, features = genes_of_interest_mouse$`10090`)+
  RotatedAxis()

FeaturePlot(merged_seurat,features = "Nfat5",split.by = "orig.ident",pt.size = 0.1)

FeaturePlot(merged_seurat,features = "Epcam",pt.size = 0.1)
FeaturePlot(merged_seurat,features = "Ptprc",pt.size = 0.1)
FeaturePlot(merged_seurat,features = "Pdgfrb",pt.size = 0.1)
FeaturePlot(merged_seurat,features = "Lpo",pt.size = 0.1)
FeaturePlot(merged_seurat,features = "Mki67",pt.size = 0.1)

###美化我们自己的图形

library(Seurat)
library(patchwork)
library(ggplot2)
library(cowplot)

# 定义要绘制的基因列表
target_genes <- c("Tsc22d2", "Phlda1", "Nfat5")

# 检查这些基因是否在数据中可用
available_genes <- target_genes[target_genes %in% rownames(merged_seurat)]
missing_genes <- setdiff(target_genes, rownames(merged_seurat))
# 第一部分：每个基因分别在HOX和NOX两组中的基因表达分布图
# 生成一个存储多个ggplot对象的列表
p_merge <- list()

# 计算每个基因在所有细胞中的表达范围，用于统一颜色标度
expression_limits <- list()
for(gene in available_genes) {
  expr_data <- FetchData(merged_seurat, vars = gene)[,1]
  expression_limits[[gene]] <- c(0, quantile(expr_data[expr_data > 0], 0.95)) # 使用95%分位数作为上限
}

# 循环绘制每个基因在两个样本中的图
plot_index <- 1
for(gene in available_genes) {
  # 打印当前基因
  print(gene)
  
  # 获取该基因的表达范围
  gene_limits <- expression_limits[[gene]]
  
  # 绘制HOX样本的图
  hox_cells <- WhichCells(merged_seurat, expression = orig.ident == "HOX")
  p_merge[[plot_index]] <- FeaturePlot(merged_seurat, 
                                       cells = hox_cells,
                                       features = gene, 
                                       order = TRUE, 
                                       pt.size = 0.4) +
    scale_color_gradientn(
      colors = c("#fef9ef", "#f9dbb4", "#f1b782", "#e67c54", "#d13f2d", "#a2292d", "#7c2124"),
      limits = gene_limits,
      name = "Expression"
    ) +
    theme_void() +
    labs(title = paste(gene, "(HOX)")) +  # 添加基因名和样本标识作为标题
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      plot.background = element_blank(),
      legend.position = "none"  # 每个单独图表不显示图例
    )
  plot_index <- plot_index + 1
  
  # 绘制NOX样本的图
  nox_cells <- WhichCells(merged_seurat, expression = orig.ident == "NOX")
  p_merge[[plot_index]] <- FeaturePlot(merged_seurat, 
                                       cells = nox_cells,
                                       features = gene, 
                                       order = TRUE, 
                                       pt.size = 0.4) +
    scale_color_gradientn(
      colors = c("#fef9ef", "#f9dbb4", "#f1b782", "#e67c54", "#d13f2d", "#a2292d", "#7c2124"),
      limits = gene_limits,
      name = "Expression"
    ) +
    theme_void() +
    labs(title = paste(gene, "(NOX)")) +  # 添加基因名和样本标识作为标题
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      plot.background = element_blank(),
      legend.position = "none"  # 每个单独图表不显示图例
    )
  plot_index <- plot_index + 1
}

# 创建一个单独的图例（使用第一个基因的第一个样本）
legend_plot <- FeaturePlot(merged_seurat, 
                           cells = WhichCells(merged_seurat, expression = orig.ident == "HOX"),
                           features = available_genes[1], 
                           order = TRUE) +
  scale_color_gradientn(
    colors = c("#fef9ef", "#f9dbb4", "#f1b782", "#e67c54", "#d13f2d", "#a2292d", "#7c2124"),
    limits = expression_limits[[available_genes[1]]],
    name = "Expression"
  ) +
  theme(legend.position = "right")

# 提取图例
legend <- cowplot::get_legend(legend_plot)

# 使用 patchwork 的 wrap_plots() 函数组合图表
# 每个基因占一行，每行有两个子图（HOX和NOX）
plot_grid <- wrap_plots(p_merge, ncol = 2)  # 2列布局，每个基因一行

# 添加图例到网格
combined_plot <- cowplot::plot_grid(
  plot_grid, 
  legend,
  rel_widths = c(1, 0.15)  # 调整相对宽度
)

# 显示合并后的图
print(combined_plot)

# 保存基因表达图
ggsave("target_genes_expression_by_sample.pdf", combined_plot, width = 7.9, height = 16)

# 第二部分：为每个基因创建并排的HOX和NOX表达图（更紧凑的布局）
# 生成一个存储多个ggplot对象的列表（每个基因一个包含两个样本的图）
compact_plots <- list()

for(i in 1:length(available_genes)) {
  gene <- available_genes[i]
  
  # 获取该基因的表达范围
  gene_limits <- expression_limits[[gene]]
  
  # 绘制HOX样本的图
  hox_cells <- WhichCells(merged_seurat, expression = orig.ident == "HOX")
  p_hox <- FeaturePlot(merged_seurat, 
                       cells = hox_cells,
                       features = gene, 
                       order = TRUE, 
                       pt.size = 0.4) +
    scale_color_gradientn(
      colors = c("#fef9ef", "#f9dbb4", "#f1b782", "#e67c54", "#d13f2d", "#a2292d", "#7c2124"),
      limits = gene_limits,
      name = "Expression"
    ) +
    theme_void() +
    labs(title = "HOX") +
    theme(
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
      legend.position = "none"
    )
  
  # 绘制NOX样本的图
  nox_cells <- WhichCells(merged_seurat, expression = orig.ident == "NOX")
  p_nox <- FeaturePlot(merged_seurat, 
                       cells = nox_cells,
                       features = gene, 
                       order = TRUE, 
                       pt.size = 0.4) +
    scale_color_gradientn(
      colors = c("#fef9ef", "#f9dbb4", "#f1b782", "#e67c54", "#d13f2d", "#a2292d", "#7c2124"),
      limits = gene_limits,
      name = "Expression"
    ) +
    theme_void() +
    labs(title = "NOX") +
    theme(
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
      legend.position = "none"
    )
  
  # 合并两个样本的图
  compact_plots[[i]] <- (p_hox | p_nox) + 
    plot_annotation(title = gene, theme = theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold")))
}

# 组合所有基因的图
compact_combined <- wrap_plots(compact_plots, ncol = 1)  # 每个基因一行

# 添加图例
compact_legend_plot <- FeaturePlot(merged_seurat, 
                                   cells = WhichCells(merged_seurat, expression = orig.ident == "HOX"),
                                   features = available_genes[1], 
                                   order = TRUE) +
  scale_color_gradientn(
    colors = c("#fef9ef", "#f9dbb4", "#f1b782", "#e67c54", "#d13f2d", "#a2292d", "#7c2124"),
    limits = expression_limits[[available_genes[1]]],
    name = "Expression"
  ) +
  theme(legend.position = "right")

compact_legend <- cowplot::get_legend(compact_legend_plot)

# 最终组合图
final_compact_plot <- cowplot::plot_grid(
  compact_combined, 
  compact_legend,
  rel_widths = c(1, 0.1)
)

# 显示紧凑布局的图
print(final_compact_plot)

# 保存紧凑布局的图
ggsave("target_genes_expression_compact.pdf", final_compact_plot, width = 8.5, height = 16)


# 第三部分：更详细的小提琴图，按细胞类型分面
# 提取表达数据
exp_data <- FetchData(merged_seurat, vars = c(available_genes, "celltype_big", "orig.ident"))

# 转换数据为长格式
library(tidyr)
exp_data_long <- pivot_longer(exp_data, 
                              cols = all_of(available_genes),
                              names_to = "Gene",
                              values_to = "Expression")

# 创建分面小提琴图
faceted_violin <- ggplot(exp_data_long, aes(x = orig.ident, y = Expression, fill = orig.ident)) +
  geom_violin(scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.1, fill = "white", alpha = 0.7, outlier.shape = NA) +
  facet_grid(Gene ~ celltype, scales = "free_y") +
  scale_fill_manual(values = c("HOX" = "#E69F00", "NOX" = "#56B4E9")) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 10),
    legend.position = "bottom"
  ) +
  labs(
    x = "Sample",
    y = "Expression Level",
    title = "Gene Expression by Cell Type and Sample",
    fill = "Sample")
ggsave("target_genes_faceted_violin.pdf", faceted_violin, width = 17.8, height = 5.9)

faceted_violin_big <- ggplot(exp_data_long, aes(x = orig.ident, y = Expression, fill = orig.ident)) +
  geom_violin(scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.1, fill = "white", alpha = 0.7, outlier.shape = NA) +
  facet_grid(Gene ~ celltype_big, scales = "free_y") +
  scale_fill_manual(values = c("HOX" = "#E69F00", "NOX" = "#56B4E9")) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 10),
    legend.position = "bottom"
  ) +
  labs(
    x = "Sample",
    y = "Expression Level",
    title = "Gene Expression by Cell Type and Sample",
    fill = "Sample")
print(faceted_violin_big)
ggsave("target_genes_faceted_violin_big.pdf", faceted_violin_big, width = 7.2, height = 5.9)

# 第四部分：点图显示平均表达和表达百分比
# 计算每个细胞类型和样本组合的平均表达和表达百分比
avg_exp <- AverageExpression(merged_seurat, 
                             features = available_genes,
                             group.by = c("celltype_big", "orig.ident"),
                             assays = "RNA")$RNA

# 重组数据以便绘图
plot_data <- data.frame()
for(gene in available_genes) {
  for(group in colnames(avg_exp)) {
    # 解析组名 (celltype_big_sample)
    parts <- strsplit(group, "_")[[1]]
    celltype_big <- paste(parts[1:(length(parts)-1)], collapse = "_")
    sample <- parts[length(parts)]
    
    # 计算表达百分比
    pct_exp <- mean(FetchData(merged_seurat, 
                              vars = gene,
                              cells = WhichCells(merged_seurat, 
                                                 expression = celltype_big == celltype_big & orig.ident == sample))[,1] > 0)
    
    plot_data <- rbind(plot_data, data.frame(
      Gene = gene,
      celltype_big = celltype_big,
      Sample = sample,
      AvgExpression = avg_exp[gene, group],
      PercentExpressed = pct_exp * 100
    ))
  }
}

# 创建点图
dot_plot <- ggplot(plot_data, aes(x = Sample, y = Gene)) +
  geom_point(aes(size = PercentExpressed, color = AvgExpression)) +
  scale_color_gradientn(values = seq(0,1,0.1),colours = c("#4575b4","#abd9e9","#e0f3f8","#ffffbf","#fdae61","#d73027","#800026"))+
  scale_size(range = c(1, 8), name = "Percent\nExpressed") +
  facet_grid(. ~ celltype_big) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 8)) +
  labs(
    x = "Sample",
    y = "Gene",
    title = "Gene Expression by Cell Type and Sample")

print(dot_plot)
ggsave("target_genes_dotplot_by_sample.pdf", dot_plot, width = 6.7, height = 3.5)

# 第五部分：统计分析和结果输出
# 执行统计检验（Wilcoxon检验）
stat_results <- data.frame()

for(gene in available_genes) {
  for(celltype_big in unique(merged_seurat$celltype_big)) {
    # 提取该细胞类型中两个样本的表达数据
    hox_cells <- WhichCells(merged_seurat, 
                            expression = celltype_big == celltype_big & orig.ident == "HOX")
    nox_cells <- WhichCells(merged_seurat, 
                            expression = celltype_big == celltype_big & orig.ident == "NOX")
    
    if(length(hox_cells) > 3 && length(nox_cells) > 3) {  # 确保有足够的细胞进行检验
      hox_expr <- FetchData(merged_seurat, vars = gene, cells = hox_cells)[,1]
      nox_expr <- FetchData(merged_seurat, vars = gene, cells = nox_cells)[,1]
      
      # 执行Wilcoxon检验
      test_result <- wilcox.test(hox_expr, nox_expr)
      
      stat_results <- rbind(stat_results, data.frame(
        Gene = gene,
        celltype_big = celltype_big,
        HOX_Mean = mean(hox_expr),
        NOX_Mean = mean(nox_expr),
        HOX_Median = median(hox_expr),
        NOX_Median = median(nox_expr),
        P_Value = test_result$p.value,
        Significant = test_result$p.value < 0.05
      ))
    }
  }
}

# 添加FDR校正
stat_results$FDR <- p.adjust(stat_results$P_Value, method = "fdr")
stat_results$Significant_FDR <- stat_results$FDR < 0.05

# 输出统计结果
print(stat_results)
write.csv(stat_results, "gene_expression_statistics.csv", row.names = FALSE)

# 创建显著性标记的热图
library(ggplot2)
sig_heatmap <- ggplot(stat_results, aes(x = celltype_big, y = Gene, fill = -log10(P_Value))) +
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
    title = "Significance of Expression Differences (HOX vs NOX)")
print(sig_heatmap)
ggsave("significance_heatmap.pdf", sig_heatmap, width = 4, height = 2.5)
