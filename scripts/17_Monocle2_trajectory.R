###############################################################################
# Pseudotime trajectory (Monocle 2, DDRTree) and trajectory-associated gene modules
#
# Source      : r.03_Cell_prop.r (lines 4190-5236)
# Input       : meso_adventitial_FB
# Output      : Fig. 4G-J; Supplementary Tables 12-14
# Note        : The Methods of the submitted manuscript state Monocle 3; the code uses monocle (v2) with DDRTree and orderCells(root_state = 2). See revision_scripts/B4_trajectory_LODO.R.
#
# Code below is unchanged from the original analysis script.
###############################################################################
# -----------------------------------------------------------------------------
# PATHS: this file keeps the absolute paths of the analysis server on which the
# study was run (/data/nas1/xiangxiaoyu_OD/1.BPD+MSCs/1.scRNA/... and
# /home/xxyloveme/study/BPD+MSCs/...). Set PROJ_DIR in config/paths.R and
# replace those prefixes with your own project directory before running.
# -----------------------------------------------------------------------------


####################3.5. 拟时序分析#############
setwd("/data/nas1/xiangxiaoyu_OD/1.BPD+MSCs/1.scRNA/04_Adventitial_FB/08_Monocle3/")
meso_adventitial_FB <- readRDS(file = "/data/nas1/xiangxiaoyu_OD/1.BPD+MSCs/1.scRNA/04_Adventitial_FB/meso_adventitial_FB.rds")
# 综合处理方法
meso_adventitial_FB$cell_subtype <- ifelse(
  # 情况1：有聚类分配且是Adventitial FB
  !is.na(meso_adventitial_FB$SCT_snn_res.0.5) & meso_adventitial_FB$celltype == "Adventitial FB",
  paste0("Adventitial_FB_Cluster", meso_adventitial_FB$SCT_snn_res.0.5),
  
  # 情况2：有聚类分配但不是Adventitial FB（保持原celltype）
  ifelse(
    !is.na(meso_adventitial_FB$SCT_snn_res.0.5) & meso_adventitial_FB$celltype != "Adventitial FB",
    meso_adventitial_FB$celltype,
    
    # 情况3：没有聚类分配（保持原celltype）
    meso_adventitial_FB$celltype
  )
)

Idents(cellchat_seurat) <- "condition"

library(Seurat)
library(ggplot2)
library(dplyr)
library(viridis)
library(monocle3)

# 查看当前库路径
.libPaths()
# 临时只使用conda环境的库路径
.libPaths("/data/home/xiangxiaoyu_OD/.conda/envs/sc_analysis/lib/R/library")
# 移除有问题的包
remove.packages("Biobase", lib = .libPaths()[2:3])
# 重新安装
BiocManager::install("Biobase")
library(Biobase)


# 将Seurat对象转换为Monocle3对象
# 安装 devtools 包（如果尚未安装）
install.packages("devtools")
# 使用 devtools 从 GitHub 安装 SeuratWrappers
library(devtools)
install_github("satijalab/seurat-wrappers")

library(SeuratWrappers)
cds <- as.cell_data_set(meso_adventitial_FB, assay = "SCT")

##细胞轨迹分析
library(ggplot2)
library(magrittr)
library(monocle)
library(Seurat)
library(tidyverse)
library(monocle)
library(Seurat)
library(tidyverse)

DefaultAssay(meso_adventitial_FB) <- "RNA"
EC_tj <- meso_adventitial_FB
EC_matrix = GetAssayData(EC_tj, slot = "count", assay = "RNA")
feature_ann = data.frame(gene_id=rownames(EC_matrix),gene_short_name=rownames(EC_matrix))
rownames(feature_ann) = rownames(EC_matrix)
EC_fd = new("AnnotatedDataFrame", data = feature_ann)
sample_ann = EC_tj@meta.data
EC_pd = new("AnnotatedDataFrame", data =sample_ann)
EC.cds = newCellDataSet(EC_matrix, phenoData =EC_pd, featureData=EC_fd, expressionFamily=negbinomial.size())
rm(EC_tj)
EC.cds = estimateSizeFactors(EC.cds)
EC.cds = estimateDispersions(EC.cds)
WGCNA::collectGarbage()
# Filtering low-quality cells 过滤低质量的细胞
# 这一操作会在fData(cds)中添加一列num_cells_expressed
mycds <- detectGenes(EC.cds, min_expr = 0.1)

print(head(fData(mycds))) 

expressed_genes <- row.names(
  subset(fData(mycds), 
         num_cells_expressed >= 10)
) # #过滤掉在小于10个细胞中表达的基因
print(head(pData(mycds)))

diff_test_res <- differentialGeneTest(
  mycds[expressed_genes,],
  fullModelFormulaStr = "~cell_subtype")  #~后面是表示对谁做差异分析的变量，理论上可以为p_data的任意列名
# ##差异表达基因作为轨迹构建的基因,差异基因的选择标准是qval<0.01,decreasing=F表示按数值增加排序
deg <- subset(diff_test_res, qval < 0.05) 
deg <- deg[order(deg$qval,decreasing=F),]
# head(deg)
# 
# ##差异基因的结果文件保存
write.table(deg,file="train.key.cells.DEG.xls",col.names=T,row.names=F,sep="\t",quote=F)

## 轨迹构建基因可视化
ordergene <- rownames(deg)
mycds1 <- setOrderingFilter(mycds, ordergene)
#这一步是很重要的，在我们得到想要的基因列表后，我们需要使用setOrderingFilter将它嵌入cds对象，后续的一系列操作都要依赖于这个list。
#setOrderingFilter之后，这些基因被储存在cds@featureData@data[["use_for_ordering"]]，可以通过table(cds@featureData@data[["use_for_ordering"]])查看

pdf("ordergenes.cds_cells.pdf")
plot_ordering_genes(mycds1)
dev.off()
png("ordergenes.cds_cells.png")
plot_ordering_genes(mycds1)
dev.off()


mycds2 <- reduceDimension(mycds1, max_components = 2,
                          method = 'DDRTree')
save(mycds2, file = "mycds_order_cells_cds.RData")


# mycds3 <- subset(mycds2, downsample = 15)
mycds3 <- orderCells(mycds2) 
mycds3 = orderCells(mycds3,root_state = 2)
save(mycds3, file = "mycds_order_cells_cds_2_final.RData")

##########绘制相应的图片
setwd("/data/nas1/xiangxiaoyu_OD/1.BPD+MSCs/1.scRNA/04_Adventitial_FB/08_Monocle3/")

library(monocle)
library(ggplot2)
library(dplyr)
library(patchwork)
library(ggridges)
library(viridis)

# 定义颜色方案
celltype_color <- c("#a4cdad","#f4c5cd","#c88688","#fde2e7","#fbf4be","#789eb1",
                    "#e0f2fe","#bbb4d5","#acceea","#f4c78c","#c9bad0","#a07a95",
                    "#d55431","#844b44","#038061","#376352")

condition_colors <- c("Term infant" = "grey", "BPD" = "#c69189", "BPD+PH" = "#e4cd87")

subtype_colors <- c("Adventitial_FB_Cluster0" = "#a4cdad", 
                    "Adventitial_FB_Cluster1" = "#f4c5cd", 
                    "Adventitial_FB_Cluster2" = "#c88688",
                    "Adventitial_FB_Cluster3" = "#fde2e7",
                    "Adventitial_FB_Cluster4" = "#789eb1",
                    "Adventitial_FB_Cluster5" = "#bbb4d5")

state_colors <- c("#fde2e7","#fbf4be","#789eb1","#e0f2fe","#bbb4d5","#acceea",
                  "#f4c78c","#c9bad0","#a07a95","#d55431","#844b44","#038061","#376352")

## 1. 拟时序轨迹图
p_trajectory <- plot_cell_trajectory(mycds3, color_by = "Pseudotime", cell_size = 1, theta = 180,
                                     size = 1, show_backbone = TRUE, show_branch_points = TRUE) +
  theme(text = element_text(size = 14))

# 保存拟时序轨迹图
ggsave("pseudotime_trajectory.pdf", p_trajectory, width = 10, height = 8, device = cairo_pdf)
ggsave("pseudotime_trajectory.png", p_trajectory, width = 10, height = 8, dpi = 300)

## 2. 细胞亚型轨迹图
p_subtype <- plot_cell_trajectory(mycds3, color_by = "cell_subtype", cell_size = 1, theta = 180,
                                  size = 1, show_backbone = TRUE, show_branch_points = TRUE) +
  scale_color_manual(values = subtype_colors) +
  theme(text = element_text(size = 14))

# 保存细胞亚型轨迹图
ggsave("subtype_trajectory.pdf", p_subtype, width = 10, height = 8, device = cairo_pdf)
ggsave("subtype_trajectory.png", p_subtype, width = 10, height = 8, dpi = 300)

## 3. 状态轨迹图
p_state <- plot_cell_trajectory(mycds3, color_by = "State", cell_size = 1, theta = 180,
                                size = 1, show_backbone = TRUE, show_branch_points = TRUE) +
  scale_color_manual(values = state_colors) +
  theme(text = element_text(size = 14))

# 保存状态轨迹图
ggsave("state_trajectory.pdf", p_state, width = 10, height = 8, device = cairo_pdf)
ggsave("state_trajectory.png", p_state, width = 10, height = 8, dpi = 300)

## 4. 条件轨迹图
p_condition <- plot_cell_trajectory(mycds3, color_by = "condition", cell_size = 1, theta = 180,
                                    size = 1, show_backbone = TRUE, show_branch_points = TRUE) +
  scale_color_manual(values = condition_colors) +
  theme(text = element_text(size = 14))

# 保存条件轨迹图
ggsave("condition_trajectory.pdf", p_condition, width = 10, height = 8, device = cairo_pdf)
ggsave("condition_trajectory.png", p_condition, width = 10, height = 8, dpi = 300)

## 5. 拟时序区间条件占比扇形图
# 提取拟时序和条件数据
pseudotime_data <- data.frame(
  pseudotime = mycds3@phenoData@data$Pseudotime,
  condition = mycds3@phenoData@data$condition
)

# 移除NA值
pseudotime_data <- pseudotime_data[!is.na(pseudotime_data$pseudotime), ]

# 创建拟时序区间
pseudotime_data$time_interval <- cut(pseudotime_data$pseudotime,
                                     breaks = c(0, 5, 10, 15, max(pseudotime_data$pseudotime)),
                                     labels = c("0-5", "5-10", "10-15", "15+"),
                                     include.lowest = TRUE)

# 计算每个区间的条件占比
condition_proportions <- pseudotime_data %>%
  group_by(time_interval, condition) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(time_interval) %>%
  mutate(proportion = count / sum(count),
         label_y = cumsum(proportion) - 0.5 * proportion)

# 创建条件占比扇形图
p_condition_pie <- ggplot(condition_proportions, aes(x = "", y = proportion, fill = condition)) +
  geom_bar(width = 1, stat = "identity", color = "black", linewidth = 0.5) +
  coord_polar("y", start = 0) +
  facet_wrap(~ time_interval, nrow = 1) +
  geom_text(aes(label = paste0(round(proportion * 100), "%")), 
            position = position_stack(vjust = 0.5),
            size = 4, color = "black", fontface = "bold") +
  scale_fill_manual(values = condition_colors) +
  theme_void() +
  theme(legend.position = "none",
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold",
                                  margin = margin(b = 10))) +
  labs(title = "Condition Distribution by Pseudotime Interval")

# 保存条件占比扇形图
ggsave("condition_pie_chart.pdf", p_condition_pie, width = 12, height = 8, device = cairo_pdf)
ggsave("condition_pie_chart.png", p_condition_pie, width = 12, height = 8, dpi = 300)

## 6. 细胞亚型占比扇形图
# 提取拟时序和细胞亚型数据
subtype_data <- data.frame(
  Pseudotime = mycds3@phenoData@data$Pseudotime,
  subtype = mycds3@phenoData@data$cell_subtype
)

# 移除NA值
subtype_data <- subtype_data[!is.na(subtype_data$Pseudotime), ]

# 创建拟时序区间
subtype_data$time_interval <- cut(subtype_data$Pseudotime,
                                  breaks = c(0, 5, 10, 15, max(subtype_data$Pseudotime)),
                                  labels = c("0-5", "5-10", "10-15", "15+"),
                                  include.lowest = TRUE)

# 计算每个区间的细胞亚型占比
subtype_proportions <- subtype_data %>%
  group_by(time_interval, subtype) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(time_interval) %>%
  mutate(proportion = count / sum(count),
         label_y = cumsum(proportion) - 0.5 * proportion)

# 创建细胞亚型占比扇形图
p_subtype_pie <- ggplot(subtype_proportions, aes(x = "", y = proportion, fill = subtype)) +
  geom_bar(width = 1, stat = "identity", color = "black", linewidth = 0.5) +
  coord_polar("y", start = 0) +
  facet_wrap(~ time_interval, nrow = 1) +
  geom_text(aes(label = paste0(round(proportion * 100), "%")), 
            position = position_stack(vjust = 0.5),
            size = 4, color = "black", fontface = "bold") +
  scale_fill_manual(values = subtype_colors) +
  theme_void() +
  theme(legend.position = "none",
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold",
                                  margin = margin(b = 10))) +
  labs(title = "Cell Subtype Distribution by Pseudotime Interval")

# 保存细胞亚型占比扇形图
ggsave("subtype_pie_chart.pdf", p_subtype_pie, width = 12, height = 8, device = cairo_pdf)
ggsave("subtype_pie_chart.png", p_subtype_pie, width = 12, height = 8, dpi = 300)

## 7. 状态中条件占比扇形图
# 提取State和condition数据
state_condition_data <- data.frame(
  state = mycds3@phenoData@data$State,
  condition = mycds3@phenoData@data$condition
)

# 移除NA值
state_condition_data <- state_condition_data[!is.na(state_condition_data$state), ]

# 计算每个State中condition的占比
state_condition_proportions <- state_condition_data %>%
  group_by(state, condition) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(state) %>%
  mutate(proportion = count / sum(count),
         label_y = cumsum(proportion) - 0.5 * proportion)

# 创建状态中条件占比扇形图
p_state_pie <- ggplot(state_condition_proportions, aes(x = "", y = proportion, fill = condition)) +
  geom_bar(width = 1, stat = "identity", color = "black", linewidth = 0.8) +
  coord_polar("y", start = 0) +
  facet_wrap(~ state, nrow = 2) +
  geom_text(aes(label = ifelse(proportion > 0.05, paste0(round(proportion * 100), "%"), "")), 
            position = position_stack(vjust = 0.5),
            size = 4, color = "black", fontface = "bold") +
  scale_fill_manual(values = condition_colors, name = "Condition") +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    strip.text = element_text(size = 12, face = "bold", margin = margin(5, 0, 5, 0)),
    plot.title = element_text(size = 16, hjust = 0.5, face = "bold",
                              margin = margin(b = 15))) +
  labs(title = "Condition Distribution by Trajectory State")

# 保存状态中条件占比扇形图
ggsave("state_condition_pie_chart.pdf", p_state_pie, width = 12, height = 10, device = cairo_pdf)
ggsave("state_condition_pie_chart.png", p_state_pie, width = 12, height = 10, dpi = 300)

## 8. 创建组合图
# 轨迹图组合
trajectory_combined <- (p_trajectory + p_subtype) / (p_state + p_condition) +
  plot_layout(guides = 'collect') +
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(face = "bold", size = 16))

# 保存组合轨迹图
ggsave("combined_trajectories.pdf", trajectory_combined, width = 16, height = 12, device = cairo_pdf)
ggsave("combined_trajectories.png", trajectory_combined, width = 16, height = 12, dpi = 300)

# 扇形图组合
pie_combined <- (p_condition_pie + p_subtype_pie) / p_state_pie +
  plot_layout(heights = c(1, 2)) +
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(face = "bold", size = 16))

# 保存组合扇形图
ggsave("combined_pie_charts.pdf", pie_combined, width = 16, height = 14, device = cairo_pdf)
ggsave("combined_pie_charts.png", pie_combined, width = 16, height = 14, dpi = 300)

## 9. 生成文件列表
cat("=== 拟时序分析图形保存完成 ===\n")
cat("生成的PDF文件:\n")
cat("- pseudotime_trajectory.pdf: 拟时序轨迹图\n")
cat("- subtype_trajectory.pdf: 细胞亚型轨迹图\n")
cat("- state_trajectory.pdf: 状态轨迹图\n")
cat("- condition_trajectory.pdf: 条件轨迹图\n")
cat("- condition_pie_chart.pdf: 条件占比扇形图\n")
cat("- subtype_pie_chart.pdf: 细胞亚型占比扇形图\n")
cat("- state_condition_pie_chart.pdf: 状态中条件占比扇形图\n")
cat("- combined_trajectories.pdf: 组合轨迹图\n")
cat("- combined_pie_charts.pdf: 组合扇形图\n")
cat("\n生成的PNG文件:\n")
cat("- pseudotime_trajectory.png: 拟时序轨迹图\n")
cat("- subtype_trajectory.png: 细胞亚型轨迹图\n")
cat("- state_trajectory.png: 状态轨迹图\n")
cat("- condition_trajectory.png: 条件轨迹图\n")
cat("- condition_pie_chart.png: 条件占比扇形图\n")
cat("- subtype_pie_chart.png: 细胞亚型占比扇形图\n")
cat("- state_condition_pie_chart.png: 状态中条件占比扇形图\n")
cat("- combined_trajectories.png: 组合轨迹图\n")
cat("- combined_pie_charts.png: 组合扇形图\n")


###轨迹相关基因分析
# 提取Term infant和BPD的细胞
term_bpd_cells <- colnames(meso_adventitial_FB)[meso_adventitial_FB$condition %in% c("Term infant", "BPD")]

# 创建子集
term_bpd_cds <- mycds3[, term_bpd_cells]

# 验证细胞数量
cat("Term infant细胞数:", sum(pData(term_bpd_cds)$condition == "Term infant"), "\n")
cat("BPD细胞数:", sum(pData(term_bpd_cds)$condition == "BPD"), "\n")

# 寻找条件特异性轨迹基因
condition_diff_genes <- differentialGeneTest(
  term_bpd_cds[expressed_genes,],
  fullModelFormulaStr = "~condition + sm.ns(Pseudotime)",
  reducedModelFormulaStr = "~sm.ns(Pseudotime)"
)

# 提取显著的条件特异性基因
condition_sig_genes <- subset(condition_diff_genes, qval < 0.05)
condition_sig_genes <- condition_sig_genes[order(condition_sig_genes$qval), ]

cat("发现", nrow(condition_sig_genes), "个Term vs BPD分化差异相关基因\n")

# 保存结果
write.table(condition_sig_genes, file = "term_bpd_trajectory_diff_genes.xls", 
            col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE)


# 选择top条件特异性基因
top_condition_genes <- rownames(condition_sig_genes)[1:min(50, nrow(condition_sig_genes))]

# 提取表达数据并按拟时序排序
library(Biobase)
expression_data <- Biobase::exprs(term_bpd_cds)[top_condition_genes, ]
pseudotime_order <- order(Biobase::pData(term_bpd_cds)$Pseudotime)
expression_data_ordered <- expression_data[, pseudotime_order]

# 创建热图注释
annotation_col <- data.frame(
  Pseudotime = Biobase::pData(term_bpd_cds)$Pseudotime[pseudotime_order],
  Condition = Biobase::pData(term_bpd_cds)$condition[pseudotime_order],
  State = Biobase::pData(term_bpd_cds)$State[pseudotime_order],
  Cell_Subtype = Biobase::pData(term_bpd_cds)$cell_subtype[pseudotime_order]
)
rownames(annotation_col) <- colnames(expression_data_ordered)

# 定义颜色
condition_colors <- c("Term infant" = "#1f77b4", "BPD" = "#ff7f0e")

# 获取细胞亚型颜色
subtype_colors <- setNames(
  scales::hue_pal()(length(unique(annotation_col$Cell_Subtype))),
  unique(annotation_col$Cell_Subtype)
)

annotation_colors <- list(
  Pseudotime = viridis::viridis(100),
  Condition = condition_colors,
  State = setNames(
    c("#fde2e7","#fbf4be","#789eb1","#e0f2fe","#bbb4d5","#acceea",
      "#f4c78c","#c9bad0","#a07a95","#d55431","#844b44","#038061","#376352"),
    unique(annotation_col$State)
  ),
  Cell_Subtype = subtype_colors
)

# 绘制条件特异性热图
p_condition_heatmap <- pheatmap::pheatmap(
  expression_data_ordered,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_colnames = FALSE,
  show_rownames = TRUE,
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  color = viridis::viridis(100),
  scale = "row",
  fontsize_row = 8,
  fontsize_col = 6,
  main = "Term vs BPD Differentiation Trajectory Genes"
)

# 保存热图
ggsave("term_bpd_trajectory_heatmap.png", p_condition_heatmap, 
       width = 16, height = 12, dpi = 300)

library(ComplexHeatmap)
library(RColorBrewer)
library(tidyverse)
library(circlize)

### 提取条件特异性轨迹基因（使用您之前分析的结果）
# condition_sig_genes 已经包含3049个Term vs BPD分化差异相关基因
condition_sig_genes <- condition_sig_genes[order(condition_sig_genes$qval), ]

# 选择top基因进行热图绘制
top_genes <- rownames(condition_sig_genes)[1:min(500, nrow(condition_sig_genes))]

### 对矩阵行列重新排序
# 提取表达矩阵并按拟时序排序细胞
matrix_ori <- Biobase::exprs(term_bpd_cds)[match(rownames(condition_sig_genes), rownames(Biobase::exprs(term_bpd_cds))), 
                                  order(Biobase::pData(term_bpd_cds)$Pseudotime)]

### 样条平滑、标准化
matrix_smooth <- t(apply(matrix_ori, 1, function(x){smooth.spline(x, df=3)$y}))
matrix_scaled <- t(apply(matrix_smooth, 1, function(x){(x-mean(x))/sd(x)}))

rownames(matrix_scaled) <- rownames(matrix_ori)
colnames(matrix_scaled) <- colnames(matrix_ori)

### 定义颜色方案
# 条件颜色
condition_cols <- c("Term infant" = "grey", "BPD" = "#c69189")

# 细胞亚型颜色（根据您数据中的cell_subtype）
subtype_levels <- unique(Biobase::pData(term_bpd_cds)$cell_subtype)
subtype_cols <- setNames(
  scales::hue_pal()(length(subtype_levels)),
  subtype_levels
)

# State颜色（如果需要）
state_levels <- unique(Biobase::pData(term_bpd_cds)$State)
state_cols <- setNames(
  viridis::viridis(length(state_levels)),
  state_levels
)

### 准备列注释信息
# 创建注释数据框
cell_anno <- data.frame(
  cell = colnames(matrix_scaled),
  Pseudotime = Biobase::pData(term_bpd_cds)[colnames(matrix_scaled), "Pseudotime"],
  Condition = Biobase::pData(term_bpd_cds)[colnames(matrix_scaled), "condition"],
  CellSubtype = Biobase::pData(term_bpd_cds)[colnames(matrix_scaled), "cell_subtype"],
  State = Biobase::pData(term_bpd_cds)[colnames(matrix_scaled), "State"]
) %>% column_to_rownames(var = "cell")

### 创建列注释
ha <- HeatmapAnnotation(
  Condition = cell_anno$Condition,
  CellSubtype = cell_anno$CellSubtype,
  State = cell_anno$State,
  Pseudotime = anno_lines(cell_anno$Pseudotime, add_points = TRUE),
  col = list(
    Condition = condition_cols,
    CellSubtype = subtype_cols,
    State = state_cols
  ),
  annotation_name_gp = gpar(fontsize = 10),
  annotation_legend_param = list(
    Condition = list(title = "Condition", title_gp = gpar(fontsize = 12)),
    CellSubtype = list(title = "Cell Subtype", title_gp = gpar(fontsize = 12)),
    State = list(title = "State", title_gp = gpar(fontsize = 12))
  )
)

### 绘制热图
p <- Heatmap(
  matrix_scaled,
  name = "z-score",
  show_row_names = F,
  show_column_names = FALSE,
  row_names_gp = gpar(fontsize = 8),
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  row_title_rot = 0,
  cluster_rows = TRUE,
  top_annotation = ha,
  cluster_row_slices = FALSE,
  cluster_columns = FALSE,
  border = TRUE,
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 12),
    labels_gp = gpar(fontsize = 10)
  )
)

print(p)

km1 <- kmeans(matrix_scaled, 4)
p2 <- Heatmap(
  matrix_scaled,
  name = "z-score",
  show_row_names = F,
  show_column_names = FALSE,
  row_names_gp = gpar(fontsize = 8),
  split = km1$cluster,
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  row_title_rot = 0,
  cluster_rows = TRUE,
  top_annotation = ha,
  cluster_row_slices = FALSE,
  cluster_columns = FALSE,
  border = TRUE,
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 12),
    labels_gp = gpar(fontsize = 10)
  )
)

print(p2)

### 保存热图
png("term_bpd_trajectory_heatmap_complex.png", width = 5.5, height = 8, units = "in", res = 600)
print(p2)
dev.off()

pdf("term_bpd_trajectory_heatmap_complex.pdf", width = 5.5, height = 8)
print(p2)
dev.off()


##GO富集分析
model <- "org.Hs.eg.db"
library(ComplexHeatmap)
library(RColorBrewer)
library(tidyverse)
library(circlize)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)

deg1 <- names(km1$cluster[km1$cluster=='1'])
full <- bitr(deg1, fromType = "SYMBOL" ,toType = c("ENSEMBL","ENTREZID") , OrgDb = model)
go1 <- enrichGO(gene = full$ENTREZID , OrgDb = model, ont = "BP", pAdjustMethod = "fdr", pvalueCutoff = 0.05, readable = TRUE)
# go1 <- head(as.data.frame(go1), 5)

deg2 <- names(km1$cluster[km1$cluster=='2'])
full <- bitr(deg2, fromType = "SYMBOL" ,toType = c("ENSEMBL","ENTREZID") , OrgDb = model)
go2 <- enrichGO(gene = full$ENTREZID , OrgDb = model, ont = "BP", pAdjustMethod = "fdr", pvalueCutoff = 0.05, readable = TRUE)
# go2 <- head(as.data.frame(go2), 5)

deg3 <- names(km1$cluster[km1$cluster=='3'])
full <- bitr(deg3, fromType = "SYMBOL" ,toType = c("ENSEMBL","ENTREZID") , OrgDb = model)
go3 <- enrichGO(gene = full$ENTREZID , OrgDb = model, ont = "BP", pAdjustMethod = "fdr", pvalueCutoff = 0.05, readable = TRUE)
# go3 <- head(as.data.frame(go3), 5)

deg4 <- names(km1$cluster[km1$cluster=='4'])
full <- bitr(deg4, fromType = "SYMBOL" ,toType = c("ENSEMBL","ENTREZID") , OrgDb = model)
go4 <- enrichGO(gene = full$ENTREZID , OrgDb = model, ont = "BP", pAdjustMethod = "fdr", pvalueCutoff = 0.05, readable = TRUE)
# go4 <- head(as.data.frame(go4), 5)


# 提取每个簇的top5 GO结果
go1_top5 <- head(as.data.frame(go1), 5)
go2_top5 <- head(as.data.frame(go2), 5)
go3_top5 <- head(as.data.frame(go3), 5)
go4_top5 <- head(as.data.frame(go4), 5)

# 添加簇标识列
go1_top5$Cluster <- "Cluster 1"
go2_top5$Cluster <- "Cluster 2"
go3_top5$Cluster <- "Cluster 3"
go4_top5$Cluster <- "Cluster 4"

# 合并所有结果
all_go <- rbind(go1_top5, go2_top5, go3_top5, go4_top5)

# 确保有足够的结果可绘制
if(nrow(all_go) > 0) {
  # 加载必要的包
  library(ggplot2)
  library(dplyr)
  library(forcats)
  
  # 方法1：分面显示四个簇的GO结果
  p1 <- ggplot(all_go, aes(x = -log10(pvalue), y = fct_reorder(Description, -log10(pvalue)))) +
    geom_col(aes(fill = Cluster), width = 0.7) +
    facet_wrap(~Cluster, scales = "free_y", ncol = 2) +
    labs(x = "-log10(p-value)", y = "GO Biological Process", 
         title = "Top 5 GO Enrichment Analysis for Each Cluster") +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 8),
      strip.text = element_text(face = "bold"),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  
  print(p1)
  
  # 方法2：将所有GO项合并显示，用颜色区分簇
  p2 <- ggplot(all_go, aes(x = -log10(pvalue), 
                           y = fct_reorder(Description, -log10(pvalue)))) +
    geom_col(aes(fill = Cluster), width = 0.7) +
    labs(x = "-log10(p-value)", y = "GO Biological Process", 
         title = "Top 5 GO Enrichment Analysis for All Clusters",
         fill = "Cluster") +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 8),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  
  print(p2)
  
  # 方法3：点图显示，同时展示pvalue和GeneRatio
  p3 <- ggplot(all_go, aes(x = GeneRatio, y = Description)) +
    geom_point(aes(size = Count, color = -log10(pvalue))) +
    facet_wrap(~Cluster, scales = "free_y", ncol = 2) +
    scale_color_gradient(low = "blue", high = "red") +
    labs(x = "Gene Ratio", y = "GO Biological Process",
         title = "Top 5 GO Enrichment Analysis - Dot Plot",
         size = "Gene Count", color = "-log10(p-value)") +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 7),
      strip.text = element_text(face = "bold"),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  
  print(p3)
  
} else {
  print("No significant GO terms found in any cluster.")
}

# 如果希望保存图片
ggsave("GO_enrichment_clusters.pdf", p1, width = 12, height = 10)
ggsave("GO_enrichment_combined.pdf", p2, width = 10, height = 8)
ggsave("GO_enrichment_dotplot.pdf", p3, width = 12, height = 10)

# 分别绘制每个簇的GO富集结果
plot_go_for_cluster <- function(go_data, cluster_name) {
  if(nrow(go_data) > 0) {
    p <- ggplot(go_data, aes(x = -log10(pvalue), 
                             y = fct_reorder(Description, -log10(pvalue)))) +
      geom_col(fill = "steelblue", width = 0.7) +
      labs(x = "-log10(p-value)", y = "",
           title = paste("Top 5 GO Terms -", cluster_name)) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold")
      )
    return(p)
  }
}

# 绘制每个簇的图
p1 <- plot_go_for_cluster(go1_top5, "Cluster 1")
p2 <- plot_go_for_cluster(go2_top5, "Cluster 2")  
p3 <- plot_go_for_cluster(go3_top5, "Cluster 3")
p4 <- plot_go_for_cluster(go4_top5, "Cluster 4")

# 使用patchwork包组合图片
library(patchwork)
combined_plot <- (p1 | p2) / (p3 | p4)
print(combined_plot)

##保存拟时序变化的结果
write.csv(data.frame(Cluster1_geneName = deg1), "deg_Cluster_1.csv", row.names = FALSE)
write.csv(data.frame(Cluster2_geneName = deg2), "deg_Cluster_2.csv", row.names = FALSE)
write.csv(data.frame(Cluster3_geneName = deg3), "deg_Cluster_3.csv", row.names = FALSE)
write.csv(data.frame(Cluster4_geneName = deg4), "deg_Cluster_4.csv", row.names = FALSE)

########对拟时序的基因基因富集分析#######
setwd("/data/nas1/xiangxiaoyu_OD/1.BPD+MSCs/1.scRNA/04_Adventitial_FB/08_Monocle3/")
# 加载必要的包
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(enrichplot)

# 读取四个cluster的差异基因
deg1 <- read.csv("deg_Cluster_1.csv")$Cluster1_geneName
deg2 <- read.csv("deg_Cluster_2.csv")$Cluster2_geneName
deg3 <- read.csv("deg_Cluster_3.csv")$Cluster3_geneName
deg4 <- read.csv("deg_Cluster_4.csv")$Cluster4_geneName


# 创建cluster列表
clusters <- list(
  Cluster1 = deg1,
  Cluster2 = deg2,
  Cluster3 = deg3,
  Cluster4 = deg4)

# 函数：进行GO和KEGG富集分析
perform_enrichment_analysis <- function(gene_list, cluster_name) {
  # 将基因符号转换为ENTREZ ID
  gene_entrez <- bitr(gene_list, 
                      fromType = "SYMBOL", 
                      toType = "ENTREZID", 
                      OrgDb = org.Hs.eg.db)
  
  if(nrow(gene_entrez) == 0) {
    return(list(GO = NULL, KEGG = NULL))
  }
  
  # GO富集分析（生物过程）
  go_enrich <- enrichGO(
    gene = gene_entrez$ENTREZID,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",  # 生物过程
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    readable = TRUE
  )
  
  # KEGG通路富集分析
  kegg_enrich <- enrichKEGG(
    gene = gene_entrez$ENTREZID,
    organism = 'hsa',
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    qvalueCutoff = 0.2
  )
  
  # 如果KEGG结果不为空，转换为可读的基因符号
  if(!is.null(kegg_enrich) && nrow(kegg_enrich) > 0) {
    kegg_enrich <- setReadable(kegg_enrich, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  }
  
  # 添加cluster信息
  if(!is.null(go_enrich) && nrow(go_enrich) > 0) {
    go_enrich@result$Cluster <- cluster_name
  }
  if(!is.null(kegg_enrich) && nrow(kegg_enrich) > 0) {
    kegg_enrich@result$Cluster <- cluster_name
  }
  
  return(list(GO = go_enrich, KEGG = kegg_enrich))
}

# 对每个cluster进行富集分析
enrichment_results <- list()
for(cluster_name in names(clusters)) {
  cat("正在分析", cluster_name, "...\n")
  enrichment_results[[cluster_name]] <- perform_enrichment_analysis(
    clusters[[cluster_name]], 
    cluster_name
  )
}

# 提取并合并GO结果
go_combined <- data.frame()
for(cluster_name in names(enrichment_results)) {
  go_result <- enrichment_results[[cluster_name]]$GO
  if(!is.null(go_result) && nrow(go_result) > 0) {
    go_df <- go_result@result
    go_combined <- rbind(go_combined, go_df)
  }
}

# 提取并合并KEGG结果
kegg_combined <- data.frame()
for(cluster_name in names(enrichment_results)) {
  kegg_result <- enrichment_results[[cluster_name]]$KEGG
  if(!is.null(kegg_result) && nrow(kegg_result) > 0) {
    kegg_df <- kegg_result@result
    kegg_combined <- rbind(kegg_combined, kegg_df)
  }
}

# 保存富集分析结果
write.csv(go_combined, "GO_enrichment_results.csv", row.names = FALSE)
write.csv(kegg_combined, "KEGG_enrichment_results.csv", row.names = FALSE)

# 函数：创建富集分析气泡图
create_enrichment_bubble_plot <- function(enrichment_data, title, top_n = 10) {
  # 如果没有数据，返回空图
  if(nrow(enrichment_data) == 0) {
    return(ggplot() + 
             geom_text(aes(x = 0.5, y = 0.5, label = "No enrichment results")) +
             theme_void())
  }
  
  # 为每个cluster选择top_n个最显著的条目
  top_pathways <- enrichment_data %>%
    group_by(Cluster) %>%
    arrange(p.adjust) %>%
    slice_head(n = top_n) %>%
    ungroup()
  
  # 如果数据仍然为空，返回空图
  if(nrow(top_pathways) == 0) {
    return(ggplot() + 
             geom_text(aes(x = 0.5, y = 0.5, label = "No significant pathways")) +
             theme_void())
  }
  
  # 创建气泡图
  p <- ggplot(top_pathways, aes(x = Cluster, y = Description)) +
    geom_point(aes(size = Count, color = -log10(p.adjust)), alpha = 0.8) +
    scale_size_continuous(
      name = "Gene Count",
      range = c(2, 8)
    ) +
    scale_color_gradientn(
      colors = c("blue", "yellow", "red"),
      name = "-log10(Adj.P)"
    ) +
    labs(
      title = title,
      x = "Cluster",
      y = "Pathway Description"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 9),
      axis.title = element_text(size = 12, face = "bold"),
      panel.grid.major = element_line(color = "grey90", size = 0.2),
      panel.grid.minor = element_blank(),
      legend.position = "right",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 8)
    )
  
  return(p)
}

# 创建GO气泡图
go_bubble_plot <- create_enrichment_bubble_plot(
  go_combined, 
  "GO Biological Process Enrichment Analysis"
)

# 创建KEGG气泡图
kegg_bubble_plot <- create_enrichment_bubble_plot(
  kegg_combined, 
  "KEGG Pathway Enrichment Analysis"
)

# 显示图形
print(go_bubble_plot)
print(kegg_bubble_plot)

# 保存图形
ggsave("GO_enrichment_bubble_plot.pdf", go_bubble_plot, width = 5.5, height = 10)
ggsave("KEGG_enrichment_bubble_plot.pdf", kegg_bubble_plot, width = 4.9, height = 10)

# 如果使用patchwork包，可以组合两个图形
if(require(patchwork)) {
  combined_plot <- go_bubble_plot / kegg_bubble_plot + 
    plot_annotation(tag_levels = 'A')
  
  print(combined_plot)
  ggsave("combined_enrichment_plots.pdf", combined_plot, width = 14, height = 16)
}

# 可选：创建更详细的可视化，显示所有cluster的top通路
create_detailed_bubble_plot <- function(enrichment_data, title, top_n = 15) {
  # 如果没有数据，返回空图
  if(nrow(enrichment_data) == 0) {
    return(ggplot() + 
             geom_text(aes(x = 0.5, y = 0.5, label = "No enrichment results")) +
             theme_void())
  }
  
  # 为每个cluster选择top_n个最显著的条目
  top_pathways <- enrichment_data %>%
    group_by(Cluster) %>%
    arrange(p.adjust) %>%
    slice_head(n = top_n) %>%
    ungroup()
  
  # 如果数据仍然为空，返回空图
  if(nrow(top_pathways) == 0) {
    return(ggplot() + 
             geom_text(aes(x = 0.5, y = 0.5, label = "No significant pathways")) +
             theme_void())
  }
  
  # 截断过长的描述
  top_pathways$Description_short <- str_trunc(top_pathways$Description, width = 50)
  
  # 创建气泡图
  p <- ggplot(top_pathways, aes(x = Cluster, y = Description_short)) +
    geom_point(aes(size = Count, color = -log10(p.adjust)), alpha = 0.8) +
    scale_size_continuous(
      name = "Gene Count",
      range = c(2, 10),
      breaks = c(5, 10, 15, 20)
    ) +
    scale_color_gradientn(
      colors = c("#4575b4", "#91bfdb", "#e0f3f8", "#fee090", "#fc8d59", "#d73027"),
      name = "-log10(Adj.P)"
    ) +
    labs(
      title = title,
      x = "Cluster",
      y = "Pathway Description"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold"),
      axis.text.y = element_text(size = 10),
      axis.title = element_text(size = 14, face = "bold"),
      panel.grid.major = element_line(color = "grey90", size = 0.2),
      panel.grid.minor = element_blank(),
      legend.position = "right",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10)
    )
  
  return(p)
}

# 创建详细的气泡图
go_detailed_plot <- create_detailed_bubble_plot(
  go_combined, 
  "GO Biological Process Enrichment Analysis",
  top_n = 15
)

kegg_detailed_plot <- create_detailed_bubble_plot(
  kegg_combined, 
  "KEGG Pathway Enrichment Analysis", 
  top_n = 15
)

# 显示和保存详细图形
print(go_detailed_plot)
print(kegg_detailed_plot)

ggsave("GO_enrichment_detailed_plot.pdf", go_detailed_plot, width = 14, height = 12)
ggsave("KEGG_enrichment_detailed_plot.pdf", kegg_detailed_plot, width = 14, height = 12)

# 输出每个cluster的富集分析统计
cat("GO富集分析统计:\n")
for(cluster_name in names(clusters)) {
  go_count <- nrow(enrichment_results[[cluster_name]]$GO@result)
  kegg_count <- nrow(enrichment_results[[cluster_name]]$KEGG@result)
  cat(cluster_name, ": GO -", go_count, "pathways, KEGG -", kegg_count, "pathways\n")
}

cat("分析完成！所有结果已保存。\n")


#####################3.6后续分析###################
# 计算BPD特异性指数 (BPD细胞数 / Term细胞数)
# 计算每个cluster的条件特异性指数
trajectory_cluster1 <- read.csv(file = "/data/nas1/xiangxiaoyu_OD/1.BPD+MSCs/1.scRNA/04_Adventitial_FB/09_Continued_analysis/deg_Cluster_1.csv")
trajectory_cluster2 <- read.csv(file = "/data/nas1/xiangxiaoyu_OD/1.BPD+MSCs/1.scRNA/04_Adventitial_FB/08_Monocle3/deg_Cluster_2.csv")
trajectory_cluster3 <- read.csv(file = "/data/nas1/xiangxiaoyu_OD/1.BPD+MSCs/1.scRNA/04_Adventitial_FB/08_Monocle3/deg_Cluster_3.csv")
trajectory_cluster4 <- read.csv(file = "/data/nas1/xiangxiaoyu_OD/1.BPD+MSCs/1.scRNA/04_Adventitial_FB/08_Monocle3/deg_Cluster_4.csv")


cell_counts <- data.frame(
  Cluster = c("Cluster0", "Cluster1", "Cluster2", "Cluster3", "Cluster4", "Cluster5"),
  Term = c(40, 90, 287, 236, 158, 1),
  BPD = c(385, 263, 29, 10, 22, 1),
  BPD_PH = c(0, 0, 0, 0, 0, 98)
)

# 计算BPD特异性指数 (BPD细胞数 / Term细胞数)
cell_counts$BPD_specificity <- cell_counts$BPD / cell_counts$Term
cell_counts$Total_cells <- cell_counts$Term + cell_counts$BPD + cell_counts$BPD_PH

print(cell_counts)
# Cluster Term BPD BPD_PH BPD_specificity Total_cells
# 1 Cluster0   40 385      0      9.62500000         425
# 2 Cluster1   90 263      0      2.92222222         353
# 3 Cluster2  287  29      0      0.10104530         316
# 4 Cluster3  236  10      0      0.04237288         246
# 5 Cluster4  158  22      0      0.13924051         180
# 6 Cluster5    1   1     98      1.00000000         100

# [1] "/data/home/xiangxiaoyu_OD/.conda/envs/sc_analysis/lib/R/library"
# [2] "/data/home/xiangxiaoyu_OD/R/x86_64-conda-linux-gnu-library/4.3" 
# [3] "/data/nas2/software/miniconda3/envs/R4.3.3/lib/R/library" 
library(ggvenn)
library(VennDiagram)

venn_data <- list(DEGs = significant_genes$gene,
                  Trajectory_Genes = unique(c(trajectory_cluster1$Cluster1_geneName,
                                              trajectory_cluster2$Cluster2_geneName,
                                              trajectory_cluster3$Cluster3_geneName,
                                              trajectory_cluster4$Cluster4_geneName)))
ggvenn(venn_data,
       show_percentage = F,
       stroke_color = "white",
       fill_color = c("#b2e7cb","#b2d4ec"),
       set_name_color = c("#4a9b83","#1d6295"))

inter <- get.venn.partitions(venn_data)
for (i in 1:nrow(inter)) inter[i,'values'] <- paste(inter[[i,'..values..']], collapse = ', ')
write.table(inter[-c(5, 6)], './09_Continued_analysis/venn_inter.txt', row.names = FALSE, sep = '\t', quote = FALSE)
ggVenn_candiated_genes <- inter$..values..[[1]]
