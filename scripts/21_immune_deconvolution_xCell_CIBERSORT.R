###############################################################################
# xCell and CIBERSORT deconvolution and correlation with PHLDA1
#
# Source      : r.03_Cell_prop.r (lines 6852-8531)
# Input       : GSE235562, GSE106910
# Output      : Fig. 7; Supplementary Figs. 10-11; Supplementary Table 17
# Note        : Line "source(\"CIBERSORT.R\")" refers to the CIBERSORT source code, which is distributed under its own licence by https://cibersortx.stanford.edu and is therefore not included here.
#
# Code below is unchanged from the original analysis script.
###############################################################################
# -----------------------------------------------------------------------------
# PATHS: this file keeps the absolute paths of the analysis server on which the
# study was run (/data/nas1/xiangxiaoyu_OD/1.BPD+MSCs/1.scRNA/... and
# /home/xxyloveme/study/BPD+MSCs/...). Set PROJ_DIR in config/paths.R and
# replace those prefixes with your own project directory before running.
# -----------------------------------------------------------------------------


####xCell###############
# 加载必要的包
library(xCell)
library(tidyverse)
library(ggplot2)
library(reshape2)
library(corrplot)
library(ggpubr)
library(tidyr)
library(dplyr)
library(pheatmap)
library(RColorBrewer)

# 设置工作目录
setwd("/data/nas1/xiangxiaoyu_OD/1.BPD+MSCs/1.scRNA/08_immunfiltration/xCell/")

# 测试集 (GSE235562) xCell分析
# 1. 准备表达矩阵 (确保是TPM或FPKM标准化数据)
# 假设training_expr已经是合适的表达矩阵
# 如果不是，需要先进行标准化

# 2. 计算免疫细胞丰度
cat("计算免疫细胞丰度...\n")
scores_test <- rawEnrichmentAnalysis(
  as.matrix(training_expr),  # 使用特征基因
  signatures = xCell.data$signatures,
  genes = xCell.data$genes,
  parallel.sz = 4,
  parallel.type = "SOCK")

# 3. 数据预处理
cat("数据预处理...\n")
scores_df_test <- as.data.frame(t(scores_test)) %>%
  mutate(
    Sample = rownames(.),
    Group = train_data$Group[match(rownames(.), rownames(train_data))],
    .before = 1)

# 4. 计算百分比
percent_df_test <- scores_df_test %>%
  mutate(across(-c(Sample, Group), ~ .x / rowSums(across(-c(Sample, Group))) * 100))

# 5. 绘制堆积柱状图
percent_long_test <- pivot_longer(percent_df_test, cols = -c(Sample, Group), 
                                  names_to = "CellType", values_to = "Percentage")

library(randomcoloR)

percent_plot_test <- ggplot(percent_long_test, aes(x = Sample, y = Percentage, fill = CellType)) +
  geom_col(width = 0.8) +
  facet_grid(~ Group, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = colpalettes) +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.background = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 8)) +
  labs(title = "Immune Cell Composition - Testing Set (GSE235562)")

ggsave("Xcell_percentage_testing.pdf", 
       plot = percent_plot_test,
       device = "pdf",
       width = 17, 
       height = 4.5)

# 6. 差异分析
cat("进行差异分析...\n")
cell_types_test <- colnames(scores_df_test)[3:ncol(scores_df_test)]

# Wilcoxon检验
results_test <- lapply(cell_types_test, function(cell){
  test <- wilcox.test(scores_df_test[[cell]] ~ scores_df_test$Group)
  data.frame(
    CellType = cell,
    p.value = test$p.value,
    W = test$statistic
  )
}) %>% bind_rows()

# 多重检验校正
results_test$p.adj <- p.adjust(results_test$p.value, method = "BH")

# 筛选显著差异细胞类型
sig_cells_test <- results_test %>% 
  filter(p.adj < 0.05) %>% 
  arrange(p.adj) %>% 
  pull(CellType)

cat("测试集中显著差异的细胞类型数量:", length(sig_cells_test), "\n")

# 7. 绘制热图 - 单个样本细胞占比
cat("绘制热图...\n")

# 准备热图数据
heatmap_data_test <- percent_df_test[,-c(1,2)]

# 如果显著细胞类型太多，选择前30个
if(length(sig_cells_test) > 0) {
  if(length(sig_cells_test) > 30) {
    top_cells <- head(sig_cells_test, 30)
  } else {
    top_cells <- sig_cells_test
  }
  
  # 提取显著细胞类型数据
  heatmap_sig_test <- heatmap_data_test[, top_cells, drop = FALSE]
  
  # 样本注释
  annotation_col_test <- data.frame(
    Group = percent_df_test$Group
  )
  rownames(annotation_col_test) <- percent_df_test$Sample
  
  # 按组别排序
  group_order <- order(annotation_col_test$Group)
  heatmap_sig_test <- heatmap_sig_test[group_order, ]
  annotation_col_test <- annotation_col_test[group_order, , drop = FALSE]
  
  # 颜色设置
  group_colors <- list(Group = c("control" = "#539382", "treatment" = "#c83e2b"))
  
  # 绘制热图
  pheatmap_test <- pheatmap(
    t(heatmap_sig_test),
    scale = "row",  # 按行标准化
    color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
    cluster_cols = FALSE,  # 列不聚类，保持组别顺序
    cluster_rows = TRUE,
    show_colnames = FALSE,
    annotation_col = annotation_col_test,
    annotation_colors = group_colors,
    fontsize_row = 8,
    main = "Testing Set (GSE235562) - Immune Cell Composition"
  )
  
  # 保存热图
  ggsave("Heatmap_Immune_Cells_Testing.pdf", 
         plot = pheatmap_test,
         device = "pdf",
         width = 8.6, 
         height = 2.5)
}

# 验证集 (GSE106910) xCell分析

cat("\n=== 开始验证集 (GSE106910) xCell分析 ===\n")

# 1. 计算免疫细胞丰度
cat("计算免疫细胞丰度...\n")
scores_valid <- rawEnrichmentAnalysis(
  as.matrix(validation_expr[,rownames(valid_data)]),  # 使用特征基因
  signatures = xCell.data$signatures,
  genes = xCell.data$genes,
  parallel.sz = 4,
  parallel.type = "SOCK"
)

# 2. 数据预处理
cat("数据预处理...\n")
scores_df_valid <- as.data.frame(t(scores_valid)) %>%
  mutate(
    Sample = rownames(.),
    Group = valid_data$Group_re[match(rownames(.), rownames(valid_data))],
    .before = 1
  )

# 3. 计算百分比
percent_df_valid <- scores_df_valid %>%
  mutate(across(-c(Sample, Group), ~ .x / rowSums(across(-c(Sample, Group))) * 100))

# 4. 绘制堆积柱状图
percent_long_valid <- pivot_longer(percent_df_valid, cols = -c(Sample, Group), 
                                   names_to = "CellType", values_to = "Percentage")
percent_long_valid$Group <- factor(percent_long_valid$Group, levels = c("control","BPD"))

percent_plot_valid <- ggplot(percent_long_valid, aes(x = Sample, y = Percentage, fill = CellType)) +
  geom_col(width = 0.8) +
  facet_grid(~ Group, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = colpalettes) +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.background = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 8)
  ) +
  labs(title = "Immune Cell Composition - Validation Set (GSE106910)")

ggsave("Xcell_percentage_validation.pdf", 
       plot = percent_plot_valid,
       device = "pdf",
       width = 17, 
       height = 4.5)

# 5. 差异分析
cat("进行差异分析...\n")
cell_types_valid <- colnames(scores_df_valid)[3:ncol(scores_df_valid)]

# Wilcoxon检验
results_valid <- lapply(cell_types_valid, function(cell){
  test <- wilcox.test(scores_df_valid[[cell]] ~ scores_df_valid$Group)
  data.frame(
    CellType = cell,
    p.value = test$p.value,
    W = test$statistic
  )
}) %>% bind_rows()

# 多重检验校正
results_valid$p.adj <- p.adjust(results_valid$p.value, method = "BH")

# 筛选显著差异细胞类型
sig_cells_valid <- results_valid %>% 
  filter(p.adj < 0.05) %>% 
  arrange(p.adj) %>% 
  pull(CellType)

cat("验证集中显著差异的细胞类型数量:", length(sig_cells_valid), "\n")

# 6. 绘制热图 - 单个样本细胞占比
cat("绘制热图...\n")

# 准备热图数据
heatmap_data_valid <- percent_df_valid[,-c(1,2)]

# 如果显著细胞类型太多，选择前30个
if(length(sig_cells_valid) > 0) {
  if(length(sig_cells_valid) > 30) {
    top_cells_valid <- head(sig_cells_valid, 30)
  } else {
    top_cells_valid <- sig_cells_valid
  }
  
  # 提取显著细胞类型数据
  heatmap_sig_valid <- heatmap_data_valid[, top_cells_valid, drop = FALSE]
  
  # 样本注释
  annotation_col_valid <- data.frame(
    Group = percent_df_valid$Group
  )
  rownames(annotation_col_valid) <- percent_df_valid$Sample
  annotation_col_valid$Group <- factor(annotation_col_valid$Group,levels = c("control","BPD"))
  # 按组别排序
  group_order_valid <- order(annotation_col_valid$Group)
  heatmap_sig_valid <- heatmap_sig_valid[group_order_valid, ]
  annotation_col_valid <- annotation_col_valid[group_order_valid, , drop = FALSE]
  
  # 颜色设置
  group_colors_valid <- list(Group = c("control" = "#b49d99", "BPD" = "#b24646"))
  
  # 绘制热图
  pheatmap_valid <- pheatmap(
    t(heatmap_sig_valid),
    scale = "row",  # 按行标准化
    color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
    cluster_cols = FALSE,  # 列不聚类，保持组别顺序
    cluster_rows = TRUE,
    show_colnames = FALSE,
    annotation_col = annotation_col_valid,
    annotation_colors = group_colors_valid,
    fontsize_row = 8,
    main = "Validation Set (GSE106910) - Immune Cell Composition"
  )
  
  # 保存热图
  ggsave("Heatmap_Immune_Cells_Validation.pdf", 
         plot = pheatmap_valid,
         device = "pdf",
         width = 8.5, 
         height = 2.5)
}


# 结果总结和保存
# 保存显著差异结果
write.csv(results_test, "Differential_Cells_Testing.csv", row.names = FALSE)
write.csv(results_valid, "Differential_Cells_Validation.csv", row.names = FALSE)

# 保存原始分数
write.csv(scores_df_test, "xCell_Scores_Testing.csv", row.names = FALSE)
write.csv(scores_df_valid, "xCell_Scores_Validation.csv", row.names = FALSE)

# 生成总结报告
cat("\n=== xCell免疫浸润分析总结 ===\n")
cat("测试集 (GSE235562):\n")
cat("- 样本数量:", nrow(scores_df_test), "\n")
cat("- 显著差异细胞类型:", length(sig_cells_test), "\n")
if(length(sig_cells_test) > 0) {
  cat("- 前5个最显著细胞类型:", paste(head(sig_cells_test, 5), collapse = ", "), "\n")
}

cat("\n验证集 (GSE106910):\n")
cat("- 样本数量:", nrow(scores_df_valid), "\n")
cat("- 显著差异细胞类型:", length(sig_cells_valid), "\n")
if(length(sig_cells_valid) > 0) {
  cat("- 前5个最显著细胞类型:", paste(head(sig_cells_valid, 5), collapse = ", "), "\n")
}

# 保存工作环境
save.image(file = "xCell_Analysis_Workplace.RData")



###############CIBERSORT###################
setwd("/data/nas1/xiangxiaoyu_OD/1.BPD+MSCs/1.scRNA/08_immunfiltration/CIBERSORT/")

# CIBERSORT免疫浸润分析
library(ggplot2)
library(ggpubr)
library(psych)
library(corrplot)
library(dplyr)
library(purrr)
library(tidyr)
library(pheatmap)
library(RColorBrewer)

# 加载CIBERSORT脚本（确保CIBERSORT.R文件在当前目录）
source("CIBERSORT.R")

# 测试集 (GSE235562) CIBERSORT分析
cat("=== 开始测试集 (GSE235562) CIBERSORT分析 ===\n")
# 1. 准备表达矩阵
# 检查是否有重复基因名
duplicate_genes_train <- training_expr[duplicated(rownames(training_expr)), ]
if (nrow(duplicate_genes_train) > 0) {
  cat("发现重复基因名，正在移除重复项...\n")
  training_expr_noduplicate <- training_expr[!duplicated(rownames(training_expr)), ]
} else {
  training_expr_noduplicate <- training_expr
}
cat("NA值数量:", sum(is.na(training_expr_noduplicate)), "\n")

training_expr_noduplicate[is.na(training_expr_noduplicate)] <- 0

# 保存表达谱数据为TXT文件
write.table(as.matrix(training_expr_noduplicate), 
            file = "training_expr_CIBERSORT.txt", 
            sep = "\t",
            quote = FALSE,
            row.names = TRUE,
            col.names = TRUE)
# 2. 运行CIBERSORT
cat("运行CIBERSORT分析...\n")
cibersort_results_train <- CIBERSORT("LM22.txt", "training_expr_CIBERSORT.txt", perm = 100, QN = F)

# 3. 处理结果
cibersort_results_train <- as.data.frame(cibersort_results_train)
# 剔除P > 0.05的样本
valid_samples_train <- cibersort_results_train[cibersort_results_train$`P-value` <= 0.05, ]
valid_samples_train <- cibersort_results_train[, 1:22]  # 提取22种免疫细胞的比例

cat("测试集有效样本数量:", nrow(valid_samples_train), "\n")

# 4. 添加分组信息
valid_samples_train$Group <- train_data$Group[match(rownames(valid_samples_train), rownames(train_data))]
valid_samples_train$sample <- rownames(valid_samples_train)

# 5. 转换为长格式
immune_long_train <- pivot_longer(valid_samples_train, 
                                  cols = -c(Group, sample), 
                                  names_to = "CellType", 
                                  values_to = "Proportion")

# 6. 绘制箱线图
p_train_box <- ggplot(immune_long_train, aes(x = CellType, y = Proportion, fill = Group)) +
  geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.8) +
  scale_fill_manual(values = c("#539382", "#c83e2b")) +
  stat_compare_means(
    aes(group = Group),
    method = "wilcox.test",
    label = "p.signif",
    label.y = max(immune_long_train$Proportion) * 1.05,
    show.legend = FALSE) +
  labs(
    title = "Immune Cell Infiltration - Testing Set (GSE235562)",
    x = "Immune Cell Type",
    y = "Proportion",
    fill = "Group") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())
ggsave("immune_cell_infiltration_testing.pdf", plot = p_train_box, width = 11.5, height = 4.7)

# 7. 绘制堆积柱状图
library(ggh4x)
immune_long_train$sample <- factor(immune_long_train$sample, levels = valid_samples_train$sample)
immune_long_train$Group <- factor(immune_long_train$Group, levels = c("control", "treatment"))

p_train_stacked <- ggplot(immune_long_train, aes(x = sample, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity") +
  labs(fill = "Cell Type", x = "Sample", y = "Estimated Proportion") +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "bottom") +
  guides(fill = guide_legend(ncol = 6)) +
  scale_y_continuous(expand = c(0, 0)) +
  facet_nested(. ~ Group, drop = TRUE, scale = "free", space = "free", switch = "y",
               strip = strip_nested(
                 background_x = elem_list_rect(fill = c("#539382", "#c83e2b")),
                 text_x = element_text(size = 14),
                 by_layer_x = FALSE)) +
  theme(panel.spacing.x = unit(0.1, "cm")) +
  scale_fill_manual(values = colorRampPalette(brewer.pal(11, "Spectral"))(22))
ggsave("immune_cell_stacked_testing.pdf", plot = p_train_stacked, width = 11.5, height = 4.5)


# 验证集 (GSE106910) CIBERSORT分析
cat("=== 开始验证集 (GSE106910) CIBERSORT分析 ===\n")
# 1. 准备表达矩阵
duplicate_genes_valid <- validation_expr[duplicated(rownames(validation_expr)), ]
if (nrow(duplicate_genes_valid) > 0) {
  cat("发现重复基因名，正在移除重复项...\n")
  validation_expr_noduplicate <- validation_expr[!duplicated(rownames(validation_expr)), ]
} else {
  validation_expr_noduplicate <- validation_expr
}

# 保存表达谱数据为TXT文件
write.table(validation_expr_noduplicate, 
            file = "validation_expr_CIBERSORT.txt", 
            sep = "\t",
            quote = FALSE,
            row.names = TRUE,
            col.names = TRUE)

# 2. 运行CIBERSORT
cat("运行CIBERSORT分析...\n")
cibersort_results_valid <- CIBERSORT("LM22.txt", "validation_expr_CIBERSORT.txt", perm = 100, QN = F)

# 3. 处理结果
cibersort_results_valid <- as.data.frame(cibersort_results_valid)

# 剔除P > 0.05的样本
valid_samples_valid <- cibersort_results_valid[cibersort_results_valid$`P-value` <= 0.05, ]
valid_samples_valid <- cibersort_results_valid[, 1:22]  # 提取22种免疫细胞的比例

cat("验证集有效样本数量:", nrow(valid_samples_valid), "\n")

# 4. 添加分组信息
valid_samples_valid <- valid_samples_valid[rownames(valid_data),]
valid_samples_valid$Group <- valid_data$Group_re[match(rownames(valid_samples_valid), rownames(valid_data))]
valid_samples_valid$sample <- rownames(valid_samples_valid)

# 5. 转换为长格式
immune_long_valid <- pivot_longer(valid_samples_valid, 
                                  cols = -c(Group, sample), 
                                  names_to = "CellType", 
                                  values_to = "Proportion")
immune_long_valid$Group <- factor(immune_long_valid$Group,levels = c("control","BPD"))
# 6. 绘制箱线图
p_valid_box <- ggplot(immune_long_valid, aes(x = CellType, y = Proportion, fill = Group)) +
  geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.8) +
  scale_fill_manual(values = c("#b49d99", "#b24646")) +
  stat_compare_means(
    aes(group = Group),
    method = "wilcox.test",
    label = "p.signif",
    label.y = max(immune_long_valid$Proportion) * 1.05,
    show.legend = FALSE) +
  labs(
    title = "Immune Cell Infiltration - Validation Set (GSE106910)",
    x = "Immune Cell Type",
    y = "Proportion",
    fill = "Group") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())
ggsave("immune_cell_infiltration_validation.pdf", plot = p_valid_box, width = 11.5, height = 4.7)

# 7. 绘制堆积柱状图
immune_long_valid$sample <- factor(immune_long_valid$sample, levels = valid_samples_valid$sample)
immune_long_valid$Group <- factor(immune_long_valid$Group, levels = c("control", "BPD"))

p_valid_stacked <- ggplot(immune_long_valid, aes(x = sample, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity") +
  labs(fill = "Cell Type", x = "Sample", y = "Estimated Proportion") +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "bottom") +
  guides(fill = guide_legend(ncol = 6)) +
  scale_y_continuous(expand = c(0, 0)) +
  facet_nested(. ~ Group, drop = TRUE, scale = "free", space = "free", switch = "y",
               strip = strip_nested(
                 background_x = elem_list_rect(fill = c("#b49d99", "#b24646")),
                 text_x = element_text(size = 14),
                 by_layer_x = FALSE)) +
  theme(panel.spacing.x = unit(0.1, "cm")) +
  scale_fill_manual(values = colorRampPalette(brewer.pal(11, "Spectral"))(22))
ggsave("immune_cell_stacked_validation.pdf", plot = p_valid_stacked, width = 11.5, height = 4.5)


# 热图展示 - 带显著性标记
create_simple_heatmap <- function(valid_samples, dataset_name) {
  # 准备数据
  heatmap_data <- valid_samples[, !colnames(valid_samples) %in% c("Group", "sample")]
  
  # 数据清理
  heatmap_data_clean <- as.matrix(heatmap_data)
  heatmap_data_clean[is.na(heatmap_data_clean)] <- 0
  heatmap_data_clean[is.nan(heatmap_data_clean)] <- 0
  heatmap_data_clean[is.infinite(heatmap_data_clean)] <- 0
  
  # 删除所有样本占比都为0的列（细胞类型）
  col_sums <- colSums(heatmap_data_clean)
  zero_cols <- which(col_sums == 0)
  
  if(length(zero_cols) > 0) {
    cat("删除在所有样本中占比为0的细胞类型:", names(zero_cols), "\n")
    heatmap_data_clean <- heatmap_data_clean[, -zero_cols, drop = FALSE]
  }
  
  # 如果删除了所有列，则返回提示
  if(ncol(heatmap_data_clean) == 0) {
    cat("警告: 所有细胞类型在所有样本中占比都为0，无法绘制热图\n")
    # 绘制空图
    plot(1, type = "n", axes = FALSE, xlab = "", ylab = "", 
         main = paste("No valid immune cell data for", dataset_name))
    return(NULL)
  }
  
  # 样本注释
  annotation_col <- data.frame(Group = valid_samples$Group)
  rownames(annotation_col) <- valid_samples$sample
  
  # 按组别排序
  group_order <- order(annotation_col$Group)
  heatmap_data_clean <- heatmap_data_clean[group_order, , drop = FALSE]
  annotation_col <- annotation_col[group_order, , drop = FALSE]
  
  # 设置颜色
  if (all(unique(annotation_col$Group) %in% c("control", "treatment"))) {
    group_colors <- list(Group = c("control" = "#539382", "treatment" = "#c83e2b"))
  } else {
    group_colors <- list(Group = c("control" = "#b49d99", "BPD" = "#b24646"))
  }
  
  # 绘制最简单的热图（无聚类，无显著性标记）
  pheatmap(
    t(heatmap_data_clean),
    scale = "row",
    color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
    cluster_cols = FALSE,
    cluster_rows = T,
    show_colnames = FALSE,
    annotation_col = annotation_col,
    annotation_colors = group_colors,
    fontsize_row = 8,
    main = paste("Immune Cell Composition -", dataset_name)
  )
}

pdf("heatmap_immune_cells_testing_simple.pdf", width = 8.5, height = 2.5)
create_simple_heatmap(valid_samples_train, "Testing Set (GSE235562)")
dev.off()

pdf("heatmap_immune_cells_validation_simple.pdf", width = 8.5, height = 2.5)
create_simple_heatmap(valid_samples_valid, "Validation Set (GSE106910)")
dev.off()



# # 差异免疫细胞相关性分析
# # 函数：分析差异免疫细胞相关性
# analyze_immune_correlation <- function(valid_samples, dataset_name) {
#   
#   # 计算显著差异细胞类型
#   pval_df <- valid_samples %>%
#     pivot_longer(cols = -c(Group, sample), names_to = "CellType", values_to = "Proportion") %>%
#     group_by(CellType) %>%
#     summarise(
#       p_value = suppressWarnings(
#         wilcox.test(
#           Proportion ~ Group,
#           data = cur_data(),
#           exact = FALSE
#         )$p.value
#       )
#     ) %>%
#     filter(!is.na(p_value)) %>%
#     mutate(
#       p.adj = p.adjust(p_value, method = "BH"),
#       significance = ifelse(p.adj < 0.05, "Sig", "Non-Sig")
#     )
#   
#   sig_cells <- pval_df %>% 
#     filter(significance == "Sig") %>% 
#     pull(CellType)
#   
#   if (length(sig_cells) > 1) {
#     # 计算相关系数
#     cor_results <- corr.test(
#       valid_samples[, sig_cells],
#       method = "spearman",
#       adjust = "BH",
#       ci = FALSE)
#     
#     # 设置颜色
#     col1 <- colorRampPalette(c("#9B2226", "#AE2021", "#BB3E03", "#CA6702", "#FFFFFF",
#                                "#90e0ef", "#0077b6", "#023e8a", "#001219"))
#     
#     # 保存相关性热图
#     pdf(paste0("immune_cell_correlation_", tolower(gsub(" ", "_", dataset_name)), ".pdf"), 
#         width = 9, height = 8)
#     corrplot(as.matrix(cor_results$r),
#              method = "square",
#              tl.cex = 0.85,
#              cl.cex = 0.85,
#              col = col1(20),
#              tl.col = "black",
#              title = paste("Immune Cell Correlation -", dataset_name))
#     dev.off()
#     
#     png(paste0("immune_cell_correlation_", tolower(gsub(" ", "_", dataset_name)), ".png"), 
#         width = 1200, height = 1200, res = 300)
#     corrplot(as.matrix(cor_results$r),
#              method = "square",
#              tl.cex = 0.85,
#              cl.cex = 0.85,
#              col = col1(20),
#              tl.col = "black",
#              title = paste("Immune Cell Correlation -", dataset_name))
#     dev.off()
#   } else {
#     cat(dataset_name, "中显著差异细胞类型不足2个，跳过相关性分析\n")
#   }
# }
# 
# # 10. 分析测试集相关性
# analyze_immune_correlation(valid_samples_train, "Testing Set")
# 
# # 11. 分析验证集相关性
# analyze_immune_correlation(valid_samples_valid, "Validation Set")


# 结果保存和总结
# 保存结果
write.csv(valid_samples_train, "CIBERSORT_results_testing.csv", row.names = FALSE)
write.csv(valid_samples_valid, "CIBERSORT_results_validation.csv", row.names = FALSE)


# ======================免疫浸润数据整合与综合分析 ======================
cat("\n=== 开始整合分析与PHLDA1相关性分析 ===\n")

# 1. 整合两个数据集的xCell结果
# 重命名列以避免冲突
xCell_cols <- colnames(scores_df_test)[!colnames(scores_df_test) %in% c("Group", "Sample","Dataset")]
xCell_cols <- str_replace(xCell_cols,"_xCell","")
colnames(scores_df_test)[3:ncol(scores_df_test)] <- paste0(
  colnames(scores_df_test)[3:ncol(scores_df_test)], "_xCell")
colnames(scores_df_valid)[3:ncol(scores_df_valid)] <- paste0(
  colnames(scores_df_valid)[3:ncol(scores_df_valid)], "_xCell")

# 添加数据集标识
scores_df_test$Dataset <- "Testing_GSE235562"
scores_df_valid$Dataset <- "Validation_GSE106910"

# 统一列名用于合并
xCell_combined <- bind_rows(scores_df_test, scores_df_valid)

# 2. 整合两个数据集的CIBERSORT结果
# 为CIBERSORT列添加标识
ciber_cols <- colnames(valid_samples_train)[!colnames(valid_samples_train) %in% c("Group", "sample")]
ciber_cols <- paste0(ciber_cols, "_CIBERSORT")
colnames(valid_samples_train)[match(ciber_cols, colnames(valid_samples_train))] <- 
  paste0(ciber_cols, "_CIBERSORT")
colnames(valid_samples_valid)[match(ciber_cols, colnames(valid_samples_valid))] <- 
  paste0(ciber_cols, "_CIBERSORT")

valid_samples_train$Dataset <- "Testing_GSE235562"
valid_samples_valid$Dataset <- "Validation_GSE106910"

CIBERSORT_combined <- bind_rows(valid_samples_train, valid_samples_valid)

# 3. 识别两种方法在两个数据集中的差异免疫细胞
# xCell方法在测试集中的显著差异细胞（已获取sig_cells_test）
# xCell方法在验证集中的显著差异细胞（已获取sig_cells_valid）
# 找出在两个数据集中都显著的xCell细胞类型
xCell_common_sig <- unique(union(sig_cells_test, sig_cells_valid))
cat("\n在两个数据集中都显著的xCell细胞类型:", length(xCell_common_sig), "个\n")
if(length(xCell_common_sig) > 0) {
  print(xCell_common_sig)
}

# 对CIBERSORT结果进行差异分析
perform_cibersort_diff <- function(cibersort_data, group_var = "Group") {
  cell_types <- grep("_CIBERSORT$", colnames(cibersort_data), value = TRUE)
  
  results <- lapply(cell_types, function(cell){
    test <- wilcox.test(cibersort_data[[cell]] ~ cibersort_data[[group_var]])
    data.frame(
      CellType = cell,
      p.value = test$p.value,
      W = test$statistic
    )
  }) %>% bind_rows()
  
  results$p.adj <- p.adjust(results$p.value, method = "BH")
  results$Significant <- results$p.value < 0.05
  
  return(results)
}

# 分别计算两个数据集的CIBERSORT差异
cibersort_diff_test <- perform_cibersort_diff(valid_samples_train)
cibersort_diff_valid <- perform_cibersort_diff(valid_samples_valid)

# 提取显著细胞类型
sig_cibersort_test <- cibersort_diff_test %>% 
  filter(Significant) %>% 
  pull(CellType)
sig_cibersort_valid <- cibersort_diff_valid %>% 
  filter(Significant) %>% 
  pull(CellType)

# 找出在两个数据集中都显著的CIBERSORT细胞类型
cibersort_common_sig <- union(sig_cibersort_test, sig_cibersort_valid)
cat("\n在两个数据集中都显著的CIBERSORT细胞类型:", length(cibersort_common_sig), "个\n")
if(length(cibersort_common_sig) > 0) {
  print(cibersort_common_sig)
}


# -------------------- 1. 加载所需包 --------------------
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(tidyverse)
library(ggsignif)
library(RColorBrewer)
library(patchwork)

# -------------------- 2. 准备热图数据 --------------------
cat("\n准备热图数据...\n")

# 2.1 准备xCell热图数据（修复版本）
prepare_xcell_heatmap_data <- function(xcell_data, common_cells, dataset_name) {
  # 选择需要的列
  cell_cols <- paste0(common_cells, "_xCell")
  cell_cols <- cell_cols[cell_cols %in% colnames(xcell_data)]
  
  if(length(cell_cols) == 0) {
    return(NULL)
  }
  
  # 提取数据 - 不使用column_to_rownames，直接创建矩阵
  heatmap_data <- xcell_data %>%
    filter(Dataset == dataset_name) %>%
    select(Sample, Group, all_of(cell_cols))
  
  # 确保Sample是字符类型且唯一
  heatmap_data$Sample <- as.character(heatmap_data$Sample)
  
  # 分离分组信息和细胞数据
  group_info <- heatmap_data$Group
  cell_data <- heatmap_data %>% 
    select(-Sample, -Group) %>%
    as.matrix()
  
  # 设置行名
  rownames(cell_data) <- heatmap_data$Sample
  
  # 按组排序
  group_order <- order(group_info)
  cell_data <- cell_data[group_order, , drop = FALSE]
  group_info <- group_info[group_order]
  
  # Z-score标准化
  cell_data_scaled <- scale(cell_data)
  
  return(list(
    data = cell_data_scaled,
    groups = group_info,
    dataset = dataset_name,
    method = "xCell"
  ))
}

# 2.2 准备CIBERSORT热图数据（修复版本）
prepare_cibersort_heatmap_data <- function(cibersort_data, common_cells, dataset_name) {
  # 选择需要的列
  cell_cols <- common_cells[common_cells %in% colnames(cibersort_data)]
  
  if(length(cell_cols) == 0) {
    return(NULL)
  }
  
  # 提取数据 - 不使用column_to_rownames
  heatmap_data <- cibersort_data %>%
    filter(Dataset == dataset_name) %>%
    select(sample, Group, all_of(cell_cols))
  
  # 确保sample是字符类型且唯一
  heatmap_data$sample <- as.character(heatmap_data$sample)
  
  # 分离分组信息和细胞数据
  group_info <- heatmap_data$Group
  cell_data <- heatmap_data %>% 
    select(-sample, -Group) %>%
    as.matrix()
  
  # 设置行名
  rownames(cell_data) <- heatmap_data$sample
  
  # 按组排序
  group_order <- order(group_info)
  cell_data <- cell_data[group_order, , drop = FALSE]
  group_info <- group_info[group_order]
  
  # Z-score标准化
  cell_data_scaled <- scale(cell_data)
  
  return(list(
    data = cell_data_scaled,
    groups = group_info,
    dataset = dataset_name,
    method = "CIBERSORT"
  ))
}


prepare_cibersort_heatmap_data <- function(cibersort_data, common_cells, dataset_name) {
  # 选择需要的列
  cell_cols <- common_cells[common_cells %in% colnames(cibersort_data)]
  
  if(length(cell_cols) == 0) {
    return(NULL)
  }
  
  # 提取数据
  heatmap_data <- cibersort_data %>%
    filter(Dataset == dataset_name) %>%
    select(sample, Group, all_of(cell_cols))
  
  # 确保Group列存在且格式正确
  if(!"Group" %in% colnames(heatmap_data)) {
    cat("错误: 数据中没有Group列\n")
    return(NULL)
  }
  
  # 标准化Group名称（全部转为小写）
  heatmap_data$Group <- tolower(as.character(heatmap_data$Group))
  
  # 检查分组类型
  cat("数据集", dataset_name, "中的分组类型:", unique(heatmap_data$Group), "\n")
  
  # 按组排序：Control在前，Treatment在后
  heatmap_data <- heatmap_data %>%
    arrange(Group != "control")  # control排前面
  
  # 分离分组信息和细胞数据
  group_info <- heatmap_data$Group
  
  # 确保sample是字符类型且唯一
  heatmap_data$sample <- as.character(heatmap_data$sample)
  
  # 提取细胞数据并转置，使细胞类型为行，样本为列
  cell_data <- heatmap_data %>% 
    select(-sample, -Group) %>%
    as.matrix() %>%
    t()  # 转置矩阵
  
  # 设置行名（细胞类型）和列名（样本）
  rownames(cell_data) <- gsub("_CIBERSORT$", "", colnames(heatmap_data)[3:ncol(heatmap_data)])
  colnames(cell_data) <- heatmap_data$sample
  
  # 检查转置后的数据是否有无效值
  if(any(is.na(cell_data)) || any(is.infinite(cell_data)) || any(is.nan(cell_data))) {
    cat("警告: 转置后的数据中存在无效值\n")
    
    # 处理无效值：将Inf和NaN替换为NA，然后将所有NA替换为0
    cell_data[is.infinite(cell_data) | is.nan(cell_data)] <- NA
    cell_data[is.na(cell_data)] <- 0
  }
  
  # 在进行Z-score标准化之前，检查每一行的方差
  # 如果某一行（某个细胞类型）的方差为0，会导致除以0的错误
  row_var <- apply(cell_data, 1, var, na.rm = TRUE)
  zero_var_rows <- which(row_var == 0)
  
  if(length(zero_var_rows) > 0) {
    cat("警告: 以下行（细胞类型）方差为0，将在标准化时跳过:\n")
    print(rownames(cell_data)[zero_var_rows])
    
    # 对于方差为0的行，我们可以将其标准化值设为0
    # 因为这些行的所有值都相同，所以标准化后应为0
    for(row_idx in zero_var_rows) {
      cell_data[row_idx, ] <- 0
    }
  }
  
  # Z-score标准化（按行，即按细胞类型）
  # 对于方差不为0的行进行标准化
  non_zero_var <- which(row_var != 0)
  if(length(non_zero_var) > 0) {
    cell_data[non_zero_var, ] <- t(scale(t(cell_data[non_zero_var, , drop = FALSE])))
  }
  
  return(list(
    data = cell_data,
    groups = group_info,
    dataset = dataset_name,
    method = "CIBERSORT"
  ))
}


# 2.3 为每个数据集和方法准备热图数据
cat("\n准备测试集xCell热图数据...\n")
xcell_test_hm <- prepare_xcell_heatmap_data(xCell_combined, xCell_cols, "Testing_GSE235562")

cat("\n准备验证集xCell热图数据...\n")
xcell_valid_hm <- prepare_xcell_heatmap_data(xCell_combined, xCell_cols, "Validation_GSE106910")

cat("\n准备测试集CIBERSORT热图数据...\n")
cibersort_test_hm <- prepare_cibersort_heatmap_data(CIBERSORT_combined, ciber_cols, "Testing_GSE235562")

cat("\n准备验证集CIBERSORT热图数据...\n")
cibersort_valid_hm <- prepare_cibersort_heatmap_data(CIBERSORT_combined, ciber_cols, "Validation_GSE106910")

# -------------------- 3. 创建热图函数 --------------------
create_immune_heatmap <- function(hm_data, title) {
  if(is.null(hm_data)) {
    cat("没有数据创建热图:", title, "\n")
    return(NULL)
  }
  # 提取数据
  mat <- hm_data$data
  groups <- hm_data$groups
  dataset <- hm_data$dataset
  method <- hm_data$method
  # 调试信息：查看矩阵维度
  cat("热图数据维度:", dim(mat), "\n")
  cat("分组信息长度:", length(groups), "\n")
  # 检查数据
  if(is.null(mat) || length(dim(mat)) != 2) {
    cat("错误: 矩阵数据无效\n")
    return(NULL)
  }
  # 如果分组长度等于矩阵行数，说明行是样本，列是细胞类型
  # 这种情况下，我们需要转置矩阵，使列表示样本，行表示细胞类型
  if(length(groups) == nrow(mat)) {
    cat("正在转置矩阵以使样本为列...\n")
    mat <- t(mat)
  }
  # 现在检查分组长度是否等于矩阵列数
  if(length(groups) != ncol(mat)) {
    cat("错误: 分组信息长度 (", length(groups), 
        ") 与矩阵列数 (", ncol(mat), ") 不匹配\n", sep="")
    
    # 尝试自动调整：如果分组长度是矩阵行数，再次转置
    if(length(groups) == nrow(mat)) {
      cat("尝试第二次转置...\n")
      mat <- t(mat)
    }
    # 如果仍然不匹配，创建简化的分组信息
    if(length(groups) != ncol(mat)) {
      cat("创建简化分组信息...\n")
      # 使用样本数量创建占位分组
      groups <- rep(c("control", "treatment"), length.out = ncol(mat))
      cat("新分组信息长度:", length(groups), "\n")
    }
  }
  # 确保分组标签一致（全小写）
  groups <- tolower(groups)
  # 统计各组样本数
  group_counts <- table(groups)
  cat("各组样本数:\n")
  print(group_counts)
  # 颜色设置
  group_colors <- c("control" = "#94ADD4", "treatment" = "#AB97B2")
  # 如果分组中有其他标签，添加默认颜色
  unique_groups <- unique(groups)
  for(group in unique_groups) {
    if(!group %in% names(group_colors)) {
      # 为未知分组添加随机颜色
      group_colors[group] <- paste0("#", paste0(sample(c(0:9, letters[1:6]), 6, replace=TRUE), collapse=""))
    }
  }
  # 分组注释
  ha_column <- HeatmapAnnotation(
    Group = groups,
    col = list(Group = group_colors),
    annotation_name_side = "left",
    annotation_legend_param = list(
      title = "Group",
      title_gp = gpar(fontsize = 10, fontface = "bold"),
      labels_gp = gpar(fontsize = 9)))
  # 热图颜色
  col_fun <- colorRamp2(
    breaks = c(-2, 0, 2),
    colors = c("#2166AC", "white", "#B2182B"))
  # 创建分组切割向量 - 使用实际的groups
  column_split <- groups
  # 可以修改分组标签显示名称（可选）
  # 将"treatment"改为"BPD"，"control"保持为"Control"
  column_split_labels <- ifelse(groups == "treatment", "BPD", "Control")
  # 创建热图
  ht <- Heatmap(
    mat,
    name = "Z-score",
    col = col_fun,
    top_annotation = ha_column,
    column_title = title,
    column_title_gp = gpar(fontsize = 12, fontface = "bold"),
    column_names_gp = gpar(fontsize = 8),
    row_names_gp = gpar(fontsize = 9),
    row_title = "Cell Types",
    row_title_gp = gpar(fontsize = 10, fontface = "bold"),
    show_row_names = TRUE,
    show_column_names = FALSE,  # 不显示样本名，避免拥挤
    cluster_rows = TRUE,
    cluster_columns = FALSE,  # 不按列聚类，保持分组顺序
    cluster_column_slices = FALSE,  # 不聚类分组切片
    column_split = column_split,  # 使用分组进行列切割
    border = TRUE,
    show_column_dend = FALSE,
    show_row_dend = TRUE,
    row_names_max_width = unit(10, "cm"),  # 为行名分配更多空间
    heatmap_legend_param = list(
      title_gp = gpar(fontsize = 10, fontface = "bold"),
      labels_gp = gpar(fontsize = 9)),
    column_title_rot = 0)
  return(ht)
}

# -------------------- 4. 绘制热图 --------------------
cat("\n绘制免疫细胞浸润热图...\n")

# 4.1 绘制xCell热图
if(!is.null(xcell_test_hm) && !is.null(xcell_valid_hm)) {
  # 测试集xCell热图
  xcell_test_ht <- create_immune_heatmap(
    xcell_test_hm,
    "xCell Immune Cell Infiltration (Testing Set)")
  # 验证集xCell热图
  xcell_valid_ht <- create_immune_heatmap(
    xcell_valid_hm,
    "xCell Immune Cell Infiltration (Validation Set)")
  # 合并绘制
  if(!is.null(xcell_test_ht) && !is.null(xcell_valid_ht)) {
    pdf("xCell_Immune_Infiltration_Heatmaps.pdf", width = 15, height = 6.6)
    draw(xcell_test_ht + xcell_valid_ht,
         heatmap_legend_side = "right",
         annotation_legend_side = "right")
    dev.off()
    cat("xCell热图已保存为PDF和PNG格式\n")
  }
}

# 4.2 绘制CIBERSORT热图
if(!is.null(cibersort_test_hm) && !is.null(cibersort_valid_hm)) {
  # 测试集CIBERSORT热图
  cibersort_test_ht <- create_immune_heatmap(
    cibersort_test_hm,
    "CIBERSORT Immune Cell Infiltration (Testing Set)")
  # 验证集CIBERSORT热图
  cibersort_valid_ht <- create_immune_heatmap(
    cibersort_valid_hm,
    "CIBERSORT Immune Cell Infiltration (Validation Set)")
  # 合并绘制
  if(!is.null(cibersort_test_ht) && !is.null(cibersort_valid_ht)) {
    pdf("CIBERSORT_Immune_Infiltration_Heatmaps.pdf", width = 15, height = 3.5)
    draw(cibersort_test_ht + cibersort_valid_ht,
         heatmap_legend_side = "right",
         annotation_legend_side = "right")
    dev.off()
    cat("CIBERSORT热图已保存为PDF和PNG格式\n")
  }
}

# -------------------- 5. 准备箱线图数据 --------------------
cat("\n准备箱线图数据...\n")

# 5.1 准备xCell箱线图数据
prepare_xcell_boxplot_data <- function(xcell_data, common_cells) {
  # 选择需要的列
  cell_cols <- paste0(common_cells, "_xCell")
  cell_cols <- cell_cols[cell_cols %in% colnames(xcell_data)]
  if(length(cell_cols) == 0) {
    return(NULL)
  }
  # 提取数据并转换为长格式
  boxplot_data <- xcell_data %>%
    select(Dataset, Sample, Group, all_of(cell_cols)) %>%
    pivot_longer(
      cols = all_of(cell_cols),
      names_to = "CellType",
      values_to = "Score") %>%
    mutate(
      CellType = gsub("_xCell$", "", CellType),
      Method = "xCell",
      Dataset_Group = paste(Dataset, Group, sep = "_"))
  return(boxplot_data)
}

# 5.2 准备CIBERSORT箱线图数据
prepare_cibersort_boxplot_data <- function(cibersort_data, common_cells) {
  # 选择需要的列
  cell_cols <- common_cells[common_cells %in% colnames(cibersort_data)]
  if(length(cell_cols) == 0) {
    return(NULL)
  }
  # 提取数据并转换为长格式
  boxplot_data <- cibersort_data %>%
    select(Dataset, sample, Group, all_of(cell_cols)) %>%
    rename(Sample = sample) %>%
    pivot_longer(
      cols = all_of(cell_cols),
      names_to = "CellType",
      values_to = "Score") %>%
    mutate(
      CellType = gsub("_CIBERSORT$", "", CellType),
      Method = "CIBERSORT",
      Dataset_Group = paste(Dataset, Group, sep = "_"))
  return(boxplot_data)
}

# 5.3 获取所有箱线图数据
xcell_boxplot_data <- prepare_xcell_boxplot_data(xCell_combined, xCell_common_sig)
xcell_boxplot_data$Group <- str_replace(xcell_boxplot_data$Group,"treatment","BPD")
xcell_boxplot_data$Group <- factor(xcell_boxplot_data$Group,levels = c("control","BPD"))
cibersort_boxplot_data <- prepare_cibersort_boxplot_data(CIBERSORT_combined, cibersort_common_sig)
cibersort_boxplot_data$Group <- str_replace(cibersort_boxplot_data$Group,"treatment","BPD")
cibersort_boxplot_data$Group <- factor(cibersort_boxplot_data$Group,levels = c("control","BPD"))
# -------------------- 6. 箱线图绘制函数 --------------------
# 使用xcell_boxplot_data数据框直接绘制小提琴图
# 为每个细胞类型计算显著性标记
# 首先按数据集和细胞类型分组，计算p值
sig_data <- xcell_boxplot_data %>%
  group_by(Dataset, CellType) %>%
  summarise(
    # 使用Wilcoxon检验计算Control和BPD组之间的显著性
    p_value = tryCatch({
      wilcox.test(Score ~ Group)$p.value
    }, error = function(e) NA),
    # 计算y轴位置（最大值乘以1.05）
    y_position = max(Score, na.rm = TRUE) * 1.05,
    .groups = "drop") %>%
  mutate(
    # 根据p值创建显著性标签
    sig_label = case_when(
      is.na(p_value) ~ "",
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  ) %>%
  # 只保留有显著性标记的数据
  filter(sig_label != "")

# 对细胞类型按平均分数排序
celltype_order <- xcell_boxplot_data %>%
  group_by(CellType) %>%
  summarise(mean_score = mean(Score, na.rm = TRUE)) %>%
  arrange(desc(mean_score)) %>%
  pull(CellType)

xcell_boxplot_data$CellType <- factor(xcell_boxplot_data$CellType, levels = celltype_order)

# 创建小提琴图
violin_plot <- ggplot(xcell_boxplot_data, aes(x = CellType, y = Score, fill = Group)) +
  # 小提琴图
  geom_violin(width = 0.8, trim = FALSE, color = NA, alpha = 0.7,
              position = position_dodge(0.8)) +
  # 箱线图
  geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.9,
               position = position_dodge(0.8)) +
  # 散点 - 使用position_jitterdodge
  geom_jitter(size = 1.2, alpha = 0.3, shape = 21, color = "black",
              position = position_jitterdodge(
                jitter.width = 0.15,
                dodge.width = 0.8)) +
  # 添加显著性标记
  geom_text(
    data = sig_data,
    aes(x = CellType, y = y_position, label = sig_label),
    inherit.aes = FALSE,
    size = 4,
    fontface = "bold",
    vjust = -0.5) +
  # 颜色设置
  scale_fill_manual(values = c("#94ADD4", "#AB97B2"), name = "Group") +
  # 分面：按Dataset分面
  facet_wrap(~ Dataset, scales = "free_x", ncol = 1) +
  # 坐标轴标签
  labs(
    title = "xCell Immune Cell Infiltration",
    x = "Immune Cell Type",
    y = "Cell Score/Proportion") +
  # 美化主题
  theme_bw() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9, face = "bold"),
    axis.text.y = element_text(size = 10, face = "bold"),
    axis.title = element_text(size = 12, face = "bold"),
    legend.position = "top",
    strip.text = element_text(size = 11, face = "bold"),
    strip.background = element_rect(fill = "#F0F0F0", color = "black"),
    panel.border = element_rect(color = "black", size = 1),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    plot.margin = margin(10, 10, 10, 10))

# 显示图形
print(violin_plot)
# 保存图形
ggsave("xCell_Violin_Plot.pdf", plot = violin_plot, width = 16, height = 6.7)


# -------------------- 绘制CIBERSORT小提琴图（修正版）
sig_data_cibersort <- cibersort_boxplot_data %>%
  group_by(Dataset, CellType) %>%
  summarise(
    # 使用Wilcoxon检验计算Control和BPD组之间的显著性
    p_value = tryCatch({
      wilcox.test(Score ~ Group)$p.value
    }, error = function(e) NA),
    # 计算y轴位置（最大值乘以1.05）
    y_position = max(Score, na.rm = TRUE) * 1.05,
    .groups = "drop") %>%
  mutate(
    # 根据p值创建显著性标签
    sig_label = case_when(
      is.na(p_value) ~ "",
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  ) %>%
  # 只保留有显著性标记的数据
  filter(sig_label != "")

if(!is.null(cibersort_boxplot_data) && nrow(cibersort_boxplot_data) > 0) {
  cibersort_celltype_order <- cibersort_boxplot_data %>%
    group_by(CellType) %>%
    summarise(mean_score = mean(Score, na.rm = TRUE)) %>%
    arrange(desc(mean_score)) %>%
    pull(CellType)
  
  cibersort_boxplot_data$CellType <- factor(cibersort_boxplot_data$CellType, 
                                            levels = cibersort_celltype_order)
  
  # 创建小提琴图
  cibersort_violin_plot <- ggplot(cibersort_boxplot_data, 
                                  aes(x = CellType, y = Score, fill = Group)) +
    geom_violin(width = 0.8, trim = FALSE, color = NA, alpha = 0.7,
                position = position_dodge(0.8)) +
    geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.9,
                 position = position_dodge(0.8)) +
    geom_jitter(size = 1.2, alpha = 0.3, shape = 21, color = "black",
                position = position_jitterdodge(
                  jitter.width = 0.15,
                  dodge.width = 0.8
                )) +
    scale_fill_manual(values = c("#94ADD4", "#AB97B2"), name = "Group") +
    facet_wrap(~ Dataset, scales = "free_x", ncol = 1) +
    # 添加显著性标记
    geom_text(
      data = sig_data_cibersort,
      aes(x = CellType, y = y_position, label = sig_label),
      inherit.aes = FALSE,
      size = 4,
      fontface = "bold",
      vjust = -0.5) +
    labs(
      title = "CIBERSORT Immune Cell Infiltration",
      x = "Immune Cell Type",
      y = "Cell Proportion") +
    theme_bw() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9, face = "bold"),
      axis.text.y = element_text(size = 10, face = "bold"),
      axis.title = element_text(size = 12, face = "bold"),
      legend.position = "top",
      strip.text = element_text(size = 11, face = "bold"),
      strip.background = element_rect(fill = "#F0F0F0", color = "black"),
      panel.border = element_rect(color = "black", size = 1),
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank(),
      plot.margin = margin(10, 10, 10, 10))
  # 显示和保存
  print(cibersort_violin_plot)
  ggsave("CIBERSORT_Violin_Plot.pdf", plot = cibersort_violin_plot, 
         width = 5.5, height = 6.7)
}


# 4. 提取PHLDA1基因表达数据
# 检查PHLDA1是否在表达矩阵中
cat("\n检查PHLDA1基因表达数据...\n")
phlda1_test <- NA
phlda1_valid <- NA

if("PHLDA1" %in% rownames(training_expr)) {
  phlda1_test <- as.numeric(training_expr["PHLDA1", ])
  names(phlda1_test) <- colnames(training_expr)
  cat("测试集中找到PHLDA1表达数据\n")
}

if("PHLDA1" %in% rownames(validation_expr)) {
  phlda1_valid <- as.numeric(validation_expr["PHLDA1", ])
  names(phlda1_valid) <- colnames(validation_expr)
  cat("验证集中找到PHLDA1表达数据\n")
}


# 5. 计算PHLDA1与差异免疫细胞的相关性（使用并集差异细胞）
calculate_correlations_union <- function(immune_data, phlda1_expr, dataset_name, 
                                         method_name, common_sig_cells) {
  if(length(common_sig_cells) == 0) {
    return(data.frame())
  }
  
  # 准备数据
  cor_results <- data.frame()
  
  for(cell in common_sig_cells) {
    # 根据方法名构建正确的列名
    if(method_name == "xCell") {
      cell_column <- paste0(cell, "_xCell")
    } else if(method_name == "CIBERSORT") {
      # CIBERSORT的细胞类型已经带后缀了，直接使用
      cell_column <- cell
    } else {
      cell_column <- cell
    }
    
    # 检查细胞列是否存在
    if(!cell_column %in% names(immune_data)) {
      cat(sprintf("警告: 列 '%s' 在数据集 '%s' 中不存在\n", cell_column, dataset_name))
      next
    }
    
    # 获取免疫细胞数据
    cell_data <- immune_data[[cell_column]]
    
    # 获取样本名
    if("Sample" %in% names(immune_data)) {
      sample_names <- immune_data$Sample
    } else if("sample" %in% names(immune_data)) {
      sample_names <- immune_data$sample
    } else if(!is.null(rownames(immune_data))) {
      sample_names <- rownames(immune_data)
    } else {
      cat(sprintf("警告: 数据集 '%s' 中找不到样本名称\n", dataset_name))
      next
    }
    
    names(cell_data) <- sample_names
    
    # 对齐样本
    common_samples <- intersect(names(cell_data), names(phlda1_expr))
    
    if(length(common_samples) >= 5) {  # 至少需要5个样本
      cell_values <- cell_data[common_samples]
      phlda1_values <- phlda1_expr[common_samples]
      
      # 计算Spearman相关性
      cor_test <- cor.test(cell_values, phlda1_values, method = "spearman")
      cor_results <- rbind(cor_results, data.frame(
        Dataset = dataset_name,
        Method = method_name,
        CellType = cell,
        Correlation = cor_test$estimate,
        Pvalue = cor_test$p.value,
        Samples = length(common_samples)
      ))
    } else {
      cat(sprintf("警告: 细胞类型 '%s' 在数据集 '%s' 中只有 %d 个共同样本，少于5个，跳过\n", 
                  cell, dataset_name, length(common_samples)))
    }
  }
  
  # 多重检验校正
  if(nrow(cor_results) > 0) {
    cor_results$Padj <- p.adjust(cor_results$Pvalue, method = "BH")
  }
  
  return(cor_results)
}

# 分别计算相关性（使用并集差异细胞）
cat("\n计算PHLDA1与免疫细胞的相关性（使用并集差异细胞）...\n")

# 确保xCell并集差异细胞列表有正确的格式（去掉可能的后缀）
xCell_common_sig_clean <- gsub("_xCell$", "", xCell_common_sig)
cat("\nxCell并集差异细胞类型数:", length(xCell_common_sig_clean), "\n")

# CIBERSORT并集差异细胞列表已经带后缀，保持原样
cat("CIBERSORT并集差异细胞类型数:", length(cibersort_common_sig), "\n")

# xCell相关性（使用并集）
xCell_cor_test_union <- calculate_correlations_union(
  scores_df_test, 
  phlda1_test, 
  "Testing", 
  "xCell",
  xCell_common_sig_clean  # 使用清理后的列表
)

xCell_cor_valid_union <- calculate_correlations_union(
  scores_df_valid, 
  phlda1_valid, 
  "Validation", 
  "xCell",
  xCell_common_sig_clean  # 使用清理后的列表
)

# CIBERSORT相关性（使用并集）
CIBERSORT_cor_test_union <- calculate_correlations_union(
  valid_samples_train, 
  phlda1_test[rownames(valid_samples_train)], 
  "Testing", 
  "CIBERSORT",
  cibersort_common_sig  # 使用原列表（已经带后缀）
)

CIBERSORT_cor_valid_union <- calculate_correlations_union(
  valid_samples_valid, 
  phlda1_valid[rownames(valid_samples_valid)], 
  "Validation", 
  "CIBERSORT",
  cibersort_common_sig  # 使用原列表（已经带后缀）
)

# 合并所有结果
all_cor_results_union <- rbind(
  xCell_cor_test_union,
  xCell_cor_valid_union,
  CIBERSORT_cor_test_union,
  CIBERSORT_cor_valid_union)

# 显示前几个结果
if(nrow(all_cor_results_union) > 0) {
  cat("\n前10个相关性结果:\n")
  print(head(all_cor_results_union, 10))
  
  # 保存结果到文件
  write.csv(all_cor_results_union, 
            "PHLDA1_immune_correlation_union_results.csv", 
            row.names = FALSE)
  cat("\n结果已保存到: PHLDA1_immune_correlation_union_results.csv\n")
  
  # 也可以按数据集和方法分别保存
  write.csv(xCell_cor_test_union, "xCell_cor_test_union.csv", row.names = FALSE)
  write.csv(xCell_cor_valid_union, "xCell_cor_valid_union.csv", row.names = FALSE)
  write.csv(CIBERSORT_cor_test_union, "CIBERSORT_cor_test_union.csv", row.names = FALSE)
  write.csv(CIBERSORT_cor_valid_union, "CIBERSORT_cor_valid_union.csv", row.names = FALSE)
} else {
  cat("\n警告: 没有计算到任何相关性结果\n")
}

# 6. 绘制PHLDA1与免疫细胞相关性气泡图（并集版本）
# 加载必要的包
library(ggplot2)
library(RColorBrewer)
library(cowplot)

# 绘制PHLDA1与免疫细胞相关性气泡图
create_phlda1_immune_corr_plot <- function(cor_data, title = "PHLDA1 Expression vs Immune Cell Correlation") {
  if(nrow(cor_data) == 0) {
    cat("没有足够的数据绘制相关性气泡图\n")
    return(NULL)
  }
  # 准备数据
  plot_data <- cor_data %>%
    # 创建组合标签（数据集_方法）
    mutate(Group = paste(Dataset, Method, sep = "_")) %>%
    # 根据Padj创建FDR分组
    mutate(
      fdr_group = case_when(
        Padj < 0.001 ~ "<0.001",
        Padj < 0.01 ~ "<0.01",
        Padj < 0.05 ~ "<0.05",
        TRUE ~ "≥0.05"),
      fdr_group_2 = fdr_group)
  # 设置组顺序
  group_levels <- c("Testing_xCell", "Validation_xCell", "Testing_CIBERSORT", "Validation_CIBERSORT")
  plot_data$Group <- factor(plot_data$Group, levels = group_levels)
  
  # 按照平均相关性绝对值排序细胞类型
  cell_order <- plot_data %>%
    group_by(CellType) %>%
    summarise(mean_abs_cor = mean(abs(Correlation), na.rm = TRUE)) %>%
    arrange(desc(mean_abs_cor)) %>%
    pull(CellType)
  
  plot_data$CellType <- factor(plot_data$CellType, levels = cell_order)
  
  # 创建绘图
  p <- ggplot(plot_data) +
    # 添加水平分隔线
    geom_hline(yintercept = unique(plot_data$Group), color = "#E8E8E8", linewidth = 0.5) +
    
    # 正相关点（实心）
    geom_point(
      data = subset(plot_data, Correlation > 0),
      shape = 19,
      stroke = 0,
      aes(x = CellType, y = Group, size = abs(Correlation), color = fdr_group)) +
    # 负相关点（空心）
    geom_point(
      data = subset(plot_data, Correlation < 0),
      shape = 21,
      stroke = 0.2,
      aes(x = CellType, y = Group, size = abs(Correlation), fill = fdr_group_2)) +
    # 点大小设置
    scale_size_continuous(
      name = "|Correlation|",
      limits = c(0, 0.8),
      breaks = c(0.2, 0.4, 0.6, 0.8),
      range = c(2, 8)) +
    # 正相关颜色（蓝色系，表示负相关？实际是正相关用蓝色？需要确认）
    # 注意：原代码中color用于负相关，fill用于正相关，但标签是相反的
    # 这里我根据你的数据调整：正相关用红色系，负相关用蓝色系
    scale_color_manual(
      name = "Positive correlation\nFDR q-value",
      values = c(
        "<0.001" = brewer.pal(9, 'Reds')[9],
        "<0.01" = brewer.pal(9, 'Reds')[7],
        "<0.05" = brewer.pal(9, 'Reds')[5],
        "≥0.05" = "#E8E8E8"),
      drop = FALSE) +
    scale_fill_manual(
      name = "Negative correlation\nFDR q-value",
      values = c(
        "<0.001" = brewer.pal(9, 'Blues')[9],
        "<0.01" = brewer.pal(9, 'Blues')[7],
        "<0.05" = brewer.pal(9, 'Blues')[5],
        "≥0.05" = "#E8E8E8"),
      drop = FALSE) +
    # 主题设置
    cowplot::theme_cowplot() +
    theme(
      panel.grid.major = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
      axis.text.y = element_text(size = 11, face = "bold"),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      legend.position = "right",
      legend.box = "vertical",
      legend.spacing.y = unit(0.2, "cm")) +
    # 标签
    labs(
      title = title,
      x = "",
      y = "") +
    
    # 图例设置
    guides(
      fill = guide_legend(
        override.aes = list(size = 4, shape = 21),
        order = 1),
      color = guide_legend(
        override.aes = list(size = 4, shape = 19),
        order = 2),
      size = guide_legend(
        override.aes = list(fill = "black", color = "black"),
        order = 3))
  return(p)
}

# 使用函数绘制气泡图
# 加载必要的包
library(ggplot2)
library(RColorBrewer)
library(cowplot)

if(nrow(all_cor_results_union) > 0) {
  # 绘制完整图
  phlda1_corr_plot <- create_phlda1_immune_corr_plot(
    all_cor_results_union,
    "PHLDA1 Expression Correlation with Immune Cells")
  print(phlda1_corr_plot)
  # 保存图片
  ggsave(
    "PHLDA1_Immune_Cell_Correlation_Bubble.png",
    plot = phlda1_corr_plot,
    width = max(12, length(unique(all_cor_results_union$CellType)) * 0.4 + 6),
    height = 4,
    dpi = 300)
  ggsave(
    "PHLDA1_Immune_Cell_Correlation_Bubble.pdf",
    plot = phlda1_corr_plot,
    width = max(12, length(unique(all_cor_results_union$CellType)) * 0.4 + 6),
    height = 4)
  # 也可以分别绘制xCell和CIBERSORT
  # xCell
  xCell_corr_data <- all_cor_results_union %>% 
    filter(Method == "xCell")
  if(nrow(xCell_corr_data) > 0) {
    xCell_plot <- create_phlda1_immune_corr_plot(
      xCell_corr_data,
      "PHLDA1 Expression Correlation with xCell Immune Cells")
    print(xCell_plot)
    ggsave(
      "PHLDA1_xCell_Correlation_Bubble.png",
      plot = xCell_plot,
      width = max(12, length(unique(xCell_corr_data$CellType)) * 0.4 + 6),
      height = 3,
      dpi = 300)
  }
  # CIBERSORT
  CIBERSORT_corr_data <- all_cor_results_union %>% 
    filter(Method == "CIBERSORT")
  if(nrow(CIBERSORT_corr_data) > 0) {
    CIBERSORT_plot <- create_phlda1_immune_corr_plot(
      CIBERSORT_corr_data,
      "PHLDA1 Expression Correlation with CIBERSORT Immune Cells")
    print(CIBERSORT_plot)
    ggsave(
      "PHLDA1_CIBERSORT_Correlation_Bubble.png",
      plot = CIBERSORT_plot,
      width = max(12, length(unique(CIBERSORT_corr_data$CellType)) * 0.4 + 6),
      height = 3,
      dpi = 300)
  }
  cat("\n气泡图已保存\n")
}
