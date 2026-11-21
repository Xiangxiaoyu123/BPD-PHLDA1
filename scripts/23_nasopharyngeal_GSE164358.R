###############################################################################
# PHLDA1 and candidate gene expression in nasopharyngeal basal cells
#
# Source      : r.03_Cell_prop.r (lines 9803-10177)
# Input       : GSE164358
# Output      : Fig. 6I
#
# Code below is unchanged from the original analysis script.
###############################################################################
# -----------------------------------------------------------------------------
# PATHS: this file keeps the absolute paths of the analysis server on which the
# study was run (/data/nas1/xiangxiaoyu_OD/1.BPD+MSCs/1.scRNA/... and
# /home/xxyloveme/study/BPD+MSCs/...). Set PROJ_DIR in config/paths.R and
# replace those prefixes with your own project directory before running.
# -----------------------------------------------------------------------------


############bulkRNA验证 GSE164358############
setwd("/data/nas1/xiangxiaoyu_OD/1.BPD+MSCs/1.scRNA/10_bulkRNA/3.GSE164358/")

# 导入所需的包
library(ggplot2)
library(ggsignif)
library(ggpubr)
library(RColorBrewer)
library(tidyr)
library(dplyr)
library(readr)

# 读取数据
# 注意：数据前两行是复杂的表头结构
data <- read.delim("GSE164358_tpm_Genelist.txt", header = FALSE, stringsAsFactors = FALSE)

# 处理表头
# 提取第一行（Day0, Day21信息）
header_row1 <- data[1, ]
# 提取第二行（preT, Term信息）
header_row2 <- data[2, ]
# 提取第三行（样本名称）
header_row3 <- data[3, ]

# 基因表达数据从第4行开始
expression_data <- data[4:nrow(data), ]
colnames(expression_data) <- header_row3

# 设置第一列为基因名
colnames(expression_data)[1] <- "Gene"
expression_data <- expression_data[!is.na(expression_data$Gene) & expression_data$Gene != "", ]
expression_data <- expression_data[,-c(6,11,16)]
# 定义感兴趣的基因
genes_of_interest <- c("AES", "ACTG1", "TSC22D2", "HSPA1A", "PHLDA1", "SLC26A7", "NFAT5", "TLE5", "MT1A")

# 检查哪些基因在数据中可用
available_genes <- genes_of_interest[genes_of_interest %in% expression_data$Gene]
missing_genes <- setdiff(genes_of_interest, expression_data$Gene)

cat("可用的基因:", paste(available_genes, collapse = ", "), "\n")
if(length(missing_genes) > 0) {
  cat("缺失的基因:", paste(missing_genes, collapse = ", "), "\n")
}

# 如果没有任何可用基因，退出
if(length(available_genes) == 0) {
  stop("没有找到任何目标基因！")
}

# 准备绘图数据
prepare_plot_data <- function(expression_data, genes) {
  # 提取目标基因的表达数据
  target_data <- expression_data[expression_data$Gene %in% genes, ]
  
  # 转置数据，使基因作为列，样本作为行
  target_data_t <- as.data.frame(t(target_data[, -1]))
  colnames(target_data_t) <- target_data$Gene
  
  # 添加样本信息
  target_data_t$Sample <- rownames(target_data_t)
  
  # 根据样本名称添加分组信息
  # 根据样本名称的前缀确定分组
  target_data_t$Group <- NA
  
  # WH001-005: Day0-preT
  target_data_t$Group[target_data_t$Sample %in% c("WH001", "WH002", "WH003", "WH005")] <- "Day0_preT"
  # WH004, WH006-008: Day0-Term
  target_data_t$Group[target_data_t$Sample %in% c("WH004", "WH006", "WH007", "WH008")] <- "Day0_Term"
  # WH009-011, WH013: Day21-preT
  target_data_t$Group[target_data_t$Sample %in% c("WH009", "WH010", "WH011", "WH013")] <- "Day21_preT"
  # WH012, WH014-016: Day21-Term
  target_data_t$Group[target_data_t$Sample %in% c("WH012", "WH014", "WH015", "WH016")] <- "Day21_Term"
  
  # 移除NA值
  target_data_t <- target_data_t[!is.na(target_data_t$Group), ]
  
  # 转换为长格式
  plot_data_long <- pivot_longer(
    target_data_t,
    cols = all_of(genes),
    names_to = "Gene",
    values_to = "Expression"
  )
  
  # 确保表达值为数值型
  plot_data_long$Expression <- as.numeric(plot_data_long$Expression)
  
  return(plot_data_long)
}

# 准备绘图数据
plot_data <- prepare_plot_data(expression_data, available_genes)

# 查看绘图数据结构
cat("绘图数据维度:", dim(plot_data), "\n")
print(head(plot_data))

# 定义组别顺序
group_levels <- c("Day0_preT", "Day0_Term", "Day21_preT", "Day21_Term")
plot_data$Group <- factor(plot_data$Group, levels = group_levels)

# 定义颜色
group_colors <- c("#85B22E", "#5F80B4", "#E29827", "#922927")

# 为每个基因创建箱线图
for(gene in available_genes) {
  # 提取当前基因的数据
  gene_data <- plot_data[plot_data$Gene == gene, ]
  
  # 计算y轴最大值，用于设置显著性标记的位置
  y_max <- max(gene_data$Expression, na.rm = TRUE)
  y_step <- y_max * 0.1  # 每个显著性标记之间的间隔
  
  # 创建箱线图
  p <- ggplot(data = gene_data) + 
    geom_boxplot(
      mapping = aes(x = Group, y = Expression, colour = Group),
      alpha = 0.5,
      size = 1.5,
      width = 0.6
    ) + 
    geom_jitter(
      mapping = aes(x = Group, y = Expression, colour = Group),
      alpha = 0.3, 
      size = 3
    ) +
    scale_color_manual(
      limits = group_levels,
      values = group_colors
    ) +
    # 添加显著性比较
    geom_signif(
      mapping = aes(x = Group, y = Expression),
      comparisons = list(
        c("Day0_preT", "Day0_Term"),
        c("Day0_preT", "Day21_preT"),
        c("Day0_preT", "Day21_Term"),
        c("Day0_Term", "Day21_preT"),
        c("Day0_Term", "Day21_Term"),
        c("Day21_preT", "Day21_Term")
      ),
      map_signif_level = T,
      tip_length = 0.01,
      y_position = c(
        y_max + y_step * 1,
        y_max + y_step * 2,
        y_max + y_step * 3,
        y_max + y_step * 4,
        y_max + y_step * 5,
        y_max + y_step * 6
      ),
      size = 1,
      textsize = 4,
      test = "wilcox.test"  # 使用Wilcoxon检验
    ) +
    theme_classic(
      base_line_size = 1
    ) +
    labs(
      title = gene,
      x = "",
      y = "TPM Expression"
    ) +
    theme(
      plot.title = element_text(
        size = 15,
        colour = "black",
        hjust = 0.5
      ),
      axis.title.y = element_text(
        size = 15,
        color = "black",
        face = "bold",
        vjust = 1.9,
        hjust = 0.5,
        angle = 90
      ),
      legend.title = element_text(
        color = "black",
        size = 15,
        face = "bold"
      ),
      legend.text = element_text(
        color = "black",
        size = 10,
        face = "bold"
      ),
      axis.text.x = element_text(
        size = 13,
        color = "black",
        face = "bold",
        vjust = 0.5,
        hjust = 0.5,
        angle = 45  # x轴标签倾斜45度
      ),
      axis.text.y = element_text(
        size = 13,
        color = "black",
        face = "bold",
        vjust = 0.5,
        hjust = 0.5,
        angle = 0
      )
    )
  
  # 显示图形
  print(p)
  
  # 保存图形
  filename <- paste0("boxplot_", gene, ".pdf")
  ggsave(filename, p, width = 7, height = 8)
}

# 创建所有基因的组合图（分面显示）
combined_plot <- ggplot(data = plot_data) + 
  geom_boxplot(
    mapping = aes(x = Group, y = Expression, colour = Group),
    alpha = 0.5,
    size = 1,
    width = 0.6
  ) + 
  geom_jitter(
    mapping = aes(x = Group, y = Expression, colour = Group),
    alpha = 0.3, 
    size = 1.5
  ) +
  scale_color_manual(
    limits = group_levels,
    values = group_colors
  ) +
  facet_wrap(~ Gene, scales = "free_y", ncol = 3) +
  theme_classic(
    base_line_size = 1
  ) +
  labs(
    title = "Gene Expression Across Groups",
    x = "",
    y = "TPM Expression"
  ) +
  theme(
    plot.title = element_text(
      size = 15,
      colour = "black",
      hjust = 0.5
    ),
    axis.title.y = element_text(
      size = 15,
      color = "black",
      face = "bold",
      vjust = 1.9,
      hjust = 0.5,
      angle = 90
    ),
    legend.title = element_text(
      color = "black",
      size = 12,
      face = "bold"
    ),
    legend.text = element_text(
      color = "black",
      size = 10,
      face = "bold"
    ),
    axis.text.x = element_text(
      size = 10,
      color = "black",
      face = "bold",
      vjust = 0.5,
      hjust = 0.5,
      angle = 45
    ),
    axis.text.y = element_text(
      size = 10,
      color = "black",
      face = "bold",
      vjust = 0.5,
      hjust = 0.5,
      angle = 0
    ),
    strip.text = element_text(
      size = 10,
      face = "bold"
    )
  )

# 显示组合图
print(combined_plot)

# 保存组合图
ggsave("combined_gene_expression_boxplot.pdf", combined_plot, width = 12, height = 9)

# 统计分析和结果输出
# 执行统计检验（Wilcoxon检验）
stat_results <- data.frame()

for(gene in available_genes) {
  cat("处理基因:", gene, "\n")
  
  # 所有组别两两比较
  groups <- unique(plot_data$Group)
  
  for(i in 1:(length(groups)-1)) {
    for(j in (i+1):length(groups)) {
      group1 <- groups[i]
      group2 <- groups[j]
      
      cat("比较:", group1, "vs", group2, "\n")
      
      # 提取两个组的数据 - 使用 [[1]] 或 $Expression 提取数值向量
      group1_data <- plot_data[plot_data$Gene == gene & plot_data$Group == group1, "Expression"][[1]]
      group2_data <- plot_data[plot_data$Gene == gene & plot_data$Group == group2, "Expression"][[1]]
      
      # 调试信息
      cat("group1_data类型:", class(group1_data), "\n")
      cat("group2_data类型:", class(group2_data), "\n")
      cat("group1_data值:", group1_data, "\n")
      cat("group2_data值:", group2_data, "\n")
      
      # 确保数据是数值向量
      group1_data <- as.numeric(group1_data)
      group2_data <- as.numeric(group2_data)
      
      # 移除NA值
      group1_data <- group1_data[!is.na(group1_data)]
      group2_data <- group2_data[!is.na(group2_data)]
      
      # 检查是否有足够的数据进行检验
      if(length(group1_data) > 1 && length(group2_data) > 1) {
        # 执行Wilcoxon检验
        test_result <- wilcox.test(group1_data, group2_data)
        
        stat_results <- rbind(stat_results, data.frame(
          Gene = gene,
          Group1 = group1,
          Group2 = group2,
          Group1_Mean = mean(group1_data, na.rm = TRUE),
          Group2_Mean = mean(group2_data, na.rm = TRUE),
          Group1_Median = median(group1_data, na.rm = TRUE),
          Group2_Median = median(group2_data, na.rm = TRUE),
          P_Value = test_result$p.value,
          Significant = test_result$p.value < 0.05
        ))
        
        cat("检验完成，p值:", test_result$p.value, "\n")
      } else {
        cat("跳过检验 - 数据不足\n")
        cat("group1_data长度:", length(group1_data), "\n")
        cat("group2_data长度:", length(group2_data), "\n")
      }
      cat("---\n")
    }
  }
}

# 添加FDR校正
if(nrow(stat_results) > 0) {
  stat_results$FDR <- p.adjust(stat_results$P_Value, method = "fdr")
  stat_results$Significant_FDR <- stat_results$FDR < 0.05
  
  # 输出统计结果
  print(stat_results)
  write.csv(stat_results, "gene_expression_statistics_TPM.csv", row.names = FALSE)
} else {
  cat("没有足够的有效数据进行统计分析\n")
}

save.image(file = "GSE164358.bulkRNA.wprkplace.RData")

getwd()
