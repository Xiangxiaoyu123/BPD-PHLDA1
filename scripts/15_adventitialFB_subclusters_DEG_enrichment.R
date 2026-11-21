###############################################################################
# Adventitial fibroblast subclustering, differential expression and GO/KEGG enrichment
#
# Source      : r.03_Cell_prop.r (lines 1281-2976)
# Input       : meso_adventitial_FB
# Output      : Fig. 3; Supplementary Fig. 6; Supplementary Tables 9-10
#
# Code below is unchanged from the original analysis script.
###############################################################################
# -----------------------------------------------------------------------------
# PATHS: this file keeps the absolute paths of the analysis server on which the
# study was run (/data/nas1/xiangxiaoyu_OD/1.BPD+MSCs/1.scRNA/... and
# /home/xxyloveme/study/BPD+MSCs/...). Set PROJ_DIR in config/paths.R and
# replace those prefixes with your own project directory before running.
# -----------------------------------------------------------------------------


##############################3.继续分析Adventitial FB#####################
setwd("/data/nas1/xiangxiaoyu_OD/1.BPD+MSCs/1.scRNA/04_Adventitial_FB/")

#############3.1 三组差异基因的火山图#############
library(ggplot2)
library(ggrepel)
library(dplyr)

# 1. 从原始数据中提取 Adventitial FB 的所有基因（包括不显著的）
adv_fb_data_1 <- volcano_data1 %>% filter(celltype == "Adventitial FB")
adv_fb_data_2 <- volcano_data2 %>% filter(celltype == "Adventitial FB") 
adv_fb_data_3 <- volcano_data3 %>% filter(celltype == "Adventitial FB")

# 2. 创建绘制火山图的函数
create_adv_fb_volcano <- function(df, title, pval_threshold = 1e-05, logfc_threshold = 1) {
  
  # 添加变化状态列 - 包含所有基因
  df$change <- ifelse(df$p_val_adj < pval_threshold & df$avg_log2FC > logfc_threshold, "Up",
                      ifelse(df$p_val_adj < pval_threshold & df$avg_log2FC < -logfc_threshold, "Down", "Not"))
  
  # 设置颜色
  df$color <- ifelse(df$change == "Up", "#CD2626",
                     ifelse(df$change == "Down", "#000080", "#BEBEBE"))
  
  # 统计显著基因数量
  up_count <- sum(df$change == "Up", na.rm = TRUE)
  down_count <- sum(df$change == "Down", na.rm = TRUE)
  total_genes <- nrow(df)
  
  # 选择要标记的基因（最显著的前10个上调和下调基因）
  top_up <- df %>% filter(change == "Up") %>% arrange(p_val_adj) %>% head(10)
  top_down <- df %>% filter(change == "Down") %>% arrange(p_val_adj) %>% head(10)
  label_genes <- bind_rows(top_up, top_down)
  
  # 创建火山图
  p <- ggplot(df, aes(x = avg_log2FC, y = -log10(p_val_adj))) + 
    geom_point(alpha = 0.6, aes(color = change, size = pct.1)) +
    scale_color_manual(values = c("Down" = "#000080", "Not" = "#BEBEBE", "Up" = "#CD2626")) +
    
    # 添加阈值线
    geom_vline(xintercept = c(-logfc_threshold, logfc_threshold), 
               linetype = "dashed", color = "black", alpha = 0.7, size = 0.5) +
    geom_hline(yintercept = -log10(pval_threshold), 
               linetype = "dashed", color = "black", alpha = 0.7, size = 0.5) +
    
    # 添加基因标签
    geom_text_repel(
      data = label_genes,
      aes(label = gene),
      size = 3,
      box.padding = 0.5,
      max.overlaps = 20,
      segment.color = "grey50",
      segment.size = 0.2
    ) +
    
    # 坐标轴和标签
    labs(x = "Log2[FC]", y = "-Log10(Padj)", 
         title = title,
         subtitle = paste0("Total genes: ", total_genes, " | Up: ", up_count, " | Down: ", down_count)) +
    
    # 坐标轴范围
    xlim(-4, 4) +
    
    # 点的大小范围
    scale_size_continuous(range = c(0.5, 4), name = "Expression\nPercentage") +
    
    # 主题
    theme_classic(base_size = 14) +
    theme(
      axis.text = element_text(size = 12, color = "black"),
      axis.title = element_text(size = 12),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 10, color = "grey40", hjust = 0.5),
      legend.position = "right",
      panel.grid.major = element_line(color = "grey90", size = 0.2),
      panel.grid.minor = element_blank()
    )
  
  return(p)
}

# 3. 绘制三个比较组的火山图
# BPD vs Term infant
p1 <- create_adv_fb_volcano(adv_fb_data_1, 
                            "Adventitial FB: BPD vs Term Infant")

# BPD+PH vs Term infant  
p2 <- create_adv_fb_volcano(adv_fb_data_2,
                            "Adventitial FB: BPD+PH vs Term Infant")

# BPD+PH vs BPD
p3 <- create_adv_fb_volcano(adv_fb_data_3,
                            "Adventitial FB: BPD+PH vs BPD")

# 4. 显示图形
print(p1)
print(p2) 
print(p3)

# 5. 保存图形
ggsave("./01_DEGs/Adventitial_FB_BPD_vs_Term_volcano.pdf", p1, 
       width = 10, height = 8, dpi = 300)
ggsave("./01_DEGs/Adventitial_FB_BPD+PH_vs_Term_volcano.pdf", p2, 
       width = 10, height = 8, dpi = 300)
ggsave("./01_DEGs/Adventitial_FB_BPD+PH_vs_BPD_volcano.pdf", p3, 
       width = 10, height = 8, dpi = 300)


# 美化版本的火山图函数
create_beautified_adv_fb_volcano <- function(df, title, pval_threshold = 1e-05, logfc_threshold = 1) {
  
  # 添加变化状态 - 包含所有基因
  df <- df %>%
    mutate(
      change = case_when(
        p_val_adj < pval_threshold & avg_log2FC > logfc_threshold ~ "Up",
        p_val_adj < pval_threshold & avg_log2FC < -logfc_threshold ~ "Down",
        TRUE ~ "Not"
      )
    )
  
  # 统计显著基因数量
  up_count <- sum(df$change == "Up", na.rm = TRUE)
  down_count <- sum(df$change == "Down", na.rm = TRUE)
  total_genes <- nrow(df)
  
  # 选择要标记的基因
  top_up <- df %>% filter(change == "Up") %>% arrange(p_val_adj) %>% head(8)
  top_down <- df %>% filter(change == "Down") %>% arrange(p_val_adj) %>% head(8)
  label_genes <- bind_rows(top_up, top_down)
  
  # 创建火山图
  p <- ggplot(df, aes(x = avg_log2FC, y = -log10(p_val_adj))) + 
    
    # 背景点（不显著的）
    geom_point(data = filter(df, change == "Not"),
               aes(size = pct.1), color = "#B3B3B3", alpha = 0.3) +
    
    # 显著点
    geom_point(data = filter(df, change != "Not"),
               aes(size = pct.1, color = change), alpha = 0.7) +
    
    # 颜色标度
    scale_color_manual(
      values = c("Down" = "#4575B4", "Up" = "#D73027"),
      labels = c("Down", "Up"),
      name = "Regulation"
    ) +
    
    # 阈值线
    geom_vline(xintercept = c(-logfc_threshold, logfc_threshold), 
               linetype = "dashed", color = "grey40", alpha = 0.8, size = 0.6) +
    geom_hline(yintercept = -log10(pval_threshold), 
               linetype = "dashed", color = "grey40", alpha = 0.8, size = 0.6) +
    
    # 基因标签
    geom_text_repel(
      data = label_genes,
      aes(label = gene, color = change),
      size = 3.2,
      box.padding = 0.6,
      point.padding = 0.3,
      max.overlaps = 25,
      segment.color = "grey40",
      segment.size = 0.3,
      min.segment.length = 0.2,
      show.legend = FALSE
    ) +
    
    # 坐标轴和标签
    labs(
      x = "Log2 Fold Change", 
      y = "-Log10 Adjusted P-value",
      title = title,
      subtitle = paste0("Total genes: ", total_genes, " | Up: ", up_count, " | Down: ", down_count)
    ) +
    
    # 坐标轴范围
    xlim(-4, 4) +
    ylim(0, max(-log10(df$p_val_adj), na.rm = TRUE) * 1.05) +
    
    # 点的大小
    scale_size_continuous(
      range = c(1, 5), 
      name = "Expression\nPercentage",
      breaks = c(0.2, 0.4, 0.6, 0.8),
      labels = c("20%", "40%", "60%", "80%")
    ) +
    
    # 主题
    theme_minimal() +
    theme(
      axis.text = element_text(size = 11, color = "black"),
      axis.title = element_text(size = 12, face = "bold"),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, color = "grey40", hjust = 0.5, 
                                   margin = margin(b = 15)),
      legend.title = element_text(face = "bold", size = 10),
      legend.text = element_text(size = 9),
      panel.grid.major = element_line(color = "grey90", size = 0.3),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "grey70", fill = NA, size = 0.5)
    )
  
  return(p)
}

# 使用美化版本绘制
p1_beauty <- create_beautified_adv_fb_volcano(adv_fb_data_1, 
                                              "Adventitial FB: BPD vs Term Infant")
p2_beauty <- create_beautified_adv_fb_volcano(adv_fb_data_2,
                                              "Adventitial FB: BPD+PH vs Term Infant")
p3_beauty <- create_beautified_adv_fb_volcano(adv_fb_data_3,
                                              "Adventitial FB: BPD+PH vs BPD")

# 显示和保存美化版本
print(p1_beauty)
print(p2_beauty)
print(p3_beauty)

ggsave("./01_DEGs/Adventitial_FB_BPD_vs_Term_volcano_beautified.pdf", p1_beauty, 
       width = 11, height = 9, dpi = 300)
ggsave("./01_DEGs/Adventitial_FB_BPD+PH_vs_Term_volcano_beautified.pdf", p2_beauty, 
       width = 11, height = 9, dpi = 300)
ggsave("./01_DEGs/Adventitial_FB_BPD+PH_vs_BPD_volcano_beautified.pdf", p3_beauty, 
       width = 11, height = 9, dpi = 300)

#############3.2 提取Adventitial FB进行后续分析#############
# 提取Adventitial FB细胞
#绘制只有Adventitial FB的umap图
# 手动创建分面图以获得完全控制
library(ggplot2)
library(dplyr)

# 提取UMAP坐标
umap_coords <- Embeddings(meso_relabel_no_preterm, reduction = "umap")
# 提取UMAP坐标和元数据
umap_data <- Embeddings(meso_relabel_no_preterm, reduction = "umap") %>% 
  as.data.frame() %>% 
  setNames(c("UMAP1", "UMAP2"))

umap_data$celltype <- meso_relabel_no_preterm$celltype
umap_data$is_adventitial <- umap_data$celltype == "Adventitial FB"
umap_data$condition <- factor(meso_relabel_no_preterm$condition,
                              levels = c("Term infant", "BPD", "BPD+PH"))


# 创建突出显示的UMAP图
custom_highlight_plot <- ggplot(umap_data, aes(x = UMAP1, y = UMAP2)) +
  
  # 首先绘制所有灰色背景细胞
  geom_point(data = subset(umap_data, !is_adventitial),
             color = "grey90", size = 0.8, alpha = 0.6) +
  
  # 然后突出显示Adventitial FB细胞
  geom_point(data = subset(umap_data, is_adventitial),
             color = "#58B3C7", size = 1.2, alpha = 0.8) +
  
  # 添加标题和说明
  labs(title = "Adventitial FB in Mesenchymal Cell Landscape",
       x = "UMAP1", y = "UMAP2",
       caption = paste("Total cells:", nrow(umap_data), 
                       "| Adventitial FB:", sum(umap_data$is_adventitial))) +
  
  # 主题设置
  theme_minimal() +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5, 
                              margin = margin(b = 10)),
    plot.subtitle = element_text(size = 12, color = "grey40", hjust = 0.5,
                                 lineheight = 1.3, margin = margin(b = 15)),
    plot.caption = element_text(size = 10, color = "grey50", hjust = 0.5,
                                margin = margin(t = 10)),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  
  # 确保坐标轴比例正确
  coord_fixed()

print(custom_highlight_plot)
# 保存图形
ggsave("./02_UMAP/Adventitial_FB_highlight_detailed.pdf", 
       custom_highlight_plot, 
       width = 12, height = 10, dpi = 300)


###每个分面中添加细胞数量统计
# 计算每个条件下的细胞数量
cell_counts <- umap_data %>%
  group_by(condition) %>%
  summarise(
    total_cells = n(),
    adventitial_cells = sum(celltype == "Adventitial FB"),
    adventitial_percent = round(adventitial_cells / total_cells * 100, 1)
  )

# 创建带统计信息的分面图
facet_with_stats <- ggplot(umap_data, aes(x = UMAP1, y = UMAP2)) +
  # 绘制所有细胞（灰色）
  geom_point(data = subset(umap_data, celltype != "Adventitial FB"),
             color = "grey90", size = 0.5, alpha = 0.5) +
  # 突出显示Adventitial FB
  geom_point(data = subset(umap_data, celltype == "Adventitial FB"),
             color = "#58B3C7", size = 0.8, alpha = 0.8) +
  # 分面
  facet_wrap(~ condition, nrow = 1) +
  # 在每个分面中添加统计信息
  geom_text(data = cell_counts, 
            aes(x = -Inf, y = Inf, 
                label = paste0("Adventitial FB: ", adventitial_percent, "%\n",
                               "(", adventitial_cells, "/", total_cells, ")")),
            hjust = -0.1, vjust = 1.5, size = 3.5, fontface = "bold") +
  # 标签和标题
  labs(title = "Adventitial FB Distribution Across Conditions",
       subtitle = "Term infant → BPD → BPD+PH: Showing cell proportion changes",
       x = "UMAP1", y = "UMAP2") +
  # 主题
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, color = "grey40", hjust = 0.5),
    strip.text = element_text(size = 12, face = "bold"),
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  ) +
  coord_fixed()

print(facet_with_stats)

# 保存带统计信息的分面图
ggsave("./02_UMAP/Adventitial_FB_by_condition_with_stats.pdf", 
       facet_with_stats, 
       width = 15, height = 5, dpi = 300)

#######2. 细胞状态转变分析
# 提取 Adventitial FB 细胞
meso_adventitial_FB <- subset(meso_relabel_no_preterm, celltype == "Adventitial FB")
meso_adventitial_FB@meta.data <- meso_adventitial_FB@meta.data[,-c(12:15)]
# 重新运行 SCTransform（确保数据正确处理）
meso_adventitial_FB <- SCTransform(meso_adventitial_FB,
                                   variable.features.n = 2000,
                                   conserve.memory = FALSE,
                                   vars.to.regress = c("percent.mt", "percent.rb"),
                                   method = "glmGamPoi")


# 运行修正后的聚类函数
# 运行修正后的聚类函数
meso_adventitial_FB <- cluster_pca_umap_fixed(meso_adventitial_FB,
                                              dims_umap = 1:30,
                                              dims_neighbors = 1:30,
                                              k_param = 20,
                                              cluster_res = seq(from = 0.1, to = 1.0,by = 0.1),
                                              harmony = FALSE)

library(clustree)
pdf(file = "./02_UMAP/clustree.pdf",height = 6,width = 15)
clustree(meso_adventitial_FB) +
  theme(legend.position = "bottom") + 
  scale_color_manual(values = color_scanpy_viridis20)+
  scale_edge_color_continuous(low = "grey80", high = "red")
dev.off()

# 可视化聚类结果 SCT_snn_res.0.5
meso_adventitial_FB$condition <- factor(meso_adventitial_FB$condition, levels = c("Term infant","BPD","BPD+PH"))

p1 <- DimPlot(meso_adventitial_FB, group.by = "SCT_snn_res.0.5", label = TRUE) +
  ggtitle("Adventitial FB Subclusters")
p2 <- DimPlot(meso_adventitial_FB, group.by = "condition", label = TRUE) +
  ggtitle("By Condition")
p3 <- DimPlot(meso_adventitial_FB, group.by = "SCT_snn_res.0.5", label = TRUE,split.by = "condition") +
  ggtitle("Adventitial FB Subclusters by Condition")

(p1 + p2) / p3

###重新绘制进行美化

library(ggplot2)
library(patchwork)
library(RColorBrewer)
# "#01353D" "#088096" "#58B3C7" "#7AD4E4" "#B8FCFC"
# "#dadada" "#dadada" "#dadada" "#E7A19E"
# "#84D7E1" "#C71000" "#008EA0" "#8A4198" "#5A9599" "#FF6348" "#FF95A8" "#3D3B25" "#FF6F00" "#ADE2D0"
# "#1A5354" "#3F4041"
# "#2ca02c" "#1f77b4" "#ff7f0e" "#d62728" "#9467bd"

# 设置颜色
condition_colors <- c("Term infant" = "#dadada", "BPD" = "#088096", "BPD+PH" = "#E7A19E")

# 创建聚类颜色的函数 - 使用与之前一致的调色板
get_cluster_colors <- function(n) {
  if (n <= 8) {
    return(brewer.pal(n, "Set2"))
  } else if (n <= 12) {
    return(brewer.pal(n, "Set3"))
  } else {
    return(colorRampPalette(brewer.pal(8, "Set2"))(n))
  }
}

# 获取聚类数量
n_clusters <- length(unique(meso_adventitial_FB$SCT_snn_res.0.5))
cluster_colors <- get_cluster_colors(n_clusters)
#"#66C2A5" "#FC8D62" "#8DA0CB" "#E78AC3" "#A6D854" "#FFD92F"

# 1. 美化聚类图
p1 <- DimPlot(meso_adventitial_FB, 
              group.by = "SCT_snn_res.0.5", 
              label = TRUE,
              label.size = 4,
              pt.size = 0.8,
              cols = cluster_colors) +
  ggtitle("Adventitial FB Subclusters") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5, 
                              margin = margin(b = 10)),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  ) +
  labs(subtitle = paste(n_clusters, "distinct subclusters identified"))

# 2. 美化条件图
p2 <- DimPlot(meso_adventitial_FB, 
              group.by = "condition", 
              label = FALSE,  # 改为FALSE避免标签重叠
              pt.size = 0.8,
              cols = condition_colors) +
  ggtitle("Adventitial FB by Condition") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5, 
                              margin = margin(b = 10)),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  ) +
  labs(subtitle = "Distribution across disease states")

# 3. 美化分面图
p3 <- DimPlot(meso_adventitial_FB, 
              group.by = "SCT_snn_res.0.5", 
              label = TRUE,
              label.size = 3,
              pt.size = 0.6,
              split.by = "condition",
              cols = cluster_colors) +
  ggtitle("Adventitial FB Subclusters by Condition") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5, 
                              margin = margin(b = 10)),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    strip.text = element_text(size = 12, face = "bold")
  ) +
  labs(subtitle = "Subcluster distribution across disease progression")

# 4. 创建组合图
combined_plot <- (p1 + p2) / p3 + 
  plot_annotation(
    title = "Adventitial FB Subpopulation Analysis in BPD Pathogenesis",
    subtitle = paste(
      "Comprehensive characterization of", n_clusters, "Adventitial FB subclusters",
      "across Term infant, BPD, and BPD+PH conditions"
    ),
    theme = theme(
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5,
                                margin = margin(b = 5)),
      plot.subtitle = element_text(size = 14, color = "grey40", hjust = 0.5,
                                   lineheight = 1.4, margin = margin(b = 15))
    )
  )

# 显示所有图形
print(p1)
print(p2)
print(p3)
print(combined_plot)

# 保存单独图形
ggsave("./03_Cluster/Adventitial_FB_subclusters.pdf", p1, 
       width = 8, height = 6, dpi = 300)
ggsave("./03_Cluster/Adventitial_FB_by_condition.pdf", p2, 
       width = 8, height = 6, dpi = 300)
ggsave("./03_Cluster/Adventitial_FB_subclusters_by_condition.pdf", p3, 
       width = 12, height = 6, dpi = 300)
ggsave("./03_Cluster/Adventitial_FB_combined_analysis.pdf", combined_plot, 
       width = 16, height = 12, dpi = 300)

# 输出统计信息
cat("Adventitial FB 亚群分析统计:\n")
cat("- 总细胞数:", ncol(meso_adventitial_FB), "\n")
cat("- 亚群数量:", n_clusters, "\n")
cat("- 各条件细胞分布:\n")
print(table(meso_adventitial_FB$condition))
# Term infant         BPD      BPD+PH 
# 812         710          98 


# 科学论文风格版本
p1_scientific <- DimPlot(meso_adventitial_FB, 
                         group.by = "SCT_snn_res.0.5", 
                         label = TRUE,
                         label.size = 3.5,
                         pt.size = 0.7,
                         cols = cluster_colors) +
  ggtitle("a") +  # 用于论文中的面板标签
  theme_classic() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0),
    legend.position = "right",
    legend.title = element_blank(),
    axis.line = element_line(size = 0.5),
    axis.title = element_text(face = "bold", size = 12),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  ) +
  labs(x = "UMAP1", y = "UMAP2")

p2_scientific <- DimPlot(meso_adventitial_FB, 
                         group.by = "condition", 
                         label = FALSE,
                         pt.size = 0.7,
                         cols = condition_colors) +
  ggtitle("b") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0),
    legend.position = "right",
    legend.title = element_blank(),
    axis.line = element_line(size = 0.5),
    axis.title = element_text(face = "bold", size = 12),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  ) +
  labs(x = "UMAP1", y = "UMAP2")

p3_scientific <- DimPlot(meso_adventitial_FB, 
                         group.by = "SCT_snn_res.0.5", 
                         label = TRUE,
                         label.size = 2.5,
                         pt.size = 0.5,
                         split.by = "condition",
                         cols = cluster_colors) +
  ggtitle("c") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0),
    legend.position = "right",
    legend.title = element_blank(),
    axis.line = element_line(size = 0.5),
    axis.title = element_text(face = "bold", size = 10),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    strip.text = element_text(size = 11, face = "bold")
  ) +
  labs(x = "UMAP1", y = "UMAP2")

# 科学论文风格的组合图
combined_scientific <- (p1_scientific + p2_scientific) / p3_scientific + 
  plot_annotation(
    title = "Figure 3: Adventitial Fibroblast Heterogeneity in Bronchopulmonary Dysplasia",
    theme = theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
    )
  )

# 保存科学论文风格图形
ggsave("./03_Cluster/Adventitial_FB_subclusters_scientific.pdf", p1_scientific, 
       width = 6, height = 5, dpi = 1200)
ggsave("./03_Cluster/Adventitial_FB_by_condition_scientific.pdf", p2_scientific, 
       width = 6, height = 5, dpi = 1200)
ggsave("./03_Cluster/Adventitial_FB_subclusters_by_condition_scientific.pdf", p3_scientific, 
       width = 10, height = 5, dpi = 1200)
ggsave("./03_Cluster/Adventitial_FB_combined_scientific.pdf", combined_scientific, 
       width = 12, height = 10, dpi = 1200)


#############3.3 Adventitial FB:细胞占比#############
library(ggplot2)
library(dplyr)
library(scales)

# 计算亚群在不同condition中的比例
freq_table <- prop.table(table(meso_adventitial_FB$SCT_snn_res.0.5,
                               meso_adventitial_FB$condition),
                         margin = 2) * 100

# 转换为数据框
freq_df <- as.data.frame(freq_table)
names(freq_df) <- c("Cluster", "Condition", "Proportion")

# 设置顺序
condition_order <- c("Term infant", "BPD", "BPD+PH")
cluster_order <- sort(unique(meso_adventitial_FB$SCT_snn_res.0.5))

freq_df$Condition <- factor(freq_df$Condition, levels = condition_order)
freq_df$Cluster <- factor(freq_df$Cluster, levels = cluster_order)

# 获取聚类颜色（与之前保持一致）
get_cluster_colors <- function(n) {
  if (n <= 8) {
    return(brewer.pal(n, "Set2"))
  } else if (n <= 12) {
    return(brewer.pal(n, "Set3"))
  } else {
    return(colorRampPalette(brewer.pal(8, "Set2"))(n))
  }
}

cluster_colors <- get_cluster_colors(length(cluster_order))

# 创建比例图
celltype_proportions <- ggplot(freq_df, 
                               aes(x = Condition, 
                                   y = Proportion, 
                                   fill = Cluster)) +
  geom_col(position = "stack", width = 0.7) +
  scale_fill_manual(values = cluster_colors, name = "Subcluster") +
  
  # 添加百分比标签
  geom_text(aes(label = ifelse(Proportion > 5, 
                               paste0(round(Proportion, 1), "%"), 
                               "")),
            position = position_stack(vjust = 0.5),
            size = 3, color = "white", fontface = "bold") +
  
  # 坐标轴和标签
  labs(title = "Adventitial FB Subcluster Proportions by Condition",
       subtitle = "Distribution of cellular subpopulations across disease states",
       x = "Condition", 
       y = "Cell Proportion (%)") +
  
  # 主题设置
  theme_minimal() +
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5,
                              margin = margin(b = 8)),
    plot.subtitle = element_text(size = 12, color = "grey40", hjust = 0.5,
                                 margin = margin(b = 15)),
    axis.text.x = element_text(size = 12, color = "black", angle = 0),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.title.x = element_text(size = 13, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 13, face = "bold", margin = margin(r = 10)),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "grey80", fill = NA, size = 0.5)
  ) +
  
  # 坐标轴扩展
  scale_x_discrete(expand = expansion(mult = 0.1)) +
  scale_y_continuous(expand = expansion(mult = 0.02),
                     breaks = seq(0, 100, 20),
                     labels = function(x) paste0(x, "%"))

# 显示图形
print(celltype_proportions)

# 保存图形
ggsave("./04_Cell_prop/Adventitial_FB_subcluster_proportions.pdf",
       celltype_proportions,
       width = 5,
       height = 6,
       dpi = 1200,
       bg = "white")

#重新绘制
# 增强版本 - 包含更多统计信息
detailed_stats <- meso_adventitial_FB@meta.data %>%
  group_by(condition, SCT_snn_res.0.5) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  group_by(condition) %>%
  mutate(Total = sum(Count),
         Proportion = Count / Total * 100,
         Label = ifelse(Proportion >= 3, 
                        paste0(round(Proportion, 1), "%\n(n=", Count, ")"), 
                        ""))

# 创建总数标签的数据框（单独创建，避免美学映射冲突）
total_labels <- detailed_stats %>% 
  distinct(condition, Total)

# 创建增强版比例图
enhanced_proportions <- ggplot(detailed_stats, 
                               aes(x = factor(condition, levels = condition_order), 
                                   y = Proportion, 
                                   fill = factor(SCT_snn_res.0.5, levels = cluster_order))) +
  geom_col(position = "stack", width = 0.7, alpha = 0.9) +
  scale_fill_manual(values = cluster_colors, name = "Subcluster") +
  
  # 添加详细的标签
  geom_text(aes(label = Label),
            position = position_stack(vjust = 0.5),
            size = 2.8, color = "white", fontface = "bold", lineheight = 0.8) +
  
  # 添加总数标签 - 修复：使用独立的数据框并禁用美学继承
  geom_text(data = total_labels,
            aes(x = condition, y = 102, label = paste("Total:", Total)),
            inherit.aes = FALSE,  # 关键修复：不继承全局美学映射
            size = 3.5, color = "grey30", fontface = "bold") +
  
  # 坐标轴和标签
  labs(title = "Adventitial FB Subcluster Distribution in BPD Pathogenesis",
       subtitle = "Quantitative analysis of cellular heterogeneity across disease progression",
       x = "Clinical Condition", 
       y = "Cell Proportion (%)",
       caption = paste("Analysis based on", ncol(meso_adventitial_FB), 
                       "Adventitial FB cells across", length(cluster_order), "subclusters")) +
  
  # 主题设置
  theme_minimal() +
  theme(
    plot.title = element_text(size = 10, face = "bold", hjust = 0.5,
                              margin = margin(b = 8)),
    plot.subtitle = element_text(size = 8, color = "grey40", hjust = 0.5,
                                 margin = margin(b = 15)),
    plot.caption = element_text(size = 8, color = "grey50", hjust = 0.5,
                                margin = margin(t = 15)),
    axis.text.x = element_text(size = 8, color = "black", face = "bold"),
    axis.text.y = element_text(size = 8, color = "black"),
    axis.title.x = element_text(size = 8, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 8, face = "bold", margin = margin(r = 10)),
    legend.title = element_text(size = 8, face = "bold"),
    legend.text = element_text(size = 8),
    panel.grid.major.y = element_line(color = "grey90", size = 0.3),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.border = element_rect(color = "grey70", fill = NA, size = 0.8)
  ) +
  
  # 坐标轴设置
  scale_x_discrete(expand = expansion(mult = 0.1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)),
                     breaks = seq(0, 100, 20),
                     labels = function(x) paste0(x, "%")) +
  
  # 限制y轴范围以容纳总数标签
  coord_cartesian(ylim = c(0, 105))

# 显示图形
print(enhanced_proportions)

# 保存增强版图形
ggsave("./04_Cell_prop/Adventitial_FB_subcluster_proportions_enhanced.pdf",
       enhanced_proportions,
       width = 5,
       height = 6,
       dpi = 1200,
       bg = "white")


#############3.3 Adventitial FB:细胞差异基因获取#############
# 设置默认的细胞身份为聚类结果
Idents(meso_adventitial_FB) <- "SCT_snn_res.0.5"
# 设置要去除的基因
rm_genes <- c("FP236383.3", "FP671120.4")
ribo_mito_pattern <- "^RP[SL]?|^MT-"

# 寻找所有聚类的marker基因
cluster_markers <- FindAllMarkers(
  meso_adventitial_FB,
  assay = "SCT",           # 使用SCT assay
  only.pos = TRUE,         # 只找正表达的marker
  min.pct = 0.25,          # 在至少25%的细胞中表达
  logfc.threshold = 0.25,  # 对数倍变化阈值
  test.use = "wilcox"      # 使用Wilcoxon秩和检验
)

# 查看每个聚类的前几个marker基因
head(cluster_markers)

cluster_markers <- cluster_markers[c(!rownames(cluster_markers) %in% rm_genes & 
                             !grepl(ribo_mito_pattern, rownames(cluster_markers))),]
# 按聚类和p值排序
cluster_markers_sorted <- cluster_markers %>%
  group_by(cluster) %>%
  arrange(p_val_adj, .by_group = TRUE)

# 保存结果
write.csv(cluster_markers_sorted, "./05_DEGs/adventitial_FB_cluster_markers.csv", row.names = FALSE)

# 提取每个聚类的前10个marker基因
top10_markers <- cluster_markers_sorted %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) %>%
  arrange(cluster, desc(avg_log2FC))

# 查看结果
print(top10_markers, n = 50)

# 用热图可视化前5个marker基因
top5_markers <- cluster_markers_sorted %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)

all_marker_plot <- DoHeatmap(
  meso_adventitial_FB,
  features = top5_markers$gene,
  assay = "SCT",
  hjust = 1,
  vjust = 0.5,
  group.colors = cluster_colors,
  angle = 0 , group.bar.height = 0.04) + 
  labs(title = "Top 5 Marker Genes for Each Adventitial FB Subcluster")

all_marker_plot_w_theme <- all_marker_plot +
  scale_fill_gradientn(colors = rev(color_glacier),
                       name = "Expression") +
  theme(axis.title = element_blank(),
        aspect.ratio = 1,
        plot.background = element_rect(color = "transparent",
                                       fill = "transparent"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12))

ggsave("./05_DEGs/adventitial_FB_heatmap_markers.pdf",
       egg::set_panel_size(p = all_marker_plot_w_theme,
                           width = unit(10, "cm"),
                           height = unit(10, "cm")),
       dpi = 1200,
       bg = "transparent",
       width = unit(7, "cm"),
       height = unit(5, "cm"))

cat(top10_markers$gene[top10_markers$cluster == 0], sep = "\n")


#############3.4 Adventitial FB:条件特异性差异表达分析#############
library(scRNAtoolVis)
library(dplyr)
library(stringr)

# 在每个cluster内比较不同条件的差异基因
# 在每个cluster内比较不同条件的差异基因
Idents(meso_adventitial_FB) <- "SCT_snn_res.0.5"

# 设置并行计算（如果需要）
library(future.apply)
plan(multisession, workers = 4)

# 获取所有聚类
clusters <- levels(Idents(meso_adventitial_FB))

# 定义三个比较的函数
run_condition_comparisons <- function(cluster_id) {
  # 设置细胞标识为聚类
  Idents(meso_adventitial_FB) <- "SCT_snn_res.0.5"
  
  # 提取特定cluster的细胞
  cluster_cells <- subset(meso_adventitial_FB, idents = cluster_id)
  Idents(cluster_cells) <- cluster_cells$condition
  
  # 检查各条件的细胞数量
  cell_counts <- table(Idents(cluster_cells))
  cat("Cluster", cluster_id, "- Cell counts:", 
      paste(names(cell_counts), cell_counts, sep = ":", collapse = ", "), "\n")
  
  results <- list()
  
  # 比较1: BPD vs Term infant
  if(sum(cluster_cells$condition == "BPD") >= 3 && 
     sum(cluster_cells$condition == "Term infant") >= 3) {
    
    # 切换到 RNA assay 并重新标准化
    DefaultAssay(cluster_cells) <- "RNA"
    cluster_cells <- NormalizeData(cluster_cells, 
                                   normalization.method = "LogNormalize", 
                                   scale.factor = 10000)
    
    markers_bpd_vs_term <- FindMarkers(
      cluster_cells,
      ident.1 = "BPD",
      ident.2 = "Term infant",
      assay = "RNA",
      test.use = "MAST",
      logfc.threshold = 0
    )
    
    # 定义上下调基因
    markers_bpd_vs_term$up_in_bpd_down_in_term <- ifelse(
      markers_bpd_vs_term$avg_log2FC > 1 & 
        markers_bpd_vs_term$p_val_adj < 1e-05 & 
        !grepl("^RP|^MT-", rownames(markers_bpd_vs_term)), 
      TRUE, FALSE
    )
    
    markers_bpd_vs_term$up_in_term_down_in_bpd <- ifelse(
      markers_bpd_vs_term$avg_log2FC < -1 & 
        markers_bpd_vs_term$p_val_adj < 1e-05 & 
        !grepl("^RP|^MT-", rownames(markers_bpd_vs_term)), 
      TRUE, FALSE
    )
    
    results[["BPD_vs_Term"]] <- markers_bpd_vs_term
  } else {
    results[["BPD_vs_Term"]] <- data.frame()
    cat("  Skipping BPD vs Term infant - insufficient cells\n")
  }
  
  # 比较2: BPD+PH vs Term infant
  if(sum(cluster_cells$condition == "BPD+PH") >= 3 && 
     sum(cluster_cells$condition == "Term infant") >= 3) {
    
    # 确保使用 RNA assay
    DefaultAssay(cluster_cells) <- "RNA"
    cluster_cells <- NormalizeData(cluster_cells, 
                                   normalization.method = "LogNormalize", 
                                   scale.factor = 10000)
    
    markers_bpdph_vs_term <- FindMarkers(
      cluster_cells,
      ident.1 = "BPD+PH",
      ident.2 = "Term infant",
      assay = "RNA",
      test.use = "MAST",
      logfc.threshold = 0
    )
    
    # 定义上下调基因
    markers_bpdph_vs_term$up_in_bpdph_down_in_term <- ifelse(
      markers_bpdph_vs_term$avg_log2FC > 1 & 
        markers_bpdph_vs_term$p_val_adj < 1e-05 & 
        !grepl("^RP|^MT-", rownames(markers_bpdph_vs_term)), 
      TRUE, FALSE
    )
    
    markers_bpdph_vs_term$up_in_term_down_in_bpdph <- ifelse(
      markers_bpdph_vs_term$avg_log2FC < -1 & 
        markers_bpdph_vs_term$p_val_adj < 1e-05 & 
        !grepl("^RP|^MT-", rownames(markers_bpdph_vs_term)), 
      TRUE, FALSE
    )
    
    results[["BPD_PH_vs_Term"]] <- markers_bpdph_vs_term
  } else {
    results[["BPD_PH_vs_Term"]] <- data.frame()
    cat("  Skipping BPD+PH vs Term infant - insufficient cells\n")
  }
  
  # 比较3: BPD+PH vs BPD
  if(sum(cluster_cells$condition == "BPD+PH") >= 3 && 
     sum(cluster_cells$condition == "BPD") >= 3) {
    
    # 确保使用 RNA assay
    DefaultAssay(cluster_cells) <- "RNA"
    cluster_cells <- NormalizeData(cluster_cells, 
                                   normalization.method = "LogNormalize", 
                                   scale.factor = 10000)
    
    markers_bpdph_vs_bpd <- FindMarkers(
      cluster_cells,
      ident.1 = "BPD+PH",
      ident.2 = "BPD",
      assay = "RNA",
      test.use = "MAST",
      logfc.threshold = 0
    )
    
    # 定义上下调基因
    markers_bpdph_vs_bpd$up_in_bpdph_down_in_bpd <- ifelse(
      markers_bpdph_vs_bpd$avg_log2FC > 1 & 
        markers_bpdph_vs_bpd$p_val_adj < 1e-05 & 
        !grepl("^RP|^MT-", rownames(markers_bpdph_vs_bpd)), 
      TRUE, FALSE
    )
    
    markers_bpdph_vs_bpd$up_in_bpd_down_in_bpdph <- ifelse(
      markers_bpdph_vs_bpd$avg_log2FC < -1 & 
        markers_bpdph_vs_bpd$p_val_adj < 1e-05 & 
        !grepl("^RP|^MT-", rownames(markers_bpdph_vs_bpd)), 
      TRUE, FALSE
    )
    
    results[["BPD_PH_vs_BPD"]] <- markers_bpdph_vs_bpd
  } else {
    results[["BPD_PH_vs_BPD"]] <- data.frame()
    cat("  Skipping BPD+PH vs BPD - insufficient cells\n")
  }
  
  return(results)
}

# 并行运行所有比较
condition_DEGs <- future_lapply(clusters, run_condition_comparisons)
names(condition_DEGs) <- paste0("Cluster", clusters)

# 保存结果
saveRDS(condition_DEGs, "./05_DEGs/adventitial_FB_condition_specific_DEGs.rds")

# 处理结果并生成汇总统计
summarize_DEG_results <- function(deg_list) {
  summary_df <- data.frame()
  
  for(cluster_name in names(deg_list)) {
    cluster_results <- deg_list[[cluster_name]]
    
    for(comparison in names(cluster_results)) {
      deg_data <- cluster_results[[comparison]]
      
      if(nrow(deg_data) > 0) {
        # 计算上下调基因数量
        up_col <- grep("up_in.*down_in.*", colnames(deg_data), value = TRUE)[1]
        down_col <- grep("up_in.*down_in.*", colnames(deg_data), value = TRUE)[2]
        
        n_up <- sum(deg_data[[up_col]], na.rm = TRUE)
        n_down <- sum(deg_data[[down_col]], na.rm = TRUE)
        n_total <- n_up + n_down
        
        summary_df <- rbind(summary_df, data.frame(
          Cluster = cluster_name,
          Comparison = comparison,
          Total_DEGs = n_total,
          Upregulated = n_up,
          Downregulated = n_down,
          stringsAsFactors = FALSE
        ))
      } else {
        summary_df <- rbind(summary_df, data.frame(
          Cluster = cluster_name,
          Comparison = comparison,
          Total_DEGs = 0,
          Upregulated = 0,
          Downregulated = 0,
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  
  return(summary_df)
}

# 生成汇总统计
deg_summary <- summarize_DEG_results(condition_DEGs)
print(deg_summary)

# 保存汇总统计
write.csv(deg_summary, "./05_DEGs/adventitial_FB_DEG_summary.csv", row.names = FALSE)

# 提取所有显著DEGs用于后续分析
extract_significant_genes <- function(deg_list) {
  all_sig_genes <- data.frame()
  
  for(cluster_name in names(deg_list)) {
    cluster_results <- deg_list[[cluster_name]]
    
    for(comparison in names(cluster_results)) {
      deg_data <- cluster_results[[comparison]]
      
      if(nrow(deg_data) > 0) {
        # 获取上下调列名
        up_col <- grep("up_in.*down_in.*", colnames(deg_data), value = TRUE)[1]
        down_col <- grep("up_in.*down_in.*", colnames(deg_data), value = TRUE)[2]
        
        # 提取显著基因
        sig_up <- deg_data[deg_data[[up_col]] == TRUE, ]
        sig_down <- deg_data[deg_data[[down_col]] == TRUE, ]
        
        if(nrow(sig_up) > 0) {
          sig_up$direction <- "up"
          sig_up$comparison <- comparison
          sig_up$cluster <- cluster_name
          sig_up$gene <- rownames(sig_up)
          all_sig_genes <- rbind(all_sig_genes, sig_up)
        }
        
        if(nrow(sig_down) > 0) {
          sig_down$direction <- "down"
          sig_down$comparison <- comparison
          sig_down$cluster <- cluster_name
          sig_down$gene <- rownames(sig_down)
          all_sig_genes <- rbind(all_sig_genes, sig_down)
        }
      }
    }
  }
  
  return(all_sig_genes)
}

# 提取所有显著基因
significant_genes <- extract_significant_genes(condition_DEGs)

if(nrow(significant_genes) > 0) {
  write.csv(significant_genes, "./05_DEGs/adventitial_FB_significant_genes.csv", row.names = FALSE)
  cat("Found", nrow(significant_genes), "significant DEGs across all clusters and comparisons\n")
} else {
  cat("No significant DEGs found with the current thresholds\n")
}

######绘制差异基因的环状的火山图
# 使用之前定义的函数提取显著基因
# 从condition_DEGs中提取BPD_vs_Term的所有基因数据（包括不显著基因）
# 从condition_DEGs中提取BPD_vs_Term的所有基因数据（包括不显著基因）
prepare_volcano_data <- function(condition_DEGs, clusters = 0:4) {
  volcano_data <- data.frame()
  
  for(clust in clusters) {
    cluster_name <- paste0("Cluster", clust)
    
    # 检查该cluster是否有BPD_vs_Term数据
    if(!is.null(condition_DEGs[[cluster_name]]$BPD_vs_Term) && 
       nrow(condition_DEGs[[cluster_name]]$BPD_vs_Term) > 0) {
      
      deg_data <- condition_DEGs[[cluster_name]]$BPD_vs_Term
      
      # 添加必要的列
      deg_data$gene <- rownames(deg_data)
      deg_data$cluster <- as.character(clust)
      
      # 选择需要的列并重命名以匹配jjVolcano的要求
      deg_subset <- deg_data %>%
        select(
          p_val,
          avg_log2FC, 
          pct.1,
          pct.2,
          p_val_adj,
          cluster,
          gene
        )
      
      volcano_data <- rbind(volcano_data, deg_subset)
    }
  }
  
  return(volcano_data)
}

# 准备火山图数据 - 只使用cluster0-4
volcano_data <- prepare_volcano_data(condition_DEGs, clusters = 0:4)

# 检查数据
print(paste("Total genes in volcano plot:", nrow(volcano_data)))
print("Genes per cluster:")
print(table(volcano_data$cluster))

# 过滤基因（去除核糖体和线粒体基因）
volcano_data_filtered <- volcano_data %>%
  filter(!str_detect(gene, "^RP[SL]?[0-9]"),
         !str_detect(gene, "^MT-"),
         !str_detect(gene, "^RP[SL]?$"),
         !str_detect(gene, "FP236383.3"),
         !str_detect(gene, "FP671120.4"))

# 设置聚类颜色
cluster_colors_volcano <- c(
  "0" = "#66C2A5",  # 蓝色
  "1" = "#FC8D62",  # 橙色
  "2" = "#8DA0CB",  # 绿色
  "3" = "#E78AC3",  # 红色
  "4" = "#A6D854"   # 紫色
)

# 方法3：使用其他包创建环形火山图

# 绘制环形火山图
circular_volcano <- jjVolcano(
  diffData = volcano_data_filtered,
  log2FC.cutoff = 1,
  adjustP.cutoff = 1e-05,
  size = 3.5,
  pSize = 0.8,
  base_size = 12,
  legend.position = c(0.85, 0.05),
  tile.col = cluster_colors_volcano,
  topGeneN = 5,                # 每个cluster显示前3个基因
  polar = F,
  cluster.order = c("0", "1", "2", "3", "4")  # 指定聚类顺序
)

# 保存图片
ggsave(
  filename = "./05_DEGs/BPD_vs_Term_circular_volcano.pdf",
  plot = circular_volcano,height = 5,width = 9)

ggsave(
  filename = "./05_DEGs/BPD_vs_Term_circular_volcano.png",
  plot = circular_volcano,height = 5,width = 9,dpi = 1200
)

# > table(significant_genes$cluster)
# Cluster0 Cluster1 Cluster2 Cluster3 Cluster4 
# 31       89       20        1       37 



#############3.4 Adventitial FB:通路富集#############
# 对每个cluster的DEGs进行通路富集
# 基因集富集分析 (GSEA) - 使用差异基因
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(DOSE)

# 准备ORA分析函数 - 使用差异基因
run_ora_analysis_diff_genes <- function(significant_genes_data, cluster_id, background_genes = NULL) {
  # 过滤特定cluster的差异基因
  cluster_genes <- significant_genes_data %>%
    filter(cluster == paste0("Cluster", cluster_id))
  
  if(nrow(cluster_genes) == 0) {
    cat("Cluster", cluster_id, "- No significant DEGs for ORA analysis\n")
    return(NULL)
  }
  
  # 获取差异基因列表
  deg_list <- cluster_genes$gene
  
  # 去除重复基因名
  deg_list <- unique(deg_list)
  
  if(length(deg_list) < 3) {
    cat("Cluster", cluster_id, "- Too few significant genes for ORA (n =", length(deg_list), ")\n")
    return(NULL)
  }
  
  cat("Running ORA for Cluster", cluster_id, "-", length(deg_list), "significant genes\n")
  
  results <- list()
  
  # GO生物过程富集分析
  ora_go <- enrichGO(
    gene = deg_list,
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL",
    ont = "BP",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    qvalueCutoff = 0.2,
    minGSSize = 10,
    maxGSSize = 500,
    readable = FALSE
  )
  
  # KEGG通路分析
  # 首先将基因符号转换为ENTREZ ID
  gene_entrez <- bitr(deg_list, fromType = "SYMBOL", 
                      toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  
  if(nrow(gene_entrez) > 0) {
    entrez_genes <- gene_entrez$ENTREZID
    
    ora_kegg <- enrichKEGG(
      gene = entrez_genes,
      organism = "hsa",
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      qvalueCutoff = 0.2,
      minGSSize = 10,
      maxGSSize = 500
    )
  } else {
    ora_kegg <- NULL
    cat("  No gene symbols could be converted to ENTREZ ID for KEGG analysis\n")
  }
  
  # 疾病本体论分析
  ora_do <- enrichDO(
    gene = deg_list,
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    qvalueCutoff = 0.2,
    minGSSize = 5,
    maxGSSize = 500
  )
  
  # 细胞组分和分子功能GO分析
  ora_cc <- enrichGO(
    gene = deg_list,
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL",
    ont = "CC",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    qvalueCutoff = 0.2,
    minGSSize = 10,
    maxGSSize = 500
  )
  
  ora_mf <- enrichGO(
    gene = deg_list,
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL",
    ont = "MF",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    qvalueCutoff = 0.2,
    minGSSize = 10,
    maxGSSize = 500
  )
  
  results <- list(
    GO_BP = ora_go,
    GO_CC = ora_cc,
    GO_MF = ora_mf,
    KEGG = ora_kegg,
    DO = ora_do
  )
  
  return(results)
}

# 对每个有显著差异基因的cluster进行ORA分析
ora_results_diff <- list()

for(clust in 0:4) {
  ora_results_diff[[paste0("Cluster", clust)]] <- run_ora_analysis_diff_genes(significant_genes, clust)
}

# 保存ORA结果
saveRDS(ora_results_diff, "./06_GSEA/adventitial_FB_ORA_results_diff_genes.rds")

# 生成ORA结果摘要
generate_ora_summary <- function(ora_results) {
  summary_df <- data.frame()
  
  for(cluster_name in names(ora_results)) {
    cluster_results <- ora_results[[cluster_name]]
    
    if(!is.null(cluster_results)) {
      for(analysis_type in names(cluster_results)) {
        ora_obj <- cluster_results[[analysis_type]]
        
        if(!is.null(ora_obj) && nrow(ora_obj) > 0) {
          n_terms <- nrow(ora_obj)
          
          summary_df <- rbind(summary_df, data.frame(
            Cluster = cluster_name,
            Analysis = analysis_type,
            N_Terms = n_terms,
            Top_Term = ora_obj$Description[1],
            Top_P_Value = ora_obj$pvalue[1],
            Top_Genes = paste(head(strsplit(ora_obj$geneID[1], "/")[[1]], 3), collapse = ","),
            stringsAsFactors = FALSE
          ))
        } else {
          summary_df <- rbind(summary_df, data.frame(
            Cluster = cluster_name,
            Analysis = analysis_type,
            N_Terms = 0,
            Top_Term = "None",
            Top_P_Value = NA,
            Top_Genes = "",
            stringsAsFactors = FALSE
          ))
        }
      }
    } else {
      # 如果没有结果，也记录
      for(analysis_type in c("GO_BP", "GO_CC", "GO_MF", "KEGG", "DO")) {
        summary_df <- rbind(summary_df, data.frame(
          Cluster = cluster_name,
          Analysis = analysis_type,
          N_Terms = 0,
          Top_Term = "No analysis",
          Top_P_Value = NA,
          Top_Genes = "",
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  
  return(summary_df)
}

# 生成并显示摘要
ora_summary <- generate_ora_summary(ora_results_diff)
print(ora_summary)

# 保存摘要
write.csv(ora_summary, "./06_GSEA/adventitial_FB_ORA_summary_diff_genes.csv", row.names = FALSE)

# 可视化ORA结果
plot_ora_results <- function(ora_results, output_dir = "./06_GSEA/") {
  for(cluster_name in names(ora_results)) {
    cluster_results <- ora_results[[cluster_name]]
    
    if(!is.null(cluster_results)) {
      for(analysis_type in names(cluster_results)) {
        ora_obj <- cluster_results[[analysis_type]]
        
        if(!is.null(ora_obj) && nrow(ora_obj) > 0) {
          # 创建点图
          dot_plot <- dotplot(ora_obj, showCategory = 15, 
                              title = paste(cluster_name, analysis_type, "ORA - Differential Genes"))
          
          # 保存点图
          ggsave(
            filename = paste0(output_dir, "ORA_", cluster_name, "_", analysis_type, "_diff_genes_dotplot.pdf"),
            plot = dot_plot,
            width = 10,
            height = 8
          )
          
          # 对GO和KEGG结果创建网络图
          if(analysis_type %in% c("GO_BP", "GO_CC", "GO_MF", "KEGG") && nrow(ora_obj) >= 3) {
            tryCatch({
              enrich_map <- emapplot(pairwise_termsim(ora_obj), 
                                     showCategory = 15,
                                     layout = "kk")
              
              ggsave(
                filename = paste0(output_dir, "ORA_", cluster_name, "_", analysis_type, "_diff_genes_network.pdf"),
                plot = enrich_map,
                width = 10,
                height = 8
              )
            }, error = function(e) {
              cat("Could not create network plot for", cluster_name, analysis_type, ":", e$message, "\n")
            })
          }
          
          # 创建柱状图
          bar_plot <- barplot(ora_obj, showCategory = 15, 
                              title = paste(cluster_name, analysis_type, "ORA"))
          
          ggsave(
            filename = paste0(output_dir, "ORA_", cluster_name, "_", analysis_type, "_diff_genes_barplot.pdf"),
            plot = bar_plot,
            width = 10,
            height = 8
          )
        }
      }
    }
  }
}

# 绘制ORA结果图
plot_ora_results(ora_results_diff)


###重新绘制图片
# 提取Cluster1的富集分析结果并整理数据
prepare_enrichment_data <- function(ora_results, cluster_name = "Cluster1") {
  cluster_results <- ora_results[[cluster_name]]
  
  # 初始化数据框
  enrichment_data <- data.frame()
  
  # 处理每种类型的富集结果
  analysis_types <- c("GO_BP", "GO_CC", "GO_MF", "KEGG")
  type_labels <- c("BP", "CC", "MF", "KEGG")
  
  for(i in 1:length(analysis_types)) {
    analysis_type <- analysis_types[i]
    type_label <- type_labels[i]
    
    if(!is.null(cluster_results[[analysis_type]]) && 
       nrow(cluster_results[[analysis_type]]@result) > 0) {
      
      # 提取结果数据
      result_df <- cluster_results[[analysis_type]]@result
      
      # 计算GeneRatio数值
      result_df$GeneRatio_value <- sapply(strsplit(result_df$GeneRatio, "/"), 
                                          function(x) as.numeric(x[1]) / as.numeric(x[2]))
      
      # 计算Fold.Enrichment
      bg_ratio <- sapply(strsplit(result_df$BgRatio, "/"), 
                         function(x) as.numeric(x[1]) / as.numeric(x[2]))
      result_df$Fold.Enrichment <- result_df$GeneRatio_value / bg_ratio
      
      # 计算LogFDR
      result_df$LogFDR <- -log10(result_df$p.adjust)
      
      # 添加类型信息
      result_df$Type <- type_label
      
      # 选择需要的列 - 使用基础R方法而不是dplyr::select
      result_df_sub <- result_df[, c("Description", "LogFDR", "Fold.Enrichment", "GeneRatio_value", "Count", "Type")]
      colnames(result_df_sub) <- c("Term", "LogFDR", "Fold.Enrichment", "Ratio", "Count", "Type")
      
      # 按Count排序并取前5个
      result_df_top <- result_df_sub[order(-result_df_sub$Count), ]
      result_df_top <- head(result_df_top, 5)
      
      enrichment_data <- rbind(enrichment_data, result_df_top)
    }
  }
  
  return(enrichment_data)
}


# 准备绘图数据
plot.data <- prepare_enrichment_data(ora_results_diff, "Cluster1")
plot.data_cluster0 <- prepare_enrichment_data(ora_results_diff, "Cluster0")
plot.data_cluster1 <- prepare_enrichment_data(ora_results_diff, "Cluster1")
plot.data_cluster2 <- prepare_enrichment_data(ora_results_diff, "Cluster2")
plot.data_cluster4 <- prepare_enrichment_data(ora_results_diff, "Cluster4")

# 检查数据
print("Enrichment data for plotting:")
print(plot.data)

# 设置颜色和形状
color.bin <- c("#00599F", "#D01910", "#2E8B57", "#8B4513")  # 为四种类型设置不同颜色
names(color.bin) <- c("BP", "CC", "MF", "KEGG")

# 设置形状映射
shape.mapping <- c("BP" = 16, "CC" = 17, "MF" = 15, "KEGG" = 18)  # 圆形, 三角, 方形, 菱形
plot.data$Shape <- shape.mapping[plot.data$Type]

# 绘制散点图
library(ggpubr)

p <- ggscatter(plot.data_cluster4, 
               x = "LogFDR", 
               y = "Fold.Enrichment", 
               color = "Type",
               shape = "Type",
               main = "Cluster1 - Top Enriched Pathways",
               size = "Ratio", 
               label = "Term", 
               palette = color.bin,
               repel = TRUE,
               font.label = list(size = 8, color = "black")
)

# 自定义图形
p <- p + 
  scale_size(range = c(3, 8), name = "Gene Ratio") +
  scale_shape_manual(values = shape.mapping, name = "Category") +
  scale_color_manual(values = color.bin, name = "Category") +
  scale_x_continuous(limits = c(0, max(plot.data$LogFDR) * 1.1), 
                     breaks = seq(0, ceiling(max(plot.data$LogFDR)), by = 1)) +
  xlab(expression(-log[10]("Adjusted P-value"))) +
  ylab("Fold Enrichment") +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "grey50", size = 0.5) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    panel.grid.major = element_line(color = "grey90", size = 0.2),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 12, face = "bold"),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 10)
  )

# 显示图形
print(p)

# 保存图形
ggsave(
  filename = "./06_GSEA/Cluster4_enrichment_bubble_plot.pdf",
  plot = p,width = 12,height = 8,dpi = 1200)

ggsave(
  filename = "./06_GSEA/Cluster4_enrichment_bubble_plot.png",
  plot = p,width = 12,height = 8,dpi = 1200)



# 创建另一种风格的图 - 气泡图
p_bubble <- ggplot(plot.data_cluster4, aes(x = LogFDR, y = Fold.Enrichment)) +
  geom_point(aes(size = Ratio, color = Type, shape = Type), alpha = 0.8) +
  geom_text(aes(label = Term), size = 2, hjust = 0, vjust = 0, nudge_x = 0.1, check_overlap = TRUE) +
  scale_size(range = c(1, 8), name = "Gene Ratio") +
  scale_color_manual(values = color.bin, name = "Category") +
  scale_shape_manual(values = shape.mapping, name = "Category") +
  labs(
    title = "Cluster4 - Top Enriched Pathways",
    x = expression(-log[10]("Adjusted P-value")),
    y = "Fold Enrichment"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 12, face = "bold"),
    legend.position = "right"
  )

# 保存气泡图
ggsave(
  filename = "./06_GSEA/Cluster4_enrichment_bubble_plot_alt.pdf",
  plot = p_bubble,width = 12,height = 8)


# 创建详细的结果表格
detailed_results <- plot.data_cluster4 %>%
  dplyr::select(Term, Type, Fold.Enrichment, LogFDR, Count, Ratio) %>%
  arrange(Type, desc(Count))
write.csv(detailed_results, "./06_GSEA/Cluster4_enrichment_detailed_results.csv", row.names = FALSE)

# 为每个类别单独绘制图表
for(category in c("BP", "CC", "MF", "KEGG")) {
  category_data <- plot.data[plot.data$Type == category, ]
  if(nrow(category_data) > 0) {
    p_category <- ggscatter(category_data, 
                            x = "LogFDR", 
                            y = "Fold.Enrichment", 
                            color = color.bin[category],
                            main = paste("Cluster1 - Top", category, "Pathways"),
                            size = "Ratio", 
                            label = "Term", 
                            repel = TRUE,
                            font.label = list(size = 8, color = "black")) +
      scale_size(range = c(4, 15), name = "Gene Ratio") +
      scale_x_continuous(limits = c(0, max(category_data$LogFDR) * 1.1)) +
      xlab(expression(-log[10]("Adjusted P-value"))) +
      ylab("Fold Enrichment") +
      theme_bw() +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "right")
    ggsave(
      filename = paste0("./06_GSEA/Cluster1_", category, "_enrichment_plot.pdf"),
      plot = p_category,
      width = 10,
      height = 7)
  }
}

#输出详细信息
# 特别关注Cluster1的详细分析
if(!is.null(ora_results_diff$Cluster2)) {
  cat("\n=== Detailed Analysis for Cluster1 (Most DEGs) ===\n")
  
  cluster2_results <- ora_results_diff$Cluster4
  
  for(analysis_type in names(cluster2_results)) {
    ora_obj <- cluster2_results[[analysis_type]]
    
    if(!is.null(ora_obj) && nrow(ora_obj) > 0) {
      cat("\n", analysis_type, "Analysis - Top 5 Terms:\n")
      top5 <- head(ora_obj, 5)
      
      for(i in 1:nrow(top5)) {
        cat(i, ". ", top5$Description[i], 
            " (p.adjust = ", format(top5$p.adjust[i], scientific = TRUE, digits = 3), 
            ", Count = ", top5$Count[i], ")\n", sep = "")
      }
    }
  }
}

# === Detailed Analysis for Cluster1 (Most DEGs) ===
#   
#   GO_BP Analysis - Top 5 Terms:
#   1. negative regulation of growth (p.adjust = 3.69e-04, Count = 10)
# 2. cellular response to zinc ion (p.adjust = 1.64e-03, Count = 4)
# 3. cellular response to copper ion (p.adjust = 3.87e-03, Count = 4)
# 4. response to copper ion (p.adjust = 1.76e-02, Count = 4)
# 5. response to zinc ion (p.adjust = 2.6e-02, Count = 4)
# 
# GO_CC Analysis - Top 5 Terms:
#   1. ficolin-1-rich granule lumen (p.adjust = 1.57e-04, Count = 7)
# 2. ficolin-1-rich granule (p.adjust = 1.12e-03, Count = 7)
# 3. secretory granule lumen (p.adjust = 1.58e-02, Count = 7)
# 4. cytoplasmic vesicle lumen (p.adjust = 1.58e-02, Count = 7)
# 5. vesicle lumen (p.adjust = 1.58e-02, Count = 7)
# 
# GO_MF Analysis - Top 5 Terms:
#   1. unfolded protein binding (p.adjust = 3.6e-03, Count = 6)
# 2. ubiquitin protein ligase binding (p.adjust = 3.84e-02, Count = 7)
# 3. ubiquitin-like protein ligase binding (p.adjust = 3.84e-02, Count = 7)
# 4. ATP-dependent protein folding chaperone (p.adjust = 3.89e-02, Count = 3)
# 5. TFIID-class transcription factor complex binding (p.adjust = 3.89e-02, Count = 2)
# 
# KEGG Analysis - Top 5 Terms:
#   1. Mineral absorption (p.adjust = 8.86e-05, Count = 6)
# 2. Legionellosis (p.adjust = 6.8e-04, Count = 5)
# 3. Kaposi sarcoma-associated herpesvirus infection (p.adjust = 1.83e-02, Count = 6)
# 4. Prion disease (p.adjust = 1.83e-02, Count = 7)
# 5. Measles (p.adjust = 1.93e-02, Count = 5)
