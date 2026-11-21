###############################################################################
# Mesenchymal composition, marker dot plots and differential expression (three comparisons)
#
# Source      : r.03_Cell_prop.r (lines 1-1280)
# Input       : alldata_neonatal_injury.rds / annotated_meso.rds
# Output      : Fig. 2D-H; Supplementary Figs. 4-5; Supplementary Tables 6-7
#
# Code below is unchanged from the original analysis script.
###############################################################################
# -----------------------------------------------------------------------------
# PATHS: this file keeps the absolute paths of the analysis server on which the
# study was run (/data/nas1/xiangxiaoyu_OD/1.BPD+MSCs/1.scRNA/... and
# /home/xxyloveme/study/BPD+MSCs/...). Set PROJ_DIR in config/paths.R and
# replace those prefixes with your own project directory before running.
# -----------------------------------------------------------------------------


######################3.细胞类型占比###################
setwd("/data/nas1/xiangxiaoyu_OD/1.BPD+MSCs/1.scRNA/03_Cell_prop/")

source("./helper_functions/globals.R")
source("./helper_functions/libraries.R")
source("./helper_functions/cluster.R")
source("./helper_functions/colors.R")

library(Seurat)
library(ggplot2)
library(dplyr)
library(ggtext)
meso_relabel <- readRDS(file = "/data/nas1/xiangxiaoyu_OD/1.BPD+MSCs/1.scRNA/02_Cell_Clustering_and_Annotation/3_mesenchyme/annotated_meso.rds")

# 可视化细胞类型比例
# 细胞类型比例柱状图

freq_table <- prop.table(table(meso_relabel$celltype,
                               meso_relabel$condition),
                         margin = 2) * 100

freq_df <- as.data.frame(freq_table)
group_order <- c("Term infant","BPD","BPD+PH","Acute preterm injury")
celltype_order <- levels(meso_relabel$celltype)
freq_df$Var2 <- ordered(as.factor(freq_df$Var2),
                        group_order)
freq_df$Va1 <- ordered(as.factor(freq_df$Var1),
                       celltype_order)

color_safe_endo_proportionplot <- c("#2ca02c","#98df8a","#1f77b4","#aec7e8","#ff7f0e","#ffbb78","#d62728")

celltype_proportions <- ggplot(as.data.frame(freq_table),
                               aes(x=ordered(Var2, group_order),
                                   y=Freq,fill=ordered(Var1, celltype_order))) +
  geom_col() +
  scale_fill_manual(values = color_safe_endo_proportionplot,
                    name = "Condition") +
  xlab("Celltype") +
  ylab("Cell Proportion") +
  theme(legend.text = element_markdown(size = 12),
        legend.title = element_text(size = 12),
        axis.text.x  = element_text(size = 12,
                                    color = "black"),
        axis.text.y = element_text(size = 12,
                                   color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        aspect.ratio = 1.2) +
  scale_x_discrete(expand = c(.05, .05)) +
  scale_y_continuous(expand = c(.01, .01)) +
  coord_flip()

ggsave("./meso_celltype_proportions.pdf",
       celltype_proportions,
       dpi = 1200,
       bg = "transparent",
       width = unit(9, "cm"),
       height = unit(3.5, "cm"))

########################4.注释基因点图#####################
meso_marker_genes <- c("WNT2", "TCF21", # Alveolar FB
                       "CTHRC1", "POSTN", # Activated FB
                       "MFAP5", "SCARA5", # Adventitial FB
                       "PDGFRA", "FGF18", # Alveolar MyoFB
                       "WNT5A", "HHIP", # Ductal MyoFB
                       "COX4I2", "CSPG4", # Pericyte
                       "ITGA7", "NTRK3" # VSMC
)

meso_marker_plot <- DotPlot(meso_relabel,
                           features = meso_marker_genes,
                           group.by = "celltype",
                           dot.scale = 6)

meso_marker_plot_w_theme <- meso_marker_plot +
  scale_color_gradientn(colors = rev(color_glacier),
                        name = "Expression") +
  scale_y_discrete(limits = rev) +
  theme(axis.title = element_blank(),
        axis.text.x = element_text(face = "italic",
                                   angle = 45,
                                   hjust = 1,
                                   size = 12),
        aspect.ratio = 1,
        panel.background = element_rect(fill = "transparent",
                                        colour = "black"),
        plot.background = element_rect(color = "transparent",
                                       fill = "transparent"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12))

gridExtra::grid.arrange(egg::set_panel_size(p = meso_marker_plot_w_theme,
                                            width = unit(10, "cm"),
                                            height = unit(5, "cm")))

ggsave("./meso_marker_dotplot.pdf",
       egg::set_panel_size(p = meso_marker_plot_w_theme,
                           width = unit(10, "cm"),
                           height = unit(5, "cm")),
       dpi = 1200,
       bg = "transparent",
       width = unit(10, "cm"),
       height = unit(5, "cm"))

###################标志基因umap图###############

SCARA5_plot <- FeaturePlot(meso_relabel, "SCARA5") +
  featureplot_theme() +
  theme(aspect.ratio = 1,
        panel.background = element_rect(fill = "transparent",
                                        colour = "black"),
        plot.background = element_rect(color = "transparent",
                                       fill = "transparent")) +
  scale_color_gradientn(colors = rev(color_glacier),
                        name = "Expression")

gridExtra::grid.arrange(egg::set_panel_size(p = SCARA5_plot,
                                            width = unit(10, "cm"),
                                            height = unit(10, "cm")))

ggsave("./meso_SCARA5_plot.pdf",
       egg::set_panel_size(p = SCARA5_plot, width = unit(10, "cm"), height = unit(10, "cm")),
       dpi = 1200, bg = "transparent",
       width = unit(6, "cm"), height = unit(5, "cm"))

###################差异基因分析###############
###############1.差异基因upset绘制###########
Idents(meso_relabel) <- meso_relabel$condition
meso_relabel_no_preterm <- subset(meso_relabel,
                               idents = "Acute preterm injury",
                               invert = TRUE)
meso_relabel_no_preterm <- PrepSCTFindMarkers(meso_relabel_no_preterm)
celltype_for_volc_w_lineage <- c(factor(levels(meso_relabel_no_preterm$celltype)))

## BPD vs Term infant
markers_compare1 <- future.apply::future_lapply(
  celltype_for_volc_w_lineage,
  function(x) {
    # 设置细胞标识
    Idents(meso_relabel_no_preterm) <- meso_relabel_no_preterm$celltype
    meso_relabel_no_preterm_sub <- subset(meso_relabel_no_preterm, idents = x)
    Idents(meso_relabel_no_preterm_sub) <- meso_relabel_no_preterm_sub$condition
    
    # 切换到 RNA assay 并重新标准化
    DefaultAssay(meso_relabel_no_preterm_sub) <- "RNA"
    meso_relabel_no_preterm_sub <- NormalizeData(meso_relabel_no_preterm_sub, 
                                                 normalization.method = "LogNormalize", 
                                                 scale.factor = 10000)
    
    # 使用 RNA assay 进行差异表达分析
    foundmarkers <- FindMarkers(
      meso_relabel_no_preterm_sub,
      ident.1 = "BPD",
      ident.2 = "Term infant",
      assay = "RNA",
      test.use = "MAST",
      logfc.threshold = 0
    )
    
    # 定义上下调基因
    up_cat <- "up_in_bpd/down_in_term_infant"
    down_cat <- "up_in_term_infant/down_in_bpd"
    
    foundmarkers[up_cat] <- ifelse(
      foundmarkers$avg_log2FC > 1 & 
        foundmarkers$p_val_adj < 1e-05 & 
        !grepl("^RP|^MT-", rownames(foundmarkers)), 
      TRUE, FALSE
    )
    
    foundmarkers[down_cat] <- ifelse(
      foundmarkers$avg_log2FC < -1 & 
        foundmarkers$p_val_adj < 1e-05 & 
        !grepl("^RP|^MT-", rownames(foundmarkers)), 
      TRUE, FALSE
    )
    
    return(foundmarkers)
  }
)

names(markers_compare1) <- celltype_for_volc_w_lineage
wb_markers_names <- celltype_for_volc_w_lineage
wb_markers <- createWorkbook()

for (idx in seq_along(wb_markers_names)) {
  print(idx)
  addWorksheet(wb_markers,
               as.character(wb_markers_names[idx]))
  print(idx)
  writeData(wb_markers,
            as.character(wb_markers_names[idx]),
            markers_compare1[[idx]],
            colNames = TRUE,
            rowNames = TRUE)
}
saveWorkbook(wb_markers,
             file = paste0("./deg/S_excel_deg_bpd_term_infant.xlsx"),
             overwrite = TRUE)

## BPD+PH vs Term infant
markers_compare2 <- future.apply::future_lapply(
  celltype_for_volc_w_lineage,
  function(x) {
    # 设置细胞标识
    Idents(meso_relabel_no_preterm) <- meso_relabel_no_preterm$celltype
    meso_relabel_no_preterm_sub <- subset(meso_relabel_no_preterm, idents = x)
    Idents(meso_relabel_no_preterm_sub) <- meso_relabel_no_preterm_sub$condition
    
    # 检查细胞数量
    cell_counts <- table(Idents(meso_relabel_no_preterm_sub))
    bpd_ph_count <- cell_counts["BPD+PH"]
    term_count <- cell_counts["Term infant"]
    
    # 如果任一组的细胞数少于3，返回空数据框并给出警告
    if (is.na(bpd_ph_count) || bpd_ph_count < 3 || is.na(term_count) || term_count < 3) {
      warning(paste("Skipping", x, "- insufficient cells: BPD+PH =", 
                    ifelse(is.na(bpd_ph_count), 0, bpd_ph_count), 
                    "Term infant =", 
                    ifelse(is.na(term_count), 0, term_count)))
      return(data.frame())  # 返回空数据框
    }
    
    # 切换到 RNA assay 并重新标准化
    DefaultAssay(meso_relabel_no_preterm_sub) <- "RNA"
    meso_relabel_no_preterm_sub <- NormalizeData(meso_relabel_no_preterm_sub, 
                                                 normalization.method = "LogNormalize", 
                                                 scale.factor = 10000)
    
    # 使用 RNA assay 进行差异表达分析
    foundmarkers <- FindMarkers(
      meso_relabel_no_preterm_sub,
      ident.1 = "BPD+PH",
      ident.2 = "Term infant",
      assay = "RNA",
      test.use = "MAST",
      logfc.threshold = 0
    )
    
    # 定义上下调基因
    up_cat <- "up_in_bpdph/down_in_term_infant"
    down_cat <- "up_in_term_infant/down_in_bpdph"
    
    foundmarkers[up_cat] <- ifelse(
      foundmarkers$avg_log2FC > 1 & 
        foundmarkers$p_val_adj < 1e-05 & 
        !grepl("^RP|^MT-", rownames(foundmarkers)), 
      TRUE, FALSE
    )
    
    foundmarkers[down_cat] <- ifelse(
      foundmarkers$avg_log2FC < -1 & 
        foundmarkers$p_val_adj < 1e-05 & 
        !grepl("^RP|^MT-", rownames(foundmarkers)), 
      TRUE, FALSE
    )
    
    return(foundmarkers)
  }
)


names(markers_compare2) <- celltype_for_volc_w_lineage
wb_markers_names <- celltype_for_volc_w_lineage
wb_markers <- createWorkbook()

for (idx in seq_along(wb_markers_names)) {
  addWorksheet(wb_markers,
               as.character(wb_markers_names[idx]))
  writeData(wb_markers,
            wb_markers_names[idx],
            markers_compare2[[idx]],
            colNames = TRUE,
            rowNames = TRUE)
}
wb_markers
saveWorkbook(wb_markers,
             file = paste0("./deg/S_excel_deg_bpdph_term_infant.xlsx"),
             overwrite = TRUE)

##BPD+PH vs BPD
markers_compare3 <- future.apply::future_lapply(
  celltype_for_volc_w_lineage,
  function(x) {
    # 设置细胞标识
    Idents(meso_relabel_no_preterm) <- meso_relabel_no_preterm$celltype
    meso_relabel_no_preterm_sub <- subset(meso_relabel_no_preterm, idents = x)
    Idents(meso_relabel_no_preterm_sub) <- meso_relabel_no_preterm_sub$condition
    
    # 检查细胞数量
    cell_counts <- table(Idents(meso_relabel_no_preterm_sub))
    bpd_ph_count <- cell_counts["BPD+PH"]
    term_count <- cell_counts["BPD"]
    
    # 如果任一组的细胞数少于3，返回空数据框并给出警告
    if (is.na(bpd_ph_count) || bpd_ph_count < 3 || is.na(term_count) || term_count < 3) {
      warning(paste("Skipping", x, "- insufficient cells: BPD+PH =", 
                    ifelse(is.na(bpd_ph_count), 0, bpd_ph_count), 
                    "BPD =", 
                    ifelse(is.na(term_count), 0, term_count)))
      return(data.frame())  # 返回空数据框
    }
    
    # 切换到 RNA assay 并重新标准化
    DefaultAssay(meso_relabel_no_preterm_sub) <- "RNA"
    meso_relabel_no_preterm_sub <- NormalizeData(meso_relabel_no_preterm_sub, 
                                                 normalization.method = "LogNormalize", 
                                                 scale.factor = 10000)
    
    # 使用 RNA assay 进行差异表达分析
    foundmarkers <- FindMarkers(
      meso_relabel_no_preterm_sub,
      ident.1 = "BPD+PH",
      ident.2 = "BPD",
      assay = "RNA",
      test.use = "MAST",
      logfc.threshold = 0
    )
    
    # 定义上下调基因
    up_cat <- "up_in_bpdph/down_in_bpd"
    down_cat <- "up_in_bpd/down_in_bpdph"
    
    foundmarkers[up_cat] <- ifelse(
      foundmarkers$avg_log2FC > 1 & 
        foundmarkers$p_val_adj < 1e-05 & 
        !grepl("^RP|^MT-", rownames(foundmarkers)), 
      TRUE, FALSE
    )
    
    foundmarkers[down_cat] <- ifelse(
      foundmarkers$avg_log2FC < -1 & 
        foundmarkers$p_val_adj < 1e-05 & 
        !grepl("^RP|^MT-", rownames(foundmarkers)), 
      TRUE, FALSE
    )
    
    return(foundmarkers)
  }
)

names(markers_compare3) <- celltype_for_volc_w_lineage
wb_markers_names <- celltype_for_volc_w_lineage
wb_markers <- createWorkbook()


for (idx in seq_along(wb_markers_names)) {
  addWorksheet(wb_markers,
               as.character(wb_markers_names[idx]))
  writeData(wb_markers,
            as.character(wb_markers_names[idx]),
            markers_compare3[[idx]],
            colNames = TRUE,
            rowNames = TRUE)
}
saveWorkbook(wb_markers,
             file = paste0("./deg/S_excel_deg_bpdph_bpd.xlsx"),
             overwrite = TRUE)
skip_set_size_plot <- function(ups) {
  main <- ups$Main_bar
  matrix <- ups$Matrix
  ## find panel grob
  panel_idx <- grep("panel", main$layout$name, fixed = TRUE)
  ## find text grob
  text_idx <- which(
    vapply(main$grobs[[panel_idx]]$children, 
           \(x) inherits(x, "text"), 
           logical(1)))
  tG <- main$grobs[[panel_idx]]$children[[text_idx]]
  main$grobs[[panel_idx]]$children[[text_idx]] <- tG
  
  y_axis_idx <- grep("axis-l", matrix$layout$name)
  if (length(y_axis_idx) > 0) {
    matrix$grobs[y_axis_idx] <- NULL
    matrix$layout <- matrix$layout[-y_axis_idx, ]
  }
  
  # Find the panel grobs for both the main plot and the matrix plot
  main_panel_idx <- grep("panel", main$layout$name, fixed = TRUE)
  matrix_panel_idx <- grep("panel", matrix$layout$name, fixed = TRUE)
  
  # Align the widths of the panels
  main_panel <- main$grobs[[main_panel_idx]]
  matrix_panel <- matrix$grobs[[matrix_panel_idx]]
  max_width <- grid::unit.pmax(main_panel$widths, matrix_panel$widths)
  #   max_width <- sum(main_panel$widths)
  print(matrix_panel)
  main$grobs[[main_panel_idx]]$widths <- max_width
  matrix$grobs[[matrix_panel_idx]]$widths <- max_width
  
  grid.newpage()
  grid.draw(arrangeGrob(main, matrix, heights = ups$mb.ratio))
}

#绘制upset图
library(UpSetR)
library(ggplot2)

# 确保输出目录存在
dir.create("./deg", showWarnings = FALSE, recursive = TRUE)

# 为每个细胞类型绘制 UpSet 图
for(x in as.character(celltype_for_volc_w_lineage)) {
  
  # 检查该细胞类型在所有比较中是否有数据
  if (!x %in% names(markers_compare1) || 
      !x %in% names(markers_compare2) || 
      !x %in% names(markers_compare3)) {
    message(paste("Skipping", x, "- not present in all comparison results"))
    next
  }
  
  bpd_term <- markers_compare1[[x]]
  bpdph_term <- markers_compare2[[x]]
  bpdph_bpd <- markers_compare3[[x]]
  
  # 检查数据框是否为空或没有行
  if (nrow(bpd_term) == 0 || nrow(bpdph_term) == 0 || nrow(bpdph_bpd) == 0) {
    message(paste("Skipping", x, "- empty results in one or more comparisons"))
    next
  }
  
  # 安全地提取上下调基因
  bpd_term_up <- tryCatch({
    if ("up_in_bpd/down_in_term_infant" %in% colnames(bpd_term)) {
      rownames(bpd_term[bpd_term[["up_in_bpd/down_in_term_infant"]] == TRUE, ])
    } else {
      character(0)
    }
  }, error = function(e) character(0))
  
  bpd_term_down <- tryCatch({
    if ("up_in_term_infant/down_in_bpd" %in% colnames(bpd_term)) {
      rownames(bpd_term[bpd_term[["up_in_term_infant/down_in_bpd"]] == TRUE, ])
    } else {
      character(0)
    }
  }, error = function(e) character(0))
  
  bpdph_term_up <- tryCatch({
    if ("up_in_bpdph/down_in_term_infant" %in% colnames(bpdph_term)) {
      rownames(bpdph_term[bpdph_term[["up_in_bpdph/down_in_term_infant"]] == TRUE, ])
    } else {
      character(0)
    }
  }, error = function(e) character(0))
  
  bpdph_term_down <- tryCatch({
    if ("up_in_term_infant/down_in_bpdph" %in% colnames(bpdph_term)) {
      rownames(bpdph_term[bpdph_term[["up_in_term_infant/down_in_bpdph"]] == TRUE, ])
    } else {
      character(0)
    }
  }, error = function(e) character(0))
  
  bpdph_bpd_up <- tryCatch({
    if ("up_in_bpdph/down_in_bpd" %in% colnames(bpdph_bpd)) {
      rownames(bpdph_bpd[bpdph_bpd[["up_in_bpdph/down_in_bpd"]] == TRUE, ])
    } else {
      character(0)
    }
  }, error = function(e) character(0))
  
  bpdph_bpd_down <- tryCatch({
    if ("up_in_bpd/down_in_bpdph" %in% colnames(bpdph_bpd)) {
      rownames(bpdph_bpd[bpdph_bpd[["up_in_bpd/down_in_bpdph"]] == TRUE, ])
    } else {
      character(0)
    }
  }, error = function(e) character(0))
  
  # 创建 UpSet 列表
  upset_list <- list(
    "Up in BPD/Down in Term" = bpd_term_up,
    "Up in Term/Down in BPD" = bpd_term_down,
    "Up in BPD+PH/Down in Term" = bpdph_term_up,
    "Up in Term/Down in BPD+PH" = bpdph_term_down,
    "Up in BPD+PH/Down in BPD" = bpdph_bpd_up,
    "Up in BPD/Down in BPD+PH" = bpdph_bpd_down
  )
  
  # 检查是否有任何基因集是非空的
  non_empty_sets <- sapply(upset_list, length) > 0
  if (sum(non_empty_sets) < 2) {
    message(paste("Skipping", x, "- insufficient non-empty gene sets for UpSet plot"))
    next
  }
  
  # 只保留非空的集合
  upset_list_filtered <- upset_list[non_empty_sets]
  
  # 绘制 UpSet 图
  pdf_file <- paste0("./deg/Figure_upset_lineage_", gsub("/", "_", x), ".pdf")
  
  tryCatch({
    pdf(file = pdf_file, width = 12, height = 8, onefile = FALSE)
    
    upset_plot <- upset(
      fromList(upset_list_filtered),
      sets = rev(names(upset_list_filtered)),
      keep.order = TRUE,
      sets.bar.color = "#56B4E9",
      order.by = "freq",
      empty.intersections = "on",
      number.angles = 30,
      point.size = 3.5,
      line.size = 1,
      mainbar.y.label = "Gene Intersections",
      sets.x.label = "Genes per Set",
      text.scale = c(1.3, 1.3, 1, 1, 1.5, 1)
    )
    
    print(upset_plot)
    dev.off()
    
    message(paste("Successfully created UpSet plot for", x))
    
  }, error = function(e) {
    dev.off()
    message(paste("Failed to create UpSet plot for", x, ":", e$message))
  })
}

###########2.多组细胞差异基因火山图绘制###############
library(ggplot2)
library(ggthemes)
library(ggrepel)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(magrittr)
library(purrr)

prepare_volcano_data <- function(markers_list, comparison_name) {
  volcano_data <- data.frame()
  rm_genes <- c("FP236383.3", "FP671120.4")  # 定义要移除的特定基因
  
  for (cell_type in names(markers_list)) {
    if (!is.null(markers_list[[cell_type]]) && nrow(markers_list[[cell_type]]) > 0) {
      markers <- markers_list[[cell_type]]
      
      # 在添加其他列之前先移除不需要的基因
      # 1. 移除特定基因
      markers <- markers[!rownames(markers) %in% rm_genes, ]
      
      # 2. 移除以"RP"或"MT-"开头的基因（更严格的模式）
      # ^RP[SL]? 匹配以RPS或RPL开头的基因（核糖体蛋白）
      # ^MT- 匹配以MT-开头的基因（线粒体基因）
      markers <- markers[!grepl("^RP[SL]?|^MT-", rownames(markers)), ]
      
      # 如果移除后数据框为空，跳过这个细胞类型
      if (nrow(markers) == 0) {
        next
      }
      
      markers$celltype <- cell_type
      markers$comparison <- comparison_name
      markers$gene <- rownames(markers)
      
      # 添加调节方向
      markers$regulate <- case_when(
        markers$p_val_adj < 0.05 & markers$avg_log2FC >= 1 ~ 'Up',
        markers$p_val_adj < 0.05 & markers$avg_log2FC <= -1 ~ 'Down',
        TRUE ~ 'Not'
      )
      
      volcano_data <- rbind(volcano_data, markers)
    }
  }
  
  return(volcano_data)
}


########2.1 BPD vs Term infant###########
setwd("/data/nas1/xiangxiaoyu_OD/1.BPD+MSCs/1.scRNA/04_Adventitial_FB/01_DEGs")
volcano_data1 <- prepare_volcano_data(markers_compare1, "BPD_vs_Term")

# 过滤掉不显著的基因
volcano_filtered_1 <- volcano_data1 %>% 
  filter(regulate != 'Not') %>%
  arrange(celltype)
# Activated FB Adventitial FB    Alveolar FB Alveolar MyoFB   Ductal MyoFB       Pericyte           VSMC 
# 54           1646           2234            637           1524            154            760
write.csv(file = "all_cells_BPDvsTerminfant.csv",volcano_filtered_1)

# 设置细胞类型顺序（如果需要）
# 按照你想要的顺序定义细胞类型
custom_order <- c("Alveolar FB", "Activated FB", "Adventitial FB", 
                  "Alveolar MyoFB", "Ductal MyoFB", "Pericyte", "VSMC")

# 应用自定义顺序
volcano_filtered_1$celltype <- factor(volcano_filtered_1$celltype, levels = custom_order)

# 选择每个细胞类型中p值最小的5个基因进行标记
label_genes_1 <- volcano_filtered_1 %>%
  group_by(celltype) %>%
  slice_min(p_val_adj, n = 5) %>%
  ungroup()

# 计算每个细胞类型的logFC范围用于背景条
bar_ranges <- sapply(c('min', 'max'), function(x) {
  if (x == 'min') {
    sapply(levels(volcano_filtered_1$celltype), function(cell_type) {
      mm <- volcano_filtered_1 %>% filter(celltype == cell_type)
      if (nrow(mm) > 0) min(mm$avg_log2FC, na.rm = TRUE) * 1.2 else 0
    })
  } else {
    sapply(levels(volcano_filtered_1$celltype), function(cell_type) {
      mm <- volcano_filtered_1 %>% filter(celltype == cell_type)
      if (nrow(mm) > 0) max(mm$avg_log2FC, na.rm = TRUE) * 1.2 else 0
    })
  }
}) %>% as.data.frame() %>% rownames_to_column('celltype')

# 创建火山图
volcano_plot <- ggplot(volcano_filtered_1, aes(x = celltype, y = avg_log2FC)) +
  
  # 背景条
  geom_col(data = bar_ranges,
           mapping = aes(x = celltype, y = min),
           fill = "#dcdcdc", alpha = 0.6) +
  geom_col(data = bar_ranges,
           mapping = aes(x = celltype, y = max),
           fill = "#dcdcdc", alpha = 0.6) +
  
  # 点图
  geom_jitter(
    mapping = aes(size = abs(avg_log2FC),
                  alpha = abs(avg_log2FC),
                  color = regulate),
    width = 0.2) +
  scale_size_continuous(range=c(1,5))+
  # 基因标签
  geom_text_repel(
    data = label_genes_1,
    mapping = aes(label = gene),
    size = 3,
    color = 'black',
    box.padding = 0.25,
    max.overlaps = Inf) +
  
  # 颜色设置
  scale_color_manual(
    values = c(
      'Down' = "#3288bd",
      'Up' = '#d92523'),
    name = "Significant") +
  
  # 坐标轴和标签
  labs(
    y = expression(log[2]~FC),
    x = '',
    color = 'Significant',
    title = "DEG Across Cell Types: BPD vs Term Infant") +
  
  # 图例设置
  guides(
    color = guide_legend(override.aes = list(size = 3))) +
  
  # 细胞类型标签背景
  geom_tile(
    data = volcano_filtered_1 %>%
      dplyr::select(celltype) %>%
      dplyr::distinct(celltype),
    mapping = aes(x = celltype, y = 0),
    height = 0.6,
    color = colorRampPalette(brewer.pal(8, 'Accent'))(length(unique(volcano_filtered_1$celltype))),
    fill = colorRampPalette(brewer.pal(8, 'Accent'))(length(unique(volcano_filtered_1$celltype))),
    alpha = 0.6,
    show.legend = FALSE) +
  
  # 细胞类型标签文字
  geom_text(
    data = volcano_filtered_1 %>%
      dplyr::select(celltype) %>%
      dplyr::distinct(celltype),
    mapping = aes(x = celltype, y = 0, label = celltype),
    size = 4) +
  
  # 主题设置
  theme_few() +
  theme(
    axis.text = element_text(color = "#000000", size = 12),
    axis.text.x = element_blank(),
    axis.title = element_text(color = "#000000", size = 15),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.background = element_rect(color = "#969696", fill = "white"))

# 显示图形
print(volcano_plot)

# 保存图形
ggsave("./deg/Figure_BPD_vs_Term_volcano.pdf",volcano_plot,width = 14,height = 8,dpi = 1200)


# 美化版本的火山图
beautified_volcano <- ggplot(volcano_filtered_1, aes(x = celltype, y = avg_log2FC)) +
  
  # 背景条 - 使用渐变色
  geom_col(data = bar_ranges,
           mapping = aes(x = celltype, y = min),
           fill = "grey95", alpha = 0.8) +
  geom_col(data = bar_ranges,
           mapping = aes(x = celltype, y = max),
           fill = "grey95", alpha = 0.8) +
  
  # 点图 - 使用更美观的参数
  geom_jitter(
    mapping = aes(size = -log10(p_val_adj),  # 大小根据p值显著性
                  alpha = abs(avg_log2FC),   # 透明度根据FC绝对值
                  color = regulate),
    width = 0.25,
    height = 0
  ) +
  scale_size_continuous(range = c(1,10))+
  # 基因标签 - 只标记最显著的基因
  geom_text_repel(
    data = label_genes_1,
    mapping = aes(label = gene),
    size = 2.8,
    box.padding = 0.35,
    point.padding = 0.3,
    max.overlaps = 20,
    min.segment.length = 0.2,
    segment.color = "grey50",
    segment.size = 0.3,
    show.legend = FALSE
  ) +
  
  # 颜色设置 - 使用更专业的配色
  scale_color_manual(
    values = c(
      'Down' = "#3288bd",
      'Up' = '#d92523'
    ),
    name = "Regulation",
    labels = c("Downregulated", "Upregulated")
  ) +
  
  # 大小比例 - 根据p值调整点的大小范围
  scale_size_continuous(
    range = c(1, 4),
    name = "-Log10 P-value"
  ) +
  
  # 透明度比例
  scale_alpha_continuous(
    range = c(0.4, 0.9),
    guide = "none"  # 隐藏透明度图例
  ) +
  
  # 坐标轴和标签
  labs(
    y = expression("Log"[2]~"Fold Change"),
    x = '',
    title = "Differential Gene Expression: BPD vs Term Infant",
    subtitle = "Colored by regulation direction, sized by significance"
  ) +
  
  # 细胞类型标签背景 - 使用更柔和的颜色
  geom_tile(
    data = volcano_filtered_1 %>% distinct(celltype),
    mapping = aes(x = celltype, y = 0),
    height = 0.8,
    fill = "white",
    color = "grey70",
    alpha = 0.8,
    linewidth = 0.5,
    show.legend = FALSE
  ) +
  
  # 细胞类型标签文字
  geom_text(
    data = volcano_filtered_1 %>% distinct(celltype),
    mapping = aes(x = celltype, y = 0, label = celltype),
    size = 3.5,
    fontface = "bold",
    color = "black"
  ) +
  
  # 主题设置 - 使用更专业的主题
  theme_minimal() +
  theme(
    axis.text = element_text(color = "black", size = 10),
    axis.text.x = element_blank(),
    axis.title = element_text(color = "black", size = 12, face = "bold"),
    axis.title.y = element_text(margin = margin(r = 10)),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 10)),
    plot.subtitle = element_text(size = 11, color = "grey40", hjust = 0.5, margin = margin(b = 15)),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(color = "grey90", size = 0.3),
    legend.position = "right",
    legend.box = "vertical",
    legend.margin = margin(l = 10),
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 9)
  )

# 显示美化后的图形
print(beautified_volcano)

# 保存美化后的图形
ggsave("./deg/Figure_BPD_vs_Term_volcano_beautified.pdf",
  beautified_volcano, width = 16,height = 9,dpi = 1200)





########2.2 BPD+PH vs Term infant###########
# 只使用BPD+PH vs Term infant的数据
volcano_data2 <- prepare_volcano_data(markers_compare2, "BPD+PH_vs_Term") 

# 过滤掉不显著的基因
volcano_filtered_2 <- volcano_data2 %>% 
  filter(regulate != 'Not') %>%
  arrange(celltype)
# Activated FB Adventitial FB    Alveolar FB   Ductal MyoFB       Pericyte           VSMC 
# 1081           1713           4730            479           1860            363
write.csv(file = "all_cells_BPD+PHvsTerminfant.csv",volcano_filtered_2)

# 设置细胞类型顺序（如果需要）
custom_order <- c("Alveolar FB", "Activated FB", "Adventitial FB", 
                  "Alveolar MyoFB", "Ductal MyoFB", "Pericyte", "VSMC")

# 应用自定义顺序
volcano_filtered_2$celltype <- factor(volcano_filtered_2$celltype, levels = custom_order)

# 选择每个细胞类型中p值最小的5个基因进行标记
label_genes_2 <- volcano_filtered_2 %>%
  group_by(celltype) %>%
  slice_min(p_val_adj, n = 5) %>%
  ungroup()

# 计算每个细胞类型的logFC范围用于背景条
bar_ranges <- sapply(c('min', 'max'), function(x) {
  if (x == 'min') {
    sapply(levels(volcano_filtered_2$celltype), function(cell_type) {
      mm <- volcano_filtered_2 %>% filter(celltype == cell_type)
      if (nrow(mm) > 0) min(mm$avg_log2FC, na.rm = TRUE) * 1.2 else 0
    })
  } else {
    sapply(levels(volcano_filtered_2$celltype), function(cell_type) {
      mm <- volcano_filtered_2 %>% filter(celltype == cell_type)
      if (nrow(mm) > 0) max(mm$avg_log2FC, na.rm = TRUE) * 1.2 else 0
    })
  }
}) %>% as.data.frame() %>% rownames_to_column('celltype')

# 创建火山图
volcano_plot <- ggplot(volcano_filtered_2, aes(x = celltype, y = avg_log2FC)) +
  
  # 背景条
  geom_col(data = bar_ranges,
           mapping = aes(x = celltype, y = min),
           fill = "#dcdcdc", alpha = 0.6) +
  geom_col(data = bar_ranges,
           mapping = aes(x = celltype, y = max),
           fill = "#dcdcdc", alpha = 0.6) +
  
  # 点图
  geom_jitter(
    mapping = aes(size = abs(avg_log2FC),
                  alpha = abs(avg_log2FC),
                  color = regulate),
    width = 0.2) +
  scale_size_continuous(range=c(1,5))+
  # 基因标签
  geom_text_repel(
    data = label_genes_2,
    mapping = aes(label = gene),
    size = 3,
    color = 'black',
    box.padding = 0.25,
    max.overlaps = Inf) +
  
  # 颜色设置
  scale_color_manual(
    values = c(
      'Down' = "#E7A19E",
      'Up' = '#6A3D9A'),
    name = "Significant") +
  
  # 坐标轴和标签
  labs(
    y = expression(log[2]~FC),
    x = '',
    color = 'Significant',
    title = "DEG Across Cell Types: BPD vs Term Infant") +
  
  # 图例设置
  guides(
    color = guide_legend(override.aes = list(size = 3))) +
  
  # 细胞类型标签背景
  geom_tile(
    data = volcano_filtered_2 %>%
      dplyr::select(celltype) %>%
      dplyr::distinct(celltype),
    mapping = aes(x = celltype, y = 0),
    height = 0.6,
    color = colorRampPalette(brewer.pal(8, 'Accent'))(length(unique(volcano_filtered_2$celltype))),
    fill = colorRampPalette(brewer.pal(8, 'Accent'))(length(unique(volcano_filtered_2$celltype))),
    alpha = 0.6,
    show.legend = FALSE) +
  
  # 细胞类型标签文字
  geom_text(
    data = volcano_filtered_2 %>%
      dplyr::select(celltype) %>%
      dplyr::distinct(celltype),
    mapping = aes(x = celltype, y = 0, label = celltype),
    size = 4) +
  
  # 主题设置
  theme_few() +
  theme(
    axis.text = element_text(color = "#000000", size = 12),
    axis.text.x = element_blank(),
    axis.title = element_text(color = "#000000", size = 15),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.background = element_rect(color = "#969696", fill = "white"))

# 显示图形
print(volcano_plot)

# 保存图形
ggsave("./deg/Figure_BPDPH_vs_Term_volcano.pdf",volcano_plot,width = 14,height = 8,dpi = 1200)


# 美化版本的火山图
beautified_volcano <- ggplot(volcano_filtered_2, aes(x = celltype, y = avg_log2FC)) +
  
  # 背景条 - 使用渐变色
  geom_col(data = bar_ranges,
           mapping = aes(x = celltype, y = min),
           fill = "grey95", alpha = 0.8) +
  geom_col(data = bar_ranges,
           mapping = aes(x = celltype, y = max),
           fill = "grey95", alpha = 0.8) +
  
  # 点图 - 使用更美观的参数
  geom_jitter(
    mapping = aes(size = -log10(p_val_adj),  # 大小根据p值显著性
                  alpha = abs(avg_log2FC),   # 透明度根据FC绝对值
                  color = regulate),
    width = 0.25,
    height = 0
  ) +
  scale_size_continuous(range = c(1,10))+
  # 基因标签 - 只标记最显著的基因
  geom_text_repel(
    data = label_genes_2,
    mapping = aes(label = gene),
    size = 2.8,
    box.padding = 0.35,
    point.padding = 0.3,
    max.overlaps = 20,
    min.segment.length = 0.2,
    segment.color = "grey50",
    segment.size = 0.3,
    show.legend = FALSE
  ) +
  
  # 颜色设置 - 使用更专业的配色
  scale_color_manual(
    values = c(
      'Down' = "#E7A19E",
      'Up' = '#6A3D9A'),
    name = "Regulation",
    labels = c("Downregulated", "Upregulated")
  ) +
  
  # 大小比例 - 根据p值调整点的大小范围
  scale_size_continuous(
    range = c(1, 4),
    name = "-Log10 P-value"
  ) +
  
  # 透明度比例
  scale_alpha_continuous(
    range = c(0.4, 0.9),
    guide = "none"  # 隐藏透明度图例
  ) +
  
  # 坐标轴和标签
  labs(
    y = expression("Log"[2]~"Fold Change"),
    x = '',
    title = "Differential Gene Expression: BPD+PH vs Term Infant",
    subtitle = "Colored by regulation direction, sized by significance"
  ) +
  
  # 细胞类型标签背景 - 使用更柔和的颜色
  geom_tile(
    data = volcano_filtered_2 %>% distinct(celltype),
    mapping = aes(x = celltype, y = 0),
    height = 0.8,
    fill = "white",
    color = "grey70",
    alpha = 0.8,
    linewidth = 0.5,
    show.legend = FALSE
  ) +
  
  # 细胞类型标签文字
  geom_text(
    data = volcano_filtered_2 %>% distinct(celltype),
    mapping = aes(x = celltype, y = 0, label = celltype),
    size = 3.5,
    fontface = "bold",
    color = "black"
  ) +
  
  # 主题设置 - 使用更专业的主题
  theme_minimal() +
  theme(
    axis.text = element_text(color = "black", size = 10),
    axis.text.x = element_blank(),
    axis.title = element_text(color = "black", size = 12, face = "bold"),
    axis.title.y = element_text(margin = margin(r = 10)),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 10)),
    plot.subtitle = element_text(size = 11, color = "grey40", hjust = 0.5, margin = margin(b = 15)),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(color = "grey90", size = 0.3),
    legend.position = "right",
    legend.box = "vertical",
    legend.margin = margin(l = 10),
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 9)
  )

# 显示美化后的图形
print(beautified_volcano)

# 保存美化后的图形
ggsave("./deg/Figure_BPDPH_vs_Term_volcano_beautified.pdf",
       beautified_volcano, width = 16,height = 9,dpi = 1200)


########2.3 BPD+PH vs BPD###########
volcano_data3 <- prepare_volcano_data(markers_compare3, "BPD+PH_vs_BPD")

# 过滤掉不显著的基因
volcano_filtered_3 <- volcano_data3 %>% 
  filter(regulate != 'Not') %>%
  arrange(celltype)
# Activated FB Adventitial FB    Alveolar FB   Ductal MyoFB       Pericyte           VSMC 
# 538           1670           1490            210             60             29
write.csv(file = "all_cells_BPD+PHvsBPD.csv",volcano_filtered_3)

# 设置细胞类型顺序（如果需要）
# 设置细胞类型顺序（如果需要）
custom_order <- c("Alveolar FB", "Activated FB", "Adventitial FB", 
                  "Alveolar MyoFB", "Ductal MyoFB", "Pericyte", "VSMC")

# 应用自定义顺序
volcano_filtered_3$celltype <- factor(volcano_filtered_3$celltype, levels = custom_order)

# 选择每个细胞类型中p值最小的5个基因进行标记
label_genes_3 <- volcano_filtered_3 %>%
  group_by(celltype) %>%
  slice_min(p_val_adj, n = 5) %>%
  ungroup()

# 计算每个细胞类型的logFC范围用于背景条
bar_ranges <- sapply(c('min', 'max'), function(x) {
  if (x == 'min') {
    sapply(levels(volcano_filtered_3$celltype), function(cell_type) {
      mm <- volcano_filtered_3 %>% filter(celltype == cell_type)
      if (nrow(mm) > 0) min(mm$avg_log2FC, na.rm = TRUE) * 1.2 else 0
    })
  } else {
    sapply(levels(volcano_filtered_3$celltype), function(cell_type) {
      mm <- volcano_filtered_3 %>% filter(celltype == cell_type)
      if (nrow(mm) > 0) max(mm$avg_log2FC, na.rm = TRUE) * 1.2 else 0
    })
  }
}) %>% as.data.frame() %>% rownames_to_column('celltype')

# 创建火山图
volcano_plot <- ggplot(volcano_filtered_3, aes(x = celltype, y = avg_log2FC)) +
  
  # 背景条
  geom_col(data = bar_ranges,
           mapping = aes(x = celltype, y = min),
           fill = "#dcdcdc", alpha = 0.6) +
  geom_col(data = bar_ranges,
           mapping = aes(x = celltype, y = max),
           fill = "#dcdcdc", alpha = 0.6) +
  
  # 点图
  geom_jitter(
    mapping = aes(size = abs(avg_log2FC),
                  alpha = abs(avg_log2FC),
                  color = regulate),
    width = 0.2) +
  scale_size_continuous(range=c(1,5))+
  # 基因标签
  geom_text_repel(
    data = label_genes_3,
    mapping = aes(label = gene),
    size = 3,
    color = 'black',
    box.padding = 0.25,
    max.overlaps = Inf) +
  
  # 颜色设置
  scale_color_manual(
    values = c(
      'Down' = "#d92523",
      'Up' = '#6A3D9A'),
    name = "Significant") +
  
  # 坐标轴和标签
  labs(
    y = expression(log[2]~FC),
    x = '',
    color = 'Significant',
    title = "DEG Across Cell Types: BPD vs Term Infant") +
  
  # 图例设置
  guides(
    color = guide_legend(override.aes = list(size = 3))) +
  
  # 细胞类型标签背景
  geom_tile(
    data = volcano_filtered_3 %>%
      dplyr::select(celltype) %>%
      dplyr::distinct(celltype),
    mapping = aes(x = celltype, y = 0),
    height = 0.6,
    color = colorRampPalette(brewer.pal(8, 'Accent'))(length(unique(volcano_filtered_3$celltype))),
    fill = colorRampPalette(brewer.pal(8, 'Accent'))(length(unique(volcano_filtered_3$celltype))),
    alpha = 0.6,
    show.legend = FALSE) +
  
  # 细胞类型标签文字
  geom_text(
    data = volcano_filtered_3 %>%
      dplyr::select(celltype) %>%
      dplyr::distinct(celltype),
    mapping = aes(x = celltype, y = 0, label = celltype),
    size = 4) +
  
  # 主题设置
  theme_few() +
  theme(
    axis.text = element_text(color = "#000000", size = 12),
    axis.text.x = element_blank(),
    axis.title = element_text(color = "#000000", size = 15),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.background = element_rect(color = "#969696", fill = "white"))

# 显示图形
print(volcano_plot)

# 保存图形
ggsave("./deg/Figure_BPDPH_vs_BPD_volcano.pdf",volcano_plot,width = 14,height = 8,dpi = 1200)


# 美化版本的火山图
beautified_volcano <- ggplot(volcano_filtered_3, aes(x = celltype, y = avg_log2FC)) +
  
  # 背景条 - 使用渐变色
  geom_col(data = bar_ranges,
           mapping = aes(x = celltype, y = min),
           fill = "grey95", alpha = 0.8) +
  geom_col(data = bar_ranges,
           mapping = aes(x = celltype, y = max),
           fill = "grey95", alpha = 0.8) +
  
  # 点图 - 使用更美观的参数
  geom_jitter(
    mapping = aes(size = -log10(p_val_adj),  # 大小根据p值显著性
                  alpha = abs(avg_log2FC),   # 透明度根据FC绝对值
                  color = regulate),
    width = 0.25,
    height = 0
  ) +
  # 基因标签 - 只标记最显著的基因
  geom_text_repel(
    data = label_genes_3,
    mapping = aes(label = gene),
    size = 2.8,
    box.padding = 0.35,
    point.padding = 0.3,
    max.overlaps = 20,
    min.segment.length = 0.2,
    segment.color = "grey50",
    segment.size = 0.3,
    show.legend = FALSE
  ) +
  
  # 颜色设置 - 使用更专业的配色
  scale_color_manual(
    values = c(
      'Down' = "#d92523",
      'Up' = '#6A3D9A'),
    name = "Regulation",
    labels = c("Downregulated", "Upregulated")
  ) +
  
  # 大小比例 - 根据p值调整点的大小范围
  scale_size_continuous(
    range = c(1, 4),
    name = "-Log10 P-value"
  ) +
  
  # 透明度比例
  scale_alpha_continuous(
    range = c(0.4, 0.9),
    guide = "none"  # 隐藏透明度图例
  ) +
  
  # 坐标轴和标签
  labs(
    y = expression("Log"[2]~"Fold Change"),
    x = '',
    title = "Differential Gene Expression: BPD vs Term Infant",
    subtitle = "Colored by regulation direction, sized by significance"
  ) +
  
  # 细胞类型标签背景 - 使用更柔和的颜色
  geom_tile(
    data = volcano_filtered_3 %>% distinct(celltype),
    mapping = aes(x = celltype, y = 0),
    height = 0.8,
    fill = "white",
    color = "grey70",
    alpha = 0.8,
    linewidth = 0.5,
    show.legend = FALSE
  ) +
  
  # 细胞类型标签文字
  geom_text(
    data = volcano_filtered_3 %>% distinct(celltype),
    mapping = aes(x = celltype, y = 0, label = celltype),
    size = 3.5,
    fontface = "bold",
    color = "black"
  ) +
  
  # 主题设置 - 使用更专业的主题
  theme_minimal() +
  theme(
    axis.text = element_text(color = "black", size = 10),
    axis.text.x = element_blank(),
    axis.title = element_text(color = "black", size = 12, face = "bold"),
    axis.title.y = element_text(margin = margin(r = 10)),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 10)),
    plot.subtitle = element_text(size = 11, color = "grey40", hjust = 0.5, margin = margin(b = 15)),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(color = "grey90", size = 0.3),
    legend.position = "right",
    legend.box = "vertical",
    legend.margin = margin(l = 10),
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 9)
  )

# 显示美化后的图形
print(beautified_volcano)

# 保存美化后的图形
ggsave("./deg/Figure_BPDPH_vs_BPD_volcano_beautified.pdf",
       beautified_volcano, width = 16,height = 9,dpi = 1200)
