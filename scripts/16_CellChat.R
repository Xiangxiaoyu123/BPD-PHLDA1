###############################################################################
# Cell-cell communication between adventitial fibroblast subclusters and lung cell types
#
# Source      : r.03_Cell_prop.r (lines 2977-4189)
# Input       : adventitial FB + epithelial/endothelial cells
# Output      : Fig. 4A-F; Supplementary Fig. 7; Supplementary Table 11
#
# Code below is unchanged from the original analysis script.
###############################################################################
# -----------------------------------------------------------------------------
# PATHS: this file keeps the absolute paths of the analysis server on which the
# study was run (/data/nas1/xiangxiaoyu_OD/1.BPD+MSCs/1.scRNA/... and
# /home/xxyloveme/study/BPD+MSCs/...). Set PROJ_DIR in config/paths.R and
# replace those prefixes with your own project directory before running.
# -----------------------------------------------------------------------------


#############3.4 Adventitial FB:细胞通讯分析#############
library(Matrix)
library(reticulate)
library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)

alldata_neonatal_injury <- readRDS(file = "/data/nas1/xiangxiaoyu_OD/1.BPD+MSCs/1.scRNA/03_All_cells/alldata_neonatal_injury.rds")
Idents(alldata_neonatal_injury) <- "celltype"
subbedseurat <- subset(alldata_neonatal_injury, idents =  c("gCap","aCap","abCap","AT1","AT2"))

##剔除数据
Idents(subbedseurat) <- subbedseurat$condition
subbedseurat_no_preterm <- subset(subbedseurat,
                                  idents = "Acute preterm injury",
                                  invert = TRUE)

cellchat_seurat <- Reduce(function(x, y) {merge(x, y)},
                        list(subbedseurat_no_preterm, meso_adventitial_FB))
cellchat_seurat@meta.data <- cellchat_seurat@meta.data[,-c(12:15,18,19,21:24,26:29)]

# 综合处理方法
cellchat_seurat$cell_subtype <- ifelse(
  # 情况1：有聚类分配且是Adventitial FB
  !is.na(cellchat_seurat$SCT_snn_res.0.5) & cellchat_seurat$celltype == "Adventitial FB",
  paste0("Adventitial_FB_Cluster", cellchat_seurat$SCT_snn_res.0.5),
  
  # 情况2：有聚类分配但不是Adventitial FB（保持原celltype）
  ifelse(
    !is.na(cellchat_seurat$SCT_snn_res.0.5) & cellchat_seurat$celltype != "Adventitial FB",
    cellchat_seurat$celltype,
    
    # 情况3：没有聚类分配（保持原celltype）
    cellchat_seurat$celltype
  )
)

Idents(cellchat_seurat) <- "condition"

cellchat_seurat_normal <- subset(cellchat_seurat, idents = c("Term infant"))
Idents(cellchat_seurat_normal) <- "cell_subtype"
cellchat_seurat_normal <- subset(cellchat_seurat_normal, 
                                 idents = c("Adventitial_FB_Cluster5"),
                                 invert = TRUE)

cellchat_seurat_bpd <- subset(cellchat_seurat, idents = c("BPD"))
Idents(cellchat_seurat_bpd) <- "cell_subtype"
cellchat_seurat_bpd <- subset(cellchat_seurat_bpd, 
                                 idents = c("Adventitial_FB_Cluster5"),
                                 invert = TRUE)

cellchat_seurat_bpdph <- subset(cellchat_seurat, idents = c("BPD+PH"))


#
cellchat.normal <- createCellChat(object = cellchat_seurat_normal,
                                  group.by = "cell_subtype",
                                  assay = "SCT")
cellchat.bpd <- createCellChat(object = cellchat_seurat_bpd,
                               group.by = "cell_subtype",
                               assay = "SCT")
cellchat.bpdph <- createCellChat(object = cellchat_seurat_bpdph,
                                 group.by = "cell_subtype",
                                 assay = "SCT")
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)

CellChatDB.use <- CellChatDB
cellchat.normal@DB <- CellChatDB.use
cellchat.bpd@DB <- CellChatDB.use
cellchat.bpdph@DB <- CellChatDB.use

cellchat.normal <- subsetData(cellchat.normal)
cellchat.bpd <- subsetData(cellchat.bpd)
cellchat.bpdph <- subsetData(cellchat.bpdph)

future::plan("multicore", workers = 16) # do parallel

cellchat.normal <- identifyOverExpressedGenes(cellchat.normal)
cellchat.normal <- identifyOverExpressedInteractions(cellchat.normal)
cellchat.bpd <- identifyOverExpressedGenes(cellchat.bpd)
cellchat.bpd <- identifyOverExpressedInteractions(cellchat.bpd)
cellchat.bpdph <- identifyOverExpressedGenes(cellchat.bpdph)
cellchat.bpdph <- identifyOverExpressedInteractions(cellchat.bpdph)

cellchat.normal <- computeCommunProb(cellchat.normal,
                                     k.min = 5,
                                     type = "truncatedMean")
cellchat.normal <- filterCommunication(cellchat.normal,
                                       min.cells = 5)

cellchat.bpd <- computeCommunProb(cellchat.bpd,
                                  k.min = 5,
                                  type = "truncatedMean")
cellchat.bpd <- filterCommunication(cellchat.bpd,
                                    min.cells = 5)

cellchat.bpdph <- computeCommunProb(cellchat.bpdph,
                                    k.min = 5,
                                    type = "truncatedMean")
cellchat.bpdph <- filterCommunication(cellchat.bpdph,
                                      min.cells = 5)
saveRDS(cellchat.normal, "./07_CellChat/cellchat_normal_0424version.RDS")
saveRDS(cellchat.bpd, "./07_CellChat/cellchat_bpd_0424version.RDS")
saveRDS(cellchat.bpdph, "./07_CellChat/cellchat_bpd_ph_0424version.RDS")


cellchat.normal <- computeCommunProbPathway(cellchat.normal)
cellchat.bpd <- computeCommunProbPathway(cellchat.bpd)
cellchat.bpdph <- computeCommunProbPathway(cellchat.bpdph)

cellchat.normal <- aggregateNet(cellchat.normal)
cellchat.bpd <- aggregateNet(cellchat.bpd)
cellchat.bpdph <- aggregateNet(cellchat.bpdph)

options(repr.plot.width = 15,
        repr.plot.height = 8)
groupSize <- as.numeric(table(cellchat.normal@idents))
par(mfrow = c(1, 2),
    xpd = TRUE)
a <- netVisual_circle(cellchat.normal@net$count,
                      vertex.weight = groupSize,
                      weight.scale = TRUE,
                      label.edge = FALSE,
                      title.name = "Term infant: Number of interactions")
b <- netVisual_circle(cellchat.normal@net$weight,
                      vertex.weight = groupSize,
                      weight.scale = TRUE,
                      label.edge = FALSE,
                      title.name = "Term infant: Interaction weights/strength")

options(repr.plot.width = 15,
        repr.plot.height = 8)
groupSize <- as.numeric(table(cellchat.bpd@idents))
par(mfrow = c(1, 2),
    xpd = TRUE)
c <- netVisual_circle(cellchat.bpd@net$count,
                      vertex.weight = groupSize,
                      weight.scale = TRUE,
                      label.edge = FALSE,
                      title.name = "BPD: Number of interactions")
d <- netVisual_circle(cellchat.bpd@net$weight,
                      vertex.weight = groupSize,
                      weight.scale = TRUE,
                      label.edge = FALSE,
                      title.name = "BPD: Interaction weights/strength")

options(repr.plot.width = 15,
        repr.plot.height = 8)
groupSize <- as.numeric(table(cellchat.bpdph@idents))
par(mfrow = c(1, 2),
    xpd = TRUE)
e <- netVisual_circle(cellchat.bpdph@net$count,
                      vertex.weight = groupSize,
                      weight.scale = TRUE,
                      label.edge = FALSE,
                      title.name = "BPD+PH: Number of interactions")
f <- netVisual_circle(cellchat.bpdph@net$weight,
                      vertex.weight = groupSize,
                      weight.scale = TRUE,
                      label.edge = FALSE,
                      title.name = "BPD+PH: Interaction weights/strength")


# 定义我们关注的细胞类型
adventitial_fb_clusters <- paste0("Adventitial_FB_Cluster", 0:5)
epithelial_cells <- c("AT1", "AT2")
endothelial_cells <- c("gCap", "aCap", "abCap")
target_celltypes <- c(adventitial_fb_clusters, epithelial_cells, endothelial_cells)

# 整合的通讯网络可视化函数
visualize_condition_communication_integrated <- function(cellchat_obj, condition_name) {
  # 创建输出子目录
  condition_dir <- paste0("./07_CellChat/", condition_name)
  dir.create(condition_dir, showWarnings = FALSE)
  
  # 1. 总体通讯网络 - 整合两个图到一张图像
  groupSize <- as.numeric(table(cellchat_obj@idents))
  
  # 设置图形参数
  options(repr.plot.width = 15, repr.plot.height = 8)
  par(mfrow = c(1, 2), xpd = TRUE)
  
  # 通讯数量图
  netVisual_circle(cellchat_obj@net$count,
                   vertex.weight = groupSize,
                   weight.scale = TRUE,
                   label.edge = FALSE,
                   title.name = paste(condition_name, ": Number of interactions"))
  
  # 通讯权重图
  netVisual_circle(cellchat_obj@net$weight,
                   vertex.weight = groupSize,
                   weight.scale = TRUE,
                   label.edge = FALSE,
                   title.name = paste(condition_name, ": Interaction weights/strength"))
  
  # 保存整合的图像
  dev.copy(pdf, paste0(condition_dir, "/integrated_network_plots.pdf"), width = 15, height = 8)
  dev.off()
  
  # 重置图形参数
  par(mfrow = c(1, 1))
  
  # 2. Adventitial FB与上皮/内皮细胞的通讯气泡图
  # FB到上皮
  if(length(intersect(adventitial_fb_clusters, levels(cellchat_obj@idents))) > 0 &&
     length(intersect(epithelial_cells, levels(cellchat_obj@idents))) > 0) {
    netVisual_bubble(cellchat_obj, 
                     sources.use = adventitial_fb_clusters,
                     targets.use = epithelial_cells,
                     remove.isolate = FALSE,
                     title.name = paste("Adventitial FB to Epithelial -", condition_name))
    
    ggsave(paste0(condition_dir, "/FB_to_epithelial_bubble.pdf"), width = 10, height = 6)
  }
  
  # FB到内皮
  if(length(intersect(adventitial_fb_clusters, levels(cellchat_obj@idents))) > 0 &&
     length(intersect(endothelial_cells, levels(cellchat_obj@idents))) > 0) {
    netVisual_bubble(cellchat_obj, 
                     sources.use = adventitial_fb_clusters,
                     targets.use = endothelial_cells,
                     remove.isolate = FALSE,
                     title.name = paste("Adventitial FB to Endothelial -", condition_name))
    
    ggsave(paste0(condition_dir, "/FB_to_endothelial_bubble.pdf"), width = 10, height = 6)
  }
  
  # 3. 通路水平分析
  pathways_show <- cellchat_obj@netP$pathways
  if(length(pathways_show) > 0) {
    # 肺相关通路
    lung_related_pathways <- c("TGFb", "WNT", "FGF", "VEGF", "BMP", "NOTCH", "EGF", "PDGF")
    
    for(pathway in lung_related_pathways) {
      if(pathway %in% pathways_show) {
        # 气泡图
        p_bubble <- netVisual_bubble(cellchat_obj, 
                                     signaling = pathway,
                                     sources.use = adventitial_fb_clusters,
                                     targets.use = c(epithelial_cells, endothelial_cells),
                                     remove.isolate = FALSE,
                                     title.name = paste(pathway, "signaling -", condition_name))
        
        ggsave(paste0(condition_dir, "/", pathway, "_bubble.pdf"), 
               p_bubble, width = 12, height = 8)
      }
    }
  }
  
  return(cellchat_obj)
}

# 对每个条件运行整合的可视化函数
visualize_condition_communication_integrated(cellchat.normal, "Term_infant")
visualize_condition_communication_integrated(cellchat.bpd, "BPD")
visualize_condition_communication_integrated(cellchat.bpdph, "BPD_PH")

# 比较对照组和BPD组的通路差异
compare_pathways_between_conditions <- function(cellchat_normal, cellchat_bpd) {
  # 获取两个条件下的显著通路
  pathways_normal <- cellchat_normal@netP$pathways
  pathways_bpd <- cellchat_bpd@netP$pathways
  
  cat("=== Pathway Comparison: Normal vs BPD ===\n\n")
  
  # 1. 共同通路
  common_pathways <- intersect(pathways_normal, pathways_bpd)
  cat("Common pathways in both conditions:", length(common_pathways), "\n")
  print(common_pathways)
  
  # 2. 对照组特有通路
  normal_specific <- setdiff(pathways_normal, pathways_bpd)
  cat("\nPathways specific to Normal condition:", length(normal_specific), "\n")
  print(normal_specific)
  
  # 3. BPD组特有通路
  bpd_specific <- setdiff(pathways_bpd, pathways_normal)
  cat("\nPathways specific to BPD condition:", length(bpd_specific), "\n")
  print(bpd_specific)
  
  # 4. 通路强度比较
  cat("\n=== Pathway Strength Comparison ===\n")
  
  # 计算每个通路的总体通讯强度
  compute_pathway_strength <- function(cellchat_obj) {
    pathways <- cellchat_obj@netP$pathways
    strength <- sapply(pathways, function(pathway) {
      sum(cellchat_obj@netP$prob[,,pathway], na.rm = TRUE)
    })
    return(data.frame(pathway = pathways, strength = strength))
  }
  
  strength_normal <- compute_pathway_strength(cellchat_normal)
  strength_bpd <- compute_pathway_strength(cellchat_bpd)
  
  # 合并结果
  strength_comparison <- merge(strength_normal, strength_bpd, 
                               by = "pathway", 
                               suffixes = c("_normal", "_bpd"),
                               all = TRUE)
  
  # 计算差异
  strength_comparison$strength_diff <- strength_comparison$strength_bpd - strength_comparison$strength_normal
  strength_comparison$abs_diff <- abs(strength_comparison$strength_diff)
  
  # 按差异大小排序
  strength_comparison <- strength_comparison[order(-strength_comparison$abs_diff), ]
  
  cat("Top pathways with largest differences:\n")
  print(head(strength_comparison, 10))
  
  return(list(
    common_pathways = common_pathways,
    normal_specific = normal_specific,
    bpd_specific = bpd_specific,
    strength_comparison = strength_comparison
  ))
}

# 执行比较分析
pathway_comparison <- compare_pathways_between_conditions(cellchat.normal, cellchat.bpd)

# === Pathway Comparison: Normal vs BPD ===
#   
#   Common pathways in both conditions: 99 
# [1] "COLLAGEN"      "LAMININ"       "FN1"           "THBS"          "TENASCIN"      "CypA"         
# [7] "MK"            "CD99"          "SLIT"          "ADGRL"         "Cholesterol"   "MIF"          
# [13] "VEGF"          "EGF"           "APP"           "SEMA3"         "ADGRG"         "JAM"          
# [19] "VISFATIN"      "CADM"          "PTPRM"         "ANGPTL"        "FLRT"          "FGF"          
# [25] "RA"            "GAS"           "GAP"           "IGFBP"         "Prostaglandin" "SEMA4"        
# [31] "ADGRE"         "PECAM1"        "UNC5"          "PROS"          "BMP"           "PDGF"         
# [37] "EDN"           "EPHA"          "NEGR"          "IL6"           "CDH"           "TWEAK"        
# [43] "WNT"           "PTPR"          "MPZ"           "CDH5"          "GDF"           "CDH1"         
# [49] "NRG"           "NOTCH"         "NECTIN"        "ESAM"          "TGFb"          "SEMA6"        
# [55] "DHEAS"         "HSPG"          "OCLN"          "PTN"           "CXCL"          "ADGRA"        
# [61] "CD46"          "ncWNT"         "GRN"           "CALCR"         "IGF"           "CNTN"         
# [67] "LIFR"          "SEMA5"         "VCAM"          "Netrin"        "NRXN"          "EPHB"         
# [73] "CLDN"          "ANGPT"         "NCAM"          "HGF"           "Testosterone"  "RELN"         
# [79] "AGRN"          "GALECTIN"      "Desmosterol"   "DHEA"          "LIGHT"         "PERIOSTIN"    
# [85] "KIT"           "NGL"           "IL1"           "SEMA7"         "DESMOSOME"     "NT"           
# [91] "Glutamate"     "CEACAM"        "DHT"           "CSF3"          "TRAIL"         "SLITRK"       
# [97] "TULP"          "APELIN"        "PCDH"         
# 
# Pathways specific to Normal condition: 8 
# [1] "MHC-II"  "CCL"     "NPNT"    "IL16"    "PLAU"    "HH"      "ACTIVIN" "MAG"    
# 
# Pathways specific to BPD condition: 6 
# [1] "MMP"      "VWF"      "THY1"     "SELL"     "CHEMERIN" "CD96"    
# 
# === Pathway Strength Comparison ===
#   Top pathways with largest differences:
#   pathway strength_normal strength_bpd strength_diff   abs_diff
# 26  COLLAGEN       76.795747   52.4417624   -24.3539843 24.3539843
# 59   LAMININ       23.499153   15.5389084    -7.9602446  7.9602446
# 103     THBS        6.610345    3.1638950    -3.4464499  3.4464499
# 42       FN1       10.246820    7.0764218    -3.1703977  3.1703977
# 29      CypA        2.623405    0.3372214    -2.2861836  2.2861836
# 65        MK        2.541242    0.5706917    -1.9705504  1.9705504
# 17      CD99        2.216852    0.6107118    -1.6061402  1.6061402
# 100 TENASCIN        2.764314    1.6509475    -1.1133664  1.1133664
# 4      ADGRG        1.046065    0.1824179    -0.8636467  0.8636467


# 可视化差异最大的通路
visualize_top_differential_pathways <- function(cellchat_normal, cellchat_bpd, pathway_comparison, top_n = 5) {
  top_pathways <- head(pathway_comparison$strength_comparison$pathway, top_n)
  
  cat("\n=== Visualizing Top", top_n, "Differential Pathways ===\n")
  
  for(pathway in top_pathways) {
    cat("Visualizing pathway:", pathway, "\n")
    
    # 创建输出目录
    pathway_dir <- paste0("./07_CellChat/differential_pathways/", pathway)
    dir.create(pathway_dir, recursive = TRUE, showWarnings = FALSE)
    
    # 1. 气泡图比较
    tryCatch({
      # Normal condition
      p1 <- netVisual_bubble(cellchat_normal, 
                             signaling = pathway,
                             sources.use = adventitial_fb_clusters,
                             targets.use = c(epithelial_cells, endothelial_cells),
                             remove.isolate = FALSE,
                             title.name = paste(pathway, "- Normal"))
      
      # BPD condition
      p2 <- netVisual_bubble(cellchat_bpd, 
                             signaling = pathway,
                             sources.use = adventitial_fb_clusters,
                             targets.use = c(epithelial_cells, endothelial_cells),
                             remove.isolate = FALSE,
                             title.name = paste(pathway, "- BPD"))
      
      # 保存图片
      ggsave(paste0(pathway_dir, "/", pathway, "_bubble_normal.pdf"), p1, width = 10, height = 6)
      ggsave(paste0(pathway_dir, "/", pathway, "_bubble_bpd.pdf"), p2, width = 10, height = 6)
    }, error = function(e) {
      cat("Error in bubble plot for", pathway, ":", e$message, "\n")
    })
    
    # 2. 热图比较
    tryCatch({
      # Normal condition
      pdf(paste0(pathway_dir, "/", pathway, "_heatmap_normal.pdf"), width = 8, height = 6)
      netVisual_heatmap(cellchat_normal, 
                        signaling = pathway,
                        color.heatmap = "Reds",
                        title.name = paste(pathway, "- Normal"))
      dev.off()
      
      # BPD condition
      pdf(paste0(pathway_dir, "/", pathway, "_heatmap_bpd.pdf"), width = 8, height = 6)
      netVisual_heatmap(cellchat_bpd, 
                        signaling = pathway,
                        color.heatmap = "Reds",
                        title.name = paste(pathway, "- BPD"))
      dev.off()
    }, error = function(e) {
      cat("Error in heatmap for", pathway, ":", e$message, "\n")
    })
  }
}

# 可视化前5个差异最大的通路
visualize_top_differential_pathways(cellchat.normal, cellchat.bpd, pathway_comparison, 5)


# 分析BPD特异性通路的生物学意义
analyze_bpd_specific_pathways <- function(bpd_specific_pathways) {
  cat("\n=== Biological Interpretation of BPD-Specific Pathways ===\n")
  
  # 定义通路类别
  pathway_categories <- list(
    fibrosis_related = c("TGFb", "BMP", "PDGF", "COLLAGEN"),
    inflammation_related = c("IL", "TNF", "COMPLEMENT", "CHEMERIN"),
    angiogenesis_related = c("VEGF", "ANGPT", "SPHINGOSINE", "EDN"),
    stress_response = c("HIF", "NRG", "PERIOSTIN"),
    development_related = c("WNT", "FGF", "NOTCH", "EPHA", "SEMA")
  )
  
  for(category in names(pathway_categories)) {
    category_pathways <- c()
    for(pattern in pathway_categories[[category]]) {
      matched <- grep(pattern, bpd_specific_pathways, value = TRUE, ignore.case = TRUE)
      category_pathways <- c(category_pathways, matched)
    }
    
    if(length(category_pathways) > 0) {
      cat("\n", toupper(category), "PATHWAYS in BPD:\n")
      cat("  -", paste(unique(category_pathways), collapse = "\n  - "), "\n")
      
      # 生物学解释
      if(category == "fibrosis_related") {
        cat("    Potential role: Fibrotic remodeling in BPD pathogenesis\n")
      } else if(category == "inflammation_related") {
        cat("    Potential role: Chronic inflammation driving disease progression\n")
      } else if(category == "angiogenesis_related") {
        cat("    Potential role: Impaired vascular development in BPD\n")
      } else if(category == "stress_response") {
        cat("    Potential role: Cellular stress responses to injury\n")
      } else if(category == "development_related") {
        cat("    Potential role: Dysregulated developmental signaling\n")
      }
    }
  }
}

# 分析BPD特异性通路的生物学意义
if(length(pathway_comparison$bpd_specific) > 0) {
  analyze_bpd_specific_pathways(pathway_comparison$bpd_specific)
}
# === Biological Interpretation of BPD-Specific Pathways ===
#   
#   INFLAMMATION_RELATED PATHWAYS in BPD:
#   - CHEMERIN 
# Potential role: Chronic inflammation driving disease progression


# 生成差异通路报告
generate_differential_pathway_report <- function(pathway_comparison) {
  report_file <- "./07_CellChat/differential_pathway_analysis_report.txt"
  
  sink(report_file)
  
  cat("DIFFERENTIAL PATHWAY ANALYSIS REPORT: Normal vs BPD\n")
  cat("===================================================\n\n")
  
  cat("SUMMARY STATISTICS:\n")
  cat("-------------------\n")
  cat("Total pathways in Normal condition:", length(pathway_comparison$common_pathways) + length(pathway_comparison$normal_specific), "\n")
  cat("Total pathways in BPD condition:", length(pathway_comparison$common_pathways) + length(pathway_comparison$bpd_specific), "\n")
  cat("Common pathways:", length(pathway_comparison$common_pathways), "\n")
  cat("Normal-specific pathways:", length(pathway_comparison$normal_specific), "\n")
  cat("BPD-specific pathways:", length(pathway_comparison$bpd_specific), "\n\n")
  
  cat("TOP 10 PATHWAYS BY STRENGTH DIFFERENCE:\n")
  cat("---------------------------------------\n")
  print(head(pathway_comparison$strength_comparison, 10))
  
  cat("\nBPD-SPECIFIC PATHWAYS:\n")
  cat("----------------------\n")
  if(length(pathway_comparison$bpd_specific) > 0) {
    cat(paste(pathway_comparison$bpd_specific, collapse = "\n"))
  } else {
    cat("No BPD-specific pathways found")
  }
  
  cat("\n\nNORMAL-SPECIFIC PATHWAYS:\n")
  cat("-------------------------\n")
  if(length(pathway_comparison$normal_specific) > 0) {
    cat(paste(pathway_comparison$normal_specific, collapse = "\n"))
  } else {
    cat("No Normal-specific pathways found")
  }
  
  cat("\n\nKEY FINDINGS:\n")
  cat("-------------\n")
  cat("1. BPD-specific pathways may represent disease-driving mechanisms\n")
  cat("2. Lost pathways in BPD may indicate impaired protective mechanisms\n")
  cat("3. Pathway strength differences suggest altered signaling intensity\n")
  cat("4. These findings provide targets for further mechanistic studies\n")
  
  sink()
  
  cat("Differential pathway report saved to:", report_file, "\n")
}

# 生成详细报告
generate_differential_pathway_report(pathway_comparison)

# 保存比较结果
saveRDS(pathway_comparison, "./07_CellChat/differential_pathway_analysis_report/pathway_comparison_normal_vs_bpd.rds")




# 创建差异通路的汇总表格
create_pathway_summary_table <- function(pathway_comparison) {
  summary_table <- pathway_comparison$strength_comparison
  
  # 添加通路类别信息
  summary_table$category <- "Other"
  
  # 分类通路
  fibrosis_pattern <- "TGFb|BMP|PDGF|COLLAGEN"
  inflammation_pattern <- "IL|TNF|COMPLEMENT|CHEMERIN"
  angiogenesis_pattern <- "VEGF|ANGPT|SPHINGOSINE|EDN"
  stress_pattern <- "HIF|NRG|PERIOSTIN"
  development_pattern <- "WNT|FGF|NOTCH|EPHA|SEMA"
  
  summary_table$category[grep(fibrosis_pattern, summary_table$pathway)] <- "Fibrosis"
  summary_table$category[grep(inflammation_pattern, summary_table$pathway)] <- "Inflammation"
  summary_table$category[grep(angiogenesis_pattern, summary_table$pathway)] <- "Angiogenesis"
  summary_table$category[grep(stress_pattern, summary_table$pathway)] <- "Stress Response"
  summary_table$category[grep(development_pattern, summary_table$pathway)] <- "Development"
  
  # 保存表格
  write.csv(summary_table, "./07_CellChat/differential_pathway_analysis_report/pathway_comparison_summary_table.csv", row.names = FALSE)
  
  return(summary_table)
}

# 创建汇总表格
pathway_summary <- create_pathway_summary_table(pathway_comparison)


##制作图片
cellchat.normal <- netAnalysis_computeCentrality(cellchat.normal,
                                                 slot.name = "netP")
cellchat.bpd <- netAnalysis_computeCentrality(cellchat.bpd,
                                              slot.name = "netP")
cellchat.bpdph <- netAnalysis_computeCentrality(cellchat.bpdph,
                                                slot.name = "netP")
func <- function(pathway) {
  # Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
  pdf(file = paste0("./07_CellChat/data_0424version/",
                    "major_signaling_roles_celltype_",
                    pathway,
                    "_normal.pdf"),
      width = 8,
      height = 10)
  a <- netAnalysis_signalingRole_network(cellchat.normal,
                                         signaling = pathway,
                                         width = 8,
                                         height = 2.5,
                                         font.size = 10)
  print(a)
  dev.off()
}
lapply(mainpathways, func)

#normal
mainpathways <- cellchat.normal@netP$pathway
pdf(file = paste0("./07_CellChat/data_0424version/heatmap_of_outgoing_signals_per_celltype_", "normal.pdf"), width = 15, height = 20)
ht1 <- netAnalysis_signalingRole_heatmap(cellchat.normal, signaling = mainpathways, 
                                         pattern = "outgoing", width = 20, height = 40, 
                                         font.size = 14, font.size.title = 20)
ComplexHeatmap::draw(ht1)
dev.off()

pdf(file = paste0("./07_CellChat/data_0424version/heatmap_of_incoming_signals_per_celltype_", "normal.pdf"), width = 15, height = 20)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat.normal, signaling = mainpathways, 
                                         pattern = "incoming", width = 20, height = 40, 
                                         font.size = 14, font.size.title = 20)
ComplexHeatmap::draw(ht2)
dev.off()

#BPD
mainpathwaysbpd <- cellchat.bpd@netP$pathways
pdf(file = paste0("./07_CellChat/data_0424version/heatmap_of_outgoing_signals_per_celltype_", "bpd.pdf"), width = 15, height = 20)
ht1 <- netAnalysis_signalingRole_heatmap(cellchat.bpd, signaling = mainpathwaysbpd, 
                                         pattern = "outgoing", width = 20, height = 40, 
                                         font.size = 14, font.size.title = 20)
ComplexHeatmap::draw(ht1)
dev.off()

pdf(file = paste0("./07_CellChat/data_0424version/heatmap_of_incoming_signals_per_celltype_", "bpd.pdf"), width = 15, height = 20)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat.bpd, signaling = mainpathwaysbpd, 
                                         pattern = "incoming", width = 20, height = 40, 
                                         font.size = 14, font.size.title = 20)
ComplexHeatmap::draw(ht2)
dev.off()


setdiff(mainpathwaysbpd, mainpathways)
#[1] "MMP"      "VWF"      "THY1"     "SELL"     "CHEMERIN" "CD96" 


ht1 <- netAnalysis_signalingRole_heatmap(cellchat.bpd, signaling = c("MMP","CD96","VWF","THY1","SELL","CHEMERIN"), 
                                         pattern = "incoming", width = 10, height = 10, font.size = 14, 
                                         font.size.title = 20)

ht2 <- netAnalysis_signalingRole_heatmap(cellchat.bpd, signaling = c("MMP","CD96","VWF","THY1","SELL","CHEMERIN"), 
                                         pattern = "outgoing", width = 10, height = 10, font.size = 14, 
                                         font.size.title = 20)

#选择通路
recommended_pathways = c(
  'COLLAGEN',     # 最大差异，纤维化核心
  'LAMININ',      # 基底膜成分，与组织结构相关
  'FN1',          # 纤维连接蛋白，细胞外基质
  'THBS',         # 血管生成，与gCap相关
  'CypA',       # 血管生成因子
  'MK',          # 表皮生长因子，与发育相关
  'CD99',        # 信号导向，发育相关
  'TENASCIN'          # 上调通路，可能有特殊意义
)

ht1 <- netAnalysis_signalingRole_heatmap(cellchat.normal, signaling = recommended_pathways, 
                                         pattern = "outgoing", width = 10, height = 4, font.size = 14, 
                                         font.size.title = 20)
ht1 <- netAnalysis_signalingRole_heatmap(cellchat.normal, signaling = recommended_pathways, 
                                         pattern = "incoming", width = 10, height = 4, font.size = 14, 
                                         font.size.title = 20)

ht2 <- netAnalysis_signalingRole_heatmap(cellchat.bpd, signaling = recommended_pathways, 
                                         pattern = "outgoing", width = 10, height = 4, font.size = 14, 
                                         font.size.title = 20)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat.bpd, signaling = recommended_pathways, 
                                         pattern = "incoming", width = 10, height = 4, font.size = 14, 
                                         font.size.title = 20)

######寻找cluster与AT1之间的信号通路
adventitial_fb_clusters <- paste0("Adventitial_FB_Cluster", 0:5)
epithelial_cells <- c("AT1", "AT2")
endothelial_cells <- c("gCap", "aCap", "abCap")
target_celltypes <- c(adventitial_fb_clusters, epithelial_cells, endothelial_cells)

setwd("/data/nas1/xiangxiaoyu_OD/1.BPD+MSCs/1.scRNA/04_Adventitial_FB/07_CellChat/BPD_specific_pathways/")
p1 <- netVisual_bubble(cellchat.normal, 
                      sources.use = "Adventitial_FB_Cluster4",
                      targets.use = "AT2",
                      remove.isolate = FALSE,
                      title.name = "epithelial_cells",
                      signaling = c("THBS"))

p2 <- netVisual_bubble(cellchat.bpd, 
                 sources.use = "Adventitial_FB_Cluster2",
                 targets.use = "AT1",
                 remove.isolate = FALSE,
                 title.name = "Adventitial FB Cluster to AT1/AT2",
                 signaling =c("COLLAGEN"))
ggsave("FB_to_epithelial_bubble.pdf", width = 6, height = 14,plot = p)

##重新和并绘制
# 加载必要的包
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# 从CellChat对象中提取通信数据
extract_communication_data <- function(cellchat_obj, condition_name, sources.use, targets.use, signaling) {
  # 提取通信数据
  comm_data <- subsetCommunication(cellchat_obj, 
                                   sources.use = sources.use, 
                                   targets.use = targets.use, 
                                   signaling = signaling)
  
  if(is.null(comm_data) || nrow(comm_data) == 0) {
    return(NULL)
  }
  
  # 添加条件信息
  comm_data$Condition <- condition_name
  comm_data$Source_Target <- paste(sources.use, targets.use, sep = " -> ")
  
  return(comm_data)
}

# 提取两个条件的THBS信号数据
normal_data <- extract_communication_data(
  cellchat.normal, 
  "Normal", 
  "Adventitial_FB_Cluster4", 
  "AT2", 
  "THBS"
)

bpd_data <- extract_communication_data(
  cellchat.bpd, 
  "BPD", 
  "Adventitial_FB_Cluster2", 
  "AT1", 
  "THBS"
)

# 合并数据
combined_data <- bind_rows(normal_data, bpd_data)

# 检查数据
print(head(combined_data))

# 如果数据为空，尝试使用其他方法提取
if(is.null(combined_data) || nrow(combined_data) == 0) {
  # 尝试使用netAnalysis_computeCentrality提取数据
  normal_data <- cellchat.normal@net[["prob"]]
  bpd_data <- cellchat.bpd@net[["prob"]]
  
  # 这里需要根据您的具体数据结构进行调整
  # 您可能需要手动提取THBS信号相关的数据
}

# 创建自定义气泡图
create_custom_bubble_plot <- function(data, title = "THBS Signaling") {
  # 确保数据存在且不为空
  if(is.null(data) || nrow(data) == 0) {
    stop("没有找到通信数据")
  }
  
  # 创建气泡图
  p <- ggplot(data, aes(x = interaction_name, y = Source_Target)) +
    geom_point(aes(size = prob, color = prob), alpha = 0.8) +
    scale_size_continuous(
      name = "Interaction\nProbability",
      range = c(2, 10),
      limits = c(0, max(data$prob, na.rm = TRUE))
    ) +
    scale_color_gradientn(
      colors = c("blue", "yellow", "red"),
      name = "Interaction\nProbability"
    ) +
    facet_grid(Condition ~ ., scales = "free_y", space = "free_y") +
    labs(
      title = title,
      x = "Ligand-Receptor Pairs",
      y = "Cell Type Interactions"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      axis.text.x = element_text(angle = 90,hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      axis.title = element_text(size = 12, face = "bold"),
      strip.text = element_text(size = 12, face = "bold"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "right",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 8))
  
  return(p)
}

# 创建合并的气泡图
combined_plot <- create_custom_bubble_plot(combined_data, "THBS Signaling in Normal and BPD Conditions")

# 显示图形
print(combined_plot)

# 保存图形
ggsave("THBS_signaling_combined_plot.pdf", combined_plot, width = 12, height = 3.5)

#ADGRL 在normal AT1中没有，在bpd中AT1中存在
p1 <- netVisual_bubble(cellchat.normal, 
                       sources.use = adventitial_fb_clusters,
                       targets.use = epithelial_cells,
                       remove.isolate = FALSE,
                       title.name = "Adventitial FB Cluster to AT1/AT2",
                       signaling = c("ADGRL"))
p2 <- netVisual_bubble(cellchat.bpd, 
                       sources.use = adventitial_fb_clusters,
                       targets.use = epithelial_cells,
                       remove.isolate = FALSE,
                       title.name = "Adventitial FB Cluster to AT1/AT2",
                       signaling = c("ADGRL"))
pdf("FB_to_epithelial_bubble.pdf", width = 7, height = 6)
p1/p2
dev.off()

gg1 <- netAnalysis_signalingRole_scatter(cellchat.normal, signaling = c("ADGRL"))
gg2 <- netAnalysis_signalingRole_scatter(cellchat.bpd, signaling = c("ADGRL"))
pdf("netAnalysis_signalingRole_scatter_ADGRL.pdf", width = 12, height = 5)
gg1 + gg2
dev.off()

############细胞通讯分析继续研究ADGRL信号通路######
aggregateNet(cellchat.bpdph, signaling = c("ADGRL"), return.object = FALSE, remove.isolate = TRUE)

object.normal = cellchat.normal
object.bpd = cellchat.bpd
object.bpdph = cellchat.bpdph
color <- c("#E41A1C","#377EB8","#4DAF4A","#F29403","#F781BF","#BC9DCC","#A65628","#54B0E4","#222F75","#1B9E77")
color.use = color
color.use.bpdph = c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#F29403","#F781BF","#BC9DCC","#A65628","#54B0E4","#222F75","#1B9E77")
color.use.bpdph.graph = c("#E41A1C","#377EB8","#4DAF4A","#F29403","#F781BF","#BC9DCC","#A65628","#54B0E4","#222F75","#1B9E77","#984EA3")
slot.name = "netP"
group = NULL
weight.MinMax = NULL
dot.size = c(2, 6)
point.shape = c(21, 22, 24, 23, 25, 8, 3)
label.size = 3
dot.alpha = 0.6
x.measure = "outdeg"
y.measure = "indeg"
xlabel = "Outgoing interaction strength"
ylabel = "Incoming interaction strength"
title = NULL
font.size = 10
font.size.title = 10
do.label = TRUE
show.legend = TRUE
show.axes = TRUE


centr.normal <- slot(object.normal, slot.name)$centr
centr.bpd <- slot(object.bpd, slot.name)$centr
centr.bpdph <- slot(object.bpdph, slot.name)$centr

outgoing.normal <- matrix(0, nrow = nlevels(object.normal@idents), ncol = length(centr.normal))
incoming.normal <- matrix(0, nrow = nlevels(object.normal@idents), ncol = length(centr.normal))
dimnames(outgoing.normal) <- list(levels(object.normal@idents), names(centr.normal))
dimnames(incoming.normal) <- dimnames(outgoing.normal)

outgoing.bpd <- matrix(0, nrow = nlevels(object.bpd@idents), ncol = length(centr.bpd))
incoming.bpd <- matrix(0, nrow = nlevels(object.bpd@idents), ncol = length(centr.bpd))
dimnames(outgoing.bpd) <- list(levels(object.bpd@idents), names(centr.bpd))
dimnames(incoming.bpd) <- dimnames(outgoing.bpd)

outgoing.bpdph <- matrix(0, nrow = nlevels(object.bpdph@idents), ncol = length(centr.bpdph))
incoming.bpdph <- matrix(0, nrow = nlevels(object.bpdph@idents), ncol = length(centr.bpdph))
dimnames(outgoing.bpdph) <- list(levels(object.bpdph@idents), names(centr.bpdph))
dimnames(incoming.bpdph) <- dimnames(outgoing.bpdph)


for (i in 1:length(centr.normal)) {
  outgoing.normal[,i] <- centr.normal[[i]][[x.measure]]
  incoming.normal[,i] <- centr.normal[[i]][[y.measure]]
}
for (i in 1:length(centr.bpd)) {
  outgoing.bpd[,i] <- centr.bpd[[i]][[x.measure]]
  incoming.bpd[,i] <- centr.bpd[[i]][[y.measure]]
}
for (i in 1:length(centr.bpdph)) {
  outgoing.bpdph[,i] <- centr.bpdph[[i]][[x.measure]]
  incoming.bpdph[,i] <- centr.bpdph[[i]][[y.measure]]
}

outgoing.cells.normal <- rowSums(outgoing.normal)
incoming.cells.normal <- rowSums(incoming.normal)

outgoing.cells.bpd <- rowSums(outgoing.bpd)
incoming.cells.bpd <- rowSums(incoming.bpd)

outgoing.cells.bpdph <- rowSums(outgoing.bpdph)
incoming.cells.bpdph <- rowSums(incoming.bpdph)


num.link.normal <- aggregateNet(object.normal, signaling = NULL, return.object = FALSE, remove.isolate = FALSE)$count
num.link.normal <- rowSums(num.link.normal) + colSums(num.link.normal)-diag(num.link.normal)

num.link.bpd <- aggregateNet(object.bpd, signaling = NULL, return.object = FALSE, remove.isolate = FALSE)$count
num.link.bpd <- rowSums(num.link.bpd) + colSums(num.link.bpd)-diag(num.link.bpd)

num.link.bpdph <- aggregateNet(object.bpdph, signaling = NULL, return.object = FALSE, remove.isolate = FALSE)$count
num.link.bpdph <- rowSums(num.link.bpdph) + colSums(num.link.bpdph)-diag(num.link.bpdph)


df.normal <- data.frame(x = outgoing.cells.normal, y = incoming.cells.normal, labels = names(incoming.cells.normal),
                        Count = num.link.normal)

df.bpd <- data.frame(x = outgoing.cells.bpd, y = incoming.cells.bpd, labels = names(incoming.cells.bpd),
                     Count = num.link.bpd)

df.bpdph <- data.frame(x = outgoing.cells.bpdph, y = incoming.cells.bpdph, labels = names(incoming.cells.bpdph),
                       Count = num.link.bpdph)

df.normal$labels <- factor(df.normal$labels, levels = names(incoming.cells.normal))

df.bpd$labels <- factor(df.bpd$labels, levels = names(incoming.cells.bpd))

df.bpdph$labels <- factor(df.bpdph$labels, levels = names(incoming.cells.bpdph))

combined <- rbind(df.normal, df.bpd, df.bpdph)


gg <- ggplot() +
  geom_point(data = df.normal, aes(x, y, colour = labels), size = 6) +
  scale_colour_manual(values = color.use, drop = FALSE) +
  labs(colour = "Normal") +
  ggnewscale::new_scale("colour") +
  geom_point(data = df.bpd, aes(x, y, colour = labels), size = 6, shape=17) +
  scale_colour_manual(values = color.use, drop = FALSE) +
  labs(colour = "BPD") +
  ggnewscale::new_scale("colour") +
  geom_point(data = df.bpdph, aes(x, y, colour = labels), size = 6, shape=15) +
  scale_colour_manual(values = color.use.bpdph, drop = FALSE) +
  labs(colour = "BPD + PH")+theme_bw()

gg <- gg + CellChat_theme_opts() +
  theme(legend.key.height = grid::unit(0.15, "in"),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 12)) +
  labs(title = title, x = xlabel, y = ylabel) + theme(plot.title = element_text(size= font.size.title, face="plain"))+
  theme(axis.line.x = element_line(size = 0.25), axis.line.y = element_line(size = 0.25)) +
  scale_fill_manual(values = ggplot2::alpha(color.use.bpdph, alpha = dot.alpha), drop = FALSE) + guides(fill=FALSE) +
  geom_text(data = combined, mapping = aes(x = x, y = y + 0.5, label = labels), colour = "black", size = 4, show.legend = FALSE,segment.size = 0.2, segment.alpha = 0.5, max.overlaps = 12)
gg

ggsave("./cellchatsendersreceivers.pdf",
       dpi = 450, bg = "transparent",
       width=unit(12, "cm"), height=unit(9, "cm"))

pdf(file = paste0("./heatmap_of_all_signals_per_celltype_", "normal_ADGRL.pdf"), width = 7, height = 4)
ht1 <- netAnalysis_signalingRole_heatmap(cellchat.normal, signaling = c("TENASCIN"), 
                                         pattern = "all", width = 10, height = 2, font.size = 14, 
                                         font.size.title = 20)
ComplexHeatmap::draw(ht1)
dev.off()

pdf(file = paste0("./heatmap_of_all_signals_per_celltype_", "bpd_ADGRL.pdf"),  width = 7, height = 4)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat.bpd, signaling = c("TENASCIN"), 
                                         pattern = "all", width = 10, height = 2, 
                                         font.size = 14, font.size.title = 20)
ComplexHeatmap::draw(ht2)
dev.off()


pdf(file = paste0("./heatmap_of_all_signals_per_celltype_", "bpdph_ADGRL.pdf"), width = 7, height = 4)
ht3 <- netAnalysis_signalingRole_heatmap(cellchat.bpdph, signaling = c("TENASCIN"), 
                                         pattern = "all", width = 10, height = 2, 
                                         font.size = 14, font.size.title = 20)
ComplexHeatmap::draw(ht3)
dev.off()

# Circle plot show pathway
pathways.show <- c("ADGRL") 
par(mfrow=c(1,2))
# 设置图片保存参数
png_filename <- "cellchat_network_circle.png"
width <- 12
height <- 10
res <- 300

# 保存normal组的网络图
png(paste0("normal_", png_filename), width = width, height = height, units = "in", res = res)
netVisual_aggregate(cellchat.normal, signaling = pathways.show, layout = "circle")
dev.off()

# 保存BPD组的网络图
png(paste0("bpd_", png_filename), width = width, height = height, units = "in", res = res)
netVisual_aggregate(cellchat.bpd, signaling = pathways.show, layout = "circle")
dev.off()

# 保存BPDph组的网络图
png(paste0("bpdph_", png_filename), width = width, height = height, units = "in", res = res)
netVisual_aggregate(cellchat.bpdph, signaling = pathways.show, layout = "circle")
dev.off()

##ADGRL基因的表达
# VlnPlot(cellchat_seurat,
#         "ADGRL1",
#         group.by = "condition")

# 提取表达数据和元数据
#TENM4 TENM3 NRXN2 FLRT2 ADGRL2 ADGRL1
expression_data <- FetchData(cellchat_seurat, vars = c("NRXN2", "condition"))

# 使用ggplot2绘制小提琴图
library(ggplot2)

ggplot(expression_data, aes(x = condition, y = NRXN2, fill = condition)) +
  geom_violin(scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.1, fill = "white", alpha = 0.7, outlier.shape = NA) +
  theme_classic() +
  labs(title = "NRXN2 Expression by Condition",
       x = "Condition",
       y = "NRXN2 Expression Level") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1))


library(ggplot2)
library(ggpubr)
library(tidyr)
library(dplyr)

# 定义要绘制的基因列表
gene_list <- c("TENM4", "TENM3", "NRXN2", "FLRT2", "ADGRL2", "ADGRL1")
gene_list <- c("TNXB", "SDC4", "SDC1", "ITGAV", "ITGB6", "ITGA9","ITGB1")

# 提取表达数据和元数据
expression_data <- FetchData(cellchat_seurat, vars = c(gene_list, "condition"))

# 重命名condition为更易读的标签
expression_data$condition <- factor(expression_data$condition,
                                    labels = c("Term infant", "BPD", "BPD+PH"))


# 转换数据为长格式
expression_long <- expression_data %>%
  pivot_longer(cols = all_of(gene_list), 
               names_to = "gene", 
               values_to = "expression") %>%
  mutate(gene = factor(gene, levels = gene_list))  # 保持基因顺序


# 绘制分面箱线图
ggplot(expression_long, aes(x = condition, y = expression, fill = condition)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 0.8) +
  stat_compare_means(method = "wilcox.test", 
                     comparisons = list(c("normal", "bpd"), 
                                        c("normal", "bpdph"),
                                        c("bpd", "bpdph")),
                     label = "p.format", 
                     size = 3) +
  facet_wrap(~ gene, scales = "free_y", ncol = 3) +
  scale_fill_brewer(palette = "Set2") +
  theme_bw() +
  labs(title = "Expression of Adhesion GPCR and Related Genes",
       x = "Condition",
       y = "Expression Level") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none",
        strip.text = element_text(face = "bold"))

# 绘制小提琴图 + 箱线图
# 绘制图形
ggplot(expression_long, aes(x = condition, y = expression, fill = condition)) +
  geom_violin(alpha = 0.6, scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.2, fill = "white", alpha = 0.8, outlier.shape = NA) +
  stat_compare_means(method = "wilcox.test", 
                     ref.group = "Term infant",  # 以Term infant组为参照
                     label = "p.format",
                     size = 3) +
  facet_wrap(~ gene, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = c("Term infant" = "#66c2a5", "BPD" = "#fc8d62", "BPD+PH" = "#8da0cb")) +
  theme_minimal() +
  labs(title = "Gene Expression Patterns Across Conditions",
       subtitle = "Adhesion GPCR family and related genes",
       x = NULL,
       y = "Normalized Expression") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        legend.position = "none",
        strip.background = element_rect(fill = "gray90"),
        strip.text = element_text(face = "bold"))


########可视化Adventitial_FB_Cluster0与gCap（Normal组）以及Adventitial_FB_Cluster0与AT1（BPD组）之间的通路通信情况。####
# 1. 提取Normal组中Adventitial_FB_Cluster0与gCap的通信
net_normal_fb_gcap <- netVisual_bubble(cellchat.normal,
                                       sources.use = "Adventitial_FB_Cluster0",
                                       targets.use = "gCap",
                                       remove.isolate = FALSE,
                                       title.name = "Term infant: Adventitial_FB_Cluster0 → gCap Communication")

# 2. 提取BPD组中Adventitial_FB_Cluster0与AT1的通信
net_bpd_fb_at1 <- netVisual_bubble(cellchat.bpd,
                                   sources.use = "Adventitial_FB_Cluster0", 
                                   targets.use = "AT1",
                                   remove.isolate = FALSE,
                                   title.name = "BPD: Adventitial_FB_Cluster0 → AT1 Communication")

# 3. 保存图片
ggsave("normal_FB0_gCap_communication.png", net_normal_fb_gcap, 
       width = 4, height = 18, dpi = 300)
ggsave("bpd_FB0_AT1_communication.png", net_bpd_fb_at1, 
       width = 4, height = 18, dpi = 300)



# 提取Normal组中Adventitial_FB_Cluster0 → gCap的通信通路
fb0_gcap_normal <- subsetCommunication(cellchat.normal,
                                       sources.use = "Adventitial_FB_Cluster0",
                                       targets.use = "gCap")
fb4_AT2_normal <- subsetCommunication(cellchat.normal,
                                       sources.use = "Adventitial_FB_Cluster4",
                                       targets.use = "AT2")
fb4_AT1_normal <- subsetCommunication(cellchat.normal,
                                      sources.use = "Adventitial_FB_Cluster4",
                                      targets.use = "AT1")

# 提取BPD组中Adventitial_FB_Cluster0 → AT1的通信通路  
fb0_at1_bpd <- subsetCommunication(cellchat.bpd,
                                   sources.use = "Adventitial_FB_Cluster0",
                                   targets.use = "AT1")
fb2_at1_bpd <- subsetCommunication(cellchat.bpd,
                                   sources.use = "Adventitial_FB_Cluster2",
                                   targets.use = "AT1")

# 查看通路概况
cat("=== Normal组: FB_Cluster0 → gCap 通信通路 ===\n")
if(nrow(fb0_gcap_normal) > 0) {
  print(table(fb0_gcap_normal$pathway_name))
}

cat("\n=== BPD组: FB_Cluster0 → AT1 通信通路 ===\n")
if(nrow(fb0_at1_bpd) > 0) {
  print(table(fb0_at1_bpd$pathway_name))
}

library(dplyr)
library(ggplot2)

# 创建通路比较数据框
pathway_comparison_2 <- function(comm_df1, comm_df2, name1, name2) {
  if(nrow(comm_df1) == 0 | nrow(comm_df2) == 0) return(data.frame())
  
  pathways1 <- comm_df1 %>%
    group_by(pathway_name) %>%
    summarise(prob_sum = sum(prob),
              n_interactions = n()) %>%
    mutate(condition = name1)
  
  pathways2 <- comm_df2 %>%
    group_by(pathway_name) %>%
    summarise(prob_sum = sum(prob), 
              n_interactions = n()) %>%
    mutate(condition = name2)
  
  bind_rows(pathways1, pathways2)
}

comparison_data <- pathway_comparison_2(fb0_gcap_normal, fb0_at1_bpd,
                                      "FB0→gCap (Normal)", "FB0→AT1 (BPD)")
comparison_data_2 <- pathway_comparison_2(fb4_AT2_normal, fb2_at1_bpd,
                                        "FB4→AT2 (Normal)", "FB2→AT1 (BPD)")
comparison_data_3 <- pathway_comparison_2(fb4_AT1_normal, fb2_at1_bpd,
                                          "FB4→AT1 (Normal)", "FB2→AT1 (BPD)")

# 绘制通路比较图
if(nrow(comparison_data) > 0) {
  p_pathway <- ggplot(comparison_data_2, 
                      aes(x = pathway_name, y = prob_sum, fill = condition)) +
    geom_col(position = "dodge", alpha = 0.8) +
    coord_flip() +
    scale_fill_manual(values = c("FB4→AT2 (Normal)" = "grey", 
                                 "FB2→AT1 (BPD)" = "#fc8d62")) +
    theme_minimal() +
    labs(title = "Pathway Communication Strength Comparison",
         x = "Pathway",
         y = "Total Communication Probability",
         fill = "Interaction") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))+
    theme_bw()+
    theme(panel.grid = element_blank())
  
  print(p_pathway)
  ggsave("pathway_strength_comparison.png", p_pathway, width = 4.5, height = 7.6, dpi = 300)
  ggsave("pathway_strength_comparison.pdf", p_pathway, width = 4.5, height = 7.6)
}


# 分析具体的配体-受体机制
analyze_ligand_receptor <- function(comm_df, title) {
  if(nrow(comm_df) == 0) return(NULL)
  
  lr_analysis <- comm_df %>%
    group_by(ligand, receptor) %>%
    summarise(avg_prob = mean(prob),
              pathway_count = n_distinct(pathway_name),
              .groups = 'drop') %>%
    arrange(desc(avg_prob))
  
  cat("\n", title, " - Top Ligand-Receptor Pairs:\n")
  print(head(lr_analysis, 10))
  
  return(lr_analysis)
}

lr_normal <- analyze_ligand_receptor(fb0_gcap_normal, "Normal: FB0→gCap")
lr_bpd <- analyze_ligand_receptor(fb0_at1_bpd, "BPD: FB0→AT1")

save.image(file = "/data/nas1/xiangxiaoyu_OD/1.BPD+MSCs/1.scRNA/04_Adventitial_FB/07_CellChat/cellchat.workplace.RData")
