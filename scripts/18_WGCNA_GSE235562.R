###############################################################################
# WGCNA co-expression modules and module-trait correlation in GSE235562
#
# Source      : r.03_Cell_prop.r (lines 5237-5789)
# Input       : GSE235562 microarray matrix
# Output      : Fig. 5B; Supplementary Fig. 8; Supplementary Table 15
#
# Code below is unchanged from the original analysis script.
###############################################################################
# -----------------------------------------------------------------------------
# PATHS: this file keeps the absolute paths of the analysis server on which the
# study was run (/data/nas1/xiangxiaoyu_OD/1.BPD+MSCs/1.scRNA/... and
# /home/xxyloveme/study/BPD+MSCs/...). Set PROJ_DIR in config/paths.R and
# replace those prefixes with your own project directory before running.
# -----------------------------------------------------------------------------


#####################3.7 bulkRNA进行WGCNA分析 GSE235562数据集###################
setwd("/data/nas1/xiangxiaoyu_OD/1.BPD+MSCs/1.scRNA/04_WGCNA/1.GSE235562/")
load(file = "/data/nas1/xiangxiaoyu_OD/1.BPD+MSCs/1.scRNA/04_Adventitial_FB/venn.RData")

expr_GSE235562 <- read.table(file = "GSE235562_series_matrix.txt",header = T,sep = "\t")
colnames(expr_GSE235562)[1] <- "ID"

probe_GSE235562 = read.table(file = 'GPL20844-88004.txt',
                            sep = '\t',
                            quote = '',         
                            comment.char = '#',    # 过滤掉'GPL14951-11332.txt'文件中以‘#’开头的注释信息
                            header = T,
                            fill = T,                #  如果文件中某行的数据少于其他行，则自动添加空白域。     
                            stringsAsFactors = F)  # 字符串不改为因子

ids_GSE235562 = probe_GSE235562[probe_GSE235562$GENE_SYMBOL != '',
                              c(1,10)] # 提取探针和geneID

##下载平台文件，提取注释信息
ids_GSE235562 <- filter(ids_GSE235562, `GENE_SYMBOL` != '')

#处理一个基因对应多行表达值的情况
exprSet_GSE235562 <- expr_GSE235562 %>%
  merge(ids_GSE235562, by='ID')%>%
  dplyr::select(-ID)%>%     ## 去除多余信息
  dplyr::select(GENE_SYMBOL, everything())%>%     ## 重新排列
  mutate(rowMean = rowMeans(.[grep('GSM', names(.))]))%>%    ## 求出平均数
  arrange(desc(rowMean))%>%       ## 把表达量的平均值从大到小排序
  distinct(GENE_SYMBOL, .keep_all = T)%>%      ## symbol留下第一个
  dplyr::select(-rowMean)%>%     ## 反向选择去除rowMean这一列
  tibble::column_to_rownames(colnames(.)[1])   ## 把第一列变成行名并删除

#提取临床信息文件，准备分组
pData_GSE235562 <- read.table(file = "GSE235562_pData.txt",header = T,sep = "\t")
group_GSE235562 <- data.frame(sample = pData_GSE235562$X.Sample_geo_accession, group = pData_GSE235562$group)
table(group_GSE235562$group)
# CMS CMS+Hox40 CMS+Hox80   control     Hox40     Hox80 
# 15         8         8        16         9         8 
group_GSE235562 <- group_GSE235562[order(group_GSE235562$group),]

exprSet_GSE235562 <- exprSet_GSE235562[,group_GSE235562$sample]

##判断是否需要log2转化
qx <- as.numeric(quantile(exprSet_GSE235562, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) {
  exprSet_GSE235562[which(exprSet_GSE235562 <= 0)] <- 0
  exprSet_GSE235562 <- log2(exprSet_GSE235562 + 1)
  print("log2 transform finished")
}else{
  print("log2 transform not needed")
}

boxplot(exprSet_GSE235562,outline=FALSE,las=2,col = 'blue',xaxt = 'n',ann = F)

pdf(file = "exprSet_GSE235562.pdf",width = 6,height = 3)
boxplot(exprSet_GSE235562,outline=FALSE,las=2,col = 'blue',xaxt = 'n',ann = F)
dev.off()

png(file = "exprSet_GSE235562.png",width = 6,height = 3,units = "in",res = 300)
boxplot(exprSet_GSE235562,outline=FALSE,las=2,col = 'blue',xaxt = 'n',ann = F)
dev.off()

write.csv(exprSet_GSE235562, file = "exprSet_GSE235562.csv")
write.csv(group_GSE235562, file = "group_GSE235562.csv", row.names = F)


# 加载必要的包
library(WGCNA)
library(ggplot2)
library(dplyr)
library(corrplot)
options(stringsAsFactors = FALSE)

# 1. 数据准备和预处理
# 假设你已经有了表达矩阵 exprSet_GSE235562 和分组信息 group_GSE235562
#AC020916.1 LOC284454
#AL117190.2 MEG3
#SELENOW SEPW1
#NORAD LINC00657
#RACK1 GNB2L1
#ATP5F1E ATP5E
#ECRG4 C2orf40 
#SNHG29 LRRC75A-AS1
#TLE5 AES
#CCN2 CTGF
#ELOB TCEB2

# 定义基因名称映射关系
gene_mapping <- list(
  "SELENOW" = "SEPW1",
  "NORAD" = "LINC00657", 
  "RACK1" = "GNB2L1",
  "ATP5F1E" = "ATP5E",
  "ECRG4" = "C2orf40",
  "SNHG29" = "LRRC75A-AS1",
  "TLE5" = "AES",
  "CCN2" = "CTGF",
  "ELOB" = "TCEB2")

# 创建替换函数
replace_genes <- function(gene_list, mapping) {
  new_list <- gene_list
  for (old_name in names(mapping)) {
    if (old_name %in% new_list) {
      new_list[new_list == old_name] <- mapping[[old_name]]
      cat("已将", old_name, "替换为", mapping[[old_name]], "\n")
    } else {
      cat("警告: 未找到基因", old_name, "在候选基因列表中\n")
    }
  }
  return(new_list)
}

# 执行基因名称替换
candidate_genes <- replace_genes(ggVenn_candiated_genes, gene_mapping)

# 检查候选基因在表达矩阵中的存在情况
available_genes <- candidate_genes[candidate_genes %in% rownames(exprSet_GSE235562)]
print(paste("在表达矩阵中找到的候选基因数量:", length(available_genes)))
missing_genes <- setdiff(candidate_genes, available_genes)
# 提取候选基因的表达矩阵
expr_candidate <- exprSet_GSE235562[available_genes, ]

# 转置矩阵（WGCNA要求行是样本，列是基因）
datExpr <- t(expr_candidate)

# 2. 样本聚类检查异常样本
sampleTree <- hclust(dist(datExpr), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")

# 3. 定义临床性状（根据实验设计创建性状矩阵）
# 创建性状矩阵
datTraits <- data.frame(
  # 基础分组
  Control = as.numeric(group_GSE235562$group == "control"),
  CMS = as.numeric(group_GSE235562$group == "CMS"),
  Hox40 = as.numeric(group_GSE235562$group == "Hox40"),
  Hox80 = as.numeric(group_GSE235562$group == "Hox80"),
  CMS_Hox40 = as.numeric(group_GSE235562$group == "CMS+Hox40"),
  CMS_Hox80 = as.numeric(group_GSE235562$group == "CMS+Hox80"),
  
  # 组合性状（可能更有生物学意义）
  Any_CMS = as.numeric(group_GSE235562$group %in% c("CMS", "CMS+Hox40", "CMS+Hox80")),
  Any_Hox = as.numeric(group_GSE235562$group %in% c("Hox40", "Hox80", "CMS+Hox40", "CMS+Hox80")),
  Stress_Only = as.numeric(group_GSE235562$group == "CMS"),
  Oxygen_Only = as.numeric(group_GSE235562$group %in% c("Hox40", "Hox80")),
  Combined_Treatment = as.numeric(group_GSE235562$group %in% c("CMS+Hox40", "CMS+Hox80")))

rownames(datTraits) <- colnames(exprSet_GSE235562)

# 4. 软阈值选择
# 选择一组软阈值功率
powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))

# 调用网络拓扑分析函数
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

# 绘制软阈值选择图
par(mfrow = c(1, 2))
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R^2",
     type = "n", main = paste("Scale independence"))
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, col = "red")
abline(h = 0.90, col = "red")

plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity",
     type = "n", main = paste("Mean connectivity"))
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, col = "red")

# 选择软阈值功率
softPower <- sft$powerEstimate
if (is.na(softPower)) {
  softPower <- 6  # 如果自动选择失败，使用默认值
}
print(paste("选择的软阈值功率:", softPower))

# 5. 构建共表达网络和识别模块
net <- blockwiseModules(
  datExpr,
  power = softPower,
  TOMType = "signed",
  minModuleSize = 5,  # 由于基因数较少，设置较小的模块大小
  reassignThreshold = 0,
  mergeCutHeight = 0.25,
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs = FALSE,
  saveTOMFileBase = "GSE235562_TOM",
  verbose = 3)

# 6. 模块可视化
# 将数字标签转换为颜色
moduleColors <- labels2colors(net$colors)

# 绘制模块树状图
plotDendroAndColors(
  net$dendrograms[[1]],
  moduleColors[net$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05)

# 7. 模块-性状关联分析
# 计算模块特征基因
MEs <- net$MEs
MEs <- orderMEs(MEs)  # 根据层次聚类排序

# 计算模块与性状的相关性
moduleTraitCor <- cor(MEs, datTraits, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(datExpr))

# 可视化模块-性状相关性
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)

par(mar = c(6, 8.5, 3, 3))
labeledHeatmap(
  Matrix = moduleTraitCor,
  xLabels = names(datTraits),
  yLabels = names(MEs),
  ySymbols = names(MEs),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix,
  setStdMargins = FALSE,
  cex.text = 0.5,
  zlim = c(-1, 1),
  main = paste("Module-trait relationships"))


##重新美化
# 加载必要的包
library(ComplexHeatmap)
library(circlize)
library(scico)
# 准备数据
# moduleTraitCor <- cor(net$MEs, datTraits, use = "p")
# moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(datExpr))
# ==============================================
# ComplexHeatmap（CNS风配色，格内仅显示 r）
# ==============================================

# 连续色映射：-1,-0.5,0,0.5,1 → vik
col_fun <- colorRamp2(
  seq(-1, 1, length.out = 11),
  scico(11, palette = "vik")
)

# 左侧模块颜色：Okabe-Ito（>8个模块自动退回 scico::batlow）
okabe_ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", 
               "#0072B2", "#D55E00", "#CC79A7", "#999999")

modCols <- sub("^ME", "", rownames(moduleTraitCor))
unq_mod <- unique(modCols)

if (length(unq_mod) <= length(okabe_ito)) {
  pal_discrete <- setNames(okabe_ito[seq_along(unq_mod)], unq_mod)
} else {
  pal_discrete <- setNames(scico(length(unq_mod), palette = "batlow"), unq_mod)
}

row_anno <- rowAnnotation(
  Module = anno_simple(modCols, col = pal_discrete),
  annotation_name_side = "top",
  annotation_legend_param = list(
    title = "Module",
    title_gp = gpar(fontsize = 10, fontface = "bold"),
    labels_gp = gpar(fontsize = 8)
  )
)

# 格内仅显示 r；使用 per-cell viewport 防溢出
lab_fun <- function(j, i, x, y, width, height, fill) {
  r <- moduleTraitCor[i, j]
  pval <- moduleTraitPvalue[i, j]
  
  # 根据p值添加星号标记
  significance <- ""
  if (pval < 0.001) {
    significance <- "***"
  } else if (pval < 0.01) {
    significance <- "**"
  } else if (pval < 0.05) {
    significance <- "*"
  }
  
  grid::pushViewport(grid::viewport(x = x, y = y, width = width, height = height,
                                    just = c("center", "center"), clip = "on"))
  
  # 显示相关系数和显著性
  grid::grid.text(sprintf("%.2f%s", r, significance), 
                  x = 0.5, y = 0.5,
                  gp = grid::gpar(fontsize = 8, fontface = "bold"))
  grid::popViewport()
}

# 创建热图
ht <- row_anno + 
  Heatmap(
    moduleTraitCor,
    name = "Correlation\n(r)",
    col = col_fun,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    row_names_side = "left",
    row_names_gp = grid::gpar(fontsize = 9, fontface = "bold"),
    row_title = "Modules",
    row_title_gp = gpar(fontsize = 10, fontface = "bold"),
    column_names_rot = 45,
    column_names_gp = grid::gpar(fontsize = 9, fontface = "bold"),
    column_title = "Experimental Conditions",
    column_title_gp = gpar(fontsize = 10, fontface = "bold"),
    heatmap_legend_param = list(
      at = c(-1, -0.5, 0, 0.5, 1),
      legend_height = unit(35, "mm"),
      title_gp = gpar(fontsize = 10, fontface = "bold"),
      labels_gp = gpar(fontsize = 8)
    ),
    cell_fun = lab_fun,
    rect_gp = gpar(col = "white", lwd = 1),  # 添加白色边框
    show_heatmap_legend = TRUE
  )

# 添加标题
title_text <- "Module-Trait Relationships in GSE235562\n(Mesenchymal Stem Cells under Mechanical and Oxidative Stress)"

# 绘制并保存
pdf("WGCNA_module_trait_heatmap_complex.pdf", width = 7.67, height = 5.96)
draw(ht, 
     column_title = title_text, 
     column_title_gp = gpar(fontsize = 12, fontface = "bold"),
     heatmap_legend_side = "right")
dev.off()

# 同时保存为PNG格式（更高分辨率）
png("WGCNA_module_trait_heatmap_complex.png", width = 7.67, height = 5.96, res = 300, units = "in")
draw(ht, 
     column_title = title_text, 
     column_title_gp = gpar(fontsize = 12, fontface = "bold"),
     heatmap_legend_side = "right")
dev.off()

# 保存相关性矩阵
write.csv(moduleTraitCor, "module_trait_correlation_matrix.csv")
write.csv(moduleTraitPvalue, "module_trait_pvalue_matrix.csv")


# 8. 识别与BPD相关性状显著相关的模块
# 定义感兴趣的性状（根据BPD病理机制）
bpd_traits <- c("Any_CMS", "Any_Hox", "Combined_Treatment")

# 提取与BPD相关性状显著相关的模块
significant_modules <- data.frame()

for (trait in bpd_traits) {
  trait_cor <- moduleTraitCor[, trait]
  trait_pval <- moduleTraitPvalue[, trait]
  
  # 选择相关性高且显著的模块 (|cor| > 0.3, p < 0.05)
  sig_idx <- which(abs(trait_cor) > 0.3 & trait_pval < 0.05)
  
  if (length(sig_idx) > 0) {
    for (i in sig_idx) {
      module_name <- names(MEs)[i]
      significant_modules <- rbind(significant_modules, data.frame(
        Module = module_name,
        Trait = trait,
        Correlation = trait_cor[i],
        Pvalue = trait_pval[i],
        Genes_in_Module = sum(net$colors == as.numeric(gsub("ME", "", module_name)))
      ))
    }
  }
}

# 打印显著相关的模块
if (nrow(significant_modules) > 0) {
  print("与BPD相关性状显著相关的模块:")
  print(significant_modules)
} else {
  print("未找到与BPD相关性状显著相关的模块")
}
# [1] "与BPD相关性状显著相关的模块:"
# Module              Trait Correlation       Pvalue Genes_in_Module
# ME0     ME0            Any_Hox  -0.3657696 0.0029565604              34
# ME8     ME8            Any_Hox   0.3569482 0.0037874719               5
# ME01    ME0 Combined_Treatment  -0.3802178 0.0019412214              34

# 9. 提取关键模块的基因
#ME8的模块
if (nrow(significant_modules) > 0) {
  # 获取最显著的模块
  top_module <- significant_modules[which.max(abs(significant_modules$Correlation)), "Module"]
  module_number <- as.numeric(gsub("ME", "", top_module))
  # 提取该模块中的基因
  module_genes <- colnames(datExpr)[net$colors == module_number]
  print(paste("模块", top_module, "中的基因数量:", length(module_genes)))
  print("模块中的基因:")
  print(module_genes)
  # 保存模块基因列表
  write.csv(data.frame(Gene = module_genes), 
            paste0("WGCNA_Module_", top_module, "_Genes.csv"), 
            row.names = FALSE)
}

#
# > module_genes ME8
# [1] "FOSB"      "FTL"       "ZFPM2-AS1" "RLF"       "CDKN1A" 

##ME0的基因模块
module_genes_ME0 <- colnames(datExpr)[net$colors == "0"]
write.csv(data.frame(Gene = module_genes_ME0), 
          paste0("WGCNA_Module_", "ME0", "_Genes.csv"), 
          row.names = FALSE)
# > module_genes_ME0
# [1] "TSC22D2"  "FTH1"     "PFDN5"    "LGALS1"   "S100A11"  "SEPW1"    "JUN"      "HSPA1A"   "CCDC102B"
# [10] "FAM133B"  "INTS6"    "NFATC2"   "SLC26A7"  "RALGAPA1" "NFAT5"    "STAG1"    "TMSB4X"   "MT2A"    
# [19] "PPIA"     "ATP5E"    "MT1A"     "COMMD6"   "CFH"      "MIF"      "FMO2"     "LAMB1"    "PHLDA1"  
# [28] "CSRP1"    "AES"      "ACTG1"    "XIST"     "CRISPLD2" "ABCA8"    "ACTB"  

# 修正的基因显著性分析代码
# 选择我们感兴趣的模块：ME0 和 ME8

# 计算基因与性状的相关性
GS_Any_CMS <- as.numeric(cor(datExpr, datTraits$Any_CMS, use = "p"))
GeneSignificance <- abs(GS_Any_CMS)

# 计算模块成员关系 (Module Membership)
# 模块成员关系是基因与模块特征基因的相关性
module_membership <- cor(datExpr, net$MEs, use = "p")

# 为ME8模块绘制散点图
module_number_8 <- 8  # ME8对应的数字标签
module_genes_8 <- net$colors == module_number_8

if(sum(module_genes_8) > 0) {
  # 提取ME8模块的模块成员关系
  MM_8 <- abs(module_membership[module_genes_8, paste0("ME", module_number_8)])
  # 提取ME8模块的基因显著性
  GS_8 <- GeneSignificance[module_genes_8]
  # 绘制ME8模块的散点图
  verboseScatterplot(
    MM_8,
    GS_8,
    xlab = paste("Module Membership in ME", module_number_8),
    ylab = "Gene significance for Any_CMS",
    main = paste("ME", module_number_8, ": MM vs. GS"),
    cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, 
    col = moduleColors[module_genes_8])
  # 计算并显示相关性
  cor_8 <- cor(MM_8, GS_8, use = "p")
  print(paste("ME8模块成员与基因显著性的相关性:", round(cor_8, 3)))
}

# 为ME0模块绘制散点图
module_number_0 <- 0  # ME0对应的数字标签
module_genes_0 <- net$colors == module_number_0

if(sum(module_genes_0) > 0) {
  # 提取ME0模块的模块成员关系
  MM_0 <- abs(module_membership[module_genes_0, paste0("ME", module_number_0)])
  # 提取ME0模块的基因显著性
  GS_0 <- GeneSignificance[module_genes_0]
  # 绘制ME0模块的散点图
  verboseScatterplot(
    MM_0,
    GS_0,
    xlab = paste("Module Membership in ME", module_number_0),
    ylab = "Gene significance for Any_CMS",
    main = paste("ME", module_number_0, ": MM vs. GS"),
    cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, 
    col = moduleColors[module_genes_0])
  # 计算并显示相关性
  cor_0 <- cor(MM_0, GS_0, use = "p")
  print(paste("ME0模块成员与基因显著性的相关性:", round(cor_0, 3)))
}

# 识别hub基因（模块内连通性高的基因）
# 计算所有基因的模块内连通性
ADJ <- adjacency(datExpr, power = softPower)
connectivity <- intramodularConnectivity(ADJ, net$colors)

# 为ME8模块提取hub基因
hub_genes_8 <- connectivity[module_genes_8, ]
hub_genes_8$Gene <- rownames(hub_genes_8)
hub_genes_8 <- hub_genes_8[order(-hub_genes_8$kWithin), ]
top_hub_genes_8 <- head(hub_genes_8, 5)
print("ME8模块的top hub genes:")
print(top_hub_genes_8[, c("Gene", "kWithin")])

# 为ME0模块提取hub基因
hub_genes_0 <- connectivity[module_genes_0, ]
hub_genes_0$Gene <- rownames(hub_genes_0)
hub_genes_0 <- hub_genes_0[order(-hub_genes_0$kWithin), ]
top_hub_genes_0 <- head(hub_genes_0, 5)

print("ME0模块的top hub genes:")
print(top_hub_genes_0[, c("Gene", "kWithin")])

# 保存hub基因信息
write.csv(hub_genes_8, "ME8_Hub_Genes.csv", row.names = FALSE)
write.csv(hub_genes_0, "ME0_Hub_Genes.csv", row.names = FALSE)

# 生成模块基因的详细报告
module_report_8 <- data.frame(
  Gene = colnames(datExpr)[module_genes_8],
  Module = "ME8",
  ModuleMembership = MM_8,
  GeneSignificance = GS_8,
  Connectivity = hub_genes_8$kWithin)

module_report_0 <- data.frame(
  Gene = colnames(datExpr)[module_genes_0],
  Module = "ME0",
  ModuleMembership = MM_0,
  GeneSignificance = GS_0,
  Connectivity = hub_genes_0$kWithin)

# 合并报告并保存
all_module_report <- rbind(module_report_8, module_report_0)
write.csv(all_module_report, "ME0_ME8_Detailed_Report.csv", row.names = FALSE)

# 可视化两个模块的基因分布
library(ggplot2)

# 创建合并数据用于可视化
plot_data <- rbind(
  data.frame(Module = "ME8", MM = MM_8, GS = GS_8),
  data.frame(Module = "ME0", MM = MM_0, GS = GS_0))

p <- ggplot(plot_data, aes(x = MM, y = GS, color = Module)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(x = "Module Membership", y = "Gene Significance for Any_CMS",
       title = "Module Membership vs Gene Significance for ME0 and ME8") +
  theme_minimal() +
  facet_wrap(~ Module) +
  scale_color_manual(values = c("ME0" = "blue", "ME8" = "red"))

print(p)
ggsave("ME0_ME8_Scatter_Comparison.pdf", p, width = 10, height = 4.76)

save.image(file = "WGCNA.RData")
