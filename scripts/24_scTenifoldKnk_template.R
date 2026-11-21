###############################################################################
# Template for the scTenifoldKnk virtual knockout
#
# Source      : scTenifoldKnk1.R
# Input       : mesenchymal count matrix
# Output      : scTenifoldKnk result object
# Note        : The script that produced the published PHLDA1 result (Fig. 8A-C; Supplementary Table 18) was lost. This is the template it was adapted from; it still carries the gene and paths of the project it came from. revision_scripts/B5b_scTenifoldKnk.R re-runs the analysis for PHLDA1.
#
# Code below is unchanged from the original analysis script.
###############################################################################
# -----------------------------------------------------------------------------
# PATHS: this file keeps the absolute paths of the analysis server on which the
# study was run (/data/nas1/xiangxiaoyu_OD/1.BPD+MSCs/1.scRNA/... and
# /home/xxyloveme/study/BPD+MSCs/...). Set PROJ_DIR in config/paths.R and
# replace those prefixes with your own project directory before running.
# -----------------------------------------------------------------------------




rm(list = ls())
gc()
setwd(dir = "/data/nas1/chenyongzhi/Research/07-scTenifoldKnk/")


library(dplyr)
library(Seurat)
library(scTenifoldKnk)
library(ggplot2)
library(ggrepel) 


# nohup Rscript /data/nas1/chenyongzhi/Research/07-scTenifoldKnk/scTenifoldKnk.R > log.txt 2>&1 &

load('/data/nas1/chenyongzhi/Research/08-CMC/UMAP_singleR.Rdata')


scRNA_harmony=UMAP
table(Idents(scRNA_harmony))
table(Idents(scRNA_harmony))




# 提取count矩阵作为输入
countMatrix <- GetAssayData(scRNA_harmony, slot = "counts")



# 执行敲除分析（推荐用高变基因做）：
# 参数gKO用于指定敲除的基因；qc开头的系列参数能够选择是否对矩阵进行质控以及相应指标；nc开头的参数对应网络的构建，
# 分别包括子网络的数量、使用的PC数量、随机抽取的细胞数量，数值越大运行时间越长但是结果可能会更稳定

result <- scTenifoldKnk(countMatrix = countMatrix, 
                        gKO = 'RBMS3', #需要敲除的基因
                        # qc = TRUE,#是否进行QC
                        # qc_mtThreshold = 0.1,#mt阈值
                        # qc_minLSize = 1000,#文库阈值(细胞测到的基因总数)
                        # qc_minLSize = 10,#文库阈值(细胞测到的基因总数)
                        # nc_nNet = 10, #子网络数量
                        # nc_nCells = 500, #每个网络中随机抽取的细胞数
                        # nc_nComp = 3 #PCA 的主成分数量
)


saveRDS(result, file = "scTenifoldKnk_PCOS_RBMS3.rds")






#######################################################


# 读取
result <- readRDS("scTenifoldKnk_PCOS_CADPS2.rds")

genes <- result$diffRegulation


top_genes <- head(result$diffRegulation[order(-result$diffRegulation$FC), ], 21)[-1,]
ggplot(top_genes, aes(x=reorder(gene, FC), y=FC)) +
  geom_bar(stat='identity', fill='#8e1a7f') +
  coord_flip() +
  labs(title="Top 20 Differentially Regulated Genes",
       x="Gene", y="FC") +
  theme_minimal()


genes0.05 <- genes[genes$p.value<0.05,]   ###########341个
genes0.05adj <- genes[genes$p.adj<0.05,]   ###########165个



df <- result$diffRegulation
df$log_pval <- -log10(df$p.value)
label_genes <- subset(df, abs(Z) > 2 & p.value < 0.05)
ggplot(df, aes(x=Z, y=log_pval)) +
  geom_point(alpha=0.5) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color="red") +
  geom_vline(xintercept=c(2), linetype="dashed", color="blue") +
  geom_text_repel(data=label_genes, aes(label=gene),
                  size=3, max.overlaps=50) +
  labs(title="Z vs -log10(p-value)",
       x="Z-score", y="-log10(p-value)") +
  theme_classic()
