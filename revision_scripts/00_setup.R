###############################################################################
# 00_setup.R  —— 所有补充分析共用的环境设置（修回期间的文件名为 00_setup_R2.R）
# 每个 B*.R 脚本开头都会 source() 本文件
#
# 运行环境：服务器 R 4.3.x，Seurat 4.4.0（与原分析一致）
###############################################################################

## ---- 0. 路径（如目录不同，只需改这里） -------------------------------------
# BASE is the project directory of the revision analyses; see config/paths.R
BASE <- "/home/xxyloveme/study/BPD+MSCs"
OUT  <- file.path(BASE, "12_Revision")                 # 所有新结果统一放这里

WORKPLACE_RDATA <- file.path(BASE, "workplace.RData")
EXPR_RDATA      <- file.path(BASE, "Expression_validation+ROC+xCell+CIBERSORT.RData")

PATH_ALLDATA  <- file.path(BASE, "03_All_cells/alldata_neonatal_injury.rds")
PATH_IMM      <- file.path(BASE, "02_Cell_Clustering_and_Annotation/4_immune/annotated_imm.rds")
PATH_MONOCLE  <- file.path(BASE, "04_Adventitial_FB/08_Monocle3/mycds_order_cells_cds_2_final.RData")
PATH_GTF      <- file.path(BASE, "06_Expression_validation/gencode.v38.annotation.gtf")
PATH_TRAJ     <- c(file.path(BASE, "04_Adventitial_FB/09_Continued_analysis/deg_Cluster_1.csv"),
                   file.path(BASE, "04_Adventitial_FB/08_Monocle3/deg_Cluster_2.csv"),
                   file.path(BASE, "04_Adventitial_FB/08_Monocle3/deg_Cluster_3.csv"),
                   file.path(BASE, "04_Adventitial_FB/08_Monocle3/deg_Cluster_4.csv"))
PATH_ML_DIR   <- file.path(BASE, "05_Machine_learning")

sub_dirs <- c("B1_annotation", "B2_donor", "B3_pseudobulk", "B4_trajectory",
              "B5_sex", "B5_scTenifoldKnk", "B7_supptable", "B9_immune", "B10_datasets")
for (d in c(OUT, file.path(OUT, sub_dirs))) dir.create(d, showWarnings = FALSE, recursive = TRUE)

## ---- 1. 包 ------------------------------------------------------------------
suppressPackageStartupMessages({
  library(Seurat); library(Matrix)
  library(dplyr);  library(tidyr); library(tibble)
  library(ggplot2); library(patchwork)
})
set.seed(2026)
options(stringsAsFactors = FALSE)

## ---- 2. 读入 workplace.RData（只在对象不存在时读，避免重复耗时） ------------
if (!exists("meso_relabel") || !exists("meso_adventitial_FB")) {
  message("Loading workplace.RData ...")
  load(WORKPLACE_RDATA)
}

## ---- 3. 统一供者 / 分组信息 --------------------------------------------------
# 原流程中 orig.ident 形如 "BPD 1_7m_BPD_Severe"，dataset 为第一段 "BPD 1"，即供者编号
DONOR_ORDER <- c("Term infant 1", "Term infant 2", "BPD 1", "BPD 2",
                 "BPD+PH 1", "BPD+PH 2", "Acute preterm injury 1")
GROUP_ORDER <- c("Term infant", "BPD", "BPD+PH", "Acute preterm injury")

# GEO 样本号对照（来自 r.00_rawdata.r）
GSM_MAP <- c("Term infant 1" = "GSM8488400 (0d-norm)",
             "Term infant 2" = "GSM8488401 (20d-norm)",
             "BPD 1"         = "GSM8488396 (7mo-1)",
             "BPD 2"         = "GSM8488397 (7mo-2)",
             "BPD+PH 1"      = "GSM8488398 (5796-1)",
             "BPD+PH 2"      = "GSM8488399 (5796-2)",
             "Acute preterm injury 1" = "GSM8488395 (7466-1)")

add_donor <- function(obj) {
  md <- obj@meta.data
  donor <- if ("dataset" %in% colnames(md)) as.character(md$dataset) else
    sub("_.*$", "", as.character(md$orig.ident))
  obj$donor <- factor(donor, levels = intersect(DONOR_ORDER, unique(donor)))
  obj$group <- factor(as.character(md$condition),
                      levels = intersect(GROUP_ORDER, unique(as.character(md$condition))))
  obj
}

meso_relabel            <- add_donor(meso_relabel)
meso_relabel_no_preterm <- add_donor(meso_relabel_no_preterm)
meso_adventitial_FB     <- add_donor(meso_adventitial_FB)
meso_adventitial_FB$subcluster <- factor(meso_adventitial_FB$SCT_snn_res.0.5)

cat("\n== 供者 × 分组核对（每个供者只应对应一个分组）==\n")
print(table(meso_relabel$donor, meso_relabel$group))

## ---- 4. 小工具 ---------------------------------------------------------------
# 按供者（或任意分组）把原始 counts 求和成 pseudo-bulk
pseudobulk <- function(obj, by = "donor", assay = "RNA") {
  cts <- GetAssayData(obj, assay = assay, slot = "counts")
  f   <- droplevels(factor(obj@meta.data[[by]]))
  mm  <- sparse.model.matrix(~ 0 + f)
  colnames(mm) <- levels(f)
  as.matrix(cts %*% mm)
}

# 确保 RNA assay 有 log-normalized data（FetchData 取表达值用）
ensure_rna_norm <- function(obj) {
  DefaultAssay(obj) <- "RNA"
  NormalizeData(obj, normalization.method = "LogNormalize",
                scale.factor = 10000, verbose = FALSE)   # 与原 DEG 分析相同的标准化
}

pick_umap <- function(obj) {
  r <- Reductions(obj)
  if ("umap" %in% r) "umap" else tail(r[grepl("umap", r, ignore.case = TRUE)], 1)
}

save_both <- function(p, file_noext, w, h) {
  ggsave(paste0(file_noext, ".pdf"), p, width = w, height = h)
}

CORE_GENES <- c("AES", "TLE5", "ACTG1", "TSC22D2", "HSPA1A", "SLC26A7", "MT1A", "NFAT5")

cat("\n00_setup.R loaded. Output dir:", OUT, "\n")
