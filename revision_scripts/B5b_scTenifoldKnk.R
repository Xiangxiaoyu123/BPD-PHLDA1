###############################################################################
# B5b  scTenifoldKnk 重跑（R5）—— 回应 R2-2：PHLDA1 虚拟敲除的下游基因是否受性别影响
#
#   原始结果（11_scTenifoldKnk/）已保存，可直接作为基准，无需再复现：
#     Table1_PHLDA1_KO_ALL_Genes.csv          原始全部基因结果
#     Table1_PHLDA1_KO_Significant_Genes.csv  原始显著基因（稿件 Supp Table 13）
#     Enrichr_*.csv                           原始富集（Enrichr：GO BP 2023 / KEGG 2021 / Reactome 2022）
#
#   本脚本只需新跑两次：
#     nosex  ：同一批细胞，去掉性染色体基因
#     female ：只用女性供者的细胞
#   （all：完全复现原分析，可选，不跑也能完成比较）
#
#   R5 相对 R4 的修改
#     1) 新增 prep 步骤：只在这一步读 workplace.RData，把输入矩阵存成小文件；
#        之后的 nosex / female 任务不再读 workplace.RData，内存占用大幅下降，可放心并行
#     2) 原始显著基因直接取 Significant_Genes 表，并自动判断当时用的阈值（p.adj 还是 p.value），
#        新结果用同一阈值，保证可比
#     3) 富集分析优先用 Enrichr（与原分析同库）；服务器不能联网时退回 clusterProfiler；
#        并检查原始富集前 10 条通路在新结果中是否仍出现
#     4) 运行前检查 PHLDA1 在输入细胞中的表达比例（scTenifoldKnk 内部会去掉 <5% 细胞表达的基因）
#
#   用法（revision_code 目录下）：
#     Rscript B5b_scTenifoldKnk_R5.R prep full                            # 约几分钟，只跑一次
#     nohup Rscript B5b_scTenifoldKnk_R5.R nosex  full > knk_nosex.log  2>&1 &
#     nohup Rscript B5b_scTenifoldKnk_R5.R female full > knk_female.log 2>&1 &
#     （可选）nohup Rscript B5b_scTenifoldKnk_R5.R all full > knk_all.log 2>&1 &
#     全部结束后：Rscript B5b_scTenifoldKnk_R5.R summary full > knk_summary.log 2>&1
#   第 2 个参数 full / fast：fast 只保留 ≥10% 细胞表达的基因，明显更快（作为敏感性分析，需在 prep 时同样指定）
###############################################################################
args      <- commandArgs(trailingOnly = TRUE)
RUN       <- if (length(args) >= 1) args[1] else "summary"
GENE_MODE <- if (length(args) >= 2) args[2] else "full"
stopifnot(RUN %in% c("prep", "all", "nosex", "female", "summary"), GENE_MODE %in% c("full", "fast"))
SUF <- if (GENE_MODE == "fast") "_fast" else ""
t_start <- Sys.time()
message(sprintf("[%s] RUN = %s, GENE_MODE = %s", format(t_start), RUN, GENE_MODE))

## ---- 路径（与 00_setup.R 保持一致；若不同请修改） --------------------------------
BASE  <- "/home/xxyloveme/study/BPD+MSCs"
OUT   <- file.path(BASE, "12_Revision")
SETUP <- here::here("revision_scripts", "00_setup.R")
od    <- file.path(OUT, "B5_scTenifoldKnk"); dir.create(od, recursive = TRUE, showWarnings = FALSE)
ORIG_DIR <- file.path(BASE, "11_scTenifoldKnk")
CELLS_FROM    <- "meso_relabel_no_preterm"
FEMALE_DONORS <- c("Term infant 2", "BPD 1", "BPD 2", "BPD+PH 1", "BPD+PH 2")   # 来自 B5 性别推断
N_PARALLEL    <- 2                                                             # 同时跑几个任务
N_CORES       <- max(1, floor((parallel::detectCores() - 2) / N_PARALLEL))
f_in <- function(tag) file.path(od, paste0("input_counts_", tag, SUF, ".rds"))

## =============================================================================
## prep：读 workplace.RData，生成输入矩阵（稀疏）后退出
## =============================================================================
if (RUN == "prep") {
  source(SETUP)
  obj <- get(CELLS_FROM)
  DefaultAssay(obj) <- "RNA"
  if (inherits(obj[["RNA"]], "Assay5") && length(Layers(obj[["RNA"]], search = "counts")) > 1)
    obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
  cts <- if (inherits(obj[["RNA"]], "Assay5")) LayerData(obj[["RNA"]], layer = "counts") else
    GetAssayData(obj, assay = "RNA", slot = "counts")
  donor <- as.character(obj$donor)
  sex_genes <- readRDS(file.path(OUT, "B5_sex", "sex_genes_strict.rds"))

  make_input <- function(cells, tag, drop_genes = character(0)) {
    m <- cts[, cells]
    m <- m[!rownames(m) %in% drop_genes, ]
    if (GENE_MODE == "fast")
      m <- m[Matrix::rowMeans(m > 0) >= 0.10 | grepl("^MT-", rownames(m)) |
               rownames(m) %in% c("PHLDA1", "RPS12", "RPL29", "FAU", "RACK1", "HSPA1B", "XIST"), ]
    pct_phlda1 <- 100 * mean(m["PHLDA1", ] > 0)
    message(sprintf("[%s] %d 基因 × %d 细胞；PHLDA1 表达细胞比例 = %.1f%%（需 > 5%%）",
                    tag, nrow(m), ncol(m), pct_phlda1))
    if (pct_phlda1 <= 5) warning(tag, "：PHLDA1 表达比例 ≤ 5%，scTenifoldKnk 会把它过滤掉！")
    print(table(donor[match(cells, colnames(obj))]))
    saveRDS(m, f_in(tag))
  }
  make_input(colnames(obj), "all")
  make_input(colnames(obj), "nosex", drop_genes = sex_genes)
  make_input(colnames(obj)[donor %in% FEMALE_DONORS], "female")
  message("prep 完成，输入文件在：", od); quit(save = "no")
}

## =============================================================================
## all / nosex / female：只读输入矩阵，不读 workplace.RData
## =============================================================================
if (RUN %in% c("all", "nosex", "female")) {
  suppressPackageStartupMessages({ library(Matrix); library(scTenifoldKnk) })
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    RhpcBLASctl::blas_set_num_threads(N_CORES); RhpcBLASctl::omp_set_num_threads(N_CORES)
  }
  tag_out <- c(all = "all_genes", nosex = "no_sex_genes", female = "female_only")[[RUN]]
  f_out <- file.path(od, paste0("diffRegulation_", tag_out, SUF, ".csv"))
  if (file.exists(f_out)) { message("已有结果：", f_out); quit(save = "no") }
  stopifnot(file.exists(f_in(RUN)))
  m <- readRDS(f_in(RUN))
  message(sprintf("输入：%d 基因 × %d 细胞；每任务核数 %d", nrow(m), ncol(m), N_CORES))
  set.seed(2026)
  # 与原 scTenifoldKnk1.R 一致：除 gKO 外全部默认参数（qc = TRUE, nc_nNet = 10, nc_nCells = 500 ...）
  res <- scTenifoldKnk(countMatrix = m, gKO = "PHLDA1", nCores = N_CORES)
  saveRDS(res, file.path(od, paste0("knk_", tag_out, SUF, ".rds")))
  write.csv(res$diffRegulation, f_out, row.names = FALSE)
  message(sprintf("完成：检验基因 %d，p.adj<0.05 %d；总用时 %s", nrow(res$diffRegulation),
                  sum(res$diffRegulation$p.adj < 0.05), format(Sys.time() - t_start)))
  quit(save = "no")
}

## =============================================================================
## summary：与原始结果比较（不读 workplace.RData）
## =============================================================================
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2) })
save_both <- function(p, f, w, h) { ggsave(paste0(f, ".pdf"), p, width = w, height = h)
                                    ggsave(paste0(f, ".png"), p, width = w, height = h, dpi = 300) }
jac <- function(a, b) length(intersect(a, b)) / max(1, length(union(a, b)))

## ---- 原始结果与当时使用的显著性阈值 -------------------------------------------------
orig_all <- read.csv(file.path(ORIG_DIR, "Table1_PHLDA1_KO_ALL_Genes.csv"), check.names = FALSE)
orig_sig_tab <- read.csv(file.path(ORIG_DIR, "Table1_PHLDA1_KO_Significant_Genes.csv"), check.names = FALSE)
fix_gene_col <- function(d) {        # 统一把基因列命名为 gene（兼容 "" / X / Gene 等写法）
  colnames(d)[colnames(d) == ""] <- "X"
  g <- intersect(c("gene", "Gene", "genes", "X"), colnames(d)); g <- if (length(g)) g[1] else colnames(d)[1]
  colnames(d)[colnames(d) == g] <- "gene"; d
}
orig_all <- fix_gene_col(orig_all); orig_sig_tab <- fix_gene_col(orig_sig_tab)
cat("原始结果列名：", paste(colnames(orig_all), collapse = ", "), "\n")
n_padj <- sum(orig_all$p.adj < 0.05, na.rm = TRUE); n_pval <- sum(orig_all$p.value < 0.05, na.rm = TRUE)
SIG_COL <- if (abs(nrow(orig_sig_tab) - n_padj) <= abs(nrow(orig_sig_tab) - n_pval)) "p.adj" else "p.value"
message(sprintf("原始显著表 %d 行；全表中 p.adj<0.05 = %d，p.value<0.05 = %d → 判断原分析用 %s < 0.05",
                nrow(orig_sig_tab), n_padj, n_pval, SIG_COL))
orig_sig <- setdiff(orig_sig_tab$gene, "PHLDA1")

## ---- 读新结果 --------------------------------------------------------------------
runs <- list(original = orig_all)
for (tg in c("all_genes", "no_sex_genes", "female_only")) {
  f <- file.path(od, paste0("diffRegulation_", tg, SUF, ".csv"))
  if (file.exists(f)) runs[[tg]] <- read.csv(f) else message("尚未完成：", basename(f))
}
stopifnot(length(runs) >= 2)
sig <- lapply(runs, function(d) setdiff(d$gene[d[[SIG_COL]] < 0.05], "PHLDA1"))
sig$original <- orig_sig

## ---- 1. 关键基因 ---------------------------------------------------------------------
key <- c("RPS12", "RPL29", "FAU", "RACK1", "HSPA1B", "XIST", "RPS4Y1", "DDX3Y")
key_tab <- bind_rows(lapply(names(runs), function(nm) {
  d <- runs[[nm]]
  data.frame(run = nm, gene = key, significant = key %in% sig[[nm]],
             Z = round(d$Z[match(key, d$gene)], 2), p = signif(d[[SIG_COL]][match(key, d$gene)], 3))
}))
write.csv(key_tab, file.path(od, "S1_key_genes_all_runs.csv"), row.names = FALSE)
print(key_tab %>% select(run, gene, significant) %>% pivot_wider(names_from = run, values_from = significant))

p_key <- key_tab %>% mutate(run = factor(run, levels = names(runs)), gene = factor(gene, levels = rev(key)),
                            status = ifelse(is.na(p), "Not tested", ifelse(significant, "Significant", "NS"))) %>%
  ggplot(aes(run, gene, fill = status)) + geom_tile(colour = "white") +
  geom_text(aes(label = ifelse(is.na(Z), "", Z)), size = 3) +
  scale_fill_manual(values = c(Significant = "#00768E", NS = "grey85", "Not tested" = "white"), name = NULL) +
  labs(x = NULL, y = NULL, title = "Predicted PHLDA1 targets across runs (numbers = Z)") +
  theme_minimal(10) + theme(axis.text.y = element_text(face = "italic"), axis.text.x = element_text(angle = 30, hjust = 1))
save_both(p_key, file.path(od, "S1_key_genes_heatmap"), 6, 4)

## ---- 2. 各次运行与原结果的重叠 ------------------------------------------------------------
cmp <- bind_rows(lapply(setdiff(names(sig), "original"), function(nm) {
  m <- inner_join(runs$original[, c("gene", "Z")], runs[[nm]][, c("gene", "Z")], by = "gene", suffix = c("_orig", "_run"))
  data.frame(run = nm, n_sig_run = length(sig[[nm]]), n_sig_original = length(orig_sig),
             shared = length(intersect(sig[[nm]], orig_sig)),
             pct_original_recovered = round(100 * length(intersect(sig[[nm]], orig_sig)) / max(1, length(orig_sig)), 1),
             jaccard = round(jac(sig[[nm]], orig_sig), 3),
             genes_tested_both = nrow(m),
             spearman_Z = round(cor(m$Z_orig, m$Z_run, method = "spearman", use = "complete.obs"), 3))
}))
write.csv(cmp, file.path(od, "S2_overlap_with_original.csv"), row.names = FALSE)
print(cmp)
writeLines(sapply(names(sig), function(n) sprintf("%s (%d): %s", n, length(sig[[n]]), paste(sig[[n]], collapse = ", "))),
           file.path(od, "S2_sig_genes_per_run.txt"))

## ---- 3. 富集：与原分析相同的 Enrichr 库；不能联网则用 clusterProfiler ------------------------
dbs <- c("GO_Biological_Process_2023", "KEGG_2021_Human", "Reactome_2022")
use_enrichr <- requireNamespace("enrichR", quietly = TRUE) &&
  isTRUE(tryCatch({ enrichR::setEnrichrSite("Enrichr"); length(enrichR::listEnrichrDbs()) > 0 }, error = function(e) FALSE))
message(ifelse(use_enrichr, "使用 Enrichr（与原分析相同）", "Enrichr 不可用，改用 clusterProfiler GO BP"))

enr_new <- list()
for (nm in setdiff(names(sig), "original")) {
  g <- sig[[nm]]; if (length(g) < 5) next
  if (use_enrichr) {
    r <- enrichR::enrichr(g, dbs)
    for (db in dbs) write.csv(r[[db]], file.path(od, paste0("S3_Enrichr_", db, "_", nm, ".csv")), row.names = FALSE)
    enr_new[[nm]] <- r[["GO_Biological_Process_2023"]] %>% transmute(Term, padj = Adjusted.P.value)
  } else {
    suppressPackageStartupMessages({ library(clusterProfiler); library(org.Hs.eg.db) })
    ego <- as.data.frame(enrichGO(g, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP", pvalueCutoff = 0.05))
    write.csv(ego, file.path(od, paste0("S3_clusterProfiler_GO_BP_", nm, ".csv")), row.names = FALSE)
    enr_new[[nm]] <- ego %>% transmute(Term = Description, padj = p.adjust)
  }
}

# 原始 GO BP 前 10 条在新结果中是否仍显著
orig_go <- read.csv(file.path(ORIG_DIR, "Enrichr_GO_Biological_Process_2023.csv"), check.names = FALSE)
pcol_go <- intersect(c("Adjusted P-value", "Adjusted.P.value", "adj_pval", "P.adjust"), colnames(orig_go))[1]
orig_top <- orig_go %>% arrange(.data[[pcol_go]]) %>% head(10) %>% pull(Term)
clean <- function(x) tolower(trimws(sub("\\s*\\(GO:\\d+\\)$", "", x)))
go_cmp <- data.frame(original_top10_term = orig_top)
for (nm in names(enr_new)) {
  e <- enr_new[[nm]]
  go_cmp[[paste0("padj_", nm)]] <- signif(e$padj[match(clean(orig_top), clean(e$Term))], 3)
}
write.csv(go_cmp, file.path(od, "S3_original_top10_GO_in_new_runs.csv"), row.names = FALSE)
print(go_cmp)
# 若用 clusterProfiler，词条名称与 Enrichr 不完全相同，NA 不一定代表“不富集”，请看 S3_*.csv 中是否仍有 ribosome / translation 相关条目

message("summary 完成，结果在：", od)
