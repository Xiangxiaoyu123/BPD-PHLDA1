###############################################################################
# B4a  补救：重建 Monocle 2 拟时序对象 mycds3（原文件 mycds_order_cells_cds_2_final.RData 丢失）
#
# 步骤
#   0) 先在新旧两个目录里找有没有残留文件：
#        - mycds_order_cells_cds_2_final.RData（最终对象，找到就直接用）
#        - mycds_order_cells_cds.RData（orderCells 之前的 mycds2，找到只需重跑最后两步，几分钟）
#        - train.key.cells.DEG.xls（当时用于排序的基因，找到可保证用完全相同的基因）
#   1) 都没有：完全按原代码（r.03 第 4236-4306 行）从 meso_adventitial_FB 重建
#   2) 根节点：原代码 root_state = 2。重建后 State 编号可能变化，
#      因此改为“自动选 Term infant 占比最高的 State 作根”，并打印各 State 组成供人工核对
#   3) 与稿件 Fig. 4G 的数字核对（94% Term at 5–10；88% BPD at 15–20），判断重建是否等价
#   4) 保存新对象 + 导出 pData（以后就算 RData 丢失也能直接画 B4 第 1 部分）
#
# 环境：monocle (v2, Bioconductor) —— 原分析在 conda 环境 sc_analysis 中运行，建议仍在该环境运行
#       后台运行：nohup Rscript B4a_rebuild_monocle2.R > B4a.log 2>&1 &
###############################################################################
source(here::here("revision_scripts", "00_setup.R"))
suppressPackageStartupMessages({ library(monocle); library(Biobase) })

MONO_DIR <- file.path(BASE, "04_Adventitial_FB/08_Monocle3")
dir.create(MONO_DIR, recursive = TRUE, showWarnings = FALSE)
OUT_RDATA <- file.path(MONO_DIR, "mycds_order_cells_cds_2_final_REBUILT.RData")
N_CORES   <- max(1, parallel::detectCores() - 2)

## ---- 0. 查找残留文件（新服务器目录 + 旧 /data/nas1 目录） ----------------------
search_dirs <- unique(c(BASE, "/home/xxyloveme/study/BPD+MSCs"))
search_dirs <- search_dirs[dir.exists(search_dirs)]
find_file <- function(pattern) {
  hits <- unlist(lapply(search_dirs, list.files, pattern = pattern, recursive = TRUE, full.names = TRUE))
  unique(hits)
}
f_final <- find_file("^mycds_order_cells_cds_2_final\\.RData$")
f_mid   <- find_file("^mycds_order_cells_cds\\.RData$")
f_genes <- find_file("^train\\.key\\.cells\\.DEG\\.xls$")
cat("找到的最终对象：", if (length(f_final)) f_final else "无", "\n")
cat("找到的中间对象：", if (length(f_mid))   f_mid   else "无", "\n")
cat("找到的排序基因：", if (length(f_genes)) f_genes else "无", "\n")

meso_adventitial_FB <- readRDS("/home/xxyloveme/study/BPD+MSCs/04_Adventitial_FB/meso_adventitial_FB.rds")

## ---- 1. 得到 mycds2（DDRTree 降维后的对象） ------------------------------------
if (length(f_final)) {
  load(f_final[1])                              # 直接得到 mycds3，跳到第 3 步
  message("已找到原始最终对象，无需重建。")
} else {
  if (length(f_mid)) {
    load(f_mid[1])                              # 得到 mycds2
    message("已找到 mycds2，只重跑 orderCells。")
  } else {
    message("未找到任何残留文件，按原代码完整重建（约 0.5–2 小时）。")
    obj <- meso_adventitial_FB
    DefaultAssay(obj) <- "RNA"
    if (inherits(obj[["RNA"]], "Assay5") && length(Layers(obj[["RNA"]], search = "counts")) > 1)
      obj[["RNA"]] <- JoinLayers(obj[["RNA"]])

    # 原代码中的 cell_subtype：Adventitial_FB_Cluster0 ~ 5
    obj$cell_subtype <- ifelse(!is.na(obj$SCT_snn_res.0.5) & obj$celltype == "Adventitial FB",
                               paste0("Adventitial_FB_Cluster", obj$SCT_snn_res.0.5),
                               as.character(obj$celltype))

    mat <- GetAssayData(obj, assay = "RNA", slot = "counts")
    fd  <- new("AnnotatedDataFrame", data = data.frame(gene_id = rownames(mat), gene_short_name = rownames(mat),
                                                       row.names = rownames(mat)))
    pd  <- new("AnnotatedDataFrame", data = obj@meta.data)
    cds <- newCellDataSet(mat, phenoData = pd, featureData = fd, expressionFamily = negbinomial.size())
    rm(obj, mat); gc()
    cds <- estimateSizeFactors(cds)
    cds <- estimateDispersions(cds)
    mycds <- detectGenes(cds, min_expr = 0.1)
    expressed_genes <- rownames(subset(fData(mycds), num_cells_expressed >= 10))

    if (length(f_genes)) {                      # 用当时保存的排序基因 —— 最接近原结果
      deg <- read.delim(f_genes[1])
      ordergene <- intersect(if ("gene_short_name" %in% colnames(deg)) deg$gene_short_name else deg[[1]],
                             rownames(mycds))
      message("使用原排序基因 ", length(ordergene), " 个")
    } else {                                    # 否则按原代码重新做 differentialGeneTest
      diff_test_res <- differentialGeneTest(mycds[expressed_genes, ],
                                            fullModelFormulaStr = "~cell_subtype", cores = N_CORES)
      deg <- subset(diff_test_res, qval < 0.05)
      deg <- deg[order(deg$qval), ]
      write.table(deg, file.path(MONO_DIR, "train.key.cells.DEG_REBUILT.xls"),
                  col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
      ordergene <- rownames(deg)
      message("重新计算排序基因 ", length(ordergene), " 个")
    }
    mycds1 <- setOrderingFilter(mycds, ordergene)
    mycds2 <- reduceDimension(mycds1, max_components = 2, method = "DDRTree")
    save(mycds2, file = file.path(MONO_DIR, "mycds_order_cells_cds_REBUILT.RData"))
  }

  ## ---- 2. orderCells + 选根 ---------------------------------------------------
  mycds3 <- orderCells(mycds2)
  st_comp <- as.data.frame.matrix(table(State = pData(mycds3)$State, Condition = pData(mycds3)$condition))
  st_comp$n <- rowSums(st_comp)
  st_comp$pct_Term <- round(100 * st_comp[["Term infant"]] / st_comp$n, 1)
  cat("\n各 State 的分组构成（原代码以 State 2 为根）：\n"); print(st_comp)
  root <- rownames(st_comp)[which.max(st_comp$pct_Term * (st_comp$n >= 30))]   # 至少 30 个细胞
  message("选择 State ", root, " 作为根（Term infant 占比最高）")
  mycds3 <- orderCells(mycds3, root_state = root)
  save(mycds3, file = OUT_RDATA)
  writeLines(capture.output(print(st_comp), cat("root_state =", root, "\n")),
             file.path(MONO_DIR, "REBUILT_state_composition_and_root.txt"))
}

## ---- 3. 与稿件 Fig. 4G 核对 ----------------------------------------------------
pd3 <- pData(mycds3)
pd3$donor <- if ("dataset" %in% colnames(pd3)) as.character(pd3$dataset) else sub("_.*$", "", pd3$orig.ident)
bt <- pd3[pd3$condition %in% c("Term infant", "BPD"), ]
chk <- data.frame(
  check = c("Pseudotime range", "Term % among cells with pseudotime 5–10 (manuscript: 94%)",
            "BPD % among cells with pseudotime 15–20 (manuscript: 88%)",
            "median pseudotime Term", "median pseudotime BPD"),
  value = c(paste(round(range(pd3$Pseudotime), 2), collapse = " – "),
            round(100 * mean(bt$condition[bt$Pseudotime >= 5  & bt$Pseudotime <= 10] == "Term infant"), 1),
            round(100 * mean(bt$condition[bt$Pseudotime >= 15 & bt$Pseudotime <= 20] == "BPD"), 1),
            round(median(bt$Pseudotime[bt$condition == "Term infant"]), 2),
            round(median(bt$Pseudotime[bt$condition == "BPD"]), 2)))
print(chk)
write.csv(chk, file.path(MONO_DIR, "REBUILT_check_vs_manuscript.csv"), row.names = FALSE)

## 与原 Fig. 4H/4I 对照：各 State / 亚簇沿拟时序的分布
p_chk <- plot_cell_trajectory(mycds3, color_by = "Pseudotime", cell_size = 0.6) |
         plot_cell_trajectory(mycds3, color_by = "condition", cell_size = 0.6) |
         plot_cell_trajectory(mycds3, color_by = "State", cell_size = 0.6)
save_both(p_chk, file.path(MONO_DIR, "REBUILT_trajectory_check"), 15, 5)

## ---- 4. 导出 pData，防止再次丢失 ------------------------------------------------
write.csv(pd3[, intersect(c("donor", "condition", "orig.ident", "SCT_snn_res.0.5", "cell_subtype",
                            "Pseudotime", "State"), colnames(pd3))],
          file.path(MONO_DIR, "monocle2_pseudotime_per_cell.csv"))
message("完成。B4 中把 PATH_MONOCLE 改为：", if (length(f_final)) f_final[1] else OUT_RDATA)
