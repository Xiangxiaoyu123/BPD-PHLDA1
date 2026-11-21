###############################################################################
# B4  拟时序的供者拆分 + 留一供者（leave-one-donor-out, LODO）敏感性分析   —— 回应 R1-2
#
# 注意：原代码（r.03 第 4236-4306 行）实际使用的是 monocle (v2) 的 DDRTree + orderCells(root_state = 2)，
#       而不是 Monocle3。稿件 Methods 与参考文献需相应改为 “Monocle 2 (DDRTree)”。
###############################################################################
source(here::here("revision_scripts", "00_setup.R"))
suppressPackageStartupMessages({ library(monocle); library(Biobase); library(ggridges) })
od <- file.path(OUT, "B4_trajectory")

if (!exists("mycds3")) load(PATH_MONOCLE)       # 得到 mycds3
pd <- pData(mycds3)
pd$donor <- if ("dataset" %in% colnames(pd)) as.character(pd$dataset) else sub("_.*$", "", pd$orig.ident)
pd$cell  <- rownames(pd)
pd$condition <- factor(pd$condition, levels = c("Term infant", "BPD", "BPD+PH"))

## ---- 1. 每个供者的拟时序分布 --------------------------------------------------
p1 <- ggplot(pd, aes(Pseudotime, factor(donor, levels = rev(DONOR_ORDER)), fill = condition)) +
  geom_density_ridges(scale = 1.2, alpha = 0.8, quantile_lines = TRUE, quantiles = 2) +
  scale_fill_manual(values = c("Term infant" = "#3182bd", "BPD" = "#e6550d", "BPD+PH" = "#756bb1")) +
  labs(y = NULL, title = "Pseudotime distribution of adventitial FB per donor") + theme_classic(11)
save_both(p1, file.path(od, "1_pseudotime_ridges_per_donor"), 7, 4.5)

pData(mycds3)$donor <- pd$donor
p1b <- plot_cell_trajectory(mycds3, color_by = "donor", cell_size = 0.6) + facet_wrap(~ donor, nrow = 2)
save_both(p1b, file.path(od, "1b_trajectory_split_by_donor"), 11, 6)

med <- pd %>% group_by(donor, condition) %>%
  summarise(n = n(), median_pt = median(Pseudotime), IQR_low = quantile(Pseudotime, .25),
            IQR_high = quantile(Pseudotime, .75), .groups = "drop")
write.csv(med, file.path(od, "1c_pseudotime_summary_per_donor.csv"), row.names = FALSE)
print(med)

## ---- 2. 稿件中的早/晚期区间组成（94% Term at 5–10；88% BPD at 15–20） --------
# 若稿件中的 "interval" 不是原始 Pseudotime 数值，请把下面两个区间改成当时使用的定义
PT_EARLY <- c(5, 10); PT_LATE <- c(15, 20)
interval_comp <- function(d) {
  bt <- d %>% filter(condition %in% c("Term infant", "BPD"))
  e <- bt %>% filter(Pseudotime >= PT_EARLY[1], Pseudotime <= PT_EARLY[2])
  l <- bt %>% filter(Pseudotime >= PT_LATE[1],  Pseudotime <= PT_LATE[2])
  c(early_pct_Term = 100 * mean(e$condition == "Term infant"), early_n = nrow(e),
    late_pct_BPD   = 100 * mean(l$condition == "BPD"),         late_n  = nrow(l),
    median_pt_BPD_minus_Term = median(bt$Pseudotime[bt$condition == "BPD"]) -
                               median(bt$Pseudotime[bt$condition == "Term infant"]))
}

## ---- 3. LODO：依次去掉一个 Term/BPD 供者，重算三项指标 -------------------------
adv  <- add_donor(ensure_rna_norm(meso_adventitial_FB))
meso <- add_donor(meso_relabel_no_preterm)

## 确认两个对象都有 donor 和 group 两列，且供者名字一致
stopifnot(all(c("donor", "group") %in% colnames(meso@meta.data)),
          all(c("donor", "group") %in% colnames(adv@meta.data)))
table(meso$donor, meso$group)
table(adv$donor,  adv$group)

## v5 多 layer 时先合并，否则 FoldChange / pseudobulk 会报 GetAssayData 错误
if (length(Layers(adv[["RNA"]], search = "counts")) > 1) adv[["RNA"]] <- JoinLayers(adv[["RNA"]])

bt_donors <- c("Term infant 1", "Term infant 2", "BPD 1", "BPD 2")
stopifnot(all(bt_donors %in% unique(as.character(meso$donor))))
lodo_metrics <- function(drop = NULL) {
  keep_m <- !(as.character(meso$donor) %in% drop) & meso$group %in% c("Term infant", "BPD")
  md <- meso@meta.data[keep_m, ]
  pct <- tapply(md$celltype == "Adventitial FB", droplevels(md$group), mean) * 100

  keep_a <- !(as.character(adv$donor) %in% drop) & adv$group %in% c("Term infant", "BPD")
  sub_a  <- subset(adv, cells = colnames(adv)[keep_a])
  Idents(sub_a) <- "group"
  fc <- FoldChange(sub_a, ident.1 = "BPD", ident.2 = "Term infant",
                   features = intersect(CORE_GENES, rownames(sub_a)), assay = "RNA")
  pb  <- pseudobulk(sub_a, "donor"); cpm <- log2(t(t(pb) / colSums(pb)) * 1e6 + 1)
  g   <- ifelse(grepl("^BPD", colnames(cpm)), "BPD", "Term")
  pb_lfc <- mean(cpm["PHLDA1", g == "BPD"]) - mean(cpm["PHLDA1", g == "Term"])

  ic <- interval_comp(pd %>% filter(!donor %in% drop))
  data.frame(dropped = ifelse(is.null(drop), "none", drop),
             advFB_pct_Term = round(pct["Term infant"], 2), advFB_pct_BPD = round(pct["BPD"], 2),
             advFB_fold = round(pct["BPD"] / pct["Term infant"], 2),
             PHLDA1_cell_log2FC = round(fc["PHLDA1", "avg_log2FC"], 3),
             PHLDA1_pseudobulk_log2CPM_diff = round(pb_lfc, 3),
             t(round(ic, 2)), row.names = NULL)
}
lodo <- bind_rows(lodo_metrics(NULL), lapply(bt_donors, lodo_metrics))
write.csv(lodo, file.path(od, "3_LODO_summary.csv"), row.names = FALSE)
print(lodo)

## ---- 4. 所有 BPD × Term 供者两两配对（共 4 对）是否方向一致 --------------------
pairs <- expand.grid(BPD = c("BPD 1", "BPD 2"), Term = c("Term infant 1", "Term infant 2"),
                     stringsAsFactors = FALSE)
pair_res <- bind_rows(lapply(seq_len(nrow(pairs)), function(i) {
  d <- c(pairs$BPD[i], pairs$Term[i])
  pm <- meso@meta.data[as.character(meso$donor) %in% d, ]
  pct <- tapply(pm$celltype == "Adventitial FB", as.character(pm$donor), mean) * 100
  sa <- subset(adv, cells = colnames(adv)[as.character(adv$donor) %in% d]); Idents(sa) <- "donor"
  fc <- FoldChange(sa, ident.1 = d[1], ident.2 = d[2], features = "PHLDA1", assay = "RNA")
  data.frame(BPD = d[1], Term = d[2], advFB_pct_BPD = round(pct[d[1]], 2),
             advFB_pct_Term = round(pct[d[2]], 2), PHLDA1_log2FC = round(fc$avg_log2FC, 3),
             median_pt_diff = round(median(pd$Pseudotime[pd$donor == d[1]]) -
                                    median(pd$Pseudotime[pd$donor == d[2]]), 2))
}))
write.csv(pair_res, file.path(od, "4_pairwise_donor_comparisons.csv"), row.names = FALSE)
print(pair_res)
