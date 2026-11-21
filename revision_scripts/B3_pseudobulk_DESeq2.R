###############################################################################
# B3  Pseudo-bulk 差异表达（按供者聚合 counts，DESeq2）   —— 回应 R1-2
#   单细胞水平 MAST 把每个细胞当独立重复，n=2 vs 2 时容易夸大显著性。
#   这里对每个间充质亚群：供者内 counts 求和 -> DESeq2 (~ group)，BPD vs Term infant。
#   并与原 MAST 结果（markers_compare1）比较：方向一致率、仍显著的基因数、核心基因结果。
#
# 建议先运行 B5（性别核查）。若两组性别不平衡，看本脚本第 4 部分“剔除性染色体基因”的结果。
###############################################################################
source(here::here("revision_scripts", "00_setup.R"))
library(DESeq2)
od <- file.path(OUT, "B3_pseudobulk")

bt <- subset(meso_relabel_no_preterm, subset = group %in% c("Term infant", "BPD"))
bt$donor <- droplevels(bt$donor); bt$group <- droplevels(bt$group)
celltypes <- levels(droplevels(bt$celltype))
MIN_CELLS <- 10   # 每个供者在该亚群中至少 10 个细胞才纳入

run_pb <- function(ct, min_cells = MIN_CELLS) {
  obj <- subset(bt, subset = celltype == ct)
  n_per_donor <- table(obj$donor)
  keep_d <- names(n_per_donor)[n_per_donor >= min_cells]
  grp <- obj$group[match(keep_d, obj$donor)]
  if (sum(grp == "BPD") < 2 || sum(grp == "Term infant") < 2) {
    message(ct, ": 供者细胞数不足，跳过 (", paste(names(n_per_donor), n_per_donor, collapse = "; "), ")")
    return(NULL)
  }
  obj <- subset(obj, subset = donor %in% keep_d)
  pb  <- pseudobulk(obj, "donor")
  ## 组名去掉空格（Term infant -> Term），避免 DESeq2 的命名提示
  grp_pb <- as.character(obj$group[match(colnames(pb), obj$donor)])
  cd  <- data.frame(donor = colnames(pb),
                    group = factor(ifelse(grp_pb == "Term infant", "Term", grp_pb),
                                   levels = c("Term", "BPD")),
                    row.names = colnames(pb))
  pb  <- pb[rowSums(pb >= 10) >= 2, , drop = FALSE]
  pb  <- round(pb); storage.mode(pb) <- "integer"        # 直接转整数，避免 converting 提示
  dds <- DESeqDataSetFromMatrix(pb, cd, design = ~ group)
  dds <- DESeq(dds, quiet = TRUE)
  res <- results(dds, contrast = c("group", "BPD", "Term"))
  out <- as.data.frame(res) %>% rownames_to_column("gene") %>%
    mutate(celltype = ct, min_cells_used = min_cells)
  list(res = out, n_cells = n_per_donor,
       vst = assay(varianceStabilizingTransformation(dds, blind = TRUE)))
}

pb_list <- setNames(lapply(celltypes, run_pb), celltypes)
pb_list <- pb_list[!sapply(pb_list, is.null)]
saveRDS(pb_list, file.path(od, "pseudobulk_DESeq2_all_celltypes.rds"))

## ---- 1. 与细胞水平 MAST 结果比较 ---------------------------------------------
compare_ct <- function(ct) {
  mast <- markers_compare1[[ct]]; mast$gene <- rownames(mast)
  mast$sig_cell <- abs(mast$avg_log2FC) > 1 & mast$p_val_adj < 1e-5 & !grepl("^RP|^MT-", mast$gene)
  m <- inner_join(mast[, c("gene", "avg_log2FC", "p_val_adj", "sig_cell")],
                  pb_list[[ct]]$res[, c("gene", "log2FoldChange", "pvalue", "padj")], by = "gene")
  cell_deg <- m %>% filter(sig_cell)
  data.frame(
    celltype = ct,
    cell_level_DEG = sum(mast$sig_cell),
    tested_in_pseudobulk = nrow(cell_deg),
    same_direction_pct = round(100 * mean(sign(cell_deg$avg_log2FC) == sign(cell_deg$log2FoldChange), na.rm = TRUE), 1),
    still_sig_padj0.05 = sum(cell_deg$padj < 0.05, na.rm = TRUE),
    still_sig_pct = round(100 * mean(cell_deg$padj < 0.05, na.rm = TRUE), 1),
    pseudobulk_DEG_total_padj0.05_absLFC1 = sum(pb_list[[ct]]$res$padj < 0.05 &
                                                abs(pb_list[[ct]]$res$log2FoldChange) > 1, na.rm = TRUE),
    spearman_LFC = round(cor(m$avg_log2FC, m$log2FoldChange, method = "spearman", use = "complete.obs"), 3)
  )
}
summ <- bind_rows(lapply(names(pb_list), compare_ct))
write.csv(summ, file.path(od, "1_cell_vs_pseudobulk_summary.csv"), row.names = FALSE)
print(summ)

## ---- 2. Adventitial FB 完整结果表（Supp Table 18） ---------------------------
adv_res <- pb_list[["Adventitial FB"]]$res
mast_adv <- markers_compare1[["Adventitial FB"]] %>% rownames_to_column("gene") %>%
  select(gene, MAST_avg_log2FC = avg_log2FC, MAST_padj = p_val_adj)
supp18 <- adv_res %>% left_join(mast_adv, by = "gene") %>% arrange(padj)
openxlsx::write.xlsx(list(AdventitialFB = supp18,
                          Summary_all_celltypes = summ,
                          Cells_per_donor = as.data.frame(pb_list[["Adventitial FB"]]$n_cells)),
                     file.path(od, "Supplementary_Table_18_pseudobulk.xlsx"))

## ---- 3. 核心基因与 107 个 disease-trajectory 基因 ----------------------------
core <- supp18 %>% filter(gene %in% CORE_GENES)
cat("\n== 8 个核心基因在 Adventitial FB pseudo-bulk 中的结果 ==\n"); print(core)
write.csv(core, file.path(od, "3a_core_genes_pseudobulk_advFB.csv"), row.names = FALSE)

if (exists("ggVenn_candiated_genes")) {
  g107 <- ggVenn_candiated_genes
  r107 <- supp18 %>% filter(gene %in% g107)
  write.csv(r107, file.path(od, "3b_107genes_pseudobulk_advFB.csv"), row.names = FALSE)
  cat(sprintf("\n107 个基因中，pseudo-bulk padj<0.05 的有 %d 个，方向与 MAST 一致 %.1f%%\n",
              sum(r107$padj < 0.05, na.rm = TRUE),
              100 * mean(sign(r107$log2FoldChange) == sign(r107$MAST_avg_log2FC), na.rm = TRUE)))
}

# 每个供者的核心基因表达（VST），作图：每点一个供者
v <- pb_list[["Adventitial FB"]]$vst
vg <- intersect(CORE_GENES, rownames(v))
dfv <- as.data.frame(t(v[vg, , drop = FALSE])) %>% rownames_to_column("donor") %>%
  pivot_longer(-donor, names_to = "gene", values_to = "VST") %>%
  mutate(group = ifelse(grepl("^BPD", donor), "BPD", "Term infant"))
p3 <- ggplot(dfv, aes(group, VST, colour = donor)) + geom_point(size = 3) +
  stat_summary(aes(group = 1), fun = mean, geom = "line", colour = "grey40") +
  facet_wrap(~ gene, scales = "free_y", nrow = 2) +
  labs(x = NULL, y = "Pseudo-bulk VST expression", title = "Adventitial FB: donor-level expression") +
  theme_bw(10)
save_both(p3, file.path(od, "3c_core_genes_donor_level_VST"), 10, 5)

## ---- 4. 剔除性染色体基因后的敏感性（需先运行 B5 生成 sex_genes.rds） ------------
sexf <- file.path(OUT, "B5_sex", "sex_genes_strict.rds")
if (file.exists(sexf)) {
  sex_genes <- readRDS(sexf)
  s2 <- supp18 %>% filter(!gene %in% sex_genes)
  s2$padj_refit <- p.adjust(s2$pvalue, "BH")
  writeLines(c(sprintf("Adventitial FB pseudo-bulk DEG (padj<0.05, |LFC|>1): all genes = %d",
                       sum(supp18$padj < 0.05 & abs(supp18$log2FoldChange) > 1, na.rm = TRUE)),
               sprintf("after removing %d sex-chromosome genes = %d", length(sex_genes),
                       sum(s2$padj_refit < 0.05 & abs(s2$log2FoldChange) > 1, na.rm = TRUE))),
             file.path(od, "4_sexgene_sensitivity.txt"))
}
