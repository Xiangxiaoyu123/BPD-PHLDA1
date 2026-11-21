###############################################################################
# B5  供者性别核查 + 性染色体基因敏感性分析   —— 回应 R2-2（XIST 提示性别效应）
#   1) 用 XIST / Y 染色体基因表达判断 GSE275938 每个供者的性别
#   2) 定义性染色体基因列表（gencode v38 GTF）
#   3) 剔除性染色体基因后：107 个 disease-trajectory 基因 / ME0-ME8 候选基因有何变化
#   4) Bulk 数据（GSE106910、GSE235562）：性别推断、按性别校正后 PHLDA1 是否仍差异、
#      Fig. 8D 各严重度组的性别构成（解释 XIST 的“严重度梯度”）
#   scTenifoldKnk 的重跑见 B5b_scTenifoldKnk.R
###############################################################################
source(here::here("revision_scripts", "00_setup.R"))
od <- file.path(OUT, "B5_sex")
Y_GENES_KEY <- c("RPS4Y1", "DDX3Y", "KDM5D", "UTY", "EIF1AY", "USP9Y", "ZFY")

## ---- 1. 单细胞供者性别 --------------------------------------------------------
use_all <- file.exists(PATH_ALLDATA)
obj_for_sex <- if (use_all) add_donor(readRDS(PATH_ALLDATA)) else meso_relabel
pseudobulk <- function(obj, group_by, assay = "RNA") {
  a <- obj[[assay]]
  if (inherits(a, "Assay5") && length(Layers(a, search = "counts")) > 1)
    a <- JoinLayers(a)
  cnt <- if (inherits(a, "Assay5")) LayerData(a, layer = "counts")
  else GetAssayData(obj, assay = assay, layer = "counts")
  grp <- factor(obj[[group_by, drop = TRUE]])
  mm  <- Matrix::sparse.model.matrix(~ 0 + grp)
  colnames(mm) <- levels(grp)
  as.matrix(cnt %*% mm)        # 基因 × 供者 的 counts 矩阵
}
pb  <- pseudobulk(obj_for_sex, "donor")
cpm <- t(t(pb) / colSums(pb)) * 1e6
gx  <- intersect(c("XIST", "TSIX", Y_GENES_KEY), rownames(cpm))
sex_tab <- as.data.frame(t(round(cpm[gx, , drop = FALSE], 2))) %>% rownames_to_column("donor")
obj_for_sex[["RNA"]] <- JoinLayers(obj_for_sex[["RNA"]])
xist_frac <- tapply(GetAssayData(obj_for_sex, assay = "RNA", layer = "counts")["XIST", ] > 0,
                    obj_for_sex$donor, mean)
sex_tab$pct_cells_XIST_pos <- round(100 * xist_frac[sex_tab$donor], 1)
sex_tab$Y_sum_CPM <- rowSums(sex_tab[, intersect(Y_GENES_KEY, colnames(sex_tab)), drop = FALSE])
sex_tab$inferred_sex <- with(sex_tab, ifelse(XIST > 10 * pmax(Y_sum_CPM, 1) , "Female",
                                      ifelse(Y_sum_CPM > 10 * pmax(XIST, 1), "Male", "Check")))
sex_tab$group <- as.character(meso_relabel$group[match(sex_tab$donor, meso_relabel$donor)])
sex_tab$GEO <- GSM_MAP[sex_tab$donor]
write.csv(sex_tab, file.path(od, "1_GSE275938_donor_sex_inference.csv"), row.names = FALSE)
cat("\n== GSE275938 供者性别推断 ==\n"); print(sex_tab[, c("donor", "group", "XIST", "Y_sum_CPM", "inferred_sex")])

# 额外提示：BPD 1/2 来自 7mo-1/7mo-2，BPD+PH 1/2 来自 5796-1/5796-2，命名上像同一患儿的两个样本。
# 若同组两样本性别相同且表达高度相关，请与原文 Methods / 原作者核实它们是否为不同个体，
# 否则“每组 2 个生物学重复”的表述需改为“2 个样本”。
cor_pb <- cor(log2(cpm + 1), method = "spearman")
write.csv(round(cor_pb, 4), file.path(od, "1b_donor_pseudobulk_correlation.csv"))

## 莫兰迪色系：XIST 用暖玫瑰色，Y 染色体基因用冷色调渐变
y_genes <- intersect(Y_GENES_KEY, colnames(sex_tab))
morandi_y <- c("#5E7C8C", "#7895A2", "#93AEB5", "#A9BFB4",
               "#8FA89A", "#B7C4C9", "#C9D3CC", "#6F8F86")
gene_cols <- c(XIST = "#C47F7A", setNames(morandi_y[seq_along(y_genes)], y_genes))

p1 <- sex_tab %>% select(donor, group, XIST, all_of(y_genes)) %>%
  pivot_longer(-c(donor, group), names_to = "gene", values_to = "CPM") %>%
  mutate(gene = factor(gene, levels = c("XIST", y_genes))) %>%
  ggplot(aes(factor(donor, levels = DONOR_ORDER), log2(CPM + 1), fill = gene)) +
  geom_col(position = position_dodge(width = 0.85), width = 0.8,
           colour = "white", linewidth = 0.2) +
  scale_fill_manual(values = gene_cols, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = NULL, y = "log2(CPM+1)", title = "Sex-marker genes per donor") +
  theme_classic(10) +
  theme(axis.text.x  = element_text(angle = 45, hjust = 1, colour = "black"),
        axis.text.y  = element_text(colour = "black"),
        plot.title   = element_text(face = "bold"),
        legend.key.size = unit(0.35, "cm"),
        legend.text  = element_text(face = "italic"))

save_both(p1, file.path(od, "1c_sex_markers_per_donor"), 9, 4.5)
if (use_all) { rm(obj_for_sex); gc() }

## ---- 2. 性染色体基因列表 ------------------------------------------------------
gtf <- rtracklayer::import(PATH_GTF)
gtf <- as.data.frame(gtf[gtf$type == "gene"])
chrY_genes <- unique(gtf$gene_name[gtf$seqnames == "chrY"])
chrX_genes <- unique(gtf$gene_name[gtf$seqnames == "chrX"])
x_escape <- c("XIST", "TSIX", "JPX", "FTX", "KDM5C", "KDM6A", "DDX3X", "EIF1AX", "ZFX", "USP9X",
              "PRKX", "NLGN4X", "TXLNG", "SMC1A", "RPS4X", "EIF2S3", "CA5B", "STS", "PUDP")
sex_genes_strict <- unique(c(chrY_genes, x_escape))           # 主分析用
sex_genes_broad  <- unique(c(chrY_genes, chrX_genes))         # 补充：剔除全部 X/Y 基因
saveRDS(sex_genes_strict, file.path(od, "sex_genes_strict.rds"))
saveRDS(sex_genes_broad,  file.path(od, "sex_genes_broad.rds"))

## ---- 3. 对候选基因链条的影响 ---------------------------------------------------
traj <- lapply(PATH_TRAJ, read.csv)
traj_genes <- unique(unlist(lapply(seq_along(traj), function(i) traj[[i]][[paste0("Cluster", i, "_geneName")]])))
deg_genes  <- unique(significant_genes$gene)
g107       <- intersect(deg_genes, traj_genes)
g107_nosex <- intersect(setdiff(deg_genes, sex_genes_strict), traj_genes)
cat(sprintf("\n原 disease-trajectory 基因 = %d（稿件为 107），剔除性染色体基因后 = %d\n",
            length(g107), length(g107_nosex)))
cat("被剔除的基因：", paste(setdiff(g107, g107_nosex), collapse = ", "), "\n")
cat("核心 8 基因中位于性染色体的：", paste(intersect(CORE_GENES, sex_genes_broad), collapse = ", "), "(应为空)\n")

me0 <- file.path(PATH_ML_DIR, "ME0_Hub_Genes.csv"); me8 <- file.path(PATH_ML_DIR, "ME8_Hub_Genes.csv")
ml_input <- if (file.exists(me0)) c(read.csv(me0)$Gene, read.csv(me8)$Gene) else character(0)
sex_in_ml <- intersect(ml_input, sex_genes_strict)
writeLines(c(sprintf("disease-trajectory genes: %d; without sex genes: %d", length(g107), length(g107_nosex)),
             paste("removed:", paste(setdiff(g107, g107_nosex), collapse = ", ")),
             paste("ML input genes (ME0+ME8):", length(ml_input)),
             paste("sex genes in ML input:", paste(sex_in_ml, collapse = ", ")),
             "",
             "如果上面一行不为空（例如 XIST），请重跑 r.03 第 3.8 节机器学习：",
             "在第 5817 行 available_genes <- c(...) 之后加一行：",
             "  available_genes <- setdiff(available_genes, readRDS('<OUT>/B5_sex/sex_genes_strict.rds'))",
             "其余代码不变，比较新的四算法交集是否仍为 AES/TLE5, ACTG1, TSC22D2, HSPA1A, PHLDA1, SLC26A7, MT1A, NFAT5。"),
           file.path(od, "3_candidate_chain_without_sex_genes.txt"))
write.csv(data.frame(gene = g107_nosex), file.path(od, "3b_disease_trajectory_genes_no_sex.csv"), row.names = FALSE)

## ---- 4. Bulk 数据：性别与 PHLDA1 / XIST ---------------------------------------
if (!exists("validation_expr") && file.exists(EXPR_RDATA)) load(EXPR_RDATA)
library(limma)

infer_sex_bulk <- function(expr) {
  y <- intersect(Y_GENES_KEY, rownames(expr))
  d <- data.frame(sample = colnames(expr), XIST = as.numeric(expr["XIST", ]),
                  Yscore = colMeans(as.matrix(expr[y, , drop = FALSE])))
  # 两个指标都标准化后取差值，>0 判为女性
  d$sex_inferred <- ifelse(as.numeric(scale(d$XIST) - scale(d$Yscore)) > 0, "Female", "Male")
  d
}

# ### 4a. GSE106910（tracheal-aspirate MSC；No BPD / mild / moderate / severe）
# hit <- sapply(vi, function(col) any(grepl("\\bmale\\b|female|\\bsex\\b|gender", col, ignore.case = TRUE)))
# names(hit)[hit]
# if (any(hit)) print(head(vi[, hit, drop = FALSE]))
# 
# sex_col
# if (length(sex_col)) print(head(vi[, sex_col, drop = FALSE]))
# table(s106$sex_metadata, useNA = "ifany")
# table(s106$sex_inferred, useNA = "ifany")
# intersect(c("XIST","RPS4Y1","DDX3Y","KDM5D","UTY","EIF1AY","USP9Y","ZFY"), rownames(validation_expr))
# 
# 
# if (exists("validation_expr")) {
#   vi <- as.data.frame(validation_info)
#   sex_col <- grep("sex|gender", colnames(vi), ignore.case = TRUE, value = TRUE)
#   s106 <- infer_sex_bulk(validation_expr)
#   s106$Group <- vi$Group[match(s106$sample, vi$`!Sample_title`)]
#   if (length(sex_col)) s106$sex_metadata <- vi[[sex_col[1]]][match(s106$sample, vi$`!Sample_title`)]
#   s106 <- s106 %>% filter(!is.na(Group), Group != "N/A")
#   s106$Group <- factor(s106$Group, levels = c("No BPD", "mild", "moderate", "severe"))
#   s106$BPD   <- factor(ifelse(s106$Group == "No BPD", "control", "BPD"), levels = c("control", "BPD"))
#   s106$PHLDA1 <- as.numeric(validation_expr["PHLDA1", s106$sample])
#   s106$sex <- if (length(sex_col)) s106$sex_metadata else s106$sex_inferred
# 
#   tab_sev <- table(Severity = s106$Group, Sex = s106$sex)
#   write.csv(as.data.frame.matrix(tab_sev), file.path(od, "4a_GSE106910_sex_by_severity.csv"))
#   cat("\n== GSE106910 各严重度组性别构成（用于解释 Fig. 8D 的 XIST 梯度）==\n"); print(tab_sev)
#   print(fisher.test(tab_sev))
# 
#   m0 <- lm(PHLDA1 ~ BPD, data = s106); m1 <- lm(PHLDA1 ~ BPD + sex, data = s106)
#   m2 <- lm(PHLDA1 ~ as.numeric(Group) + sex, data = s106)          # 严重度趋势 + 性别
#   m3 <- lm(PHLDA1 ~ BPD * sex, data = s106)                        # 性别交互
#   capture.output(
#     cat("PHLDA1 ~ BPD\n"), print(summary(m0)$coefficients),
#     cat("\nPHLDA1 ~ BPD + sex\n"), print(summary(m1)$coefficients),
#     cat("\nPHLDA1 ~ severity(ordinal) + sex\n"), print(summary(m2)$coefficients),
#     cat("\nPHLDA1 ~ BPD * sex\n"), print(summary(m3)$coefficients),
#     file = file.path(od, "4b_GSE106910_PHLDA1_sex_adjusted.txt"))
#   print(summary(m1)$coefficients)
# 
#   s106$XIST <- as.numeric(validation_expr["XIST", s106$sample])
#   p4 <- ggplot(s106, aes(Group, XIST, colour = sex)) +
#     geom_boxplot(outlier.shape = NA, colour = "grey50") + geom_jitter(width = 0.15, size = 2) +
#     labs(x = NULL, title = "GSE106910: XIST by severity, coloured by sex") + theme_bw(11)
#   p5 <- ggplot(s106, aes(Group, PHLDA1, colour = sex)) +
#     geom_boxplot(outlier.shape = NA, colour = "grey50") + geom_jitter(width = 0.15, size = 2) +
#     labs(x = NULL, title = "GSE106910: PHLDA1 by severity, coloured by sex") + theme_bw(11)
#   save_both(p4 | p5, file.path(od, "4c_GSE106910_XIST_PHLDA1_by_sex"), 11, 4.5)
#   write.csv(s106, file.path(od, "4d_GSE106910_sample_table.csv"), row.names = FALSE)
# }

### 4b. GSE235562（MSC 体外 CMS/高氧处理；"training/test" 集）
if (exists("training_expr") && exists("training_info")) {
  s235 <- infer_sex_bulk(training_expr)
  s235$group <- training_info$group[match(s235$sample, rownames(training_info))]
  write.csv(s235, file.path(od, "4e_GSE235562_sex_inference.csv"), row.names = FALSE)
  cat("\n== GSE235562 处理组 × 推断性别 ==\n"); print(table(s235$group, s235$sex_inferred))
  # 若 control 与处理组的性别构成不同，XIST 被选入 WGCNA ME0 即可由此解释
}


###############################################################################
# B5c  性别相关补充图 —— 回应 R2-2（在 B5_sex_check.R 结果基础上）
#
# B5 已得到：GSE275938 中 Term infant 1 = 男，Term infant 2 / BPD 1 / BPD 2 / BPD+PH 1 / BPD+PH 2 = 女，
#            Acute preterm injury 1 = 男。也就是 BPD 组 2 例全为女性，Term 组 1 男 1 女。
#
#   图 1  供者性别判定散点图（XIST vs Y 基因），按分组着色
#   图 2  性别匹配比较（只用女性：BPD 1、BPD 2 vs Term infant 2）——关键分析
#         a) Adventitial FB 占比：全部供者 vs 仅女性供者
#         b) XIST / RPS4Y1 / 8 个核心基因的 log2FC：全部供者 vs 仅女性供者
#   图 3  GSE235562：各处理组性别构成 + XIST 按处理组（点按性别着色）
#   图 4  GSE106910：严重度 × 性别（修复 B5 第 4a 部分输出为空的问题）
###############################################################################
source(here::here("revision_scripts", "00_setup.R"))
od <- file.path(OUT, "B5_sex")

group_cols <- c("Term infant" = "#D7D7D7", "BPD" = "#00768E", "BPD+PH" = "#E59995",
                "Acute preterm injury" = "#B8A07E")
sex_cols   <- c(Female = "#C47F7A", Male = "#5E7C8C")
sex_shapes <- c(Female = 16, Male = 17)

sex_tab <- read.csv(file.path(od, "1_GSE275938_donor_sex_inference.csv"))
donor_sex <- setNames(sex_tab$inferred_sex, sex_tab$donor)

## ---- 图 1：供者性别判定 --------------------------------------------------------
p1 <- sex_tab %>% mutate(group = factor(group, levels = names(group_cols))) %>%
  ggplot(aes(log2(XIST + 1), log2(Y_sum_CPM + 1))) +
  geom_point(aes(fill = group, shape = inferred_sex), size = 4, colour = "grey20", stroke = 0.4) +
  ggrepel::geom_text_repel(aes(label = donor), size = 3, box.padding = 0.5, min.segment.length = 0) +
  scale_fill_manual(values = group_cols, name = "Group") +
  scale_shape_manual(values = c(Female = 21, Male = 24), name = "Inferred sex") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  labs(x = "XIST, log2(CPM + 1)", y = "Y-chromosome genes, log2(summed CPM + 1)",
       title = "GSE275938: donor sex inferred from expression") + theme_classic(11)
save_both(p1, file.path(od, "5a_donor_sex_scatter"), 6.5, 4.8)

## ---- 图 2：性别匹配比较（仅女性供者） ------------------------------------------
female_bt <- c("Term infant 2", "BPD 1", "BPD 2")

# 2a. Adventitial FB 占间充质比例
meso_bt <- meso_relabel_no_preterm@meta.data %>%
  filter(group %in% c("Term infant", "BPD")) %>%
  mutate(donor = as.character(donor), group = as.character(group))
prop_set <- function(d, label) d %>% group_by(group) %>%
  summarise(pct_AdvFB = 100 * mean(celltype == "Adventitial FB"), n_donors = n_distinct(donor), .groups = "drop") %>%
  mutate(comparison = label)
prop_cmp <- bind_rows(prop_set(meso_bt, "All donors (Term: 1 M + 1 F; BPD: 2 F)"),
                      prop_set(filter(meso_bt, donor %in% female_bt), "Female donors only (Term: 1 F; BPD: 2 F)"))
per_donor <- meso_bt %>% group_by(donor, group) %>%
  summarise(pct_AdvFB = 100 * mean(celltype == "Adventitial FB"), .groups = "drop") %>%
  mutate(sex = donor_sex[donor])
write.csv(prop_cmp, file.path(od, "5b_AdvFB_pct_all_vs_female_only.csv"), row.names = FALSE)
print(prop_cmp); print(per_donor)

p2a <- ggplot(per_donor, aes(factor(group, levels = c("Term infant", "BPD")), pct_AdvFB)) +
  stat_summary(aes(fill = group), fun = mean, geom = "col", width = 0.6, colour = "grey30", linewidth = 0.2) +
  geom_point(aes(shape = sex, colour = sex), size = 3.5, position = position_jitter(width = 0.08, seed = 1)) +
  ggrepel::geom_text_repel(aes(label = donor), size = 2.8, min.segment.length = 0) +
  scale_fill_manual(values = group_cols, guide = "none") +
  scale_colour_manual(values = sex_cols, name = "Sex") + scale_shape_manual(values = sex_shapes, name = "Sex") +
  labs(x = NULL, y = "Adventitial FB (% of mesenchymal cells)",
       title = "Expansion is present in the sex-matched comparison",
       subtitle = "Female Term infant 2 vs female BPD 1 and BPD 2") + theme_classic(11)
save_both(p2a, file.path(od, "5c_AdvFB_pct_by_donor_sex"), 5.5, 4.5)

# 2b. 关键基因 log2FC：全部供者 vs 仅女性供者（Adventitial FB，MAST，仅检验这些基因）
adv <- ensure_rna_norm(meso_adventitial_FB)
if (inherits(adv[["RNA"]], "Assay5") && length(Layers(adv[["RNA"]], search = "counts")) > 1)
  adv[["RNA"]] <- JoinLayers(adv[["RNA"]])
key <- intersect(c("XIST", "RPS4Y1", "DDX3Y", CORE_GENES), rownames(adv))
run_fc <- function(donors, label) {
  s <- subset(adv, cells = colnames(adv)[as.character(adv$donor) %in% donors & adv$group %in% c("BPD", "Term infant")])
  Idents(s) <- "group"
  fm <- FindMarkers(s, ident.1 = "BPD", ident.2 = "Term infant", features = key, assay = "RNA",
                    test.use = "MAST", logfc.threshold = 0, min.pct = 0)
  fm$gene <- rownames(fm); fm$comparison <- label; fm
}
fc_all <- run_fc(c("Term infant 1", "Term infant 2", "BPD 1", "BPD 2"), "All donors")
fc_fem <- run_fc(female_bt, "Female donors only")
fc_cmp <- bind_rows(fc_all, fc_fem) %>% select(comparison, gene, avg_log2FC, pct.1, pct.2, p_val_adj)
write.csv(fc_cmp, file.path(od, "5d_key_genes_log2FC_all_vs_female_only.csv"), row.names = FALSE)
print(fc_cmp %>% select(comparison, gene, avg_log2FC, p_val_adj) %>%
        pivot_wider(names_from = comparison, values_from = c(avg_log2FC, p_val_adj)))

p2b <- fc_cmp %>%
  mutate(gene = factor(gene, levels = key),
         is_sex = ifelse(gene %in% c("XIST", "RPS4Y1", "DDX3Y"), "Sex-chromosome genes", "Core genes")) %>%
  ggplot(aes(gene, avg_log2FC, fill = comparison)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75, colour = "grey30", linewidth = 0.2) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  facet_grid(~ is_sex, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = c("All donors" = "#9BB0BC", "Female donors only" = "#C47F7A"), name = NULL) +
  labs(x = NULL, y = "log2FC, BPD vs Term infant (adventitial FB)",
       title = "Sex-chromosome signals disappear in the sex-matched comparison; core genes persist") +
  theme_classic(10) + theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),
                            legend.position = "top")
save_both(p2b, file.path(od, "5e_key_genes_log2FC_all_vs_female_only"), 8, 4.5)

## ---- 图 3：GSE235562 性别构成 + XIST -----------------------------------------
s235 <- read.csv(file.path(od, "4e_GSE235562_sex_inference.csv"))
grp_lv <- c("control", "CMS", "Hox40", "Hox80", "CMS+Hox40", "CMS+Hox80")
s235$group <- factor(s235$group, levels = grp_lv)
tab235 <- table(Treatment = s235$group, Sex = s235$sex_inferred)
print(tab235); print(fisher.test(tab235))
write.csv(as.data.frame.matrix(tab235), file.path(od, "5f_GSE235562_sex_by_treatment.csv"))

p3a <- as.data.frame(tab235) %>% ggplot(aes(Treatment, Freq, fill = Sex)) +
  geom_col(width = 0.7, colour = "white") + scale_fill_manual(values = sex_cols) +
  labs(x = NULL, y = "Samples", title = sprintf("GSE235562 (Fisher P = %.2f)", fisher.test(tab235)$p.value)) +
  theme_classic(10) + theme(axis.text.x = element_text(angle = 30, hjust = 1))
p3b <- ggplot(s235, aes(group, XIST)) +
  geom_boxplot(outlier.shape = NA, fill = "grey95", colour = "grey50") +
  geom_jitter(aes(colour = sex_inferred), width = 0.15, size = 1.8) +
  scale_colour_manual(values = sex_cols, name = "Sex") +
  labs(x = NULL, y = "XIST expression", title = "XIST tracks donor sex, not treatment") +
  theme_classic(10) + theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_both(p3a | p3b, file.path(od, "5g_GSE235562_sex_and_XIST"), 10, 4)


## ---- 图 4：GSE106910 严重度 × 性别（元数据无性别 -> 表达推断；按婴儿汇总） --------
if (!exists("validation_expr")) load(EXPR_RDATA)
vi <- as.data.frame(validation_info)
sex_genes_avail <- intersect(c("XIST", "RPS4Y1", "DDX3Y", "KDM5D", "UTY", "EIF1AY", "USP9Y", "ZFY"),
                             rownames(validation_expr))
cat("\n[诊断] GSE106910 可用的性别基因：", paste(sex_genes_avail, collapse = ", "), "\n")
stopifnot("XIST" %in% sex_genes_avail)

s106 <- data.frame(sample = colnames(validation_expr)) %>%
  mutate(title   = vi$`!Sample_title`[match(sample, vi$`!Sample_title`)],
         Group   = vi$Group[match(sample, vi$`!Sample_title`)],
         subject = toupper(sub("\\s*rep\\d+$", "", title)),            # "exs2 rep1" -> "EXS2"
         XIST    = as.numeric(validation_expr["XIST", sample]),
         PHLDA1  = as.numeric(validation_expr["PHLDA1", sample]))
y <- setdiff(sex_genes_avail, "XIST")
s106$Yscore <- if (length(y)) colMeans(as.matrix(validation_expr[y, s106$sample, drop = FALSE])) else NA

# 样本水平：XIST（和 Y 基因）标准化后差值做 kmeans 两类
sc <- as.numeric(scale(s106$XIST))
if (length(y)) sc <- sc - as.numeric(scale(s106$Yscore))
km <- kmeans(sc, centers = 2, nstart = 50)
fem <- which.max(tapply(sc, km$cluster, mean))
s106$sex_sample <- ifelse(km$cluster == fem, "Female", "Male")

# 婴儿水平：重复样本性别是否一致（质控）+ 多数决
subj <- s106 %>% filter(!is.na(Group)) %>% group_by(subject, Group) %>%
  summarise(n_rep = n(), n_female = sum(sex_sample == "Female"),
            sex = ifelse(n_female > n_rep / 2, "Female", "Male"),
            consistent = n_female %in% c(0, n_rep),
            XIST = mean(XIST), Yscore = mean(Yscore), PHLDA1 = mean(PHLDA1), .groups = "drop")
print(subj, n = Inf)
if (!all(subj$consistent)) warning("有婴儿的重复样本推断性别不一致，请检查：",
                                   paste(subj$subject[!subj$consistent], collapse = ", "))
write.csv(s106, file.path(od, "4d_GSE106910_sample_table.csv"), row.names = FALSE)
write.csv(subj, file.path(od, "4d2_GSE106910_subject_table.csv"), row.names = FALSE)

# 只保留有 BPD 分级的 15 名婴儿
subj_b <- subj %>% filter(Group != "N/A") %>%
  mutate(Group = factor(Group, levels = c("No BPD", "mild", "moderate", "severe")),
         BPD = factor(ifelse(Group == "No BPD", "control", "BPD"), levels = c("control", "BPD")))
tab_sev <- table(Severity = subj_b$Group, Sex = subj_b$sex)
print(tab_sev)
write.csv(as.data.frame.matrix(tab_sev), file.path(od, "4a_GSE106910_sex_by_severity_subjects.csv"))

capture.output(
  cat("Fisher test, severity x sex (subjects)\n"), print(fisher.test(tab_sev)),
  cat("\nPHLDA1 ~ BPD + sex (subject means)\n"),
  print(summary(lm(PHLDA1 ~ BPD + sex, data = subj_b))$coefficients),
  cat("\nPHLDA1 ~ severity(ordinal) + sex (subject means)\n"),
  print(summary(lm(PHLDA1 ~ as.numeric(Group) + sex, data = subj_b))$coefficients),
  cat("\nXIST ~ severity(ordinal) + sex (subject means)\n"),
  print(summary(lm(XIST ~ as.numeric(Group) + sex, data = subj_b))$coefficients),
  file = file.path(od, "4b_GSE106910_subject_level_sex_adjusted.txt"))

# 图：样本水平性别判定（同一婴儿用线连起来）+ 婴儿水平 XIST / PHLDA1
p4a <- s106 %>% filter(!is.na(Group), Group != "N/A") %>%
  ggplot(aes(XIST, Yscore, colour = sex_sample)) + geom_point(size = 2) +
  geom_line(aes(group = subject), colour = "grey70", linewidth = 0.3) +
  ggrepel::geom_text_repel(data = subj_b, aes(XIST, Yscore, label = subject), colour = "grey30",
                           size = 2.6, inherit.aes = FALSE) +
  scale_colour_manual(values = sex_cols, name = "Inferred sex") +
  labs(x = "XIST", y = "Y-chromosome genes (mean)", title = "GSE106910: sex inferred per sample") +
  theme_classic(10)
p4b <- as.data.frame(tab_sev) %>% ggplot(aes(Severity, Freq, fill = Sex)) +
  geom_col(width = 0.7, colour = "white") + scale_fill_manual(values = sex_cols) +
  labs(x = NULL, y = "Infants", title = "Sex by BPD severity (n = 15 infants)") + theme_classic(10)
p4c <- ggplot(subj_b, aes(Group, XIST)) +
  geom_boxplot(outlier.shape = NA, fill = "grey95", colour = "grey50") +
  geom_jitter(aes(colour = sex), width = 0.12, size = 2.5) + scale_colour_manual(values = sex_cols) +
  labs(x = NULL, title = "XIST (infant mean)") + theme_classic(10)
p4d <- ggplot(subj_b, aes(Group, PHLDA1)) +
  geom_boxplot(outlier.shape = NA, fill = "grey95", colour = "grey50") +
  geom_jitter(aes(colour = sex), width = 0.12, size = 2.5) + scale_colour_manual(values = sex_cols) +
  labs(x = NULL, title = "PHLDA1 (infant mean)") + theme_classic(10)
save_both((p4a | p4b) / (p4c | p4d), file.path(od, "5h_GSE106910_sex_severity_XIST_PHLDA1"), 11, 8)

