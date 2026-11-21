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

## ---- 图 4：GSE106910 严重度 × 性别（修复版） -----------------------------------
if (!exists("validation_expr")) load(EXPR_RDATA)
vi <- as.data.frame(validation_info)
cat("\n[诊断] GSE106910 中可用的性别基因：",
    paste(intersect(c("XIST", "RPS4Y1", "DDX3Y", "KDM5D", "UTY", "EIF1AY", "USP9Y", "ZFY"),
                    rownames(validation_expr)), collapse = ", "), "\n")
sex_cols_meta <- names(vi)[sapply(vi, function(x) any(grepl("\\bmale\\b|female", x, ignore.case = TRUE)))]
cat("[诊断] 元数据中含 male/female 的列：", paste(sex_cols_meta, collapse = ", "), "\n")

s106 <- data.frame(sample = colnames(validation_expr))
s106$Group <- vi$Group[match(s106$sample, vi$`!Sample_title`)]
cat("[诊断] Group 匹配情况：\n"); print(table(s106$Group, useNA = "ifany"))

# 性别：元数据优先；否则 XIST（+Y 基因）表达推断
if (length(sex_cols_meta)) {
  raw <- vi[[sex_cols_meta[1]]][match(s106$sample, vi$`!Sample_title`)]
  s106$sex <- ifelse(grepl("female", raw, ignore.case = TRUE), "Female",
                     ifelse(grepl("male", raw, ignore.case = TRUE), "Male", NA))
} else s106$sex <- NA
if (all(is.na(s106$sex)) && "XIST" %in% rownames(validation_expr)) {
  y <- intersect(c("RPS4Y1", "DDX3Y", "KDM5D", "UTY", "EIF1AY", "USP9Y", "ZFY"), rownames(validation_expr))
  sc <- as.numeric(scale(as.numeric(validation_expr["XIST", ])))
  if (length(y)) sc <- sc - as.numeric(scale(colMeans(as.matrix(validation_expr[y, , drop = FALSE]))))
  km <- kmeans(sc, centers = 2, nstart = 20)
  fem <- which.max(tapply(sc, km$cluster, mean))
  s106$sex <- ifelse(km$cluster == fem, "Female", "Male")
  message("GSE106910 性别由表达推断（kmeans 两类）")
}
s106 <- s106 %>% filter(!is.na(Group), Group != "N/A") %>%
  mutate(Group = factor(Group, levels = c("No BPD", "mild", "moderate", "severe")),
         XIST = as.numeric(validation_expr["XIST", sample]),
         PHLDA1 = as.numeric(validation_expr["PHLDA1", sample]))
tab_sev <- table(Severity = s106$Group, Sex = s106$sex)
print(tab_sev)
write.csv(as.data.frame.matrix(tab_sev), file.path(od, "4a_GSE106910_sex_by_severity.csv"))
write.csv(s106, file.path(od, "4d_GSE106910_sample_table.csv"), row.names = FALSE)

if (length(unique(na.omit(s106$sex))) == 2) {
  m1 <- lm(PHLDA1 ~ I(Group != "No BPD") + sex, data = s106)
  m2 <- lm(PHLDA1 ~ as.numeric(Group) + sex, data = s106)
  capture.output(print(fisher.test(tab_sev)), print(summary(m1)$coefficients), print(summary(m2)$coefficients),
                 file = file.path(od, "4b_GSE106910_PHLDA1_sex_adjusted.txt"))
  p4a <- as.data.frame(tab_sev) %>% ggplot(aes(Severity, Freq, fill = Sex)) +
    geom_col(width = 0.7, colour = "white") + scale_fill_manual(values = sex_cols) +
    labs(x = NULL, y = "Samples", title = "GSE106910: sex by BPD severity") + theme_classic(10)
  p4b <- ggplot(s106, aes(Group, XIST)) + geom_boxplot(outlier.shape = NA, fill = "grey95", colour = "grey50") +
    geom_jitter(aes(colour = sex), width = 0.15, size = 2) + scale_colour_manual(values = sex_cols) +
    labs(x = NULL, title = "XIST by severity") + theme_classic(10)
  p4c <- ggplot(s106, aes(Group, PHLDA1)) + geom_boxplot(outlier.shape = NA, fill = "grey95", colour = "grey50") +
    geom_jitter(aes(colour = sex), width = 0.15, size = 2) + scale_colour_manual(values = sex_cols) +
    labs(x = NULL, title = "PHLDA1 by severity") + theme_classic(10)
  save_both(p4a | p4b | p4c, file.path(od, "5h_GSE106910_sex_severity_XIST_PHLDA1"), 13, 4)
}
