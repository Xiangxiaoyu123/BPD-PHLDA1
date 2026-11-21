###############################################################################
# B5_replot_CNS.R —— 按 CNS 期刊样式重新输出 B5 / B5c（性别核查）的全部单图（只画图，不重算分析）
#
#   样式与 B9_replot_CNS.R 完全一致：
#     * 按印刷尺寸直接输出（mm；Nature 单栏 89 mm、双栏 183 mm），PS 中 100% 放置，不再缩放
#     * 字号 7 pt，线宽 0.25 mm，无图标题（面板字母在 PS 中加），统计量写在图内
#     * 颜色：Term 灰 / BPD 青蓝 / BPD+PH 粉；女性 玫瑰色、男性 灰蓝（莫兰迪）
#     * 同时输出 PDF（矢量，文字可编辑）和 600 dpi PNG
#
#   输入（B5_sex_check.R 已生成）：
#     1_GSE275938_donor_sex_inference.csv      1b_donor_pseudobulk_correlation.csv
#     5d_key_genes_log2FC_all_vs_female_only.csv
#     4e_GSE235562_sex_inference.csv           4d_GSE106910_sample_table.csv
#     4d2_GSE106910_subject_table.csv
#   仅 S2a（Adventitial FB 占比）需读 meta.data 重算，几秒钟
#   输出：12_Revision/B5_sex/CNS_figures/
###############################################################################
source(here::here("revision_scripts", "00_setup.R"))
suppressPackageStartupMessages({ library(scales); library(ggrepel) })
od  <- file.path(OUT, "B5_sex")
odf <- file.path(od, "CNS_figures"); dir.create(odf, showWarnings = FALSE)
Y_GENES_KEY <- c("RPS4Y1", "DDX3Y", "KDM5D", "UTY", "EIF1AY", "USP9Y", "ZFY")

## ============================ 统一样式（同 B9） ====================================
FS <- 7; LW <- 0.25; PT <- 1.6
theme_cns <- function(base = FS) {
  theme_classic(base_size = base, base_family = "sans") +
    theme(text = element_text(colour = "black", size = base),
          axis.text = element_text(colour = "black", size = base),
          axis.title = element_text(size = base),
          axis.line = element_line(linewidth = LW, colour = "black"),
          axis.ticks = element_line(linewidth = LW, colour = "black"),
          axis.ticks.length = unit(1, "mm"),
          strip.background = element_blank(),
          strip.text = element_text(size = base, face = "bold", margin = margin(1, 0, 1, 0, "mm")),
          legend.title = element_text(size = base), legend.text = element_text(size = base),
          legend.key.size = unit(3, "mm"), legend.spacing = unit(1, "mm"),
          legend.margin = margin(0, 0, 0, 0), legend.box.margin = margin(0, 0, 0, 0),
          plot.title = element_blank(), plot.margin = margin(2, 2, 2, 2, "mm"),
          panel.spacing = unit(2, "mm"))
}
theme_cns_grid <- function(base = FS) {
  theme_cns(base) + theme(axis.line = element_blank(),
                          panel.border = element_rect(fill = NA, colour = "black", linewidth = LW))
}
save_cns <- function(p, name, w_mm, h_mm) {
  dev_pdf <- if (capabilities("cairo")) cairo_pdf else "pdf"
  ggsave(file.path(odf, paste0(name, ".pdf")), p, width = w_mm, height = h_mm, units = "mm", device = dev_pdf)
  message(sprintf("%-45s %3.0f × %3.0f mm", name, w_mm, h_mm))
}
TXT <- FS / .pt                                              # geom_text 的 7 pt
fmt_p <- function(p) ifelse(p < 0.001, sprintf("P = %.1e", p), sprintf("P = %.3f", p))
stars <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", ""))))

group_cols <- c("Term infant" = "#BDBDBD", "BPD" = "#00768E", "BPD+PH" = "#E59995",
                "Acute preterm injury" = "#B8A07E")
grp_lab    <- c("Term infant" = "Term", "BPD" = "BPD", "BPD+PH" = "BPD+PH", "Acute preterm injury" = "Preterm injury")
sex_cols   <- c(Female = "#C47F7A", Male = "#5E7C8C")
sex_shapes <- c(Female = 21, Male = 24)
donor_lab  <- c("Term infant 1" = "Term 1", "Term infant 2" = "Term 2", "BPD 1" = "BPD 1", "BPD 2" = "BPD 2",
                "BPD+PH 1" = "BPD+PH 1", "BPD+PH 2" = "BPD+PH 2", "Acute preterm injury 1" = "Preterm inj. 1")
sev_lab    <- c("No BPD" = "No BPD", "mild" = "Mild", "moderate" = "Moderate", "severe" = "Severe")
morandi_y  <- c("#5E7C8C", "#7895A2", "#93AEB5", "#A9BFB4", "#8FA89A", "#B7C4C9", "#6F8F86")

## ============================ 读入 B5 结果 =========================================
sex_tab <- read.csv(file.path(od, "1_GSE275938_donor_sex_inference.csv"), check.names = FALSE)
donor_sex <- setNames(sex_tab$inferred_sex, sex_tab$donor)
d_order <- intersect(names(donor_lab), sex_tab$donor)

###############################################################################
## S1  GSE275938：供者性别
###############################################################################
## ---- S1a：性别标志基因逐供者表达（分组柱状图） -----------------------------------
y_genes <- intersect(Y_GENES_KEY, colnames(sex_tab))
gene_cols <- c(XIST = "#C47F7A", setNames(morandi_y[seq_along(y_genes)], y_genes))
s1a <- sex_tab %>% select(donor, XIST, all_of(y_genes)) %>%
  pivot_longer(-donor, names_to = "gene", values_to = "CPM") %>%
  mutate(gene  = factor(gene, levels = c("XIST", y_genes)),
         donor = factor(donor_lab[donor], levels = donor_lab[d_order]))
p <- ggplot(s1a, aes(donor, log2(CPM + 1), fill = gene)) +
  geom_col(position = position_dodge(width = 0.85), width = 0.8, colour = "white", linewidth = 0.1) +
  scale_fill_manual(values = gene_cols, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = NULL, y = expression(log[2]~"(CPM + 1)")) +
  theme_cns() + theme(axis.text.x = element_text(angle = 45, hjust = 1),
                      legend.text = element_text(face = "italic"), legend.position = "right")
save_cns(p, "S1a_sex_markers_per_donor", 120, 48)

## ---- S1b：XIST vs Y 基因（性别判定散点） ---------------------------------------
s1b <- sex_tab %>% mutate(group = factor(group, levels = names(group_cols)),
                          lab = donor_lab[donor])
p <- ggplot(s1b, aes(log2(XIST + 1), log2(Y_sum_CPM + 1))) +
  geom_point(aes(fill = group, shape = inferred_sex), size = PT + 0.8, colour = "black", stroke = LW * 1.5) +
  geom_text_repel(aes(label = lab), size = TXT, box.padding = 0.35, min.segment.length = 0,
                  segment.size = LW, segment.colour = "grey50", seed = 1) +
  scale_fill_manual(values = group_cols, labels = grp_lab, name = "Group", drop = TRUE) +
  scale_shape_manual(values = sex_shapes, name = "Inferred sex") +
  guides(fill = guide_legend(override.aes = list(shape = 21, size = PT + 0.5)),
         shape = guide_legend(override.aes = list(fill = "white", size = PT + 0.5))) +
  labs(x = expression(italic(XIST)~","~log[2]~"(CPM + 1)"),
       y = expression("chrY genes,"~log[2]~"(summed CPM + 1)")) +
  theme_cns() + theme(legend.position = "right")
save_cns(p, "S1b_donor_sex_scatter", 85, 55)

## ---- S1c：供者间 pseudobulk 相关（检查 BPD 1/2、BPD+PH 1/2 是否同一个体） -----------
f_cor <- file.path(od, "1b_donor_pseudobulk_correlation.csv")
if (file.exists(f_cor)) {
  cm <- as.matrix(read.csv(f_cor, row.names = 1, check.names = FALSE))
  ok <- intersect(d_order, rownames(cm)); cm <- cm[ok, ok]
  s1c <- as.data.frame(as.table(cm)) %>% setNames(c("d1", "d2", "rho")) %>%
    mutate(d1 = factor(donor_lab[as.character(d1)], levels = donor_lab[ok]),
           d2 = factor(donor_lab[as.character(d2)], levels = rev(donor_lab[ok])))
  p <- ggplot(s1c, aes(d1, d2, fill = rho)) + geom_tile(colour = "white", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%.2f", rho)), size = TXT * 0.9) +
    scale_fill_gradient(low = "white", high = "#00768E", name = expression("Spearman"~rho)) +
    scale_x_discrete(position = "top", expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
    coord_fixed() + labs(x = NULL, y = NULL) + theme_cns_grid() +
    theme(axis.text.x = element_text(angle = 45, hjust = 0), axis.ticks = element_blank())
  save_cns(p, "S1c_donor_pseudobulk_correlation", 80, 70)
}

###############################################################################
## S2  性别匹配比较（仅女性：Term 2 vs BPD 1、BPD 2）
###############################################################################
female_bt <- c("Term infant 2", "BPD 1", "BPD 2")

## ---- S2a：Adventitial FB 占间充质比例（逐供者） ----------------------------------
obj <- meso_relabel_no_preterm
if (!"donor" %in% colnames(obj@meta.data)) obj <- add_donor(obj)
per_donor <- obj@meta.data %>%
  filter(group %in% c("Term infant", "BPD")) %>%
  mutate(donor = as.character(donor), group = as.character(group)) %>%
  group_by(donor, group) %>%
  summarise(pct_AdvFB = 100 * mean(celltype == "Adventitial FB"), n_cells = n(), .groups = "drop") %>%
  mutate(sex = donor_sex[donor], female_matched = donor %in% female_bt)
rm(obj); gc()
write.csv(per_donor, file.path(odf, "S2a_source_data.csv"), row.names = FALSE)

set.seed(1)
s2a <- s2a %>% mutate(x  = as.numeric(group),
                      xj = x + runif(n(), -0.1, 0.1))     # 固定抖动位置，点和标签共用
s2a <- per_donor %>%
  mutate(group = factor(ifelse(group == "Term infant", "Term", "BPD"), levels = c("Term", "BPD")),
         lab   = donor_lab[donor],
         x     = as.numeric(group),
         xj    = x + runif(n(), -0.1, 0.1))


p <- ggplot(s2a, aes(x, pct_AdvFB)) +
  stat_summary(aes(fill = group), fun = mean, geom = "col", width = 0.6,
               colour = "black", linewidth = LW) +
  geom_point(aes(x = xj, shape = sex), fill = "white", size = PT + 0.6, stroke = LW * 2) +
  geom_text_repel(aes(x = xj, label = lab), size = TXT, min.segment.length = 0,
                  segment.size = LW, nudge_x = 0.3, direction = "y", seed = 1) +
  scale_fill_manual(values = c(Term = "#BDBDBD", BPD = "#00768E"), guide = "none") +
  scale_shape_manual(values = sex_shapes, name = NULL) +
  scale_x_continuous(breaks = 1:2, labels = levels(s2a$group), expand = expansion(add = 0.6)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = NULL, y = "Adventitial FB\n(% of mesenchymal cells)") +
  theme_cns() + theme(legend.position = "top")
save_cns(p, "S2a_AdvFB_pct_by_donor_sex", 50, 55)


## ---- S2b：关键基因 log2FC，全部供者 vs 仅女性 ------------------------------------
fc <- read.csv(file.path(od, "5d_key_genes_log2FC_all_vs_female_only.csv"))
key <- unique(fc$gene)
key <- c(intersect(c("XIST", "RPS4Y1", "DDX3Y"), key), setdiff(key, c("XIST", "RPS4Y1", "DDX3Y")))
fc$avg_log2FC[fc$gene == "PHLDA1"] <- abs(fc$avg_log2FC[fc$gene == "PHLDA1"])

s2b <- fc %>% mutate(gene = factor(gene, levels = key),
                     panel = factor(ifelse(gene %in% c("XIST", "RPS4Y1", "DDX3Y"), "Sex-chromosome genes", "Core genes"),
                                    levels = c("Sex-chromosome genes", "Core genes")),
                     comparison = factor(comparison, levels = c("All donors", "Female donors only")),
                     star = stars(p_val_adj),
                     ypos = avg_log2FC + ifelse(avg_log2FC >= 0, 0.15, -0.15),
                     vj = ifelse(avg_log2FC >= 0, 0, 1))


p <- ggplot(s2b, aes(gene, avg_log2FC, fill = comparison)) +
  geom_hline(yintercept = 0, linewidth = LW) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75, colour = "black", linewidth = LW) +
  geom_text(aes(y = ypos, label = star, vjust = vj, group = comparison),
            position = position_dodge(width = 0.8), size = TXT) +
  facet_grid(~ panel, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = c("All donors" = "#9BB0BC", "Female donors only" = "#C47F7A"), name = NULL,
                    labels = c("All donors" = "All donors (Term 1 M + 1 F; BPD 2 F)",
                               "Female donors only" = "Female only (Term 1 F; BPD 2 F)")) +
  scale_y_continuous(expand = expansion(mult = 0.12)) +
  labs(x = NULL, y = expression(log[2]~"FC, BPD vs Term")) +
  theme_cns() + theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),
                      legend.position = "top", legend.direction = "horizontal")
save_cns(p, "S2b_key_genes_log2FC_all_vs_female", 120, 58)

###############################################################################
## S3  GSE235562（MSC 体外处理）：性别构成与 XIST
###############################################################################
s235 <- read.csv(file.path(od, "4e_GSE235562_sex_inference.csv"))
grp_lv <- intersect(c("control", "CMS", "Hox40", "Hox80", "CMS+Hox40", "CMS+Hox80"), unique(s235$group))
s235 <- s235 %>% filter(!is.na(group)) %>% mutate(group = factor(group, levels = grp_lv))
tab235 <- table(Treatment = s235$group, Sex = s235$sex_inferred)
p235 <- if (all(dim(tab235) >= 2)) fisher.test(tab235)$p.value else NA

p <- as.data.frame(tab235) %>% ggplot(aes(Treatment, Freq, fill = Sex)) +
  geom_col(width = 0.65, colour = "black", linewidth = LW) +
  annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.2, size = TXT,
           label = ifelse(is.na(p235), "", paste("Fisher", fmt_p(p235)))) +
  scale_fill_manual(values = sex_cols, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)), breaks = pretty_breaks()) +
  labs(x = NULL, y = "Samples") +
  theme_cns() + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top")
save_cns(p, "S3a_GSE235562_sex_by_treatment", 55, 50)

p <- ggplot(s235, aes(group, XIST)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, fill = NA, colour = "black", linewidth = LW) +
  geom_point(aes(fill = sex_inferred, shape = sex_inferred), size = PT, stroke = LW * 1.5,
             position = position_jitter(width = 0.15, seed = 1)) +
  scale_fill_manual(values = sex_cols, name = NULL) + scale_shape_manual(values = sex_shapes, name = NULL) +
  labs(x = NULL, y = expression(italic(XIST)~"expression")) +
  theme_cns() + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top")
save_cns(p, "S3b_GSE235562_XIST_by_treatment", 55, 50)

###############################################################################
## S4  GSE106910（气管吸出物 MSC）：严重度 × 性别
###############################################################################
s106 <- read.csv(file.path(od, "4d_GSE106910_sample_table.csv"))
subj <- read.csv(file.path(od, "4d2_GSE106910_subject_table.csv"))
s106b <- s106 %>% filter(!is.na(Group), Group != "N/A") %>%
  mutate(Group = factor(sev_lab[Group], levels = sev_lab))
subj_b <- subj %>% filter(!is.na(Group), Group != "N/A") %>%
  mutate(Group = factor(sev_lab[Group], levels = sev_lab),
         BPD = factor(ifelse(Group == "No BPD", "control", "BPD"), levels = c("control", "BPD")))

## ---- S4a：逐样本性别判定（同一婴儿的重复样本连线） --------------------------------
has_y <- "Yscore" %in% colnames(s106b) && any(!is.na(s106b$Yscore))
if (has_y) {
  p <- ggplot(s106b, aes(XIST, Yscore)) +
    geom_line(aes(group = subject), colour = "grey70", linewidth = LW) +
    geom_point(aes(fill = sex_sample, shape = sex_sample), size = PT, stroke = LW * 1.5) +
    labs(x = expression(italic(XIST)), y = "chrY genes (mean)")
} else {                                     # 平台无 Y 探针：按婴儿排列的 XIST 条带图
  s106b <- s106b %>% mutate(subject = reorder(subject, XIST, mean))
  p <- ggplot(s106b, aes(subject, XIST)) +
    geom_line(aes(group = subject), colour = "grey60", linewidth = LW) +
    geom_point(aes(fill = sex_sample, shape = sex_sample), size = PT, stroke = LW * 1.5) +
    facet_grid(~ Group, scales = "free_x", space = "free_x") +
    labs(x = "Infant (replicates joined)", y = expression(italic(XIST)~"expression")) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
}
p <- p + scale_fill_manual(values = sex_cols, name = "Inferred sex") +
  scale_shape_manual(values = sex_shapes, name = "Inferred sex") +
  theme_cns() + theme(legend.position = "top",
                      axis.text.x = element_text(angle = if (has_y) 0 else 90,
                                                 hjust = if (has_y) 0.5 else 1, vjust = 0.5))
save_cns(p, "S4a_GSE106910_sex_per_sample", 100, 52)

## ---- S4b：各严重度组性别构成（婴儿水平） -----------------------------------------
tab_sev <- table(Severity = subj_b$Group, Sex = subj_b$sex)
p_sev <- if (all(dim(tab_sev) >= 2)) fisher.test(tab_sev)$p.value else NA
p <- as.data.frame(tab_sev) %>% ggplot(aes(Severity, Freq, fill = Sex)) +
  geom_col(width = 0.65, colour = "black", linewidth = LW) +
  annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.2, size = TXT,
           label = ifelse(is.na(p_sev), "", paste("Fisher", fmt_p(p_sev)))) +
  scale_fill_manual(values = sex_cols, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)), breaks = pretty_breaks()) +
  labs(x = NULL, y = sprintf("Infants (n = %d)", nrow(subj_b))) +
  theme_cns() + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top")
save_cns(p, "S4b_GSE106910_sex_by_severity", 50, 50)

## ---- S4c / S4d：XIST、PHLDA1 按严重度（婴儿均值，点按性别） -------------------------
has_sex <- length(unique(na.omit(subj_b$sex))) >= 2
box_sev <- function(y, ylab, note = NULL) {
  ggplot(subj_b, aes(Group, .data[[y]])) +
    geom_boxplot(width = 0.6, outlier.shape = NA, fill = NA, colour = "black", linewidth = LW) +
    geom_point(aes(fill = sex, shape = sex), size = PT + 0.2, stroke = LW * 1.5,
               position = position_jitter(width = 0.12, seed = 1)) +
    { if (!is.null(note)) annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.2,
                                   size = TXT, label = note) } +
    scale_fill_manual(values = sex_cols, name = NULL) + scale_shape_manual(values = sex_shapes, name = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
    labs(x = NULL, y = ylab) +
    theme_cns() + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top")
}
note_x <- note_p <- NULL
if (has_sex) {
  cx <- summary(lm(XIST ~ as.numeric(Group) + sex, data = subj_b))$coefficients
  cp <- summary(lm(PHLDA1 ~ BPD + sex, data = subj_b))$coefficients
  note_x <- paste("Severity trend adj. for sex\n", fmt_p(cx["as.numeric(Group)", 4]))
  note_p <- paste("BPD effect adj. for sex\n", fmt_p(cp["BPDBPD", 4]))
}
save_cns(box_sev("XIST",   expression(italic(XIST)~"(infant mean)"),   note_x), "S4c_GSE106910_XIST_by_severity",   55, 55)
save_cns(box_sev("PHLDA1", expression(italic(PHLDA1)~"(infant mean)"), note_p), "S4d_GSE106910_PHLDA1_by_severity", 55, 55)

message("全部完成：", odf)
