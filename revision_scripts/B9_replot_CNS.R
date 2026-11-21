###############################################################################
# B9_replot_CNS.R —— 按 CNS 期刊样式重新输出 B9 / B9b / B9c 的全部单图（只画图，不重算分析）
#
#   样式原则（Nature 图表规范）：
#     * 直接按印刷尺寸输出（单位 mm；Nature 单栏 89 mm、双栏 183 mm），在 PS 中 100% 放置，不要再缩放
#     * 字号 7 pt（坐标轴标题、刻度、图例统一），线宽 0.25 mm，无图标题（面板字母在 PS 中加）
#     * 颜色统一：Term = 灰 #BDBDBD，BPD = 青蓝 #00768E；供者同色系深浅
#     * 同时输出 PDF（矢量，文字可编辑）和 600 dpi PNG
#
#   输入：已保存的 csv 与 CellChat rds；C2（CD8 状态小提琴图）和 C（AdvFB 侧基因）需从对象重算，耗时几分钟
#   输出：12_Revision/B9_immune/7_CNS_figures/
###############################################################################
source(here::here("revision_scripts", "00_setup.R"))
suppressPackageStartupMessages({ library(CellChat); library(scales) })
od  <- file.path(OUT, "B9_immune")
od5 <- file.path(od, "5_focus_AdvFB_CD8")
od6 <- file.path(od, "6_immune_evidence")
odf <- file.path(od, "7_CNS_figures"); dir.create(odf, showWarnings = FALSE)

## ============================ 统一样式 ============================================
FS <- 7                                                     # 字号 pt
LW <- 0.25                                                  # 线宽 mm
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
theme_cns_grid <- function(base = FS) {                      # 热图 / 气泡图用（带边框、无轴线）
  theme_cns(base) + theme(axis.line = element_blank(),
                          panel.border = element_rect(fill = NA, colour = "black", linewidth = LW))
}
save_cns <- function(p, name, w_mm, h_mm) {
  dev_pdf <- if (capabilities("cairo")) cairo_pdf else "pdf"
  ggsave(file.path(odf, paste0(name, ".pdf")), p, width = w_mm, height = h_mm, units = "mm", device = dev_pdf)
  ggsave(file.path(odf, paste0(name, ".png")), p, width = w_mm, height = h_mm, units = "mm", dpi = 600, bg = "white")
  message(sprintf("%-45s %3.0f × %3.0f mm", name, w_mm, h_mm))
}
PT <- 1.6                                                    # 点大小（ggplot size）
cond_cols  <- c("Term" = "#BDBDBD", "BPD" = "#00768E")
group_cols <- c("Term infant" = "#BDBDBD", "BPD" = "#00768E", "BPD+PH" = "#E59995")
donor_cols <- c("Term infant 1" = "#9E9E9E", "Term infant 2" = "#616161",
                "BPD 1" = "#00768E", "BPD 2" = "#5DB3C6", "BPD+PH 1" = "#E59995", "BPD+PH 2" = "#C25550")
donor_lab  <- c("Term infant 1" = "Term 1", "Term infant 2" = "Term 2", "BPD 1" = "BPD 1", "BPD 2" = "BPD 2",
                "BPD+PH 1" = "BPD+PH 1", "BPD+PH 2" = "BPD+PH 2")
grp_lab    <- c("Term infant" = "Term", "BPD" = "BPD", "BPD+PH" = "BPD+PH")
fix_sub <- function(x) gsub("Alveolar macrophage", "Alv. macrophage", gsub(" Cell$", " cells", x))

## ================= C1a：免疫亚群逐供者比例（9 个小面板） ===========================
show_subsets <- c("CD8 T (all)", "CD8 naive", "CD4 Th1-like", "CD4 memory", "NK Cell",
                  "Alveolar macrophage", "Monocyte", "cDC", "B Cell")
imm_comp <- read.csv(file.path(od, "1a_immune_composition_per_donor.csv"))
t_comp   <- read.csv(file.path(od, "2a_T_subsets_per_donor.csv"))
comp_long <- bind_rows(
  imm_comp %>% transmute(donor, group, subset = celltype, pct_of_all),
  t_comp %>% filter(Tsub != "Unassigned") %>% transmute(donor, group, subset = Tsub, pct_of_all),
  t_comp %>% filter(grepl("^CD8", Tsub)) %>% group_by(donor, group) %>%
    summarise(pct_of_all = sum(pct_of_all), .groups = "drop") %>% mutate(subset = "CD8 T (all)"))
c1 <- comp_long %>% filter(subset %in% show_subsets, group %in% c("Term infant", "BPD")) %>%
  mutate(group = factor(grp_lab[group], levels = c("Term", "BPD")),
         subset = factor(fix_sub(subset), levels = fix_sub(show_subsets)),
         donor = factor(donor_lab[donor], levels = donor_lab))
p <- ggplot(c1, aes(group, pct_of_all)) +
  stat_summary(aes(fill = group), fun = mean, geom = "col", width = 0.65, colour = "black", linewidth = LW) +
  geom_point(aes(colour = donor), size = PT, position = position_jitter(width = 0.12, seed = 1)) +
  facet_wrap(~ subset, scales = "free_y", nrow = 2) +
  scale_fill_manual(values = cond_cols, guide = "none") +
  scale_colour_manual(values = setNames(donor_cols, donor_lab[names(donor_cols)]), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(x = NULL, y = "% of all lung cells") + theme_cns() + theme(legend.position = "right")
save_cns(p, "C1a_immune_subsets_per_donor", 150, 62)

## ================= C1b：log2FC 棒棒糖图 ===========================================
fc <- read.csv(file.path(od6, "C1_immune_subsets_BPD_vs_Term_scRNA.csv")) %>%
  mutate(subset = factor(fix_sub(subset), levels = rev(fix_sub(show_subsets))))
p <- ggplot(fc, aes(log2FC, subset)) +
  geom_vline(xintercept = 0, linewidth = LW, colour = "grey40") +
  geom_segment(aes(x = 0, xend = log2FC, yend = subset), linewidth = 0.4, colour = "grey60") +
  geom_point(aes(fill = donor_consistent), shape = 21, size = 2.2, stroke = LW * 1.5) +
  scale_fill_manual(values = c(`TRUE` = "#00768E", `FALSE` = "white"),
                    labels = c(`TRUE` = "Consistent across donors", `FALSE` = "Overlapping donors"), name = NULL) +
  labs(x = expression(log[2]~"(BPD / Term)"), y = NULL) +
  theme_cns() + theme(legend.position = "bottom", legend.direction = "vertical")
save_cns(p, "C1b_immune_subsets_log2FC", 60, 62)

## ================= C2：CD8+ T 状态评分（需重算单细胞评分） =========================
imm <- add_donor(readRDS(PATH_IMM))
tc  <- subset(imm, subset = celltype %in% c("T Cell", "NKT Cell") & condition %in% c("Term infant", "BPD"))
rm(imm); gc()
DefaultAssay(tc) <- "RNA"
if (inherits(tc[["RNA"]], "Assay5")) tc[["RNA"]] <- JoinLayers(tc[["RNA"]])
tc  <- NormalizeData(tc, assay = "RNA", verbose = FALSE)
cnt <- GetAssayData(tc, assay = "RNA", layer = "counts")
pos <- function(x) if (x %in% rownames(cnt)) cnt[x, ] > 0 else rep(FALSE, ncol(cnt))
tc$Tclass <- ifelse((pos("CD8A") | pos("CD8B")) & !pos("CD4"), "CD8 T",
                    ifelse(pos("CD4") & !(pos("CD8A") | pos("CD8B")), "CD4 T", "Other"))
cd8 <- subset(tc, subset = Tclass == "CD8 T")
modules <- list(Naive        = c("CCR7", "SELL", "LEF1", "TCF7", "IL7R"),
                Cytotoxic    = c("GZMB", "GZMA", "GZMK", "PRF1", "NKG7", "GNLY", "CCL5"),
                `IFN response` = c("ISG15", "IFI6", "IFIT1", "IFIT3", "MX1", "STAT1"),
                Exhaustion   = c("PDCD1", "LAG3", "HAVCR2", "TIGIT", "TOX", "CTLA4"),
                Activation   = c("CD69", "FOS", "JUN", "NR4A1", "HLA-DRA", "CD38"))
modules <- lapply(modules, intersect, rownames(cd8)); modules <- modules[sapply(modules, length) >= 2]
cd8 <- AddModuleScore(cd8, features = modules, name = "MOD_", assay = "RNA",
                      ctrl = min(100, max(10, floor(ncol(cd8) / 10))), seed = 2026)
sc <- cd8@meta.data[, paste0("MOD_", seq_along(modules))]; colnames(sc) <- names(modules)
sc <- cbind(cd8@meta.data[, c("donor", "group")], sc) %>%
  mutate(donor = factor(donor_lab[as.character(donor)], levels = donor_lab),
         group = factor(grp_lab[as.character(group)], levels = c("Term", "BPD"))) %>%
  pivot_longer(all_of(names(modules)), names_to = "module", values_to = "score") %>%
  mutate(module = factor(module, levels = names(modules)))
sc_donor <- sc %>% group_by(module, donor, group) %>% summarise(score = mean(score), .groups = "drop")
p <- ggplot(sc, aes(group, score)) +
  geom_violin(aes(fill = group), scale = "width", width = 0.8, colour = "black", linewidth = LW, trim = TRUE) +
  geom_point(data = sc_donor, aes(colour = donor), size = PT, position = position_jitter(width = 0.1, seed = 1)) +
  facet_wrap(~ module, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = cond_cols, guide = "none") +
  scale_colour_manual(values = setNames(c("#000000", "#555555", "#003F4D", "#1F6F80"), c("Term 1", "Term 2", "BPD 1", "BPD 2")),
                      name = "Donor mean") +
  labs(x = NULL, y = "Module score") + theme_cns()
save_cns(p, "C2_CD8_state_module_scores", 150, 48)

## ================= C3：供者水平热图 ===============================================
dm <- read.csv(file.path(od6, "C3_donor_level_matrix.csv")) %>%
  mutate(donor = factor(donor_lab[donor], levels = donor_lab),
         subset = fix_sub(subset),
         subset = factor(subset, levels = rev(c("Adventitial FB (% mesenchyme)", "PHLDA1 in adventitial FB",
                                                fix_sub(show_subsets)))))
p <- ggplot(dm, aes(donor, subset, fill = z)) + geom_tile(colour = "white", linewidth = 0.3) +
  geom_text(aes(label = ifelse(abs(value) >= 10, sprintf("%.0f", value), sprintf("%.2g", value))), size = FS / .pt * 0.9) +
  scale_fill_gradient2(low = "#3B6E8F", mid = "white", high = "#B2182B", name = "Row\nz-score") +
  geom_vline(xintercept = c(2.5, 4.5), linewidth = LW, colour = "black") +
  scale_x_discrete(position = "top", expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
  labs(x = NULL, y = NULL) + theme_cns_grid() +
  theme(axis.text.x = element_text(angle = 45, hjust = 0), axis.ticks = element_blank())
save_cns(p, "C3_donor_level_heatmap", 95, 70)

## ================= C4：相关方向对照 ===============================================
ct <- read.csv(file.path(od6, "C4_PHLDA1_immune_correlation_scRNA_vs_deconv.csv")) %>%
  filter(!is.na(corr_with_PHLDA1)) %>%
  mutate(deconv_sign = ifelse(corr_with_PHLDA1 > 0, "Positive", "Negative"),
         agree = ifelse(same_direction, "Same direction", "Opposite direction"))
p <- ggplot(ct, aes(rho_scRNA, reorder(subset, rho_scRNA))) +
  geom_vline(xintercept = 0, linewidth = LW, colour = "grey40") +
  geom_point(aes(fill = deconv_sign, shape = agree), size = 2.4, stroke = LW * 1.5) +
  scale_fill_manual(values = c(Positive = "#B2182B", Negative = "#3B6E8F"), name = "Deconvolution") +
  scale_shape_manual(values = c("Same direction" = 21, "Opposite direction" = 24), name = NULL) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  scale_x_continuous(limits = c(-1, 1), breaks = c(-1, -0.5, 0, 0.5, 1)) +
  labs(x = expression("Spearman"~rho~"(donor level, n = 6)"), y = NULL) +
  theme_cns() + theme(legend.position = "right")
save_cns(p, "C4_PHLDA1_immune_direction", 85, 35)

## ================= CellChat：读取对象 ==============================================
cc_term <- readRDS(file.path(od, "cellchat_advFB_immune_Term_infant.rds"))
cc_bpd  <- readRDS(file.path(od, "cellchat_advFB_immune_BPD.rds"))
grp <- intersect(levels(cc_term@idents), levels(cc_bpd@idents))

## ================= 4a：总通讯数量 / 强度 ===========================================
tot <- data.frame(cond = factor(c("Term", "BPD"), levels = c("Term", "BPD")),
                  count  = c(sum(cc_term@net$count[grp, grp]),  sum(cc_bpd@net$count[grp, grp])),
                  weight = c(sum(cc_term@net$weight[grp, grp]), sum(cc_bpd@net$weight[grp, grp])))
bar1 <- function(y, ylab, digits) ggplot(tot, aes(cond, .data[[y]], fill = cond)) +
  geom_col(width = 0.65, colour = "black", linewidth = LW) +
  geom_text(aes(label = formatC(.data[[y]], format = "f", digits = digits, big.mark = ",")), vjust = -0.4, size = FS / .pt) +
  scale_fill_manual(values = cond_cols, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)), labels = label_comma()) +
  labs(x = NULL, y = ylab) + theme_cns()
save_cns(bar1("count",  "Number of inferred interactions", 0), "4a_cellchat_total_count",  35, 45)
save_cns(bar1("weight", "Interaction strength",            1), "4a_cellchat_total_weight", 35, 45)

## ================= A1：差异最大的通路（哑铃图） =====================================
top_pw <- read.csv(file.path(od5, "A1_top_pathways.csv")) %>%
  mutate(direction = factor(sub("AdvFB", "Adv. FB", direction), levels = c("Adv. FB -> T", "T -> Adv. FB")))
dd <- top_pw %>% pivot_longer(c(Term, BPD), names_to = "cond", values_to = "prob") %>%
  mutate(cond = factor(cond, levels = c("Term", "BPD")),
         pathway = reorder(paste(pathway, direction, sep = "___"), delta))
p <- ggplot(dd, aes(prob, pathway)) +
  geom_line(aes(group = pathway), colour = "grey70", linewidth = 0.4) +
  geom_point(aes(colour = cond), size = PT + 0.3) +
  facet_wrap(~ direction, scales = "free") +
  scale_y_discrete(labels = function(x) sub("___.*$", "", x)) +
  scale_colour_manual(values = cond_cols, name = NULL) +
  labs(x = "Communication probability", y = NULL) +
  theme_cns() + theme(legend.position = "top", panel.grid.major.y = element_line(colour = "grey92", linewidth = LW))
save_cns(p, "A1_top_pathways_dumbbell", 110, 62)

## ================= A1b：各 AdvFB 亚群的贡献（热图） ================================
sum_prob <- function(P, src, tgt) {
  s <- intersect(src, dimnames(P)[[1]]); t <- intersect(tgt, dimnames(P)[[2]])
  if (!length(s) || !length(t)) return(setNames(numeric(0), character(0)))
  apply(P[s, t, , drop = FALSE], 3, sum)
}
ADV <- grep("^AdvFB_C", grp, value = TRUE); TG <- intersect("T Cell", grp)
sub_delta <- bind_rows(lapply(c("AdvFB -> T", "T -> AdvFB"), function(dn) bind_rows(lapply(ADV, function(a) {
  src <- if (dn == "AdvFB -> T") a else TG; tgt <- if (dn == "AdvFB -> T") TG else a
  x <- sum_prob(cc_term@netP$prob, src, tgt); y <- sum_prob(cc_bpd@netP$prob, src, tgt)
  pw <- top_pw$pathway[as.character(top_pw$direction) == sub("AdvFB", "Adv. FB", dn)]
  data.frame(direction = sub("AdvFB", "Adv. FB", dn), subcluster = sub("AdvFB_", "", a), pathway = pw,
             delta = unname(ifelse(is.na(y[pw]), 0, y[pw]) - ifelse(is.na(x[pw]), 0, x[pw])))
}))))
lim <- max(abs(sub_delta$delta))
p <- ggplot(sub_delta, aes(subcluster, pathway, fill = delta)) + geom_tile(colour = "white", linewidth = 0.3) +
  facet_wrap(~ direction, scales = "free_y") +
  scale_fill_gradient2(low = "#7F7F7F", mid = "white", high = "#00768E", limits = c(-lim, lim),
                       name = "BPD - Term") +
  scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
  labs(x = "Adventitial FB subcluster", y = NULL) + theme_cns_grid() + theme(axis.ticks = element_blank())
save_cns(p, "A1b_top_pathways_by_subcluster", 100, 60)

## ================= A2：差异最大的 L–R 对（气泡图） ==================================
top_lr <- read.csv(file.path(od5, "A2_top_LR.csv")) %>%
  mutate(direction = factor(sub("AdvFB", "Adv. FB", direction), levels = c("Adv. FB -> T", "T -> Adv. FB")))
dl <- top_lr %>% pivot_longer(c(Term, BPD), names_to = "cond", values_to = "prob") %>%
  mutate(cond = factor(cond, levels = c("Term", "BPD")),
         label = reorder(paste0(interaction_name_2, "___", direction), delta))
p <- ggplot(dl, aes(cond, label)) +
  geom_point(aes(size = prob, colour = prob)) +
  facet_wrap(~ direction, scales = "free_y") +
  scale_y_discrete(labels = function(x) sub("___.*$", "", x)) +
  scale_colour_gradient(low = "grey85", high = "#00768E", name = "Prob.") +
  scale_size_continuous(range = c(0.3, 3.5), name = "Prob.") +
  guides(colour = guide_legend(), size = guide_legend()) +
  labs(x = NULL, y = NULL) + theme_cns_grid() +
  theme(panel.grid.major = element_line(colour = "grey92", linewidth = LW))
save_cns(p, "A2_top_LR_bubble", 120, 70)

## ================= B：T 细胞侧基因 CD8 vs CD4（点图） ===============================
dot <- read.csv(file.path(od5, "B_Tside_genes_CD8_vs_CD4.csv")) %>%
  mutate(column = factor(gsub("\n", " ", column),
                         levels = c("CD8 T Term", "CD8 T BPD", "CD4 T Term", "CD4 T BPD")),
         direction = sub("AdvFB", "Adv. FB", direction))
ng <- length(unique(dot$gene))
p <- ggplot(dot, aes(column, gene)) +
  geom_point(aes(size = pct, colour = avg_scaled)) +
  facet_grid(direction ~ ., scales = "free_y", space = "free_y") +
  scale_colour_gradient2(low = "#3B6E8F", mid = "grey90", high = "#B2182B", name = "Scaled\nexpr.") +
  scale_size_continuous(range = c(0.2, 3.5), name = "% expr.") +
  labs(x = NULL, y = NULL) + theme_cns_grid() +
  theme(axis.text.y = element_text(face = "italic"), axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text.y = element_text(angle = 0))
save_cns(p, "B_Tside_genes_CD8_vs_CD4", 85, 20 + 3.3 * ng)

## ================= C：AdvFB 侧基因（点图，需重算） ==================================
adv <- ensure_rna_norm(meso_adventitial_FB)
if (inherits(adv[["RNA"]], "Assay5")) adv[["RNA"]] <- JoinLayers(adv[["RNA"]])
split_genes <- function(x) unique(unlist(strsplit(x, "_")))
fb_genes <- unique(unlist(lapply(seq_len(nrow(top_lr)), function(i)
  split_genes(if (grepl("^Adv", top_lr$direction[i])) top_lr$ligand[i] else top_lr$receptor[i]))))
fb_genes <- intersect(fb_genes, rownames(adv))
mdA <- adv@meta.data %>% filter(condition %in% c("Term infant", "BPD"), subcluster != "5")
dotA <- cbind(cond = grp_lab[as.character(mdA$condition)], FetchData(adv, vars = fb_genes, cells = rownames(mdA))) %>%
  pivot_longer(-cond, names_to = "gene", values_to = "expr") %>%
  group_by(gene, cond) %>% summarise(pct = 100 * mean(expr > 0), avg = mean(expm1(expr)), .groups = "drop") %>%
  group_by(gene) %>% mutate(avg_scaled = as.numeric(scale(log1p(avg)))) %>% ungroup() %>%
  mutate(cond = factor(cond, levels = c("Term", "BPD")))
p <- ggplot(dotA, aes(cond, gene)) + geom_point(aes(size = pct, colour = avg_scaled)) +
  scale_colour_gradient2(low = "#3B6E8F", mid = "grey90", high = "#B2182B", name = "Scaled\nexpr.") +
  scale_size_continuous(range = c(0.2, 3.5), name = "% expr.") +
  labs(x = NULL, y = NULL) + theme_cns_grid() + theme(axis.text.y = element_text(face = "italic"))
save_cns(p, "C_AdvFB_side_genes", 50, 18 + 3.3 * length(fb_genes))

message("全部完成：", odf)
