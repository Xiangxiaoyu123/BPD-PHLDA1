###############################################################################
# B9c  单细胞层面的免疫证据补充图（回应 R2-3）—— 不重跑 CellChat，只用已有结果 + 表达矩阵
#
#  R2-3 原意：稿件声称“免疫–间充质互作微环境”，但只有趋势，没有空间分析；
#             要么改写表述，要么做组织染色（如 CD8+ T）。
#  单细胞数据不能证明空间毗邻，但可以把“反卷积推断”升级为“组织单细胞实测”，并为 IF 实验提供依据：
#
#   C1  免疫亚群比例在 BPD vs Term 中的变化（逐供者），与 xCell/CIBERSORT 的方向逐项对照
#   C2  CD8+ T 细胞自身的状态改变（naive / 细胞毒 / 干扰素应答 / 耗竭 模块评分），Term vs BPD，逐供者
#   C3  供者水平整合热图：AdvFB 比例、AdvFB 中 PHLDA1、各免疫亚群比例（一眼看出共变）
#   C4  供者水平：AdvFB-PHLDA1 与各免疫亚群比例的相关方向 vs 反卷积中报告的相关方向
#
#  输入：B9 第 1–2 部分已输出的 1a_immune_composition_per_donor.csv、2a_T_subsets_per_donor.csv；
#        B2 输出的 2b_adventitialFB_pct_per_donor.csv；PATH_IMM（免疫细胞对象）；meso_adventitial_FB
###############################################################################
source(here::here("revision_scripts", "00_setup.R"))
od  <- file.path(OUT, "B9_immune")
od6 <- file.path(od, "6_immune_evidence"); dir.create(od6, showWarnings = FALSE)

group_cols <- c("Term infant" = "#D7D7D7", "BPD" = "#00768E", "BPD+PH" = "#E59995")
donor_cols <- c("Term infant 1" = "#C4C4C4", "Term infant 2" = "#8F8F8F",
                "BPD 1" = "#00768E", "BPD 2" = "#5DB3C6", "BPD+PH 1" = "#E59995", "BPD+PH 2" = "#C25550")
BT_DONORS <- c("Term infant 1", "Term infant 2", "BPD 1", "BPD 2")

## ---- 反卷积结果（来自稿件 Results / Supplementary Table 12，请按实际结果核对或补充） -----
##  corr_with_PHLDA1：反卷积中该细胞与 PHLDA1 的相关方向（+1 / -1）
##  BPD_vs_control  ：反卷积中 BPD 相对对照的变化方向（+1 / -1；未知填 NA）
deconv <- tibble::tribble(
  ~deconv_celltype,        ~scRNA_subset,          ~corr_with_PHLDA1, ~BPD_vs_control, ~source,
  "CD8+ naive T cells",    "CD8 naive",             1,                NA,              "xCell (r = 0.414)",
  "Th1 cells",             "CD4 Th1-like",          1,                NA,              "xCell (r = 0.430)",
  "CD4+ memory T cells",   "CD4 memory",           -1,                NA,              "xCell (r = -0.674)",
  "NK cells activated",    "NK Cell",              NA,                NA,              "CIBERSORT",
  "Dendritic cells activated", "cDC",              NA,                NA,              "CIBERSORT",
  "Macrophages",           "Alveolar macrophage",  NA,                NA,              "xCell/CIBERSORT"
)

## ---- 读入已有的比例表 ---------------------------------------------------------------
imm_comp <- read.csv(file.path(od, "1a_immune_composition_per_donor.csv"))
t_comp   <- read.csv(file.path(od, "2a_T_subsets_per_donor.csv"))
adv_pct  <- read.csv(file.path(OUT, "B2_donor", "2b_adventitialFB_pct_per_donor.csv"))

comp_long <- bind_rows(
  imm_comp %>% transmute(donor, group, subset = celltype, pct_of_all),
  t_comp   %>% filter(!Tsub %in% c("Unassigned")) %>% transmute(donor, group, subset = Tsub, pct_of_all),
  t_comp   %>% filter(grepl("^CD8", Tsub)) %>% group_by(donor, group) %>%
    summarise(pct_of_all = sum(pct_of_all), .groups = "drop") %>% mutate(subset = "CD8 T (all)")) %>%
  filter(donor %in% names(donor_cols))

## ---- C1. BPD vs Term 的变化（逐供者）+ 与反卷积方向对照 -----------------------------------
show_subsets <- c("CD8 T (all)", "CD8 naive", "CD4 Th1-like", "CD4 memory", "NK Cell",
                  "Alveolar macrophage", "Monocyte", "cDC", "B Cell")
c1 <- comp_long %>% filter(subset %in% show_subsets, donor %in% BT_DONORS)
fc <- c1 %>% group_by(subset) %>%
  summarise(mean_Term = mean(pct_of_all[group == "Term infant"]), mean_BPD = mean(pct_of_all[group == "BPD"]),
            min_BPD = min(pct_of_all[group == "BPD"]), max_BPD = max(pct_of_all[group == "BPD"]),
            min_Term = min(pct_of_all[group == "Term infant"]), max_Term = max(pct_of_all[group == "Term infant"]),
            .groups = "drop") %>%
  mutate(log2FC = log2((mean_BPD + 0.01) / (mean_Term + 0.01)),
         donor_consistent = (min_BPD > max_Term) | (max_BPD < min_Term),
         sc_direction = sign(log2FC)) %>%
  left_join(deconv %>% select(subset = scRNA_subset, BPD_vs_control, corr_with_PHLDA1, source), by = "subset") %>%
  mutate(agree_with_deconv = ifelse(is.na(BPD_vs_control), NA, sc_direction == BPD_vs_control),
         subset = factor(subset, levels = rev(show_subsets)))
write.csv(fc, file.path(od6, "C1_immune_subsets_BPD_vs_Term_scRNA.csv"), row.names = FALSE)
print(fc)

p_c1a <- c1 %>% mutate(group = factor(group, levels = c("Term infant", "BPD")),
                       subset = factor(subset, levels = show_subsets)) %>%
  ggplot(aes(group, pct_of_all)) +
  stat_summary(aes(fill = group), fun = mean, geom = "col", width = 0.6, colour = "grey30", linewidth = 0.2) +
  geom_point(aes(colour = donor), size = 2.4, position = position_jitter(width = 0.08, seed = 1)) +
  facet_wrap(~ subset, scales = "free_y", nrow = 2) +
  scale_fill_manual(values = group_cols, guide = "none") + scale_colour_manual(values = donor_cols, name = "Donor") +
  labs(x = NULL, y = "% of all lung cells", title = "Immune subsets measured directly in lung scRNA-seq (GSE275938)") +
  theme_classic(9) + theme(strip.background = element_blank(), strip.text = element_text(face = "bold"))
save_both(p_c1a, file.path(od6, "C1a_immune_subsets_per_donor"), 11, 5)

p_c1b <- ggplot(fc, aes(log2FC, subset)) +
  geom_vline(xintercept = 0, colour = "grey50") +
  geom_segment(aes(x = 0, xend = log2FC, yend = subset), colour = "grey70") +
  geom_point(aes(fill = donor_consistent), shape = 21, size = 3.5) +
  scale_fill_manual(values = c(`TRUE` = "#00768E", `FALSE` = "white"),
                    labels = c(`TRUE` = "Both BPD donors beyond both Term donors", `FALSE` = "Overlapping donors"),
                    name = NULL) +
  labs(x = "log2 (mean BPD / mean Term), % of all cells", y = NULL,
       title = "Change of immune subsets in BPD (scRNA-seq)") +
  theme_bw(10) + theme(legend.position = "bottom", panel.grid.minor = element_blank())
save_both(p_c1b, file.path(od6, "C1b_immune_subsets_log2FC_lollipop"), 6.5, 4.5)

## ---- C2. CD8+ T 细胞状态（模块评分） --------------------------------------------------
imm <- add_donor(readRDS(PATH_IMM))
tc  <- subset(imm, subset = celltype %in% c("T Cell", "NKT Cell") & condition %in% c("Term infant", "BPD"))
rm(imm); gc()
if (inherits(tc[["RNA"]], "Assay5")) tc[["RNA"]] <- JoinLayers(tc[["RNA"]])
tc  <- NormalizeData(tc, verbose = FALSE)
cnt <- GetAssayData(tc, assay = "RNA", layer = "counts")
pos <- function(x) if (x %in% rownames(cnt)) cnt[x, ] > 0 else rep(FALSE, ncol(cnt))
cd8 <- subset(tc, cells = colnames(tc)[(pos("CD8A") | pos("CD8B")) & !pos("CD4")])
cat("\nCD8+ T 细胞数（供者）：\n"); print(table(cd8$donor))

modules <- list(Naive        = c("CCR7", "SELL", "LEF1", "TCF7", "IL7R"),
                Cytotoxic    = c("GZMB", "GZMA", "GZMK", "PRF1", "NKG7", "GNLY", "CCL5"),
                IFN_response = c("ISG15", "IFI6", "IFIT1", "IFIT3", "MX1", "STAT1"),
                Exhaustion   = c("PDCD1", "LAG3", "HAVCR2", "TIGIT", "TOX", "CTLA4"),
                Activation   = c("CD69", "FOS", "JUN", "NR4A1", "HLA-DRA", "CD38"))
modules <- lapply(modules, intersect, rownames(cd8))
cd8 <- AddModuleScore(cd8, features = modules, name = "MOD_", assay = "RNA", seed = 2026)
sc  <- cd8@meta.data[, paste0("MOD_", seq_along(modules))]; colnames(sc) <- names(modules)
sc  <- cbind(cd8@meta.data[, c("donor", "group")], sc) %>%
  mutate(donor = as.character(donor), group = factor(as.character(group), levels = c("Term infant", "BPD"))) %>%
  pivot_longer(all_of(names(modules)), names_to = "module", values_to = "score")
sc_donor <- sc %>% group_by(module, donor, group) %>% summarise(score = mean(score), n = n(), .groups = "drop")
sc_test  <- sc_donor %>% group_by(module) %>%
  summarise(mean_Term = mean(score[group == "Term infant"]), mean_BPD = mean(score[group == "BPD"]),
            donor_consistent = min(score[group == "BPD"]) > max(score[group == "Term infant"]) |
                               max(score[group == "BPD"]) < min(score[group == "Term infant"]), .groups = "drop")
write.csv(sc_donor, file.path(od6, "C2_CD8_state_scores_per_donor.csv"), row.names = FALSE)
write.csv(sc_test,  file.path(od6, "C2_CD8_state_scores_summary.csv"),   row.names = FALSE)
print(sc_test)

p_c2 <- ggplot(sc, aes(group, score)) +
  geom_violin(aes(fill = group), scale = "width", colour = "grey30", linewidth = 0.2, alpha = 0.8) +
  geom_point(data = sc_donor, aes(colour = donor), size = 2.6, position = position_jitter(width = 0.08, seed = 1)) +
  facet_wrap(~ module, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = group_cols, guide = "none") + scale_colour_manual(values = donor_cols, name = "Donor mean") +
  labs(x = NULL, y = "Module score", title = "Transcriptional state of CD8⁺ T cells (cells: violins; donors: points)") +
  theme_classic(9) + theme(strip.background = element_blank(), strip.text = element_text(face = "bold"))
save_both(p_c2, file.path(od6, "C2_CD8_state_module_scores"), 11, 3.8)

## ---- C3. 供者水平整合热图 ------------------------------------------------------------
pbA <- pseudobulk(meso_adventitial_FB, "donor")
phl <- log2(pbA["PHLDA1", ] / colSums(pbA) * 1e6 + 1)
donor_mat <- comp_long %>% filter(subset %in% show_subsets) %>%
  select(donor, subset, value = pct_of_all) %>%
  bind_rows(adv_pct %>% transmute(donor, subset = "Adventitial FB (% mesenchyme)", value = pct_adv),
            tibble(donor = names(phl), subset = "PHLDA1 in adventitial FB", value = unname(phl))) %>%
  filter(donor %in% names(donor_cols)) %>%
  group_by(subset) %>% mutate(z = as.numeric(scale(value))) %>% ungroup() %>%
  mutate(donor = factor(donor, levels = names(donor_cols)),
         subset = factor(subset, levels = rev(c("Adventitial FB (% mesenchyme)", "PHLDA1 in adventitial FB", show_subsets))))
write.csv(donor_mat, file.path(od6, "C3_donor_level_matrix.csv"), row.names = FALSE)
p_c3 <- ggplot(donor_mat, aes(donor, subset, fill = z)) + geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.2g", value)), size = 2.6) +
  scale_fill_gradient2(low = "#3B6E8F", mid = "white", high = "#B2182B", name = "Row z-score") +
  geom_vline(xintercept = c(2.5, 4.5), colour = "grey30") +
  labs(x = NULL, y = NULL, title = "Donor-level overview (numbers = raw values)") +
  theme_minimal(9) + theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())
save_both(p_c3, file.path(od6, "C3_donor_level_heatmap"), 7.5, 5)

## ---- C4. 供者水平相关方向 vs 反卷积 ----------------------------------------------------
cor_tab <- comp_long %>% filter(subset %in% deconv$scRNA_subset, donor %in% c(BT_DONORS, "BPD+PH 1", "BPD+PH 2")) %>%
  mutate(PHLDA1 = phl[donor]) %>% group_by(subset) %>%
  summarise(n_donors = n(), rho_scRNA = suppressWarnings(cor(PHLDA1, pct_of_all, method = "spearman")),
            .groups = "drop") %>%
  left_join(deconv %>% select(subset = scRNA_subset, deconv_celltype, corr_with_PHLDA1, source), by = "subset") %>%
  mutate(same_direction = ifelse(is.na(corr_with_PHLDA1), NA, sign(rho_scRNA) == corr_with_PHLDA1))
write.csv(cor_tab, file.path(od6, "C4_PHLDA1_immune_correlation_scRNA_vs_deconv.csv"), row.names = FALSE)
print(cor_tab)
p_c4 <- cor_tab %>% filter(!is.na(corr_with_PHLDA1)) %>%
  mutate(deconv_sign = ifelse(corr_with_PHLDA1 > 0, "Positive in deconvolution", "Negative in deconvolution")) %>%
  ggplot(aes(rho_scRNA, reorder(subset, rho_scRNA))) +
  geom_vline(xintercept = 0, colour = "grey50") +
  geom_point(aes(fill = deconv_sign, shape = same_direction), size = 4) +
  scale_fill_manual(values = c("Positive in deconvolution" = "#B2182B", "Negative in deconvolution" = "#3B6E8F"), name = NULL) +
  scale_shape_manual(values = c(`TRUE` = 21, `FALSE` = 24), labels = c(`TRUE` = "Same direction", `FALSE` = "Opposite"),
                     name = NULL) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  xlim(-1, 1) +
  labs(x = "Spearman ρ (donor level, n = 6): PHLDA1 in adventitial FB vs subset %", y = NULL,
       title = "Direction of association: scRNA-seq vs deconvolution") +
  theme_bw(10) + theme(legend.position = "bottom")
save_both(p_c4, file.path(od6, "C4_PHLDA1_immune_direction_scRNA_vs_deconv"), 7, 3.8)
message("完成，结果在：", od6)
