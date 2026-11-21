###############################################################################
# B2  供者水平的细胞组成与可视化   —— 回应 R1-2        （R2：统一配色）
#   1) UMAP 按供者分面（整体间充质 + Adventitial FB）
#   2) 每个供者的 7 种间充质亚群比例、6 个 Adventitial FB 亚簇比例
#   3) propeller 比例检验（每组 n=2，仅作描述）
#   4) 每个供者的 PHLDA1 表达（单细胞小提琴 + 供者 pseudo-bulk）
#
# 配色（与正文 Fig. 2/3 一致，取自原图图例）：
#   分组  Term infant = 浅灰 #D7D7D7，BPD = 青蓝 #00768E，BPD+PH = 粉 #E59995
#   供者  同组两个供者用同一色系的深/浅两色
#   细胞类型、Adventitial FB 亚簇：与正文图例相同
###############################################################################
source(here::here("revision_scripts", "00_setup.R"))
od <- file.path(OUT, "B2_donor")
dir.create(od, recursive = TRUE, showWarnings = FALSE)

## ---- 0. 配色 -----------------------------------------------------------------
group_cols <- c("Term infant"          = "#D7D7D7",
                "BPD"                  = "#00768E",
                "BPD+PH"               = "#E59995",
                "Acute preterm injury" = "#B8A07E")

donor_cols <- c("Term infant 1"          = "#C4C4C4",   # 浅灰
                "Term infant 2"          = "#8F8F8F",   # 中灰（纯 #D7D7D7 在白底上太浅，作点/线时略加深）
                "BPD 1"                  = "#00768E",   # 青蓝
                "BPD 2"                  = "#5DB3C6",   # 浅青蓝
                "BPD+PH 1"               = "#E59995",   # 粉
                "BPD+PH 2"               = "#C25550",   # 深粉红
                "Acute preterm injury 1" = "#B8A07E")

celltype_cols <- c("Alveolar FB"    = "#2CA02C",
                   "Activated FB"   = "#98DF8A",
                   "Adventitial FB" = "#1F77B4",
                   "Alveolar MyoFB" = "#AEC7E8",
                   "Ductal MyoFB"   = "#FF7F0E",
                   "Pericyte"       = "#FFBB78",
                   "VSMC"           = "#D62728")

subcluster_cols <- c("0" = "#66C2A5", "1" = "#FC8D62", "2" = "#8DA0CB",
                     "3" = "#E78AC3", "4" = "#A6D854", "5" = "#FFD92F")

# 分面条带按分组上色（需要 ggh4x；没有则退回默认灰色条带）
group_strips <- function(groups) {
  if (!requireNamespace("ggh4x", quietly = TRUE)) return(NULL)
  g <- intersect(names(group_cols), groups)
  txt <- ifelse(g == "BPD", "white", "black")
  ggh4x::strip_themed(background_x = ggh4x::elem_list_rect(fill = group_cols[g], colour = NA),
                      text_x = ggh4x::elem_list_text(colour = txt, face = "bold"))
}
facet_by_group <- function(groups) {
  st <- group_strips(groups)
  if (is.null(st)) facet_grid(~ group, scales = "free_x", space = "free_x")
  else ggh4x::facet_grid2(~ group, scales = "free_x", space = "free_x", strip = st)
}

## 固定因子顺序
lvl_ct <- names(celltype_cols)
for (nm in c("meso_relabel", "meso_relabel_no_preterm")) {
  o <- get(nm); o$celltype <- factor(as.character(o$celltype), levels = lvl_ct); assign(nm, o)
}
meso_adventitial_FB$subcluster <- factor(as.character(meso_adventitial_FB$subcluster),
                                         levels = names(subcluster_cols))

## ---- 1. UMAP 按供者分面 ------------------------------------------------------
u1 <- DimPlot(meso_relabel_no_preterm, reduction = pick_umap(meso_relabel_no_preterm),
              group.by = "celltype", split.by = "donor", ncol = 6, pt.size = 0.2,
              cols = celltype_cols) +
  ggtitle("Mesenchymal cells, split by donor")
save_both(u1, file.path(od, "1a_UMAP_mesenchyme_by_donor"), 14, 3)

u1b <- DimPlot(meso_relabel_no_preterm, reduction = pick_umap(meso_relabel_no_preterm),
               group.by = "donor", cols = donor_cols, pt.size = 0.2, shuffle = TRUE, seed = 2026) +
  ggtitle("Mesenchymal cells coloured by donor")
save_both(u1b, file.path(od, "1b_UMAP_mesenchyme_donor_colour"), 7, 5.5)

u1c <- DimPlot(meso_relabel_no_preterm, reduction = pick_umap(meso_relabel_no_preterm),
               group.by = "group", cols = group_cols, pt.size = 0.2, shuffle = TRUE, seed = 2026) +
  ggtitle("Mesenchymal cells coloured by group")
save_both(u1c, file.path(od, "1b2_UMAP_mesenchyme_group_colour"), 7, 5.5)

u2 <- DimPlot(meso_adventitial_FB, reduction = pick_umap(meso_adventitial_FB),
              group.by = "subcluster", split.by = "donor", ncol = 6, label = TRUE, pt.size = 0.4,
              cols = subcluster_cols) +
  ggtitle("Adventitial FB subclusters, split by donor")
save_both(u2, file.path(od, "1c_UMAP_advFB_subclusters_by_donor"), 14, 3)

u2b <- DimPlot(meso_adventitial_FB, reduction = pick_umap(meso_adventitial_FB),
               group.by = "donor", cols = donor_cols, pt.size = 0.2, shuffle = TRUE, seed = 2026) +
  ggtitle("Adventitial FB coloured by donor")
save_both(u2b, file.path(od, "2b_UMAP_advFB_subclusters_donor_colour"), 7, 5.5)


## ---- 2. 每个供者的组成 -------------------------------------------------------
donor2group <- distinct(meso_relabel@meta.data[, c("donor", "group")]) %>%
  mutate(donor = as.character(donor), group = as.character(group))
comp_meso <- droplevels(meso_relabel_no_preterm@meta.data) %>%
  count(donor, celltype, .drop = FALSE) %>%
  group_by(donor) %>% mutate(pct = 100 * n / sum(n)) %>% ungroup() %>%
  mutate(donor = as.character(donor)) %>%
  left_join(donor2group, by = "donor") %>%
  mutate(donor = factor(donor, levels = DONOR_ORDER),
         group = factor(group, levels = names(group_cols)))
write.csv(comp_meso, file.path(od, "2a_mesenchymal_composition_per_donor.csv"), row.names = FALSE)

adv_tab <- comp_meso %>% filter(celltype == "Adventitial FB") %>%
  select(donor, group, n_adv = n, pct_adv = pct)
pooled <- meso_relabel_no_preterm@meta.data %>% group_by(group) %>%
  summarise(pooled_pct_adv = 100 * mean(celltype == "Adventitial FB"), .groups = "drop")
cat("\n== Adventitial FB 占间充质比例：逐供者 ==\n"); print(adv_tab)
cat("\n== 合并计算（稿件原数字 44.74/11.60/2.97 应与此一致）==\n"); print(pooled)
write.csv(adv_tab, file.path(od, "2b_adventitialFB_pct_per_donor.csv"), row.names = FALSE)

p2a <- ggplot(comp_meso, aes(donor, pct, fill = celltype)) + geom_col(width = 0.8) +
  facet_by_group(levels(droplevels(comp_meso$group))) +
  scale_fill_manual(values = celltype_cols) +
  labs(x = NULL, y = "% of mesenchymal cells", fill = NULL) + theme_classic(11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_both(p2a, file.path(od, "2c_mesenchymal_composition_stacked_by_donor"), 9, 5)

# 每点一个供者；横线 = 组均值；点用供者色，组用分组色的半透明底条
grp_mean <- comp_meso %>% group_by(celltype, group) %>% summarise(mean_pct = mean(pct), .groups = "drop")
p2b <- ggplot(comp_meso, aes(group, pct)) +
  geom_col(data = grp_mean, aes(y = mean_pct, fill = group), width = 0.6, alpha = 0.35) +
  geom_point(aes(colour = donor), size = 2.6, position = position_dodge(width = 0.4)) +
  facet_wrap(~ celltype, scales = "free_y", nrow = 2) +
  scale_fill_manual(values = group_cols, guide = "none") +
  scale_colour_manual(values = donor_cols, name = "Donor") +
  labs(x = NULL, y = "% of mesenchymal cells") + theme_bw(10) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), panel.grid.minor = element_blank())
save_both(p2b, file.path(od, "2d_mesenchymal_composition_dotplot_per_donor"), 11, 6)

# 只看 Adventitial FB（建议作为主图/补充图的单独面板）
p2b_adv <- comp_meso %>% filter(celltype == "Adventitial FB") %>%
  ggplot(aes(group, pct)) +
  geom_col(data = filter(grp_mean, celltype == "Adventitial FB"),
           aes(y = mean_pct, fill = group), width = 0.6, colour = "grey30", linewidth = 0.2) +
  geom_point(aes(colour = donor), size = 3, position = position_dodge(width = 0.35)) +
  scale_fill_manual(values = group_cols, guide = "none") +
  scale_colour_manual(values = donor_cols, name = "Donor") +
  labs(x = NULL, y = "Adventitial FB (% of mesenchymal cells)") + theme_classic(11)
save_both(p2b_adv, file.path(od, "2d2_AdvFB_pct_per_donor"), 5, 4)

comp_adv <- droplevels(meso_adventitial_FB@meta.data) %>%
  count(donor, subcluster, .drop = FALSE) %>%
  group_by(donor) %>% mutate(pct = 100 * n / sum(n)) %>% ungroup() %>%
  mutate(donor = as.character(donor)) %>%
  left_join(donor2group, by = "donor") %>%
  mutate(donor = factor(donor, levels = DONOR_ORDER),
         group = factor(group, levels = names(group_cols)))
write.csv(comp_adv, file.path(od, "2e_advFB_subcluster_composition_per_donor.csv"), row.names = FALSE)
p2c <- ggplot(comp_adv, aes(donor, pct, fill = subcluster)) + geom_col(width = 0.8) +
  facet_by_group(levels(droplevels(comp_adv$group))) +
  scale_fill_manual(values = subcluster_cols) +
  labs(x = NULL, y = "% of adventitial FB", fill = "Cluster") + theme_classic(11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_both(p2c, file.path(od, "2f_advFB_subcluster_stacked_by_donor"), 8, 5)

## ---- 3. propeller 比例检验（BPD vs Term，n = 2 vs 2） -------------------------
if (!requireNamespace("speckle", quietly = TRUE)) BiocManager::install("speckle")
library(speckle)
bt <- meso_relabel_no_preterm@meta.data %>% filter(group %in% c("Term infant", "BPD")) %>% droplevels()
prop_meso <- propeller(clusters = bt$celltype, sample = bt$donor, group = bt$group, transform = "asin")
write.csv(prop_meso, file.path(od, "3a_propeller_mesenchyme_BPD_vs_Term.csv"))
print(prop_meso)

bt2 <- meso_adventitial_FB@meta.data %>% filter(group %in% c("Term infant", "BPD")) %>% droplevels()
prop_adv <- propeller(clusters = bt2$subcluster, sample = bt2$donor, group = bt2$group, transform = "asin")
write.csv(prop_adv, file.path(od, "3b_propeller_advFB_subclusters_BPD_vs_Term.csv"))
print(prop_adv)
# 说明：每组仅 2 个供者，P 值功效极低；文中应以“两例 BPD 均高于两例 Term”这类逐例描述为主

# 逐例一致性：所有 BPD 供者是否都高于所有 Term 供者
consistent <- comp_meso %>% filter(group %in% c("Term infant", "BPD")) %>%
  group_by(celltype) %>%
  summarise(min_BPD = min(pct[group == "BPD"]), max_BPD = max(pct[group == "BPD"]),
            min_Term = min(pct[group == "Term infant"]), max_Term = max(pct[group == "Term infant"]),
            all_BPD_gt_all_Term = min_BPD > max_Term,
            all_BPD_lt_all_Term = max_BPD < min_Term, .groups = "drop")
write.csv(consistent, file.path(od, "3c_donor_level_consistency.csv"), row.names = FALSE)
print(consistent)

## ---- 4. PHLDA1 逐供者表达 ----------------------------------------------------
adv <- ensure_rna_norm(meso_adventitial_FB)
if (inherits(adv[["RNA"]], "Assay5") && length(Layers(adv[["RNA"]], search = "counts")) > 1)
  adv[["RNA"]] <- JoinLayers(adv[["RNA"]])                     # Seurat v5 多 layer 时先合并
ex  <- FetchData(adv, vars = c("PHLDA1", "donor", "group", "subcluster"))
ex$donor <- factor(as.character(ex$donor), levels = DONOR_ORDER)
ex$group <- factor(as.character(ex$group), levels = names(group_cols))
p4a <- ggplot(ex, aes(donor, PHLDA1, fill = donor)) +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.8, colour = "grey30", linewidth = 0.2) +
  stat_summary(fun = mean, geom = "point", size = 2.5, colour = "black") +
  facet_by_group(levels(droplevels(ex$group))) +
  scale_fill_manual(values = donor_cols, guide = "none") +
  labs(x = NULL, y = "PHLDA1 (log-normalized)", title = "PHLDA1 in adventitial FB, per donor") +
  theme_classic(11) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_both(p4a, file.path(od, "4a_PHLDA1_violin_per_donor"), 7, 4.5)

pb  <- pseudobulk(meso_adventitial_FB, "donor")
cpm <- t(t(pb) / colSums(pb)) * 1e6
phl <- data.frame(donor = colnames(cpm), PHLDA1_log2CPM = log2(cpm["PHLDA1", ] + 1),
                  pct_cells_expr = tapply(ex$PHLDA1 > 0, as.character(ex$donor), mean)[colnames(cpm)] * 100,
                  n_cells = as.vector(table(as.character(meso_adventitial_FB$donor))[colnames(cpm)]))
phl$group <- donor2group$group[match(phl$donor, donor2group$donor)]
write.csv(phl, file.path(od, "4b_PHLDA1_pseudobulk_per_donor.csv"), row.names = FALSE)
print(phl)

p4b <- phl %>% mutate(group = factor(group, levels = names(group_cols)),
                      donor = factor(donor, levels = DONOR_ORDER)) %>%
  ggplot(aes(group, PHLDA1_log2CPM)) +
  stat_summary(aes(fill = group), fun = mean, geom = "col", width = 0.6, colour = "grey30", linewidth = 0.2) +
  geom_point(aes(colour = donor), size = 3, position = position_dodge(width = 0.35)) +
  scale_fill_manual(values = group_cols, guide = "none") +
  scale_colour_manual(values = donor_cols, name = "Donor") +
  labs(x = NULL, y = "PHLDA1, pseudo-bulk log2(CPM + 1)", title = "Adventitial FB, donor level") +
  theme_classic(11)
save_both(p4b, file.path(od, "4c_PHLDA1_pseudobulk_per_donor"), 5.5, 4)
