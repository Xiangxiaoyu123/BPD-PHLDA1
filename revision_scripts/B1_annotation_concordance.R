###############################################################################
# B1  与原文（Shirazi et al., Nat Commun 2025）细胞注释的对应分析   —— 回应 R1-1
#
# 背景（读代码后确认）：
#   本文的整体/间充质注释是直接运行原作者公开代码（SucreLab/HumanNeonatal）得到的，
#   并沿用了原代码中 "cluster 编号 -> 细胞类型" 的映射（如 labels[idxes %in% c(3,12)] <- "Adventitial FB"）。
#   由于软件版本/随机性不同，重跑后的 cluster 编号不一定与原作者一致，
#   所以必须证明：①每个注释仍然由正确的 marker 支持；②与原作者发布的注释一致（若能拿到）。
#
# 原文补充数据（Supplementary_tables.zip 解压后放到 ORIG_DIR）：
#   MOESM3 = Supplementary Data 1：BPD vs Term infant 的各细胞类型/谱系 DEG（MAST）
#   MOESM4 = Supplementary Data 2：BPD+PH vs Term infant
#   MOESM5 = Supplementary Data 3：BPD+PH vs BPD
#   MOESM6 = Supplementary Data 4：33 种细胞类型的 hallmark 基因（注释所用 marker）
#   MOESM7 = Supplementary Data 5：每种细胞类型 vs 其余所有细胞的 marker DEG
#
# A 部分：用原文 hallmark 基因（MOESM6）和原文 marker DEG（MOESM7）给每个细胞打分，检查注释一致性
# B 部分：与原作者发布注释逐细胞比对（需要 cell metadata，目前拿不到则自动跳过）
# C 部分：本文 DEG 与原文 Supplementary Data 1–3 逐细胞类型比较（重叠、方向、核心基因）
# D 部分：给稿件填空用的数字
###############################################################################
source(here::here("revision_scripts", "00_setup.R"))
library(openxlsx)

BASE <- "/home/xxyloveme/study/BPD+MSCs"
OUT  <- file.path(BASE, "12_Revision")
od   <- file.path(OUT, "B1_annotation")
dir.create(od, recursive = TRUE, showWarnings = FALSE)

ORIG_DIR <- file.path(BASE, "Supplementary_tables_original")
f_orig <- function(n) file.path(ORIG_DIR, sprintf("41467_2025_60371_MOESM%d_ESM.xlsx", n))

## 检查文件是否都在（应全部为 TRUE）
file.exists(f_orig(3:7))

## 查看 MOESM6 中工作表的确切名称
getSheetNames(f_orig(6))

MESO_TYPES <- c("Alveolar FB", "Activated FB", "Adventitial FB", "Alveolar MyoFB",
                "Ductal MyoFB", "Pericyte", "VSMC")

## =============================================================================
## A. Marker 打分一致性：每个细胞按原作者 marker 集合打分 -> 取最高分 -> 与注释比较
## =============================================================================
## A1. 原文 Supplementary Data 4（MOESM6）中的 hallmark 基因；去掉标注 "(neg)" 的阴性 marker
hm <- read.xlsx(f_orig(6), sheet = "Hallmark Gene List")
hm <- tidyr::fill(hm, Lineage)              # 原表只在每个谱系第一行写了 Lineage
hm <- hm[hm$Lineage == "Mesenchymal", ]
meso_markers <- setNames(lapply(strsplit(hm$Hallmark.Genes, ","), function(x) {
  x <- trimws(x); x[!grepl("^\\(neg\\)", x)] }), trimws(hm$Cell.Type))
meso_markers <- lapply(meso_markers, intersect, rownames(meso_relabel))

meso <- ensure_rna_norm(meso_relabel)

# 查看当前有哪些 layer
Layers(meso[["RNA"]])

# 合并 layer
meso[["RNA"]] <- JoinLayers(meso[["RNA"]])

# 合并后应该只剩 counts 和 data（可能还有 scale.data）
Layers(meso[["RNA"]])

# 重新计算打分
meso <- AddModuleScore(meso, features = meso_markers, name = "MS_",
                       assay = "RNA", seed = 2026)


sc_cols <- paste0("MS_", seq_along(meso_markers))
scores  <- meso@meta.data[, sc_cols]; colnames(scores) <- names(meso_markers)
meso$marker_call <- names(meso_markers)[max.col(as.matrix(scores), ties.method = "first")]

conf <- table(Annotation = meso$celltype, MarkerCall = meso$marker_call)
conf_pct <- round(prop.table(conf, 1) * 100, 1)
agree <- mean(as.character(meso$celltype) == meso$marker_call) * 100
ari_marker <- mclust::adjustedRandIndex(as.character(meso$celltype), meso$marker_call)
cat(sprintf("\n[A] 注释与 marker 打分一致率 = %.1f%%, ARI = %.3f\n", agree, ari_marker))
write.csv(as.data.frame.matrix(conf),     file.path(od, "A_confusion_counts.csv"))
write.csv(as.data.frame.matrix(conf_pct), file.path(od, "A_confusion_rowpct.csv"))

# 每个 louvain cluster 的 marker 打分均值：检查 "原代码 cluster->label 映射" 在重跑后是否仍正确
clu_col <- if ("louvain" %in% colnames(meso@meta.data)) "louvain" else "seurat_clusters"
cl_tab <- cbind(meso@meta.data[, c(clu_col, "celltype")], scores) %>%
  group_by(across(all_of(c(clu_col, "celltype")))) %>%
  summarise(n_cells = n(), across(all_of(names(meso_markers)), mean), .groups = "drop")
cl_tab$best_marker_set <- names(meso_markers)[max.col(as.matrix(cl_tab[, names(meso_markers)]))]
cl_tab$label_supported <- cl_tab$best_marker_set == cl_tab$celltype
write.csv(cl_tab, file.path(od, "A_cluster_label_check.csv"), row.names = FALSE)
cat("\n[A] 以下 cluster 的标签与 marker 打分不一致（需人工复核，若为空则全部一致）：\n")
print(cl_tab %>% filter(!label_supported) %>% select(all_of(clu_col), celltype, best_marker_set, n_cells))

p_conf <- as.data.frame(conf_pct) %>%
  ggplot(aes(MarkerCall, Annotation, fill = Freq)) + geom_tile(color = "white") +
  geom_text(aes(label = ifelse(Freq >= 1, Freq, "")), size = 3) +
  scale_fill_gradient(low = "white", high = "#2166AC", name = "% of row") +
  labs(x = "Highest marker-set score", y = "Cell annotation",
       title = sprintf("Agreement %.1f%%, ARI %.2f", agree, ari_marker)) +
  theme_minimal(11) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_both(p_conf, file.path(od, "A_marker_concordance_heatmap"), 7, 5.5)

## A2. 原文 Supplementary Data 5（MOESM7）：每种间充质细胞 vs 其他所有细胞的 marker；
##     取 avg_log2FC > 1、padj < 1e-5、pct.1 - pct.2 最大的前 30 个作为第二套独立 marker
top_m7 <- lapply(MESO_TYPES, function(ct) {
  d <- read.xlsx(f_orig(7), sheet = ct, rowNames = TRUE)
  d$diff <- d$pct.1 - d$pct.2
  d <- d[d$avg_log2FC > 1 & d$p_val_adj < 1e-5, ]
  head(rownames(d)[order(-d$diff)], 30)
})
names(top_m7) <- MESO_TYPES
top_m7 <- lapply(top_m7, intersect, rownames(meso))
meso <- AddModuleScore(meso, features = top_m7, name = "M7_", assay = "RNA", seed = 2026)
sc7 <- meso@meta.data[, paste0("M7_", seq_along(top_m7))]; colnames(sc7) <- names(top_m7)
meso$marker_call_m7 <- names(top_m7)[max.col(as.matrix(sc7), ties.method = "first")]
agree_m7 <- mean(as.character(meso$celltype) == meso$marker_call_m7) * 100
ari_m7   <- mclust::adjustedRandIndex(as.character(meso$celltype), meso$marker_call_m7)
cat(sprintf("[A2] 原文 MOESM7 marker：一致率 = %.1f%%, ARI = %.3f\n", agree_m7, ari_m7))
write.csv(round(prop.table(table(Annotation = meso$celltype, MOESM7_call = meso$marker_call_m7), 1) * 100, 1),
          file.path(od, "A2_confusion_rowpct_MOESM7_markers.csv"))

# 已发表的 Fig. 2E 式点图（原作者 Fig.2A 同款 marker）——放入补充图作为注释证据
dp <- DotPlot(meso, assay = "RNA", group.by = "celltype",
              features = c("WNT2","TCF21","CTHRC1","POSTN","MFAP5","SCARA5","PI16",
                           "PDGFRA","FGF18","WNT5A","HHIP","COX4I2","CSPG4","ITGA7","NTRK3")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_both(dp, file.path(od, "A_marker_dotplot"), 9, 4.5)

## =============================================================================
## B. 与原作者发布的细胞注释逐细胞比对（GSE275938_cell_metadata.csv.gz）   —— R1-1 的直接证据
##    原表列：id（如 ACTGCAATCAATGTCG-1_1_1_1）, celltype, celltype_lineage, age, dataset（= 本文 donor）, ...
##    匹配键 = dataset + 16 bp 细胞条形码（不依赖 merge 时添加的 _1_1 等后缀）
##    B1  间充质：本文 7 个亚群 vs 原文注释（交叉表、ARI、Sankey、热图）
##    B2  全部细胞：本文 lineage / celltype vs 原文（若 alldata 可读）
##    B3  用原文自己的注释计算每个供者的间充质组成 —— 独立验证 Adventitial FB 在 BPD 中扩增
##    B4  未匹配细胞（只在本文 / 只在原文）按细胞类型统计
## =============================================================================
ORIG_META <- file.path(BASE, "GSE275938_cell_metadata.csv.gz")      # <- 文件实际位置
orig_all  <- read.csv(ORIG_META, check.names = FALSE)
core_bc   <- function(x) sub("^.*?([ACGTN]{16}).*$", "\\1", x, perl = TRUE)
orig_all$key <- paste(orig_all$dataset, core_bc(orig_all$id), sep = "|")
stopifnot(!anyDuplicated(orig_all$key))
cat(sprintf("\n[B] 原文注释：%d 个细胞，%d 种细胞类型\n", nrow(orig_all), length(unique(orig_all$celltype))))

MESO_ORDER <- MESO_TYPES
match_tab <- function(our_cells, our_donor, our_label, orig_df) {
  ours <- data.frame(cell = our_cells, donor = as.character(our_donor),
                     ours = as.character(our_label))
  ours$key <- paste(ours$donor, core_bc(ours$cell), sep = "|")
  stopifnot(!anyDuplicated(ours$key))
  list(ours = ours,
       m = inner_join(ours, orig_df[, c("key", "celltype", "celltype_lineage")], by = "key"))
}

## ---- B1. 间充质 ---------------------------------------------------------------
orig_meso <- orig_all[orig_all$celltype_lineage == "Mesenchymal", ]
mt <- match_tab(colnames(meso), meso$donor, meso$celltype, orig_all)   # 与原文全部细胞匹配，看本文间充质在原文被注释成什么
oursM <- mt$ours; mM <- mt$m
cat(sprintf("[B1] 本文间充质 %d 个细胞，其中 %d (%.1f%%) 在原文注释中找到\n",
            nrow(oursM), nrow(mM), 100 * nrow(mM) / nrow(oursM)))
cat(sprintf("[B1] 原文间充质 %d 个细胞，其中 %d (%.1f%%) 在本文间充质中找到\n",
            nrow(orig_meso), sum(orig_meso$key %in% oursM$key),
            100 * mean(orig_meso$key %in% oursM$key)))

mM$ours     <- factor(mM$ours, levels = MESO_ORDER)
mM$original <- factor(mM$celltype, levels = c(MESO_ORDER, setdiff(sort(unique(mM$celltype)), MESO_ORDER)))
tabB <- table(This_study = mM$ours, Shirazi = mM$original)
tabB <- tabB[, colSums(tabB) > 0, drop = FALSE]
row_pct <- round(prop.table(tabB, 1) * 100, 1)
col_pct <- round(prop.table(tabB, 2) * 100, 1)

mm <- mM[mM$celltype_lineage == "Mesenchymal", ]
ari_meso   <- mclust::adjustedRandIndex(as.character(mm$ours), as.character(mm$original))
agree_meso <- mean(as.character(mm$ours) == as.character(mm$original)) * 100
per_type <- mm %>% mutate(celltype = as.character(ours)) %>% group_by(celltype) %>%
  summarise(n_this_study = n(), same_label_pct = round(100 * mean(as.character(ours) == as.character(original)), 1),
            main_original_label = names(which.max(table(original))), .groups = "drop") %>%
  left_join(mm %>% mutate(celltype = as.character(original)) %>% group_by(celltype) %>%
              summarise(n_original = n(),
                        recovered_pct = round(100 * mean(as.character(ours) == as.character(original)), 1),
                        .groups = "drop"), by = "celltype")
cat(sprintf("[B1] 间充质细胞（双方都判为间充质）：标签一致 %.1f%%，ARI = %.3f\n", agree_meso, ari_meso))
print(per_type)

openxlsx::write.xlsx(list(counts = as.data.frame.matrix(tabB),
                          row_pct_this_study = as.data.frame.matrix(row_pct),
                          col_pct_original = as.data.frame.matrix(col_pct),
                          per_celltype = per_type,
                          summary = data.frame(
                            metric = c("this-study mesenchymal cells", "matched to original",
                                       "original mesenchymal cells", "original mesenchymal found here",
                                       "label agreement (%)", "ARI"),
                            value = c(nrow(oursM), nrow(mM), nrow(orig_meso),
                                      sum(orig_meso$key %in% oursM$key), round(agree_meso, 1), round(ari_meso, 3)))),
                     file.path(od, "B1_mesenchyme_this_study_vs_Shirazi.xlsx"), rowNames = TRUE)

# 热图（行 = 本文，列 = 原文；数值 = 行百分比）
p_tab <- as.data.frame(row_pct) %>%
  ggplot(aes(Shirazi, This_study, fill = Freq)) + geom_tile(colour = "white") +
  geom_text(aes(label = ifelse(Freq >= 1, Freq, "")), size = 3) +
  scale_fill_gradient(low = "white", high = "#2166AC", name = "% of row", limits = c(0, 100)) +
  scale_y_discrete(limits = rev) +
  labs(x = "Annotation in Shirazi et al.", y = "Annotation in this study",
       title = sprintf("Mesenchymal cells: %.1f%% identical labels (ARI = %.2f)", agree_meso, ari_meso)) +
  theme_minimal(11) + theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())
save_both(p_tab, file.path(od, "B1_heatmap_this_study_vs_Shirazi"), 7.5, 5.5)

# Sankey
library(ggalluvial)
library(ggrepel)
library(dplyr)

sk <- as.data.frame(tabB) %>% filter(Freq > 0) %>%
  mutate(This_study = as.character(This_study), Shirazi = as.character(Shirazi))

## 1. 合并占比 <0.5% 的小类别为 "Other"
tot <- sum(sk$Freq)
small_l <- sk %>% group_by(This_study) %>% summarise(n = sum(Freq)) %>% filter(n / tot < 0.005) %>% pull(This_study)
small_r <- sk %>% group_by(Shirazi)    %>% summarise(n = sum(Freq)) %>% filter(n / tot < 0.005) %>% pull(Shirazi)
sk <- sk %>%
  mutate(This_study = ifelse(This_study %in% small_l, "Other", This_study),
         Shirazi    = ifelse(Shirazi    %in% small_r, "Other", Shirazi)) %>%
  group_by(This_study, Shirazi) %>% summarise(Freq = sum(Freq), .groups = "drop")

## 2. 统一两侧顺序；标记一致 / 不一致
lv <- c(intersect(MESO_TYPES, union(sk$This_study, sk$Shirazi)),
        setdiff(union(sk$This_study, sk$Shirazi), c(MESO_TYPES, "Other")),
        if ("Other" %in% c(sk$This_study, sk$Shirazi)) "Other")
sk <- sk %>%
  mutate(This_study = factor(This_study, levels = lv),
         Shirazi    = factor(Shirazi,    levels = lv),
         flow = ifelse(as.character(This_study) == as.character(Shirazi), "Concordant", "Discordant"))
pct_conc <- 100 * sum(sk$Freq[sk$flow == "Concordant"]) / sum(sk$Freq)

## 3. 配色：细胞类型 + 流带类型共用一个 fill 标尺
ct_cols <- setNames(c("#4E79A7", "#F28E2B", "#59A14F", "#76B7B2",
                      "#EDC948", "#B07AA1", "#9C755F", "#FF9DA7",
                      "#BAB0AC", "#86BCB6")[seq_along(lv)], lv)
if ("Other" %in% lv) ct_cols["Other"] <- "grey60"
fill_all <- c(ct_cols, "Concordant" = "grey85", "Discordant" = "#E8825F")

w <- 0.06   # 色柱宽度

p_sk <- ggplot(sk, aes(axis1 = This_study, axis2 = Shirazi, y = Freq)) +
  geom_alluvium(aes(fill = flow), width = w, alpha = 0.75,
                curve_type = "sigmoid", colour = NA) +
  geom_stratum(aes(fill = after_stat(stratum)), width = w,
               colour = "white", linewidth = 0.8) +
  ## 左侧标签
  geom_text_repel(aes(label = ifelse(after_stat(x) == 1, as.character(after_stat(stratum)), NA)),
                  stat = "stratum", size = 3.2, direction = "y", hjust = 1,
                  nudge_x = -0.08, segment.size = 0.3, segment.colour = "grey50",
                  min.segment.length = 0, na.rm = TRUE) +
  ## 右侧标签
  geom_text_repel(aes(label = ifelse(after_stat(x) == 2, as.character(after_stat(stratum)), NA)),
                  stat = "stratum", size = 3.2, direction = "y", hjust = 0,
                  nudge_x = 0.08, segment.size = 0.3, segment.colour = "grey50",
                  min.segment.length = 0, na.rm = TRUE) +
  scale_fill_manual(values = fill_all,
                    breaks = c("Concordant", "Discordant"),
                    labels = c("Same annotation", "Different annotation"),
                    name = NULL) +
  scale_x_discrete(limits = c("This study", "Shirazi et al."),
                   expand = c(0.35, 0.35), position = "top") +
  labs(y = NULL,
       title = sprintf("ARI = %.2f   |   %.1f%% cells concordant", ari_meso, pct_conc)) +
  theme_void(11) +
  theme(axis.text.x = element_text(face = "bold", size = 11, margin = margin(b = 4)),
        plot.title = element_text(hjust = 0.5, size = 11),
        legend.position = "bottom",
        legend.key.size = unit(0.4, "cm"),
        plot.margin = margin(10, 10, 10, 10))

suppressWarnings(save_both(p_sk, file.path(od, "B1_sankey_this_study_vs_Shirazi"), 6.5, 7))
suppressWarnings(print(p_sk))



## 如未安装：devtools::install_github("davidsjoberg/ggsankey")
library(ggsankey)
library(dplyr)
library(tidyr)
library(ggplot2)

## ---------- 1. 数据整理 ----------
sk <- as.data.frame(tabB) %>% filter(Freq > 0) %>%
  mutate(This_study = as.character(This_study), Shirazi = as.character(Shirazi))

## 两侧占比 <0.5% 的小类别合并为 "Other"
tot <- sum(sk$Freq)
small_of <- function(col) sk %>% group_by(.data[[col]]) %>% summarise(n = sum(Freq)) %>%
  filter(n / tot < 0.005) %>% pull(1)
sl <- small_of("This_study"); sr <- small_of("Shirazi")
sk <- sk %>%
  mutate(This_study = ifelse(This_study %in% sl, "Other", This_study),
         Shirazi    = ifelse(Shirazi    %in% sr, "Other", Shirazi)) %>%
  group_by(This_study, Shirazi) %>% summarise(Freq = sum(Freq), .groups = "drop")

## 统一节点顺序（MESO_TYPES 在前，Other 最后）
lv <- c(intersect(MESO_TYPES, union(sk$This_study, sk$Shirazi)),
        setdiff(union(sk$This_study, sk$Shirazi), c(MESO_TYPES, "Other")),
        if ("Other" %in% c(sk$This_study, sk$Shirazi)) "Other")

## 展开为每个细胞一行，转成 ggsankey 的长格式；x 用数值，方便把标签挪到两侧
df <- sk %>% uncount(Freq) %>%
  make_long(This_study, Shirazi) %>%
  mutate(x      = as.numeric(x),
         next_x = as.numeric(next_x),
         node      = factor(node,      levels = rev(lv)),   # rev：第一个类别在最上方
         next_node = factor(next_node, levels = rev(lv)))

## ---------- 2. 配色（柔和色系，与参考图风格一致） ----------
pal <- c("#944D55", "#D05C5F", "#EAA475", "#96D1E2", "#4E74B4",
         "#7FA67A", "#B39CC8", "#C9A27E", "#E8C468", "#7BA7A0")
cols <- setNames(pal[seq_along(lv)], lv)
if ("Other" %in% lv) cols["Other"] <- "grey80"

## ---------- 3. 标签图层：节点位置不变，x 向外平移 ----------
SPACE <- round(0.02 * tot)     # 节点间隙，按总细胞数的 2% 设置
lab_l <- df %>% filter(x == 1) %>% mutate(x = 1 - 0.09, next_x = NA, next_node = NA)
lab_r <- df %>% filter(x == 2) %>% mutate(x = 2 + 0.09, next_x = NA, next_node = NA)

## ---------- 4. 绘图 ----------
p_sk <- ggplot(df, aes(x = x, next_x = next_x, node = node,
                       next_node = next_node, fill = node)) +
  geom_sankey(width = 0.08, smooth = 6, space = SPACE,
              flow.alpha = 0.45, node.color = "transparent", na.rm = TRUE) +
  geom_sankey_text(data = lab_l, aes(x = x, next_x = next_x, node = node,
                                     next_node = next_node, label = node),
                   inherit.aes = FALSE, space = SPACE, hjust = 1,
                   size = 3.3, colour = "black") +
  geom_sankey_text(data = lab_r, aes(x = x, next_x = next_x, node = node,
                                     next_node = next_node, label = node),
                   inherit.aes = FALSE, space = SPACE, hjust = 0,
                   size = 3.3, colour = "black") +
  scale_fill_manual(values = cols) +
  scale_x_continuous(breaks = c(1, 2), labels = c("This study", "Shirazi et al."),
                     limits = c(0.45, 2.55), position = "top") +
  labs(title = sprintf("ARI = %.2f", ari_meso)) +
  theme_sankey(base_size = 12) +
  theme(axis.title.x = element_blank(),
        axis.text.x  = element_text(colour = "black", size = 12, face = "bold"),
        plot.title   = element_text(hjust = 0.5, size = 11, colour = "grey30"),
        legend.position = "none",
        plot.margin = margin(5, 5, 5, 5))

save_both(p_sk, file.path(od, "B1_sankey_this_study_vs_Shirazi"), 7, 6)
p_sk

## ---- B2. 全部细胞（lineage 与 celltype 层面） ----------------------------------
# if (file.exists(PATH_ALLDATA)) {
#   alldB <- add_donor(readRDS(PATH_ALLDATA))
#   our_lin <- if ("bulk_celltype" %in% colnames(alldB@meta.data)) alldB$bulk_celltype else alldB$celltype_lineage
#   mtA <- match_tab(colnames(alldB), alldB$donor, alldB$celltype, orig_all)
#   mA_all <- mtA$m
#   mA_all$our_lineage <- as.character(our_lin[match(mA_all$cell, colnames(alldB))])
#   rm(alldB); gc()
#   lin_tab <- table(This_study = mA_all$our_lineage, Shirazi = mA_all$celltype_lineage)
#   ari_lin <- mclust::adjustedRandIndex(mA_all$our_lineage, mA_all$celltype_lineage)
#   ari_ct  <- mclust::adjustedRandIndex(mA_all$ours, mA_all$celltype)
#   sum_all <- data.frame(
#     metric = c("cells this study", "cells original", "matched",
#                "lineage agreement (%)", "lineage ARI", "cell-type agreement (%)", "cell-type ARI"),
#     value = c(nrow(mtA$ours), nrow(orig_all), nrow(mA_all),
#               round(100 * mean(mA_all$our_lineage == mA_all$celltype_lineage), 1), round(ari_lin, 3),
#               round(100 * mean(tolower(mA_all$ours) == tolower(mA_all$celltype)), 1), round(ari_ct, 3)))
#   print(sum_all); print(lin_tab)
#   openxlsx::write.xlsx(list(summary = sum_all, lineage_counts = as.data.frame.matrix(lin_tab),
#                             celltype_counts = as.data.frame.matrix(table(mA_all$ours, mA_all$celltype))),
#                        file.path(od, "B2_all_cells_this_study_vs_Shirazi.xlsx"), rowNames = TRUE)
# }
# 


## ---- B2. 全部细胞（lineage 与 celltype 层面）：表格 + 热图 + Sankey -------------
library(ggsankey)   # 如未安装：devtools::install_github("davidsjoberg/ggsankey")
library(tidyr)

LIN_ORDER <- c("Epithelial", "Endothelial", "Mesenchymal", "Immune")
CT_ORDER  <- c("AT1", "AT2", "RASC", "Secretory -3A1, -3A2", "Secretory MUC5B", "Basal", "Multiciliated",
               "gCap", "aCap", "abCap", "Arterial EC", "Pulmonary venous EC", "Systemic venous EC", "Lymphatic",
               "Alveolar FB", "Activated FB", "Adventitial FB", "Alveolar MyoFB", "Ductal MyoFB", "Pericyte", "VSMC",
               "T Cell", "NK Cell", "NKT Cell", "Alveolar Macrophage", "Monocyte", "pDC", "cDC",
               "B Cell", "Mast cell", "Plasma cell", "Basophil", "Neutrophil")
LIN_COLS <- c(Epithelial = "#E69F00", Endothelial = "#D55E00", Mesenchymal = "#009E73", Immune = "#0072B2")
# 统一大小写（本文 "Alveolar macrophage" vs 原文 "Alveolar Macrophage"）
harmonize_ct <- function(x) { i <- match(tolower(x), tolower(CT_ORDER)); ifelse(is.na(i), x, CT_ORDER[i]) }

# 通用热图：行 = 本文，列 = 原文，数值 = 行百分比
plot_conc_heat <- function(tab, title, xlab = "Annotation in Shirazi et al.", ylab = "Annotation in this study",
                           text_size = 3, min_label = 1) {
  pct <- as.data.frame(round(prop.table(tab, 1) * 100, 1))
  colnames(pct) <- c("This_study", "Shirazi", "Freq")
  ggplot(pct, aes(Shirazi, This_study, fill = Freq)) + geom_tile(colour = "white", linewidth = 0.3) +
    geom_text(aes(label = ifelse(Freq >= min_label, Freq, "")), size = text_size) +
    scale_fill_gradient(low = "white", high = "#2166AC", name = "% of row", limits = c(0, 100)) +
    scale_y_discrete(limits = rev) + coord_fixed() +
    labs(x = xlab, y = ylab, title = title) +
    theme_minimal(10) + theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())
}

# 通用 Sankey（ggsankey 版）：左 = 本文，右 = 原文；按 lineage 着色（同一 lineage 内为深浅渐变）
# 去掉细胞数很少的流向以免太乱；max_rows：细胞过多时按比例抽稀后再展开，避免内存不足
plot_conc_sankey <- function(df, left, right, fill_by, title, min_cells = 10, text_size = 2.6,
                             max_rows = 2e5) {
  ## 1. 统计流向，并去掉细胞数很少的流向
  sk <- df %>% count(L = .data[[left]], R = .data[[right]], name = "Freq") %>%
    filter(Freq >= min_cells, !is.na(L), !is.na(R)) %>%
    mutate(L = as.character(L), R = as.character(R))
  
  ## 2. 节点顺序：沿用因子水平（LIN_ORDER / CT_ORDER）
  lv0 <- if (is.factor(df[[left]])) levels(df[[left]]) else sort(unique(c(sk$L, sk$R)))
  lv  <- c(intersect(lv0, union(sk$L, sk$R)), setdiff(union(sk$L, sk$R), lv0))
  
  ## 3. 每个节点所属的 lineage：左侧取 fill_by 的众数；右侧用 CT_ORDER 映射
  ct2lin  <- setNames(rep(LIN_ORDER, c(7, 7, 7, 12)), CT_ORDER)
  lin_map <- df %>% filter(!is.na(.data[[left]]), !is.na(.data[[fill_by]])) %>%
    count(N = as.character(.data[[left]]), Lin = as.character(.data[[fill_by]])) %>%
    group_by(N) %>% slice_max(n, n = 1, with_ties = FALSE) %>% ungroup()
  node_lin <- setNames(lin_map$Lin, lin_map$N)
  miss <- setdiff(lv, names(node_lin))
  node_lin[miss] <- ifelse(miss %in% LIN_ORDER, miss, ct2lin[miss])
  node_lin <- node_lin[lv]; node_lin[is.na(node_lin)] <- "Unknown"
  
  ## 4. 节点颜色：同一 lineage 内由深到浅
  cols <- character(0)
  for (l in unique(node_lin)) {
    nn   <- lv[node_lin == l]
    base <- if (l %in% names(LIN_COLS)) LIN_COLS[[l]] else "grey60"
    cols[nn] <- if (length(nn) == 1) base else
      colorRampPalette(c(base, "white"))(length(nn) + 2)[seq_along(nn)]
  }
  
  ## 5. 展开为每个细胞一行（细胞过多时按比例抽稀）
  tot <- sum(sk$Freq); k <- max(1, ceiling(tot / max_rows))
  sk$n_draw <- pmax(1, round(sk$Freq / k))
  long <- sk %>% select(L, R, n_draw) %>% uncount(n_draw) %>%
    make_long(L, R) %>%
    mutate(x = as.numeric(x), next_x = as.numeric(next_x),
           node      = factor(node,      levels = rev(lv)),
           next_node = factor(next_node, levels = rev(lv)))
  
  ## 6. 节点间隙：所有间隙加起来约占总高度的 25%
  n_draw <- sum(sk$n_draw)
  SPACE  <- max(1, round(0.25 * n_draw / max(1, length(lv))))
  
  ## 7. 标签图层：节点位置不变，x 向外平移
  off   <- 0.08
  lab_l <- long %>% filter(x == 1) %>% mutate(x = 1 - off, next_x = NA, next_node = NA)
  lab_r <- long %>% filter(x == 2) %>% mutate(x = 2 + off, next_x = NA, next_node = NA)
  
  ## 8. 图例用的占位点（不可见，只用来生成 lineage 图例）
  lin_present <- intersect(LIN_ORDER, unique(node_lin))
  leg_df <- data.frame(Lineage = factor(lin_present, levels = LIN_ORDER), x = 1.5, y = 0)
  
  cap <- sprintf("Flows with < %d cells are not drawn", min_cells)
  if (k > 1) cap <- paste0(cap, sprintf("; flows drawn at 1:%d scale", k))
  
  ggplot(long, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = node)) +
    geom_sankey(width = 0.06, smooth = 6, space = SPACE,
                flow.alpha = 0.45, node.color = "transparent", na.rm = TRUE) +
    geom_sankey_text(data = lab_l, aes(x = x, next_x = next_x, node = node,
                                       next_node = next_node, label = node),
                     inherit.aes = FALSE, space = SPACE, hjust = 1,
                     size = text_size, colour = "black") +
    geom_sankey_text(data = lab_r, aes(x = x, next_x = next_x, node = node,
                                       next_node = next_node, label = node),
                     inherit.aes = FALSE, space = SPACE, hjust = 0,
                     size = text_size, colour = "black") +
    geom_point(data = leg_df, aes(x = x, y = y, colour = Lineage),
               inherit.aes = FALSE, alpha = 0) +
    scale_fill_manual(values = cols, guide = "none") +
    scale_colour_manual(values = LIN_COLS, name = "Lineage (this study)", drop = TRUE) +
    guides(colour = guide_legend(override.aes = list(alpha = 1, size = 4, shape = 15))) +
    scale_x_continuous(breaks = c(1, 2), labels = c("This study", "Shirazi et al."),
                       limits = c(0.4, 2.6), position = "top") +
    labs(title = title, caption = cap) +
    theme_sankey(base_size = 10) +
    theme(axis.title.x = element_blank(),
          axis.text.x  = element_text(colour = "black", size = 11, face = "bold"),
          plot.title   = element_text(hjust = 0.5, size = 11),
          plot.caption = element_text(size = 8, colour = "grey40"),
          legend.position = "bottom",
          plot.margin = margin(5, 5, 5, 5))
}

if (file.exists(PATH_ALLDATA)) {
  alldB <- add_donor(readRDS(PATH_ALLDATA))
  our_lin <- if ("bulk_celltype" %in% colnames(alldB@meta.data)) alldB$bulk_celltype else alldB$celltype_lineage
  mtA <- match_tab(colnames(alldB), alldB$donor, alldB$celltype, orig_all)
  mA_all <- mtA$m
  mA_all$our_lineage <- as.character(our_lin[match(mA_all$cell, colnames(alldB))])
  n_ours_all <- nrow(mtA$ours)
  rm(alldB); gc()
  
  mA_all <- mA_all %>% mutate(
    our_ct  = factor(harmonize_ct(ours),     levels = CT_ORDER),
    orig_ct = factor(harmonize_ct(celltype), levels = CT_ORDER),
    our_lin = factor(our_lineage,      levels = LIN_ORDER),
    orig_lin = factor(celltype_lineage, levels = LIN_ORDER))
  unk <- unique(c(mA_all$ours[is.na(mA_all$our_ct)], mA_all$celltype[is.na(mA_all$orig_ct)]))
  if (length(unk)) message("[B2] 以下标签不在 CT_ORDER 中，未纳入 cell-type 图表：", paste(unk, collapse = ", "))
  
  ## 汇总数字
  lin_tab <- table(This_study = mA_all$our_lin, Shirazi = mA_all$orig_lin)
  ct_tab  <- table(This_study = mA_all$our_ct,  Shirazi = mA_all$orig_ct)
  ct_tab  <- ct_tab[rowSums(ct_tab) > 0, colSums(ct_tab) > 0, drop = FALSE]
  ari_lin <- mclust::adjustedRandIndex(as.character(mA_all$our_lin), as.character(mA_all$orig_lin))
  ari_ct  <- mclust::adjustedRandIndex(as.character(mA_all$our_ct),  as.character(mA_all$orig_ct))
  agree_lin <- 100 * mean(as.character(mA_all$our_lin) == as.character(mA_all$orig_lin))
  agree_ct  <- 100 * mean(as.character(mA_all$our_ct)  == as.character(mA_all$orig_ct))
  sum_all <- data.frame(
    metric = c("cells this study", "cells original", "matched cells", "matched (% of this study)",
               "lineage agreement (%)", "lineage ARI", "cell-type agreement (%)", "cell-type ARI"),
    value = c(n_ours_all, nrow(orig_all), nrow(mA_all), round(100 * nrow(mA_all) / n_ours_all, 1),
              round(agree_lin, 1), round(ari_lin, 3), round(agree_ct, 1), round(ari_ct, 3)))
  print(sum_all); print(lin_tab)
  
  # 每个细胞类型：本文 -> 原文同名比例（precision 视角）与 原文 -> 本文（recall 视角）
  per_ct_all <- full_join(
    mA_all %>% group_by(celltype = as.character(our_ct)) %>%
      summarise(lineage = as.character(first(our_lin)), n_this_study = n(),
                pct_same_in_original = round(100 * mean(as.character(our_ct) == as.character(orig_ct)), 1),
                main_original_label = names(which.max(table(as.character(orig_ct)))), .groups = "drop"),
    mA_all %>% group_by(celltype = as.character(orig_ct)) %>%
      summarise(n_original = n(),
                pct_recovered_in_this_study = round(100 * mean(as.character(our_ct) == as.character(orig_ct)), 1),
                .groups = "drop"),
    by = "celltype") %>% arrange(match(celltype, CT_ORDER))
  # 每个 lineage 内部的 cell-type ARI
  per_lin_ari <- mA_all %>% filter(our_lin == orig_lin) %>% group_by(lineage = our_lin) %>%
    summarise(n = n(), agreement_pct = round(100 * mean(as.character(our_ct) == as.character(orig_ct)), 1),
              ARI = round(mclust::adjustedRandIndex(as.character(our_ct), as.character(orig_ct)), 3), .groups = "drop")
  print(per_lin_ari)
  
  openxlsx::write.xlsx(list(summary = sum_all, per_lineage = per_lin_ari, per_celltype = per_ct_all,
                            lineage_counts = as.data.frame.matrix(lin_tab),
                            lineage_row_pct = as.data.frame.matrix(round(prop.table(lin_tab, 1) * 100, 1)),
                            celltype_counts = as.data.frame.matrix(ct_tab),
                            celltype_row_pct = as.data.frame.matrix(round(prop.table(ct_tab, 1) * 100, 1))),
                       file.path(od, "B2_all_cells_this_study_vs_Shirazi.xlsx"), rowNames = TRUE)
  
  ## 图 1：lineage 热图（4 × 4）
  p_lin_h <- plot_conc_heat(lin_tab, sprintf("All cells, lineage: %.1f%% identical (ARI = %.2f)", agree_lin, ari_lin),
                            text_size = 3.5)
  save_both(p_lin_h, file.path(od, "B2a_heatmap_lineage_all_cells"), 5.5, 4.5)
  
  ## 图 2：cell-type 热图（33 × 33），用虚线框出四个 lineage
  p_ct_h <- plot_conc_heat(ct_tab, sprintf("All cells, cell type: %.1f%% identical (ARI = %.2f)", agree_ct, ari_ct),
                           text_size = 1.9, min_label = 5)
  lin_of <- function(v) setNames(rep(LIN_ORDER, c(7, 7, 7, 12)), CT_ORDER)[v]
  xs <- colnames(ct_tab); ys <- rev(rownames(ct_tab))
  boxes <- bind_rows(lapply(LIN_ORDER, function(l) {
    xi <- which(lin_of(xs) == l); yi <- which(lin_of(ys) == l)
    if (!length(xi) || !length(yi)) return(NULL)
    data.frame(xmin = min(xi) - 0.5, xmax = max(xi) + 0.5, ymin = min(yi) - 0.5, ymax = max(yi) + 0.5)
  }))
  p_ct_h <- p_ct_h + geom_rect(data = boxes, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                               inherit.aes = FALSE, fill = NA, colour = "grey20", linetype = 2, linewidth = 0.3)
  save_both(p_ct_h, file.path(od, "B2b_heatmap_celltype_all_cells"), 11, 10)
  
  ## 图 3：lineage Sankey
  p_lin_s <- plot_conc_sankey(mutate(mA_all, L = our_lin, R = orig_lin), "L", "R", "our_lin",
                              sprintf("All cells, lineage (ARI = %.2f)", ari_lin), text_size = 3.2)
  save_both(p_lin_s, file.path(od, "B2c_sankey_lineage_all_cells"), 6, 6)
  
  ## 图 4：cell-type Sankey（33 种，按 lineage 着色）
  p_ct_s <- plot_conc_sankey(mutate(mA_all, L = our_ct, R = orig_ct), "L", "R", "our_lin",
                             sprintf("All cells, cell type (ARI = %.2f)", ari_ct), min_cells = 10, text_size = 2.3)
  save_both(p_ct_s, file.path(od, "B2d_sankey_celltype_all_cells"), 8, 14)
  
  ## 图 5：四个 lineage 各一张 cell-type 热图（便于放进补充图分面）
  for (l in LIN_ORDER) {
    sub_l <- mA_all %>% filter(our_lin == l, orig_lin == l) %>% droplevels()
    if (!nrow(sub_l)) next
    t_l <- table(This_study = sub_l$our_ct, Shirazi = sub_l$orig_ct)
    a_l <- per_lin_ari$ARI[per_lin_ari$lineage == l]; g_l <- per_lin_ari$agreement_pct[per_lin_ari$lineage == l]
    p_l <- plot_conc_heat(t_l, sprintf("%s: %.1f%% identical (ARI = %.2f)", l, g_l, a_l), text_size = 2.8)
    w <- 2 + 0.45 * ncol(t_l); save_both(p_l, file.path(od, paste0("B2e_heatmap_", l)), w, w * 0.85)
  }
}


## ---- B3. 原文注释下每个供者的间充质组成：独立验证 Adventitial FB 扩增 ----------
comp_orig <- orig_meso %>% count(donor = dataset, celltype) %>%
  group_by(donor) %>% mutate(pct = 100 * n / sum(n)) %>% ungroup() %>% mutate(annotation = "Shirazi et al.")
comp_ours <- meso@meta.data %>% count(donor = as.character(donor), celltype = as.character(celltype)) %>%
  group_by(donor) %>% mutate(pct = 100 * n / sum(n)) %>% ungroup() %>% mutate(annotation = "This study")
comp_both <- bind_rows(comp_orig, comp_ours) %>%
  mutate(group = case_when(grepl("^Term", donor) ~ "Term infant", grepl("^BPD\\+PH", donor) ~ "BPD+PH",
                           grepl("^BPD", donor) ~ "BPD", TRUE ~ "Acute preterm injury"),
         donor = factor(donor, levels = DONOR_ORDER), celltype = factor(celltype, levels = MESO_ORDER))
write.csv(comp_both, file.path(od, "B3_mesenchymal_composition_per_donor_both_annotations.csv"), row.names = FALSE)

adv_both <- comp_both %>% filter(celltype == "Adventitial FB") %>%
  select(annotation, donor, group, n, pct) %>% arrange(annotation, donor)
cat("\n[B3] 每个供者 Adventitial FB 占间充质比例（两套注释）：\n"); print(adv_both, n = Inf)
pooled_both <- bind_rows(
  orig_meso %>% mutate(group = case_when(grepl("^Term", dataset) ~ "Term infant", grepl("^BPD\\+PH", dataset) ~ "BPD+PH",
                                         grepl("^BPD", dataset) ~ "BPD", TRUE ~ "Acute preterm injury"),
                       adv = celltype == "Adventitial FB", annotation = "Shirazi et al.") %>%
    select(annotation, group, adv),
  meso@meta.data %>% transmute(annotation = "This study", group = as.character(group),
                               adv = celltype == "Adventitial FB")) %>%
  group_by(annotation, group) %>% summarise(pooled_pct_AdvFB = round(100 * mean(adv), 2), n_meso = n(), .groups = "drop")
print(pooled_both)
write.csv(pooled_both, file.path(od, "B3_pooled_AdvFB_pct_both_annotations.csv"), row.names = FALSE)

p_adv <- adv_both %>% filter(group != "Acute preterm injury") %>%
  mutate(group = factor(group, levels = c("Term infant", "BPD", "BPD+PH"))) %>%
  ggplot(aes(group, pct, colour = annotation, shape = annotation)) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  stat_summary(aes(group = annotation), fun = mean, geom = "crossbar", width = 0.35,
               linewidth = 0.3, position = position_dodge(width = 0.5)) +
  scale_colour_manual(values = c("Shirazi et al." = "#2166AC", "This study" = "#D6604D")) +
  labs(x = NULL, y = "Adventitial FB (% of mesenchymal cells)",
       title = "Each point = one donor") + theme_bw(11)
save_both(p_adv, file.path(od, "B3_AdvFB_pct_per_donor_both_annotations"), 6, 4.5)

if (requireNamespace("speckle", quietly = TRUE)) {
  bt_o <- orig_meso %>% filter(grepl("^Term|^BPD [12]$", dataset)) %>%
    mutate(group = ifelse(grepl("^Term", dataset), "Term infant", "BPD"))
  prop_o <- speckle::propeller(clusters = bt_o$celltype, sample = bt_o$dataset,
                               group = bt_o$group, transform = "asin")
  write.csv(prop_o, file.path(od, "B3_propeller_original_annotation_BPD_vs_Term.csv"))
  print(prop_o)
}

## ---- B4. 未匹配细胞 -----------------------------------------------------------
only_ours <- oursM[!oursM$key %in% orig_all$key, ]
only_orig <- orig_meso[!orig_meso$key %in% oursM$key, ]
in_orig_other_lineage <- mM[mM$celltype_lineage != "Mesenchymal", ]
unmatched <- list(
  only_this_study_by_type = as.data.frame(table(celltype = only_ours$ours)),
  only_original_by_type   = as.data.frame(table(celltype = only_orig$celltype)),
  this_study_meso_but_original_other_lineage = as.data.frame(
    table(this_study = in_orig_other_lineage$ours, original = in_orig_other_lineage$celltype)) %>% filter(Freq > 0))
openxlsx::write.xlsx(unmatched, file.path(od, "B4_unmatched_cells.xlsx"))
cat(sprintf("\n[B4] 仅在本文间充质：%d 个；仅在原文间充质：%d 个；本文判为间充质而原文判为其他谱系：%d 个\n",
            nrow(only_ours), nrow(only_orig), nrow(in_orig_other_lineage)))

## =============================================================================
## C. 本文 DEG 与原文 Supplementary Data 1–3 比较
##    三组对比分别对应：markers_compare1 ↔ MOESM3；markers_compare2 ↔ MOESM4；markers_compare3 ↔ MOESM5
##    两边都按原文图例的同一标准判定显著：|avg_log2FC| > 1 且 p_val_adj < 1e-5，去除 ^RP / ^MT- 基因
## =============================================================================
read_deg <- function(file, sheet) {
  d <- read.xlsx(file, sheet = sheet, rowNames = TRUE, check.names = FALSE)
  d$gene <- rownames(d); d
}
call_sig <- function(d) {
  ok <- d$p_val_adj < 1e-5 & !grepl("^RP|^MT-", d$gene)
  ifelse(ok & d$avg_log2FC > 1, "up", ifelse(ok & d$avg_log2FC < -1, "down", "ns"))
}
contrasts <- list(
  "BPD_vs_Term"      = list(ours = markers_compare1, file = f_orig(3)),
  "BPDPH_vs_Term"    = list(ours = markers_compare2, file = f_orig(4)),
  "BPDPH_vs_BPD"     = list(ours = markers_compare3, file = f_orig(5)))
KEY_GENES <- unique(c(CORE_GENES, "XIST", "RPS4Y1", "DDX3Y", "SCARA5", "MFAP5"))
table(meso$celltype, meso$group)["Alveolar MyoFB", ]

need_cols <- c("avg_log2FC", "p_val_adj")

C_summary <- list(); C_key <- list(); C_genes <- list()
for (cn in names(contrasts)) {
  ours_list <- contrasts[[cn]]$ours
  orig_sheets <- getSheetNames(contrasts[[cn]]$file)
  for (ct in intersect(c(MESO_TYPES, "Mesenchyme"), intersect(names(ours_list), orig_sheets))) {
    o <- read_deg(contrasts[[cn]]$file, ct)
    u <- ours_list[[ct]]
    
    ## ---- 空表 / 缺列检查 ----
    bad_ours <- is.null(u) || nrow(u) == 0 || !all(need_cols %in% colnames(u))
    bad_orig <- nrow(o) == 0 || !all(need_cols %in% colnames(o))
    if (bad_ours || bad_orig) {
      message(sprintf("跳过 %s | %s：%s 无可用 DEG 结果", cn, ct,
                      ifelse(bad_ours, "ours", "原文")))
      C_summary[[length(C_summary) + 1]] <- data.frame(
        contrast = cn, celltype = ct,
        DEG_ours_total = NA, DEG_orig_total = NA, genes_in_both_tables = NA,
        DEG_ours_tested = NA, DEG_orig_tested = NA, overlap_same_direction = NA,
        pct_of_orig_recovered = NA, pct_of_ours_in_orig = NA,
        opposite_direction = NA, hypergeom_p = NA, spearman_log2FC = NA,
        note = ifelse(bad_ours, "ours: empty (too few cells?)", "orig: empty"))
      next
    }
    
    o$call <- call_sig(o)
    u$gene <- rownames(u); u$call <- call_sig(u)
    m <- inner_join(u[, c("gene", "avg_log2FC", "p_val_adj", "call")],
                    o[, c("gene", "avg_log2FC", "p_val_adj", "call")],
                    by = "gene", suffix = c("_ours", "_orig"))
    a <- m$gene[m$call_ours != "ns"]; b <- m$gene[m$call_orig != "ns"]
    k <- sum(m$call_ours != "ns" & m$call_ours == m$call_orig)
    hyp <- phyper(k - 1, length(b), nrow(m) - length(b), length(a), lower.tail = FALSE)
    C_summary[[length(C_summary) + 1]] <- data.frame(
      contrast = cn, celltype = ct,
      DEG_ours_total = sum(u$call != "ns"), DEG_orig_total = sum(o$call != "ns"),
      genes_in_both_tables = nrow(m),
      DEG_ours_tested = length(a), DEG_orig_tested = length(b),
      overlap_same_direction = k,
      pct_of_orig_recovered = round(100 * k / max(1, length(b)), 1),
      pct_of_ours_in_orig   = round(100 * k / max(1, length(a)), 1),
      opposite_direction = sum(m$call_ours != "ns" & m$call_orig != "ns" & m$call_ours != m$call_orig),
      hypergeom_p = signif(hyp, 3),
      spearman_log2FC = round(cor(m$avg_log2FC_ours, m$avg_log2FC_orig, method = "spearman"), 3),
      note = "")
    kg <- full_join(u[u$gene %in% KEY_GENES, c("gene", "avg_log2FC", "p_val_adj", "call")],
                    o[o$gene %in% KEY_GENES, c("gene", "avg_log2FC", "p_val_adj", "call")],
                    by = "gene", suffix = c("_ours", "_orig"))
    C_key[[length(C_key) + 1]] <- cbind(contrast = cn, celltype = ct, kg)
    C_genes[[paste(cn, ct, sep = "|")]] <- m
  }
}

C_summary_df <- dplyr::bind_rows(C_summary)
C_summary_df

C_summary <- bind_rows(C_summary); C_key <- bind_rows(C_key)
openxlsx::write.xlsx(list(summary = C_summary, key_genes = C_key,
                          AdvFB_BPD_vs_Term_merged = C_genes[["BPD_vs_Term|Adventitial FB"]]),
                     file.path(od, "C_DEG_comparison_with_original_SupplData1-3.xlsx"))
cat("\n[C] 本文 vs 原文 DEG 比较（同一判定标准）：\n"); print(C_summary)
cat("\n[C] 关键基因（Adventitial FB, BPD vs Term）：\n")
print(C_key %>% filter(contrast == "BPD_vs_Term", celltype == "Adventitial FB"))

# Adventitial FB, BPD vs Term：本文 vs 原文 log2FC 散点图
mA <- C_genes[["BPD_vs_Term|Adventitial FB"]]
mA$label <- ifelse(mA$gene %in% KEY_GENES, mA$gene, NA)
mA$sig_group <- factor(
  dplyr::case_when(
    mA$call_orig != "ns" & mA$call_ours != "ns" ~ "both",
    mA$call_orig != "ns"                         ~ "Shirazi et al",
    mA$call_ours != "ns"                         ~ "this study only",
    TRUE                                         ~ "neither"),
  levels = c("neither", "this study only", "Shirazi et al", "both"))
mA <- mA[order(mA$sig_group), ]          # 灰点在最下层，both 在最上层

fill_cols <- c("neither"         = "grey85",
               "this study only" = "#F4A582",   # 浅橙
               "Shirazi et al"   = "#2166AC",   # 深蓝
               "both"            = "#1B9E3E")   # 绿色
line_cols <- c("neither"         = "grey55",    # 灰点用灰色细边
               "this study only" = "black",
               "Shirazi et al"   = "black",
               "both"            = "black")

pC <- ggplot(mA, aes(avg_log2FC_orig, avg_log2FC_ours)) +
  geom_hline(yintercept = 0, colour = "grey30", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey30", linewidth = 0.3) +
  geom_abline(linetype = 2, colour = "grey50", linewidth = 0.3) +
  geom_point(aes(fill = sig_group, colour = sig_group),
             shape = 21, size = 1.1, stroke = 0.15, alpha = 0.75) +
  scale_fill_manual(values = fill_cols, name = "Significant in", drop = FALSE) +
  scale_colour_manual(values = line_cols, guide = "none") +
  guides(fill = guide_legend(override.aes = list(size = 3, alpha = 1,
                                                 colour = line_cols))) +
  labs(x = "log2FC, Shirazi et al. Suppl. Data 1", y = "log2FC, this study",
       title = "Adventitial FB: BPD vs Term infant") +
  theme_bw(11) +
  theme(panel.grid.minor = element_blank())

save_both(pC, file.path(od, "C_AdvFB_log2FC_ours_vs_original"), 7, 5.5)


## C2. 诊断：若本文 Adventitial FB DEG 数远多于原文，检查 FindMarkers 所用 assay 是否不同
##     原作者在 SCTransform 后的合并对象上直接 FindMarkers（默认 assay = SCT, slot = data）。
##     这里用重跑的合并对象分别以 SCT 与 RNA 重算一次，看哪一种能复现原文的 DEG 数。

RERUN_C2 <- TRUE
if (RERUN_C2 && file.exists(PATH_ALLDATA)) {
  alld <- readRDS(PATH_ALLDATA)
  advA <- subset(alld, subset = celltype == "Adventitial FB" & condition %in% c("BPD", "Term infant"))
  rm(alld); gc()
  Idents(advA) <- "condition"
  
  ## RNA assay 若是 v5 多 layer，先合并
  if (length(Layers(advA[["RNA"]], search = "counts")) > 1)
    advA[["RNA"]] <- JoinLayers(advA[["RNA"]])
  
  ## SCT：在子集上重新统一测序深度
  advA_prep <- PrepSCTFindMarkers(advA, assay = "SCT", verbose = FALSE)
  
  run_fm <- function(obj, as, tag, ...) {
    fm <- FindMarkers(obj, ident.1 = "BPD", ident.2 = "Term infant", assay = as,
                      test.use = "MAST", logfc.threshold = 0, ...)
    fm$gene <- rownames(fm); fm$call <- call_sig(fm)
    write.csv(fm, file.path(od, paste0("C2_AdvFB_BPD_vs_Term_rerun_", tag, ".csv")))
    data.frame(assay = tag,
               n_cells_BPD  = sum(obj$condition == "BPD"),
               n_cells_Term = sum(obj$condition == "Term infant"),
               DEG = sum(fm$call != "ns"),
               PHLDA1_log2FC = fm["PHLDA1", "avg_log2FC"],
               PHLDA1_padj   = fm["PHLDA1", "p_val_adj"])
  }
  
  ## 把 SCT 的 data 层复制成普通 assay，绕开 SCTAssay + MAST 的兼容问题
  sct_as_plain <- function(obj, new = "SCTd") {
    obj[[new]] <- CreateAssayObject(data = GetAssayData(obj, assay = "SCT", layer = "data"))
    obj
  }
  
  advA_prep <- PrepSCTFindMarkers(advA, assay = "SCT", verbose = FALSE)
  
  diag <- bind_rows(
    run_fm(sct_as_plain(advA_prep), "SCTd", "SCT_prep"),
    run_fm(sct_as_plain(advA),      "SCTd", "SCT_noRecorrect"),
    run_fm(NormalizeData(advA, assay = "RNA", verbose = FALSE), "RNA", "RNA")
  )
  rm(advA_prep); gc()
  
  ## 以下与原代码相同
  diag <- rbind(diag, data.frame(assay = "markers_compare1 (manuscript)", n_cells_BPD = NA, n_cells_Term = NA,
                                 DEG = sum(call_sig(transform(markers_compare1[["Adventitial FB"]],
                                                              gene = rownames(markers_compare1[["Adventitial FB"]]))) != "ns"),
                                 PHLDA1_log2FC = markers_compare1[["Adventitial FB"]]["PHLDA1", "avg_log2FC"],
                                 PHLDA1_padj   = markers_compare1[["Adventitial FB"]]["PHLDA1", "p_val_adj"]))
  oA <- read_deg(f_orig(3), "Adventitial FB")
  diag <- rbind(diag, data.frame(assay = "Shirazi Suppl. Data 1", n_cells_BPD = NA, n_cells_Term = NA,
                                 DEG = sum(call_sig(oA) != "ns"),
                                 PHLDA1_log2FC = oA["PHLDA1", "avg_log2FC"],
                                 PHLDA1_padj   = oA["PHLDA1", "p_val_adj"]))
  write.csv(diag, file.path(od, "C2_AdvFB_DEG_count_diagnosis.csv"), row.names = FALSE)
  print(diag)
}


## =============================================================================
## D. 稿件填空用的数字
## =============================================================================
packageVersion("Seurat")
alld <- readRDS(PATH_ALLDATA)
adv  <- subset(alld, subset = celltype == "Adventitial FB" & condition %in% c("BPD", "Term infant"))
x <- GetAssayData(adv, assay = "SCT", slot = "data")[c("RPS4Y1", "PHLDA1", "XIST", "TSC22D2"), ]
b <- adv$condition == "BPD"
fc_v4 <- log2(rowMeans(expm1(x[, b])) + 1) - log2(rowMeans(expm1(x[, !b])) + 1)
fc_v5 <- log2((rowSums(expm1(x[, b])) + 1) / sum(b)) - log2((rowSums(expm1(x[, !b])) + 1) / sum(!b))
round(cbind(fc_v4, fc_v5), 2)


fill <- list()
if (file.exists(PATH_ALLDATA)) {
  alld <- readRDS(PATH_ALLDATA); alld <- add_donor(alld)
  fill$total_cells_after_QC <- ncol(alld)
  write.csv(as.data.frame.matrix(table(alld$donor, alld$bulk_celltype)),
            file.path(od, "D_cells_per_donor_lineage.csv"))
  rm(alld); gc()
}
fill$mesenchymal_cells <- ncol(meso_relabel)
fill$mesenchymal_per_donor <- table(meso_relabel$donor)
fill$marker_agreement_pct <- round(agree, 1)
fill$marker_ARI <- round(ari_marker, 3)
fill$marker_agreement_MOESM7_pct <- round(agree_m7, 1)
fill$original_annotation_agreement_meso_pct <- round(agree_meso, 1)
fill$original_annotation_ARI_meso <- round(ari_meso, 3)
fill$AdvFB_pct_both_annotations <- pooled_both
fill$DEG_comparison <- C_summary %>% filter(contrast == "BPD_vs_Term")
capture.output(print(fill), file = file.path(od, "D_numbers_for_manuscript.txt"))
print(fill)
# > print(fill)
# $mesenchymal_cells
# [1] 13051
# 
# $mesenchymal_per_donor
# 
# Term infant 1          Term infant 2                  BPD 1                  BPD 2               BPD+PH 1 
# 3830                   3168                    797                    790                   2291 
# BPD+PH 2 Acute preterm injury 1 
# 1007                   1168 
# 
# $marker_agreement_pct
# [1] 81.7
# 
# $marker_ARI
# [1] 0.597
# 
# $marker_agreement_MOESM7_pct
# [1] 81.8
# 
# $DEG_comparison
# contrast       celltype DEG_ours_total DEG_orig_total genes_in_both_tables DEG_ours_tested DEG_orig_tested
# 1 BPD_vs_Term    Alveolar FB           1443             81                 6726            1411              81
# 2 BPD_vs_Term   Activated FB             11              0                 7137              11               0
# 3 BPD_vs_Term Adventitial FB           1197             51                 6826            1080              51
# 4 BPD_vs_Term Alveolar MyoFB            290             37                 7567             289              37
# 5 BPD_vs_Term   Ductal MyoFB            996             42                 6696             979              42
# 6 BPD_vs_Term       Pericyte             40             95                 7769              40              95
# 7 BPD_vs_Term           VSMC            262             56                 7404             261              56
# overlap_same_direction pct_of_orig_recovered pct_of_ours_in_orig opposite_direction hypergeom_p spearman_log2FC note
# 1                     61                  75.3                 4.3                  0    8.21e-26           0.830     
# 2                      0                   0.0                 0.0                  0    1.00e+00           0.299     
# 3                     50                  98.0                 4.6                  0    1.51e-39           0.923     
# 4                     36                  97.3                12.5                  0    3.57e-51           0.896     
# 5                     42                 100.0                 4.3                  0    3.95e-36           0.932     
# 6                     30                  31.6                75.0                  0    2.01e-51           0.787     
# 7                     46                  82.1                17.6                  0    6.67e-59           0.876  