###############################################################################
# B9（建议）GSE275938 中真实免疫细胞的供者水平比例 + Adventitial FB–免疫细胞通讯   —— 支持 R2-3
#   反卷积（xCell/CIBERSORT）是在 MSC bulk 数据上做的，审稿人认为证据弱。
#   这里直接在肺组织单细胞数据中：
#     1) 逐供者计算各免疫细胞类型比例，以及 CD8+ / CD4+、naive / memory / Th1 样 T 细胞比例
#     2) 每个供者 Adventitial FB 的 PHLDA1 水平与 CD8+ T 细胞比例的关系（n=4–6，仅描述）
#     3) CellChat：Adventitial FB 亚簇 + 免疫细胞（Term vs BPD）
###############################################################################
source(here::here("revision_scripts", "00_setup.R"))
od <- file.path(OUT, "B9_immune")

imm  <- add_donor(readRDS(PATH_IMM))
alld <- add_donor(readRDS(PATH_ALLDATA))
tot_cells <- table(alld$donor)

## ---- 1. 免疫细胞类型比例（占全部细胞 & 占免疫细胞） -----------------------------
imm_comp <- droplevels(imm@meta.data) %>% count(donor, celltype, .drop = FALSE) %>%
  group_by(donor) %>% mutate(pct_of_immune = 100 * n / sum(n)) %>% ungroup() %>%
  mutate(pct_of_all = 100 * n / as.numeric(tot_cells[as.character(donor)]),
         group = as.character(alld$group[match(donor, alld$donor)]))
write.csv(imm_comp, file.path(od, "1a_immune_composition_per_donor.csv"), row.names = FALSE)

p1 <- imm_comp %>% filter(group != "Acute preterm injury") %>%
  ggplot(aes(group, pct_of_all, colour = donor)) + geom_point(size = 2.5) +
  facet_wrap(~ celltype, scales = "free_y", nrow = 3) +
  labs(x = NULL, y = "% of all lung cells") + theme_bw(9) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_both(p1, file.path(od, "1b_immune_composition_dotplot"), 11, 7)

## ---- 2. T 细胞亚型（marker 门控） ---------------------------------------------
tc <- subset(imm, subset = celltype %in% c("T Cell", "NKT Cell"))
tc[["RNA"]] <- JoinLayers(tc[["RNA"]])
cnt <- GetAssayData(tc, assay = "RNA", layer = "counts")
g <- function(x) if (x %in% rownames(cnt)) cnt[x, ] > 0 else rep(FALSE, ncol(cnt))
cd8 <- (g("CD8A") | g("CD8B")) & !g("CD4")
cd4 <- g("CD4") & !(g("CD8A") | g("CD8B"))
naive <- (g("CCR7") | g("SELL")) & (g("LEF1") | g("TCF7"))
th1   <- cd4 & (g("TBX21") | g("IFNG") | g("CXCR3"))
mem   <- (g("IL7R") | g("CD44") | g("S100A4")) & !naive
tc$Tsub <- case_when(cd8 & naive ~ "CD8 naive", cd8 ~ "CD8 other",
                     th1 ~ "CD4 Th1-like", cd4 & mem ~ "CD4 memory",
                     cd4 & naive ~ "CD4 naive", cd4 ~ "CD4 other", TRUE ~ "Unassigned")
t_comp <- droplevels(tc@meta.data) %>% count(donor, Tsub, .drop = FALSE) %>%
  group_by(donor) %>% mutate(pct_of_T = 100 * n / sum(n)) %>% ungroup() %>%
  mutate(pct_of_all = 100 * n / as.numeric(tot_cells[as.character(donor)]),
         group = as.character(alld$group[match(donor, alld$donor)]))
write.csv(t_comp, file.path(od, "2a_T_subsets_per_donor.csv"), row.names = FALSE)
print(t_comp %>% filter(Tsub %in% c("CD8 naive", "CD4 Th1-like", "CD4 memory")))
# 与反卷积结果对照：稿件报告 PHLDA1 与 CD8+ naive、Th1 正相关，与 CD4+ memory 负相关

## ---- 3. 供者水平：Adventitial FB PHLDA1 vs CD8 比例 ----------------------------
pbA <- pseudobulk(meso_adventitial_FB, "donor")
phl <- log2(pbA["PHLDA1", ] / colSums(pbA) * 1e6 + 1)
dd <- t_comp %>% filter(Tsub == "CD8 naive") %>% select(donor, group, pct_CD8naive = pct_of_all) %>%
  mutate(PHLDA1_advFB_log2CPM = phl[as.character(donor)]) %>% filter(!is.na(PHLDA1_advFB_log2CPM))
write.csv(dd, file.path(od, "3_PHLDA1_vs_CD8naive_per_donor.csv"), row.names = FALSE)
p3 <- ggplot(dd, aes(PHLDA1_advFB_log2CPM, pct_CD8naive, colour = group, label = donor)) +
  geom_point(size = 3) + geom_text(vjust = -0.8, size = 3, show.legend = FALSE) +
  labs(title = "Donor level (n = 6): descriptive only") + theme_bw(11)
save_both(p3, file.path(od, "3_PHLDA1_vs_CD8naive_per_donor"), 6, 4.5)

## ---- 4. CellChat：Adventitial FB 亚簇 ↔ 免疫细胞 -------------------------------
library(CellChat)
adv <- meso_adventitial_FB
adv$cc_label <- paste0("AdvFB_C", adv$subcluster)
imm$cc_label <- as.character(imm$celltype)
DefaultAssay(adv) <- "RNA"; DefaultAssay(imm) <- "RNA"

## 1. 预处理阶段用单线程
future::plan("sequential")

cc_obj <- merge(DietSeurat(adv, assays = "RNA"), DietSeurat(imm, assays = "RNA"))
cc_obj <- NormalizeData(cc_obj, verbose = FALSE)
cc_obj[["RNA"]] <- JoinLayers(cc_obj[["RNA"]])
Layers(cc_obj[["RNA"]])

## 2. 跑 CellChat 前再开并行，并放宽数据大小上限
options(future.globals.maxSize = 8 * 1024^3)   # 8 GB
future::plan("multisession", workers = 8)
run_cc <- function(cond) {
  o <- subset(cc_obj, subset = condition == cond & cc_label != "AdvFB_C5")
  keep <- names(which(table(o$cc_label) >= 10)); o <- subset(o, subset = cc_label %in% keep)
  cc <- createCellChat(o, group.by = "cc_label", assay = "RNA")
  cc@DB <- CellChatDB.human
  cc <- subsetData(cc) %>% identifyOverExpressedGenes() %>% identifyOverExpressedInteractions()
  cc <- computeCommunProb(cc, type = "truncatedMean", trim = 0.1) %>% filterCommunication(min.cells = 10)
  cc <- computeCommunProbPathway(cc) %>% aggregateNet()
  saveRDS(cc, file.path(od, paste0("cellchat_advFB_immune_", gsub(" ", "_", cond), ".rds")))
  cc
}
## 3. 跑完恢复单线程，避免影响后面的 Seurat 步骤
cc_term <- run_cc("Term infant"); cc_bpd <- run_cc("BPD")

future::plan("sequential")


grp <- intersect(levels(cc_term@idents), levels(cc_bpd@idents))

str(lapply(cc_term@net,  function(x) dim(x)))
str(lapply(cc_term@netP, function(x) if (is.list(x)) "list" else dim(x)))

sub_cc <- function(cc, groups) {
  groups <- intersect(groups, levels(cc@idents))
  idx <- match(groups, levels(cc@idents))
  
  ## 1. 通讯网络：三维数组取前两维，二维矩阵取行列
  cc@net$prob  <- cc@net$prob [idx, idx, , drop = FALSE]
  cc@net$pval  <- cc@net$pval [idx, idx, , drop = FALSE]
  cc@net$count  <- cc@net$count [idx, idx, drop = FALSE]
  cc@net$weight <- cc@net$weight[idx, idx, drop = FALSE]
  
  ## 2. 细胞层面：表达矩阵、meta、idents 同步取子集
  keep <- cc@idents %in% groups
  if (ncol(cc@data.signaling) == length(keep)) cc@data.signaling <- cc@data.signaling[, keep, drop = FALSE]
  if (length(cc@data) && ncol(cc@data) == length(keep)) cc@data <- cc@data[, keep, drop = FALSE]
  if (nrow(cc@meta) == length(keep)) cc@meta <- cc@meta[keep, , drop = FALSE]
  cc@idents <- factor(as.character(cc@idents[keep]), levels = groups)
  
  ## 3. 重新计算通路层面结果、网络汇总和中心性
  cc@netP <- list()
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  cc <- netAnalysis_computeCentrality(cc, slot.name = "netP")
  cc
}

cc_term <- sub_cc(cc_term, grp)
cc_bpd  <- sub_cc(cc_bpd,  grp)

## 检查
levels(cc_term@idents); dim(cc_term@net$prob); length(cc_term@netP$pathways)
levels(cc_bpd@idents);  dim(cc_bpd@net$prob);  length(cc_bpd@netP$pathways)

cc_list <- list(Term = cc_term, BPD = cc_bpd)
cc_m <- mergeCellChat(cc_list, add.names = names(cc_list))

## ---------- 1. 通讯数量 / 强度柱状图 ----------
pdf(file.path(od, "4a_cellchat_compare_bar.pdf"), width = 8, height = 4)
print(compareInteractions(cc_m, show.legend = FALSE, group = c(1, 2)) +
        compareInteractions(cc_m, show.legend = FALSE, group = c(1, 2), measure = "weight"))
dev.off()

## ---------- 2. 差异通讯强度热图 ----------
pdf(file.path(od, "4b_cellchat_diff_heatmap.pdf"), width = 8, height = 7)
ComplexHeatmap::draw(netVisual_heatmap(cc_m, measure = "weight"))
dev.off()

## ---------- 3. 气泡图：按行列数自动设定页面大小 ----------
src <- grep("^AdvFB", grp); tgt <- which(!grepl("^AdvFB", grp))

save_bubble <- function(p, file, base_w = 3, base_h = 2) {
  d <- p$data
  n_col <- length(unique(d$source.target))
  n_row <- length(unique(d$interaction_name_2))
  w <- base_w + 0.22 * n_col; h <- base_h + 0.15 * n_row
  p <- p + theme(axis.text.y = element_text(size = 7),
                 axis.text.x = element_text(size = 7))
  pdf(file, width = min(w, 40), height = min(h, 60))
  print(p); dev.off()
  message(basename(file), ": ", n_row, " L-R pairs x ", n_col, " columns -> ",
          round(w, 1), " x ", round(h, 1), " in")
}

## 3a. 每个 AdvFB 亚群一页（AdvFB_Cx -> 各种免疫细胞）
for (s in src) {
  p <- netVisual_bubble(cc_m, sources.use = s, targets.use = tgt,
                        comparison = c(1, 2), angle.x = 45, remove.isolate = TRUE)
  save_bubble(p, file.path(od, paste0("4c_bubble_", grp[s], "_to_immune.pdf")))
}

## 3b. 只画变化的通讯：BPD 中上调 / 下调（所有 AdvFB 亚群 -> 免疫细胞）
p_up <- netVisual_bubble(cc_m, sources.use = src, targets.use = tgt, comparison = c(1, 2),
                         max.dataset = 2, title.name = "Increased signaling in BPD",
                         angle.x = 45, remove.isolate = TRUE)
p_dn <- netVisual_bubble(cc_m, sources.use = src, targets.use = tgt, comparison = c(1, 2),
                         max.dataset = 1, title.name = "Decreased signaling in BPD",
                         angle.x = 45, remove.isolate = TRUE)
save_bubble(p_up, file.path(od, "4d_bubble_AdvFB_to_immune_up_in_BPD.pdf"))
save_bubble(p_dn, file.path(od, "4d_bubble_AdvFB_to_immune_down_in_BPD.pdf"))


###############################################################################
# B9b  聚焦 Adventitial FB ↔ CD8+ T：BPD vs Term 差异最大的通路（不重跑 CellChat）
#
#   直接读取 B9 已保存的结果：
#     cellchat_advFB_immune_Term_infant.rds / cellchat_advFB_immune_BPD.rds
#   其中 T 细胞是 "T Cell" 整体（当时未拆 CD8），因此分两步：
#     A. 在 CellChat 结果中只取 AdvFB_C0–C4 ↔ T Cell，按 BPD − Term 排序，取前 N 条通路 / L-R 对
#     B. 用基因表达确认这些 L-R 在 T 细胞一侧的基因是否由 CD8+ T 表达（CD8 vs CD4，Term vs BPD）
#   全程只做矩阵取值和画图，几分钟内完成
###############################################################################
source(here::here("revision_scripts", "00_setup.R"))
library(CellChat)
od  <- file.path(OUT, "B9_immune")
od5 <- file.path(od, "5_focus_AdvFB_CD8"); dir.create(od5, showWarnings = FALSE)

TOP_N    <- 10                       # 每个方向展示前 N 条
T_GROUPS <- c("T Cell")              # 如需把 NKT 一起算，改为 c("T Cell", "NKT Cell")
cols_cond <- c(Term = "#9E9E9E", BPD = "#00768E")

## ---- 0. 读入已保存的 CellChat 对象 ----------------------------------------------
cc_term <- readRDS(file.path(od, "cellchat_advFB_immune_Term_infant.rds"))
cc_bpd  <- readRDS(file.path(od, "cellchat_advFB_immune_BPD.rds"))
cat("Term 中的细胞群：", paste(levels(cc_term@idents), collapse = ", "), "\n")
cat("BPD  中的细胞群：", paste(levels(cc_bpd@idents),  collapse = ", "), "\n")
ADV <- intersect(grep("^AdvFB_C", levels(cc_term@idents), value = TRUE),
                 grep("^AdvFB_C", levels(cc_bpd@idents),  value = TRUE))
T_GROUPS <- intersect(T_GROUPS, intersect(levels(cc_term@idents), levels(cc_bpd@idents)))
stopifnot(length(ADV) > 0, length(T_GROUPS) > 0)
cat("纳入的 AdvFB 亚群：", paste(ADV, collapse = ", "), "；T 细胞群：", paste(T_GROUPS, collapse = ", "), "\n")

## 从 3 维数组 [source, target, pathway/LR] 中取指定细胞群并求和
sum_prob <- function(P, src, tgt) {
  s <- intersect(src, dimnames(P)[[1]]); t <- intersect(tgt, dimnames(P)[[2]])
  if (!length(s) || !length(t)) return(setNames(numeric(0), character(0)))
  apply(P[s, t, , drop = FALSE], 3, sum)
}
min_p <- function(P, src, tgt) {
  s <- intersect(src, dimnames(P)[[1]]); t <- intersect(tgt, dimnames(P)[[2]])
  apply(P[s, t, , drop = FALSE], 3, min)
}
dirs <- list("AdvFB -> T" = list(src = ADV, tgt = T_GROUPS),
             "T -> AdvFB" = list(src = T_GROUPS, tgt = ADV))

## ---- A1. 通路层面 ---------------------------------------------------------------
pw_tab <- bind_rows(lapply(names(dirs), function(dn) {
  d <- dirs[[dn]]
  a <- sum_prob(cc_term@netP$prob, d$src, d$tgt); b <- sum_prob(cc_bpd@netP$prob, d$src, d$tgt)
  pw <- union(names(a), names(b))
  data.frame(direction = dn, pathway = pw, Term = unname(a[pw]), BPD = unname(b[pw]))
})) %>%
  mutate(across(c(Term, BPD), ~ replace(.x, is.na(.x), 0))) %>%
  filter(Term > 0 | BPD > 0) %>%
  mutate(delta = BPD - Term,
         change = case_when(Term == 0 ~ "BPD only", BPD == 0 ~ "Term only",
                            delta > 0 ~ "Up in BPD", TRUE ~ "Down in BPD")) %>%
  arrange(direction, desc(abs(delta)))
write.csv(pw_tab, file.path(od5, "A1_pathways_AdvFB_T_BPD_vs_Term_all.csv"), row.names = FALSE)
top_pw <- pw_tab %>% group_by(direction) %>% slice_max(abs(delta), n = TOP_N, with_ties = FALSE) %>% ungroup()
cat("\n== 差异最大的通路 ==\n"); print(top_pw, n = Inf)
write.csv(top_pw, file.path(od5, "A1_top_pathways.csv"), row.names = FALSE)

p_pw <- top_pw %>%
  pivot_longer(c(Term, BPD), names_to = "condition", values_to = "prob") %>%
  ggplot(aes(prob, reorder(pathway, delta))) +
  geom_line(aes(group = pathway), colour = "grey75", linewidth = 0.6) +
  geom_point(aes(colour = condition), size = 3) +
  facet_wrap(~ direction, scales = "free_y") +
  scale_colour_manual(values = cols_cond, name = NULL) +
  labs(x = "Communication probability (summed over L–R pairs and AdvFB subclusters)", y = NULL,
       title = "Adventitial FB ↔ T cells: pathways most changed in BPD vs Term") +
  theme_bw(10) + theme(panel.grid.minor = element_blank(), legend.position = "top")
save_both(p_pw, file.path(od5, "A1_top_pathways_dumbbell"), 10, 2.5 + 0.3 * TOP_N)

# 前 N 条通路在各 AdvFB 亚群中的变化（BPD − Term），看是哪个亚群驱动
sub_delta <- bind_rows(lapply(names(dirs), function(dn) {
  bind_rows(lapply(ADV, function(a) {
    d <- dirs[[dn]]
    src <- if (dn == "AdvFB -> T") a else d$src; tgt <- if (dn == "AdvFB -> T") d$tgt else a
    x <- sum_prob(cc_term@netP$prob, src, tgt); y <- sum_prob(cc_bpd@netP$prob, src, tgt)
    pw <- top_pw$pathway[top_pw$direction == dn]
    data.frame(direction = dn, subcluster = a, pathway = pw,
               delta = unname(ifelse(is.na(y[pw]), 0, y[pw]) - ifelse(is.na(x[pw]), 0, x[pw])))
  }))
}))
lim <- max(abs(sub_delta$delta))
p_sub <- ggplot(sub_delta, aes(subcluster, pathway, fill = delta)) + geom_tile(colour = "white") +
  facet_wrap(~ direction, scales = "free_y") +
  scale_fill_gradient2(low = "#9E9E9E", mid = "white", high = "#00768E", limits = c(-lim, lim),
                       name = "BPD − Term") +
  labs(x = NULL, y = NULL, title = "Contribution of each adventitial FB subcluster") +
  theme_minimal(10) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_both(p_sub, file.path(od5, "A1b_top_pathways_by_subcluster"), 9, 2 + 0.3 * TOP_N)

## ---- A2. 配体–受体对层面（至少一组 P < 0.05） ------------------------------------
db <- CellChatDB.human$interaction[, c("interaction_name", "interaction_name_2", "pathway_name", "ligand", "receptor")]
lr_tab <- bind_rows(lapply(names(dirs), function(dn) {
  d <- dirs[[dn]]
  a <- sum_prob(cc_term@net$prob, d$src, d$tgt); b <- sum_prob(cc_bpd@net$prob, d$src, d$tgt)
  pa <- min_p(cc_term@net$pval, d$src, d$tgt);  pb <- min_p(cc_bpd@net$pval, d$src, d$tgt)
  lr <- union(names(a), names(b))
  data.frame(direction = dn, interaction_name = lr,
             Term = unname(a[lr]), BPD = unname(b[lr]), p_Term = unname(pa[lr]), p_BPD = unname(pb[lr]))
})) %>%
  mutate(across(c(Term, BPD), ~ replace(.x, is.na(.x), 0)),
         across(c(p_Term, p_BPD), ~ replace(.x, is.na(.x), 1))) %>%
  filter((Term > 0 & p_Term < 0.05) | (BPD > 0 & p_BPD < 0.05)) %>%
  mutate(delta = BPD - Term) %>%
  left_join(db, by = "interaction_name") %>%
  arrange(direction, desc(abs(delta)))
write.csv(lr_tab, file.path(od5, "A2_LR_AdvFB_T_BPD_vs_Term_all.csv"), row.names = FALSE)
top_lr <- lr_tab %>% group_by(direction) %>% slice_max(abs(delta), n = TOP_N, with_ties = FALSE) %>% ungroup()
write.csv(top_lr, file.path(od5, "A2_top_LR.csv"), row.names = FALSE)
print(top_lr %>% select(direction, interaction_name_2, pathway_name, Term, BPD, delta), n = Inf)

p_lr <- top_lr %>%
  pivot_longer(c(Term, BPD), names_to = "condition", values_to = "prob") %>%
  mutate(condition = factor(condition, levels = c("Term", "BPD")),
         label = reorder(paste0(interaction_name_2, "  [", pathway_name, "]"), delta)) %>%
  ggplot(aes(condition, label)) +
  geom_point(aes(size = prob, colour = prob)) +
  facet_wrap(~ direction, scales = "free_y") +
  scale_colour_gradient(low = "grey85", high = "#00768E", name = "Prob.") +
  scale_size_continuous(range = c(0.5, 6), name = "Prob.") +
  labs(x = NULL, y = NULL, title = sprintf("Top %d L–R pairs changed in BPD (AdvFB ↔ T cells)", TOP_N)) +
  theme_bw(9) + theme(panel.grid.minor = element_blank())
save_both(p_lr, file.path(od5, "A2_top_LR_bubble"), 10, 2.5 + 0.35 * TOP_N)

## ---- B. CD8 特异性：T 细胞一侧的基因在 CD8+ vs CD4+ T 中的表达 ----------------------
imm <- add_donor(readRDS(PATH_IMM))
tc  <- subset(imm, subset = celltype %in% c("T Cell", "NKT Cell"))
if (inherits(tc[["RNA"]], "Assay5")) tc[["RNA"]] <- JoinLayers(tc[["RNA"]])
tc  <- NormalizeData(tc, verbose = FALSE)
cnt <- GetAssayData(tc, assay = "RNA", layer = "counts")
pos <- function(x) if (x %in% rownames(cnt)) cnt[x, ] > 0 else rep(FALSE, ncol(cnt))
tc$Tclass <- ifelse((pos("CD8A") | pos("CD8B")) & !pos("CD4"), "CD8 T",
                    ifelse(pos("CD4") & !(pos("CD8A") | pos("CD8B")), "CD4 T", "Other T/NKT"))
cat("\n== T 细胞分类 × 分组 ==\n"); print(table(tc$Tclass, tc$condition))

# 每对 L-R 中属于 T 细胞一侧的基因：AdvFB -> T 取受体；T -> AdvFB 取配体（复合物按 "_" 拆开）
split_genes <- function(x) unique(unlist(strsplit(x, "_")))
t_side <- bind_rows(lapply(seq_len(nrow(top_lr)), function(i) {
  g <- if (top_lr$direction[i] == "AdvFB -> T") top_lr$receptor[i] else top_lr$ligand[i]
  data.frame(pair = top_lr$interaction_name_2[i], direction = top_lr$direction[i],
             role = ifelse(top_lr$direction[i] == "AdvFB -> T", "receptor", "ligand"), gene = split_genes(g))
})) %>% filter(gene %in% rownames(tc)) %>% distinct()

md  <- tc@meta.data %>% filter(condition %in% c("Term infant", "BPD"), Tclass != "Other T/NKT")
ex  <- FetchData(tc, vars = unique(t_side$gene), cells = rownames(md))
dot <- cbind(md[, c("Tclass", "condition")], ex) %>%
  pivot_longer(-c(Tclass, condition), names_to = "gene", values_to = "expr") %>%
  group_by(gene, Tclass, condition) %>%
  summarise(pct = 100 * mean(expr > 0), avg = mean(expm1(expr)), n_cells = n(), .groups = "drop") %>%
  group_by(gene) %>% mutate(avg_scaled = as.numeric(scale(log1p(avg)))) %>% ungroup() %>%
  left_join(t_side %>% select(gene, pair, direction) %>% group_by(gene) %>%
              summarise(pairs = paste(unique(pair), collapse = "; "), direction = first(direction)), by = "gene") %>%
  mutate(column = factor(paste(Tclass, ifelse(condition == "BPD", "BPD", "Term"), sep = "\n"),
                         levels = c("CD8 T\nTerm", "CD8 T\nBPD", "CD4 T\nTerm", "CD4 T\nBPD")))
write.csv(dot, file.path(od5, "B_Tside_genes_CD8_vs_CD4.csv"), row.names = FALSE)

p_cd8 <- ggplot(dot, aes(column, gene)) +
  geom_point(aes(size = pct, colour = avg_scaled)) +
  facet_grid(direction ~ ., scales = "free_y", space = "free_y") +
  scale_colour_gradient2(low = "#3B6E8F", mid = "grey90", high = "#B2182B", name = "Scaled\nexpression") +
  scale_size_continuous(range = c(0.3, 6), name = "% expressing") +
  labs(x = NULL, y = NULL, title = "T-cell-side genes of the top L–R pairs: CD8⁺ vs CD4⁺ T cells") +
  theme_bw(9) + theme(axis.text.y = element_text(face = "italic"), panel.grid.minor = element_blank())
save_both(p_cd8, file.path(od5, "B_Tside_genes_CD8_vs_CD4_dotplot"), 7, 2 + 0.28 * length(unique(dot$gene)))

# 每对 L-R 是否“CD8 主导”：T 侧基因在 CD8 中的表达比例 ≥ CD4（Term 与 BPD 合并）
cd8_flag <- dot %>% group_by(gene, Tclass) %>% summarise(pct = weighted.mean(pct, n_cells), .groups = "drop") %>%
  pivot_wider(names_from = Tclass, values_from = pct) %>%
  mutate(CD8_enriched = `CD8 T` >= `CD4 T`)
top_lr_cd8 <- t_side %>% left_join(cd8_flag, by = "gene") %>%
  group_by(direction, pair) %>% summarise(T_side_genes = paste(gene, collapse = "/"),
                                          pct_CD8 = round(mean(`CD8 T`), 1), pct_CD4 = round(mean(`CD4 T`), 1),
                                          CD8_enriched = all(CD8_enriched), .groups = "drop") %>%
  left_join(top_lr %>% select(direction, pair = interaction_name_2, pathway_name, Term, BPD, delta),
            by = c("direction", "pair")) %>% arrange(direction, desc(abs(delta)))
write.csv(top_lr_cd8, file.path(od5, "B_top_LR_with_CD8_specificity.csv"), row.names = FALSE)
cat("\n== 差异最大的 L-R 对及其 CD8 特异性 ==\n"); print(top_lr_cd8, n = Inf)

## ---- C. AdvFB 一侧的基因（Term vs BPD），确认变化由成纤维细胞表达驱动 ----------------
adv <- ensure_rna_norm(meso_adventitial_FB)
if (inherits(adv[["RNA"]], "Assay5")) adv[["RNA"]] <- JoinLayers(adv[["RNA"]])
fb_side <- bind_rows(lapply(seq_len(nrow(top_lr)), function(i) {
  g <- if (top_lr$direction[i] == "AdvFB -> T") top_lr$ligand[i] else top_lr$receptor[i]
  data.frame(direction = top_lr$direction[i], gene = split_genes(g))
})) %>% filter(gene %in% rownames(adv)) %>% distinct()
mdA <- adv@meta.data %>% filter(condition %in% c("Term infant", "BPD"), subcluster != "5")
exA <- FetchData(adv, vars = unique(fb_side$gene), cells = rownames(mdA))
dotA <- cbind(mdA[, "condition", drop = FALSE], exA) %>%
  pivot_longer(-condition, names_to = "gene", values_to = "expr") %>%
  group_by(gene, condition) %>% summarise(pct = 100 * mean(expr > 0), avg = mean(expm1(expr)), .groups = "drop") %>%
  group_by(gene) %>% mutate(avg_scaled = as.numeric(scale(log1p(avg)))) %>% ungroup() %>%
  mutate(condition = factor(ifelse(condition == "BPD", "BPD", "Term"), levels = c("Term", "BPD")))
p_fb <- ggplot(dotA, aes(condition, gene)) + geom_point(aes(size = pct, colour = avg_scaled)) +
  scale_colour_gradient2(low = "#3B6E8F", mid = "grey90", high = "#B2182B", name = "Scaled\nexpression") +
  scale_size_continuous(range = c(0.3, 6), name = "% expressing") +
  labs(x = NULL, y = NULL, title = "Adventitial-FB-side genes") +
  theme_bw(9) + theme(axis.text.y = element_text(face = "italic"))
save_both(p_fb, file.path(od5, "C_AdvFB_side_genes_dotplot"), 4.5, 2 + 0.28 * length(unique(dotA$gene)))
message("完成，结果在：", od5)


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

DefaultAssay(tc) <- "RNA"                                    # 明确切换到 RNA
if (inherits(tc[["RNA"]], "Assay5")) tc[["RNA"]] <- JoinLayers(tc[["RNA"]])
tc  <- NormalizeData(tc, assay = "RNA", verbose = FALSE)     # 明确指定 assay
Layers(tc[["RNA"]])                                          # 应包含 counts 和 data

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
print(sapply(modules, length))                               # 每个模块实际可用的基因数
modules <- modules[sapply(modules, length) >= 2]             # 去掉基因太少的模块，避免报错

## 细胞较少时缩小对照基因数，避免 AddModuleScore 的分箱出错
n_ctrl <- min(100, max(10, floor(ncol(cd8) / 10)))
cd8 <- AddModuleScore(cd8, features = modules, name = "MOD_", assay = "RNA",
                      ctrl = n_ctrl, seed = 2026)

## 把 MOD_1、MOD_2… 改成模块名，方便后面使用
old <- paste0("MOD_", seq_along(modules))
colnames(cd8@meta.data)[match(old, colnames(cd8@meta.data))] <- paste0("MOD_", names(modules))
head(cd8@meta.data[, paste0("MOD_", names(modules))])

sc <- cd8@meta.data[, paste0("MOD_", names(modules))]; colnames(sc) <- names(modules)

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

## 先确认这个标签在哪一列
unique(donor_mat$subset)
unique(donor_mat$donor)

plot_df <- subset(donor_mat, subset != "PHLDA1 in adventitial FB")
## 或者用 dplyr：plot_df <- dplyr::filter(donor_mat, subset != "PHLDA1 in adventitial FB")
nrow(donor_mat); nrow(plot_df)    # 应该正好少了一行（每个供者一个格子）

p_c3 <- ggplot(plot_df, aes(donor, subset, fill = z)) + geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.2g", value)), size = 2.6) +
  scale_fill_gradient2(low = "#3B6E8F", mid = "white", high = "#B2182B", name = "Row z-score") +
  geom_vline(xintercept = c(2.5, 4.5), colour = "grey30") +
  labs(x = NULL, y = NULL, title = "Donor-level overview (numbers = raw values)") +
  theme_minimal(9) + theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())

# p_c3 <- ggplot(donor_mat[-which(donor_mat$donor == "PHLDA1 in adventitial FB"),], aes(donor, subset, fill = z)) + geom_tile(colour = "white") +
#   geom_text(aes(label = sprintf("%.2g", value)), size = 2.6) +
#   scale_fill_gradient2(low = "#3B6E8F", mid = "white", high = "#B2182B", name = "Row z-score") +
#   geom_vline(xintercept = c(2.5, 4.5), colour = "grey30") +
#   labs(x = NULL, y = NULL, title = "Donor-level overview (numbers = raw values)") +
#   theme_minimal(9) + theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())
save_both(p_c3, file.path(od6, "C3_donor_level_heatmap"), 7.5, 4.2)

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
plot_df <- subset(dm, subset != "PHLDA1 in adventitial FB")
p <- ggplot(plot_df, aes(donor, subset, fill = z)) + geom_tile(colour = "white", linewidth = 0.3) +
  geom_text(aes(label = ifelse(abs(value) >= 10, sprintf("%.0f", value), sprintf("%.2g", value))), size = FS / .pt * 0.9) +
  scale_fill_gradient2(low = "#3B6E8F", mid = "white", high = "#B2182B", name = "Row\nz-score") +
  geom_vline(xintercept = c(2.5, 4.5), linewidth = LW, colour = "black") +
  scale_x_discrete(position = "top", expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
  labs(x = NULL, y = NULL) + theme_cns_grid() +
  theme(axis.text.x = element_text(angle = 45, hjust = 0), axis.ticks = element_blank())
save_cns(p, "C3_donor_level_heatmap", 95, 60)

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

