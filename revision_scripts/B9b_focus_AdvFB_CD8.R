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
