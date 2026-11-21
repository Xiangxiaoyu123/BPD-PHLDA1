###############################################################################
# B10  验证数据集核查与重新报告   —— 回应编辑“准确报告结果”的要求
#
# 读代码后发现（r.03 第 6429-6461 行）：
#   * “Test set” 实际是 GSE235562 = 早产儿气管吸出物 MSC 体外接受 CMS（cyclic mechanical stretch，
#     周期性机械牵张）和/或高氧（Hox40/Hox80）处理 vs 对照，并不是 BPD vs 对照；
#     而且它就是 WGCNA 和机器学习的训练数据 —— 在其上算的 AUC 属于“训练集内”结果。
#   * “Validation set” 才是 GSE106910（气管吸出物 MSC；No BPD / mild / moderate / severe）。
#   * Fig. 8D 的严重度分层也来自 GSE106910。
#   稿件中 “GSE235562 = chronic-intermittent-hypoxia / shear-stress”，
#   以及 “GSE106910 的 test 与 validation 两个子集” 的描述都需要改正。
#
# 本脚本：
#   1) 输出两个数据集的样本构成表（写入 Methods 用）
#   2) 只在真正独立的 GSE106910 上重新报告 8 个基因的 BPD vs No BPD 差异与 AUC（含 95% CI）
#   3) PHLDA1 的严重度趋势检验 + 性别校正的 logistic 回归 AUC
###############################################################################
source(here::here("revision_scripts", "00_setup.R"))
library(pROC)
od <- file.path(OUT, "B10_datasets")
if (!exists("validation_expr")) load(EXPR_RDATA)

genes8 <- c("TLE5", "ACTG1", "TSC22D2", "HSPA1A", "PHLDA1", "SLC26A7", "MT1A", "NFAT5")

## ---- 1. 样本构成 ---------------------------------------------------------------
vi <- as.data.frame(validation_info)
v_samples <- data.frame(sample = colnames(validation_expr),
                        Group = vi$Group[match(colnames(validation_expr), vi$`!Sample_title`)])
comp <- list(
  GSE235562_training = as.data.frame(table(training_info$group)),
  GSE106910_independent = as.data.frame(table(v_samples$Group, useNA = "ifany")))
print(comp)
capture.output(print(comp), file = file.path(od, "1_dataset_composition.txt"))

## ---- 2. GSE106910：BPD vs No BPD（排除 N/A） ----------------------------------
v <- v_samples %>% filter(!is.na(Group), Group != "N/A")
v$BPD <- factor(ifelse(v$Group == "No BPD", "control", "BPD"), levels = c("control", "BPD"))
v$sev <- as.numeric(factor(v$Group, levels = c("No BPD", "mild", "moderate", "severe"))) - 1
X <- t(as.matrix(validation_expr[intersect(genes8, rownames(validation_expr)), v$sample]))

res <- bind_rows(lapply(colnames(X), function(g) {
  r  <- roc(v$BPD, X[, g], levels = c("control", "BPD"), direction = "<", quiet = TRUE)
  ci <- ci.auc(r, method = "delong")
  data.frame(gene = g,
             mean_control = mean(X[v$BPD == "control", g]), mean_BPD = mean(X[v$BPD == "BPD", g]),
             wilcox_p = wilcox.test(X[, g] ~ v$BPD)$p.value,
             AUC = as.numeric(auc(r)), AUC_low = ci[1], AUC_high = ci[3],
             spearman_rho_severity = cor(X[, g], v$sev, method = "spearman"),
             spearman_p_severity = cor.test(X[, g], v$sev, method = "spearman", exact = FALSE)$p.value)
}))
res$wilcox_padj <- p.adjust(res$wilcox_p, "BH")
write.csv(res, file.path(od, "2_GSE106910_8genes_BPD_vs_noBPD.csv"), row.names = FALSE)
print(res)
cat(sprintf("\nGSE106910: control n = %d, BPD n = %d\n", sum(v$BPD == "control"), sum(v$BPD == "BPD")))

## ---- 3. PHLDA1：性别校正 + 严重度趋势 ----------------------------------------
sexf <- file.path(OUT, "B5_sex", "4d_GSE106910_sample_table.csv")
if (file.exists(sexf)) {
  s <- read.csv(sexf)
  v$sex <- s$sex[match(v$sample, s$sample)]
  v$PHLDA1 <- X[, "PHLDA1"]
  fit0 <- glm(BPD ~ sex, data = v, family = binomial)
  fit1 <- glm(BPD ~ PHLDA1 + sex, data = v, family = binomial)
  r0 <- roc(v$BPD, fitted(fit0), quiet = TRUE); r1 <- roc(v$BPD, fitted(fit1), quiet = TRUE)
  capture.output(
    print(summary(fit1)$coefficients),
    cat(sprintf("\nAUC sex only = %.3f; AUC PHLDA1 + sex = %.3f; DeLong p = %.3g\n",
                auc(r0), auc(r1), roc.test(r0, r1)$p.value)),
    file = file.path(od, "3_PHLDA1_sex_adjusted_logistic.txt"))
  print(summary(fit1)$coefficients)
}
jt <- cor.test(X[, "PHLDA1"], v$sev, method = "kendall", exact = FALSE)   # 有序趋势（Kendall tau）
cat(sprintf("PHLDA1 严重度趋势：Kendall tau = %.3f, P = %.3g\n", jt$estimate, jt$p.value))

p <- data.frame(Group = factor(v$Group, levels = c("No BPD", "mild", "moderate", "severe")),
                PHLDA1 = X[, "PHLDA1"]) %>%
  ggplot(aes(Group, PHLDA1)) + geom_boxplot(outlier.shape = NA, fill = "grey90") +
  geom_jitter(width = 0.15, size = 2) +
  labs(x = NULL, y = "PHLDA1 (VST)",
       title = sprintf("GSE106910  Kendall tau = %.2f, P = %.3g", jt$estimate, jt$p.value)) + theme_bw(11)
save_both(p, file.path(od, "3_PHLDA1_by_severity_GSE106910"), 5, 4)
