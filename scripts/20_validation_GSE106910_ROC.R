###############################################################################
# Expression of the candidate genes in the test and validation sets and ROC analysis
#
# Source      : r.03_Cell_prop.r (lines 6299-6851)
# Input       : GSE235562, GSE106910
# Output      : Fig. 6A-D; Supplementary Fig. 9
#
# Code below is unchanged from the original analysis script.
###############################################################################
# -----------------------------------------------------------------------------
# PATHS: this file keeps the absolute paths of the analysis server on which the
# study was run (/data/nas1/xiangxiaoyu_OD/1.BPD+MSCs/1.scRNA/... and
# /home/xxyloveme/study/BPD+MSCs/...). Set PROJ_DIR in config/paths.R and
# replace those prefixes with your own project directory before running.
# -----------------------------------------------------------------------------


###############3.9 bulkRNA基因的验证情况###################
setwd("/data/nas1/xiangxiaoyu_OD/1.BPD+MSCs/1.scRNA/06_Expression_validation/")
rm(list = ls())

##########验证集2.GSE106910#######
# 加载必要的包
library(data.table)
library(dplyr)
library(stringr)
library(purrr)
library(readxl)

expr_GSE106910 <- fread("GSE106910_combined_counts_geo.txt",
                        header = TRUE,
                        sep = "\t",
                        quote = "",
                        check.names = FALSE)
colnames(expr_GSE106910)[1] <- "gene_id"
# 查看原始数据维度
cat("原始表达矩阵维度:", dim(expr_GSE106910), "\n")

# 读取GTF文件并创建映射
gtf_file <- rtracklayer::import("/data/nas1/xiangxiaoyu_OD/1.BPD+MSCs/1.scRNA/06_Expression_validation/gencode.v38.annotation.gtf")
humanGTF <- as.data.frame(gtf_file)

# 创建更准确的基因映射
# 方法1: 使用transcript_id映射
gene_mapping <- humanGTF %>% 
  dplyr::select(gene_name, gene_id) %>%
  distinct() %>%
  na.omit()
gene_mapping$gene_id <- sub("\\..*", "", gene_mapping$gene_id)

# 查看映射关系
cat("GTF中的transcript数量:", nrow(gene_mapping), "\n")
cat("表达矩阵中的transcript数量:", length(unique(expr_GSE106910$id)), "\n")

# 合并表达矩阵和基因名称
exprSet_GSE106910 <- inner_join(expr_GSE106910, gene_mapping, by = "gene_id")
cat("合并后的数据维度:", dim(exprSet_GSE106910), "\n")

# 1. 提取样本列（除了gene_id和gene_name）
sample_cols <- setdiff(colnames(exprSet_GSE106910), c("gene_id", "gene_name"))

# 2. 处理重复基因名 - 对相同基因的count值求和
exprSet_aggregated <- exprSet_GSE106910 %>%
  group_by(gene_name) %>%
  summarise(across(all_of(sample_cols), sum, na.rm = TRUE)) %>%
  ungroup()

# 3. 设置基因名为行名
exprSet_final <- as.data.frame(exprSet_aggregated)
rownames(exprSet_final) <- exprSet_final$gene_name
exprSet_final$gene_name <- NULL

#5. Count值转换选项
# 方法1: Counts Per Million (CPM) 标准化
cpm_matrix <- apply(exprSet_final, 2, function(x) {
  (x / sum(x)) * 1e6
})

# 方法2: log2(CPM + 1) 转换（更稳定）
log2cpm_matrix <- log2(cpm_matrix + 1)

# 方法3: 方差稳定转换 (VST) - 需要DESeq2
if (!require(DESeq2)) {
  BiocManager::install("DESeq2")
  library(DESeq2)
}

# 创建DESeqDataSet对象进行VST转换
tryCatch({
  # 创建样本信息（这里需要根据实际情况调整）
  sample_info <- data.frame(
    sample = colnames(exprSet_final),
    group = sapply(strsplit(colnames(exprSet_final), " "), function(x) x[1]))
  rownames(sample_info) <- sample_info$sample
  
  # 创建DESeqDataSet
  dds <- DESeqDataSetFromMatrix(
    countData = round(exprSet_final),  # 确保是整数
    colData = sample_info,
    design = ~ group)
  
  # 进行VST转换
  vsd <- vst(dds, blind = FALSE)
  vst_matrix <- assay(vsd)
  
  cat("VST转换完成\n")
}, error = function(e) {
  cat("VST转换失败:", e$message, "\n")
  vst_matrix <- NULL
})

# 保存原始count矩阵
write.csv(exprSet_final, "GSE106910_raw_counts_matrix.csv")
# 保存CPM矩阵
write.csv(cpm_matrix, "GSE106910_CPM_matrix.csv")
# 保存log2(CPM+1)矩阵
write.csv(log2cpm_matrix, "GSE106910_log2CPM_matrix.csv")
# 保存VST矩阵（如果成功）
write.csv(vst_matrix, "GSE106910_VST_matrix.csv")

##判断是否需要log2转化
qx <- as.numeric(quantile(vst_matrix, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) {
  vst_matrix[which(vst_matrix <= 0)] <- 0
  vst_matrix <- log2(vst_matrix + 1)
  print("log2 transform finished")
}else{
  print("log2 transform not needed")
}

pdf(file = "exprSet_GSE106910.pdf",width = 6,height = 3)
boxplot(vst_matrix,outline=FALSE,las=2,col = 'blue',xaxt = 'n',ann = F)
dev.off()

#.表达验证 -----------------------------------------------------------------
# 加载必要的包
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)
library(tidyr)
detach("package:DESeq2", unload=TRUE)
library(gghalves)
#导入数据
training_expr <- read.csv(file = "./exprSet_GSE235562.csv",row.names = 1)
training_info <- read.csv(file = "./group_GSE235562.csv",row.names = 1)

validation_expr <- as.data.frame(vst_matrix)
validation_info <- fread("GSE106910_pData.txt",
                        header = TRUE,
                        sep = "\t",
                        quote = "",
                        check.names = FALSE)
feature_genes <- read.csv(file = "common_genes.csv")
feature_genes <- feature_genes$geneNames
feature_genes_2 <- c("ACTG1","TLE5","HSPA1A","MT1A","NFAT5","PHLDA1","SLC26A7","TSC22D2")
## 步骤1：准备数据
# 检查基因是否在数据集中
available_genes_train <- feature_genes[feature_genes %in% rownames(training_expr)]
available_genes_valid <- feature_genes_2[feature_genes_2 %in% rownames(validation_expr)]

cat("训练集中可用的特征基因:", length(available_genes_train), "\n")
#训练集中可用的特征基因: 8 
cat("验证集中可用的特征基因:", length(available_genes_valid), "\n")
#验证集中可用的特征基因: 8

# 提取特征基因表达数据
train_data <- as.data.frame(t(training_expr[available_genes_train, ]))
valid_data <- as.data.frame(t(validation_expr[available_genes_valid, ]))

# 添加分组信息
train_data$Group <- training_info$group[match(rownames(train_data), rownames(training_info))]
train_data$Group <- ifelse(train_data$Group == "control", "control", "treatment")

valid_data$Group <- validation_info$Group[match(rownames(valid_data), validation_info$`!Sample_title`)]
valid_data <- valid_data[-which(valid_data$Group == "N/A"),]
valid_data$Group_re <- ifelse(valid_data$Group == "No BPD", "control", "BPD")

# 转换为长格式用于绘图
train_long <- train_data %>%
  pivot_longer(cols = -Group, names_to = "Gene", values_to = "Expression") %>%
  mutate(Dataset = "Training")

valid_long <- valid_data %>%
  pivot_longer(cols = -c(Group,Group_re), names_to = "Gene", values_to = "Expression") %>%
  mutate(Dataset = "Validation")
levels(valid_long$Group) <- c("No BPD","mild","moderate","severe")
valid_long$Group <- factor(valid_long$Group,levels = levels(valid_long$Group))
valid_long$Group_re <- factor(valid_long$Group_re,levels = c("control", "BPD"))
valid_long$Gene <- factor(valid_long$Gene, levels = feature_genes_2)


##绘制测试集的数据
# 对每个基因进行统计学检验
stat_test <- train_long %>%
  group_by(Gene) %>%
  rstatix::wilcox_test(Expression ~ Group) %>%
  rstatix::adjust_pvalue(method = "BH") %>%
  rstatix::add_xy_position(x = "Group", dodge = 0.8)
# 保存统计检验结果
library(writexl)
write_xlsx(as.data.frame(stat_test), "statistical_test_results.xlsx")
p_train_genes <- ggplot(data = train_long,
                        aes(x=1, y=Expression, fill = Group, color = Group))+
  geom_half_violin(position=position_nudge(x=-0.08),
                   side='l',adjust=1.2,trim=F, alpha = 0.4)+
  geom_point(aes(x = 1.15, y = Expression, color = Group),
             position = position_jitter(width = 0.05),
             size = 1)+
  geom_boxplot(outlier.shape = NA,width = 0.15,color="black")+
  facet_wrap(~ Gene, scales = "free_y",ncol = 8)+
  stat_compare_means(method = "wilcox.test", label = "p.signif")+
  scale_fill_manual(values = c("#539382","#c83e2b"))+
  scale_color_manual(values = c("#539382","#c83e2b"))+
  labs(x=NULL,y="Expression Level",fill=NULL, color=NULL)+
  theme_bw()+
  theme(legend.background = element_blank(),
        axis.text = element_text(color = 'black',size=9),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size=13),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
ggsave("Feature_Genes_Expression_inTesting.pdf", p_train_genes, width = 15, height = 2.6)




##绘制验证集的数据
# 对每个基因进行统计学检验
library(rstatix)
library(ggpubr)
stat_test_validation <- valid_long %>%
  group_by(Gene) %>%
  wilcox_test(Expression ~ Group, 
              ref.group = "No BPD",  # 以"No BPD"为参考组
              p.adjust.method = "BH") %>%
  add_xy_position(x = "Group", step.increase = 0.1)

stat_test_validation_re <- valid_long %>%
  group_by(Gene) %>%
  wilcox_test(Expression ~ Group_re, 
              ref.group = "control",  # 以"No BPD"为参考组
              p.adjust.method = "BH") %>%
  add_xy_position(x = "Group", step.increase = 0.1)


# 保存统计检验结果
library(writexl)
write_xlsx(as.data.frame(stat_test_validation), "statistical_valid_results.xlsx")
write_xlsx(as.data.frame(stat_test_validation_re), "statistical_valid_results_re.xlsx")

valid_long_restart <- valid_long[which(valid_long$Gene == "PHLDA1"),]
write.csv(valid_long_restart,"valid_long_restart.csv")
write.csv(valid_long,"valid_long.csv")

valid_long_restart <- read.csv(file = "/data/nas1/xiangxiaoyu_OD/1.BPD+MSCs/1.scRNA/06_Expression_validation/valid_long.csv")
levels(valid_long_restart$Group) <- c("No BPD","mild","moderate","severe")
valid_long_restart$Group <- factor(valid_long_restart$Group,levels = levels(valid_long_restart$Group))
valid_long_restart$Group_re <- factor(valid_long_restart$Group_re,levels = c("control", "BPD"))
valid_long_restart$Gene <- factor(valid_long_restart$Gene, levels = feature_genes_2)



p_valid_genes <- ggplot(data = valid_long,
                        aes(x=1, y=Expression, fill = Group, color = Group))+
  geom_half_violin(position=position_nudge(x=-0.08),
                   side='l',adjust=1.2,trim=F, alpha = 0.4)+
  geom_point(aes(x = 1.15, y = Expression, color = Group),
             position = position_jitter(width = 0.05),
             size = 1)+
  geom_boxplot(outlier.shape = NA,width = 0.15,color="black")+
  stat_compare_means(
    method = "wilcox.test",
    ref.group = "No BPD",        # 指定参考组
    label = "p.signif",hide.ns = F,
    size = 1)+
  facet_wrap(~ Gene, scales = "free_y",ncol = 8)+
  scale_fill_manual(values = c("#5f6f8a","#e7c659","#ea9c9b","#005495"))+
  scale_color_manual(values = c("#5f6f8a","#e7c659","#ea9c9b","#005495"))+
  labs(x=NULL,y="Expression Level",fill=NULL, color=NULL)+
  theme_bw()+
  theme(legend.background = element_blank(),
        axis.text = element_text(color = 'black',size=9),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size=13),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
ggsave("Feature_Genes_Expression_inValid.pdf", p_valid_genes, width = 15, height = 2.6)

p_valid_genes_PHLDA1 <- ggplot(data = valid_long_restart,
                        aes(x=1, y=Expression, fill = Group, color = Group))+
  geom_half_violin(position=position_nudge(x=-0.08),
                   side='l',adjust=1.2,trim=F, alpha = 0.4)+
  geom_point(aes(x = 1.15, y = Expression, color = Group),
             position = position_jitter(width = 0.05),
             size = 1)+
  geom_boxplot(outlier.shape = NA,width = 0.15,color="black")+
  stat_compare_means(
    method = "wilcox.test",
    ref.group = "No BPD",        # 指定参考组
    label = "p.signif",hide.ns = F,
    size = 1)+
  facet_wrap(~ Gene, scales = "free_y",ncol = 8)+
  scale_fill_manual(values = c("#5f6f8a","#e7c659","#ea9c9b","#005495"))+
  scale_color_manual(values = c("#5f6f8a","#e7c659","#ea9c9b","#005495"))+
  labs(x=NULL,y="Expression Level",fill=NULL, color=NULL)+
  theme_bw()+
  theme(legend.background = element_blank(),
        axis.text = element_text(color = 'black',size=9),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size=13),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
ggsave("Feature_Genes_Expression_inValid_PHLDA1.pdf", p_valid_genes_PHLDA1, width = 15, height = 2.6)

c("#ebd57c","#e18748","#ed4437","#fe9e37","#b24646","#96873b","#b49d99","#b37557","#fc9a9a")


p_valid_genes_re <- ggplot(data = valid_long,
                           aes(x=1, y=Expression, fill = Group_re, color = Group_re))+
  geom_half_violin(position=position_nudge(x=-0.08),
                   side='l',adjust=1.2,trim=F, alpha = 0.4)+
  geom_point(aes(x = 1.15, y = Expression, color = Group_re),
             position = position_jitter(width = 0.05),
             size = 1)+
  geom_boxplot(outlier.shape = NA,width = 0.15,color="black")+
  facet_wrap(~ Gene, scales = "free_y",ncol = 8)+
  stat_compare_means(method = "wilcox.test", label = "p.signif")+
  scale_fill_manual(values = c("#b49d99","#b24646"))+
  scale_color_manual(values = c("#b49d99","#b24646"))+
  labs(x=NULL,y="Expression Level",fill=NULL, color=NULL)+
  theme_bw()+
  theme(legend.background = element_blank(),
        axis.text = element_text(color = 'black',size=9),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size=13),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

p_valid_genes_re_PHLDA1 <- ggplot(data = valid_long_restart,
                           aes(x=1, y=Expression, fill = Group_re, color = Group_re))+
  geom_half_violin(position=position_nudge(x=-0.08),
                   side='l',adjust=1.2,trim=F, alpha = 0.4)+
  geom_point(aes(x = 1.15, y = Expression, color = Group_re),
             position = position_jitter(width = 0.05),
             size = 1)+
  geom_boxplot(outlier.shape = NA,width = 0.15,color="black")+
  facet_wrap(~ Gene, scales = "free_y",ncol = 8)+
  stat_compare_means(method = "wilcox.test", label = "p.signif")+
  scale_fill_manual(values = c("#b49d99","#b24646"))+
  scale_color_manual(values = c("#b49d99","#b24646"))+
  labs(x=NULL,y="Expression Level",fill=NULL, color=NULL)+
  theme_bw()+
  theme(legend.background = element_blank(),
        axis.text = element_text(color = 'black',size=9),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size=13),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# 保存所有特征基因的图形
ggsave("Feature_Genes_Expression_inValid_re.pdf", p_valid_genes_re, width = 15, height = 2.6)
ggsave("Feature_Genes_Expression_inValid_re_PHLDA1.pdf", p_valid_genes_re_PHLDA1, width = 15, height = 2.6)
save.image(file = "/data/nas1/xiangxiaoyu_OD/1.BPD+MSCs/1.scRNA/Expression_validation+ROC+xCell+CIBERSORT.RData")


#.ROC验证 -----------------------------------------------------------------
setwd("/data/nas1/xiangxiaoyu_OD/1.BPD+MSCs/1.scRNA/07_ROC/")
# 加载必要的包
library(pROC)
library(ggplot2)
library(dplyr)

# 绘制测试集（训练集）的ROC曲线

# 准备测试集数据
test_data <- train_data
test_data$Group_numeric <- ifelse(test_data$Group == "control", 0, 1)

# 创建ROC数据框
roc_data_test <- data.frame(Group = test_data$Group_numeric)

# 添加每个基因的表达值
for(gene in available_genes_train) {
  roc_data_test[[gene]] <- test_data[[gene]]
}

# 顶刊配色方案
top_journal_colors <- c("#ebd57c","#e18748","#ed4437","#fe9e37","#b24646","#96873b","#b49d99","#b37557","#fc9a9a")

# 如果基因超过8个，使用颜色渐变
if (length(available_genes_train) > length(top_journal_colors)) {
  bioCol <- colorRampPalette(top_journal_colors)(length(available_genes_train))
} else {
  bioCol <- top_journal_colors[1:length(available_genes_train)]
}

# 初始化绘图
pdf("ROC_Curve_Testing_Set.pdf", width = 8, height = 8)
plot(NULL, xlim = c(1, 0), ylim = c(0, 1), 
     xlab = "False Positive Rate (1 - Specificity)",
     ylab = "True Positive Rate (Sensitivity)", 
     main = "ROC Curves - Testing Set (GSE235562)",
     las = 1, bty = "l", cex.lab = 1.2, cex.main = 1.3)
grid(col = "gray90", lty = 3)
segments(x0 = 1, y0 = 0, x1 = 0, y1 = 1, col = "gray50", lty = 2, lwd = 1.5)
# 循环绘制各基因ROC曲线，并计算AUC
aucText <- c()
for (i in 1:length(available_genes_train)) {
  gene <- available_genes_train[i]
  roc_obj <- roc(roc_data_test$Group, roc_data_test[[gene]])
  auc_val <- auc(roc_obj)
  aucText <- c(aucText, paste0(gene, " (AUC=", sprintf("%.3f", auc_val), ")"))
  
  if (i == 1) {
    plot(roc_obj, col = bioCol[i], lwd = 2.5, add = TRUE)
  } else {
    lines(roc_obj$specificities, roc_obj$sensitivities, col = bioCol[i], lwd = 2.5)
  }
}

# 添加图例
legend("bottomright", legend = aucText, col = bioCol, lwd = 2.5, 
       bty = "n", cex = 0.9, title = "Genes")
dev.off()

# 保存AUC结果
auc_results_test <- data.frame(
  Gene = available_genes_train,
  AUC = sapply(available_genes_train, function(gene) {
    auc(roc(roc_data_test$Group, roc_data_test[[gene]]))
  })
)
write.csv(auc_results_test, "AUC_Results_Testing_Set.csv", row.names = FALSE)


# 绘制验证集的ROC曲线

# 准备验证集数据（使用二分类Group_re）
valid_data_roc <- valid_data
valid_data_roc$Group_numeric <- ifelse(valid_data_roc$Group_re == "control", 0, 1)

# 创建ROC数据框
roc_data_valid <- data.frame(Group = valid_data_roc$Group_numeric)

# 添加每个基因的表达值
for(gene in available_genes_valid) {
  roc_data_valid[[gene]] <- valid_data_roc[[gene]]
}

# 初始化验证集ROC图
pdf("ROC_Curve_Validation_Set.pdf", width = 8, height = 8)
plot(NULL, xlim = c(1, 0), ylim = c(0, 1), 
     xlab = "False Positive Rate (1 - Specificity)",
     ylab = "True Positive Rate (Sensitivity)", 
     main = "ROC Curves - Validation Set (GSE106910)",
     las = 1, bty = "l", cex.lab = 1.2, cex.main = 1.3)
grid(col = "gray90", lty = 3)
segments(x0 = 1, y0 = 0, x1 = 0, y1 = 1, col = "gray50", lty = 2, lwd = 1.5)

# 循环绘制各基因ROC曲线，并计算AUC
aucText_valid <- c()
for (i in 1:length(available_genes_valid)) {
  gene <- available_genes_valid[i]
  roc_obj <- roc(roc_data_valid$Group, roc_data_valid[[gene]])
  auc_val <- auc(roc_obj)
  aucText_valid <- c(aucText_valid, paste0(gene, " (AUC=", sprintf("%.3f", auc_val), ")"))
  
  if (i == 1) {
    plot(roc_obj, col = bioCol[i], lwd = 2.5, add = TRUE)
  } else {
    lines(roc_obj$specificities, roc_obj$sensitivities, col = bioCol[i], lwd = 2.5)
  }
}

# 添加图例
legend("bottomright", legend = aucText_valid, col = bioCol, lwd = 2.5, 
       bty = "n", cex = 0.9, title = "Genes")
dev.off()

# 保存验证集AUC结果
auc_results_valid <- data.frame(
  Gene = available_genes_valid,
  AUC = sapply(available_genes_valid, function(gene) {
    auc(roc(roc_data_valid$Group, roc_data_valid[[gene]]))
  })
)
write.csv(auc_results_valid, "AUC_Results_Validation_Set.csv", row.names = FALSE)


# 合并AUC结果并生成总结报告
# 合并两个数据集的AUC结果
auc_summary <- full_join(auc_results_test, auc_results_valid, 
                         by = "Gene", suffix = c("_Test", "_Valid"))

# 计算平均AUC
auc_summary$Average_AUC <- rowMeans(auc_summary[, c("AUC_Test", "AUC_Valid")], na.rm = TRUE)

# 按平均AUC排序
auc_summary <- auc_summary[order(-auc_summary$Average_AUC), ]

# 保存总结报告
write.csv(auc_summary, "AUC_Summary_Both_Datasets.csv", row.names = FALSE)

# 生成AUC总结报告
cat("=== ROC曲线分析总结 ===\n\n")
cat("测试集 (GSE235562):\n")
for(i in 1:nrow(auc_results_test)) {
  cat(sprintf("%-10s: AUC = %.3f\n", 
              auc_results_test$Gene[i], auc_results_test$AUC[i]))
}

cat("\n验证集 (GSE106910):\n")
for(i in 1:nrow(auc_results_valid)) {
  cat(sprintf("%-10s: AUC = %.3f\n", 
              auc_results_valid$Gene[i], auc_results_valid$AUC[i]))
}

cat("\nAUC排名 (按平均值):\n")
for(i in 1:nrow(auc_summary)) {
  cat(sprintf("%-10s: 测试集=%.3f, 验证集=%.3f, 平均=%.3f\n", 
              auc_summary$Gene[i], 
              auc_summary$AUC_Test[i],
              auc_summary$AUC_Valid[i],
              auc_summary$Average_AUC[i]))
}


# 可选：使用ggplot2绘制更美观的ROC曲线
if(require(ggplot2) & require(plotROC)) {
  # 准备测试集数据（长格式）
  test_long_roc <- roc_data_test %>%
    pivot_longer(cols = -Group, names_to = "Gene", values_to = "Expression")
  p_test_roc <- ggplot(test_long_roc, aes(d = Group, m = Expression, color = Gene)) +
    geom_roc(n.cuts = 0, labels = FALSE) +
    style_roc(theme = theme_bw) +
    scale_color_manual(values = bioCol) +
    labs(title = "ROC Curves - Testing Set (GSE235562)",
         color = "Genes") +
    theme(legend.position = "right")
  
  ggsave("ROC_Curve_Testing_Set_ggplot.pdf", p_test_roc, width = 8, height = 6)
  
  # 准备验证集数据（长格式）
  valid_long_roc <- roc_data_valid %>%
    pivot_longer(cols = -Group, names_to = "Gene", values_to = "Expression")
  
  p_valid_roc <- ggplot(valid_long_roc, aes(d = Group, m = Expression, color = Gene)) +
    geom_roc(n.cuts = 0, labels = FALSE) +
    style_roc(theme = theme_bw) +
    scale_color_manual(values = bioCol) +
    labs(title = "ROC Curves - Validation Set (GSE106910)",
         color = "Genes") +
    theme(legend.position = "right")
  
  ggsave("ROC_Curve_Validation_Set_ggplot.pdf", p_valid_roc, width = 8, height = 6)
}


save.image("/data/nas1/xiangxiaoyu_OD/1.BPD+MSCs/1.scRNA/06_Expression_validation/Expression_validation+ROC.RData")

#.免疫浸润分析 -----------------------------------------------------------------
