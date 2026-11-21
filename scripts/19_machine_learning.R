###############################################################################
# LASSO, Random Forest, SVM-RFE and XGBoost screening of candidate genes
#
# Source      : r.03_Cell_prop.r (lines 5790-6298)
# Input       : disease-trajectory genes + GSE235562
# Output      : Fig. 5C-H; Supplementary Table 16
#
# Code below is unchanged from the original analysis script.
###############################################################################
# -----------------------------------------------------------------------------
# PATHS: this file keeps the absolute paths of the analysis server on which the
# study was run (/data/nas1/xiangxiaoyu_OD/1.BPD+MSCs/1.scRNA/... and
# /home/xxyloveme/study/BPD+MSCs/...). Set PROJ_DIR in config/paths.R and
# replace those prefixes with your own project directory before running.
# -----------------------------------------------------------------------------


############################3.8 bulkRNA进行机器学习分析 GSE235562数据集###################

# 多组别机器学习分析代码
setwd("/data/nas1/xiangxiaoyu_OD/1.BPD+MSCs/1.scRNA/05_Machine_learning/")

# 加载必要的包
library(randomForest)
library(caret)
library(glmnet)
library(e1071)
library(pROC)
library(ggplot2)
library(patchwork)
library(xgboost)
library(Matrix)
library(PRROC)
library(dplyr)
library(purrr)
library(ggvenn)

# 读取数据
training_set <- read.csv("exprSet_GSE235562.csv", row.names = 1)
sample_info <- read.csv("group_GSE235562.csv", row.names = 1)

# 读取候选基因
available_genes_ME08 <- read.csv("ME8_Hub_Genes.csv", row.names = 1)
available_genes_ME00 <- read.csv("ME0_Hub_Genes.csv", row.names = 1)
available_genes <- c(available_genes_ME00$Gene, available_genes_ME08$Gene)

# 准备数据
candidate_expr <- as.data.frame(t(training_set[available_genes, ]))
candidate_expr$Group <- sample_info$group
candidate_expr$Group <- as.factor(candidate_expr$Group)
colnames(candidate_expr)[colnames(candidate_expr) == "ZFPM2-AS1"] <- "ZFPM2_AS1"

# 1. 随机森林 (多分类) ------------
# 检查是否有缺失值
cat("是否有缺失值:\n")
cat("Group 变量:", sum(is.na(candidate_expr$Group)), "\n")
cat("基因表达数据:", sum(is.na(candidate_expr[, -ncol(candidate_expr)])), "\n")
cat("\n=== 随机森林分析 (多分类) ===\n")
set.seed(123456)

# 创建二分类变量
candidate_expr_binary <- candidate_expr
candidate_expr_binary$Group <- ifelse(candidate_expr_binary$Group == "control", "control", "treatment")
candidate_expr_binary$Group <- as.factor(candidate_expr_binary$Group)

cat("二分类后的Group分布:\n")
print(table(candidate_expr_binary$Group))

# 测试二分类随机森林
set.seed(123456)
ntree_values <- c(10, 20, 30, 40, 50, 100, 150, 200)  # 先用较少的树测试

tryCatch({
  oob_errors_binary <- numeric(length(ntree_values))
  
  for (i in seq_along(ntree_values)) {
    rf_temp <- randomForest(Group ~ ., 
                            data = candidate_expr_binary, 
                            importance = TRUE, 
                            ntree = ntree_values[i])
    oob_errors_binary[i] <- tail(rf_temp$err.rate[, "OOB"], 1)
    cat("ntree =", ntree_values[i], "OOB error =", oob_errors_binary[i], "\n")
  }
  
  cat("二分类随机森林成功!\n")
  
}, error = function(e) {
  cat("二分类也失败:", e$message, "\n")
})

# 绘制OOB误差收敛图
error_data <- data.frame(Trees = ntree_values, OOB_Error = oob_errors_binary)
optimal_ntree <- error_data$Trees[which.min(error_data$OOB_Error)]

p_rf_error_convergence <- ggplot(error_data, aes(x = Trees, y = OOB_Error)) +
  geom_line(color = "blue", size = 1) +
  geom_point(color = "blue", size = 2) +
  geom_vline(xintercept = optimal_ntree, linetype = "dashed", color = "red", size = 1) +
  geom_text(aes(x = max(Trees), y = max(OOB_Error), 
                label = paste("Optimal ntree =", optimal_ntree)), 
            hjust = 1.1, vjust = 1, color = "red") +
  labs(title = "Random Forest - OOB Error Convergence (Multi-class)",
       x = "Number of Trees", y = "OOB Error Rate") +
  theme_classic()

# 使用最优ntree构建最终模型
rf_model <- randomForest(Group ~ ., 
                         data = candidate_expr_binary, 
                         importance = TRUE, 
                         ntree = 100,
                         proximity = TRUE)

print(rf_model)

# 特征重要性 (多分类版本)
rf_importance <- importance(rf_model)
rf_importance <- as.data.frame(rf_importance)
rf_importance$Gene <- rownames(rf_importance)

# 计算平均重要性（跨所有类别的MeanDecreaseAccuracy和MeanDecreaseGini）
rf_importance$MeanDecreaseAccuracy_Avg <- rowMeans(rf_importance[, grep("MeanDecreaseAccuracy", names(rf_importance))])
rf_importance$MeanDecreaseGini_Avg <- rf_importance$MeanDecreaseGini  # Gini已经是平均的

rf_importance <- rf_importance[order(-rf_importance$MeanDecreaseAccuracy), ]
write.csv(rf_importance, "RF_all_genes_importance_multiclass.csv", row.names = FALSE)

# 特征重要性条形图
p_rf_importance <- ggplot(rf_importance, 
                          aes(x = reorder(Gene,MeanDecreaseGini), 
                              y = MeanDecreaseGini,
                              fill = MeanDecreaseGini)) +
  geom_bar(stat = "identity", alpha = 0.8) +
  scale_fill_gradient(low = "lightblue", high = "darkblue", name = "Importance") +
  geom_hline(yintercept = 0, col = 'grey',linetype = 'dotted')+
  coord_flip() +
  labs(title = "Random Forest - Feature Importance",
       x = "Gene", y = "MeanDecreaseGini") +
  theme_classic()

# 保存图形
ggsave("RF_error_convergence_multiclass.pdf", p_rf_error_convergence, width = 4.7, height = 4)
ggsave("RF_importance_multiclass.pdf", p_rf_importance, width = 3.7, height = 7)


# 2. SVM-RFE (多分类) ------------
cat("\n=== SVM-RFE分析 (多分类) ===\n")
set.seed(123456)

X <- candidate_expr[, -ncol(candidate_expr)]
Y <- candidate_expr_binary$Group

# 多分类SVM-RFE
ctrl <- trainControl(method = "cv", number = 5, verboseIter = FALSE)
svm_profile <- rfe(X, Y,
                   sizes = c(1:length(available_genes)),  # 限制特征数量以避免计算时间过长
                   rfeControl = rfeControl(functions = caretFuncs,
                                           number = 5),
                   method = "svmRadial",
                   trControl = ctrl,
                   metric = "Accuracy")  # 使用准确率作为评估指标
print(svm_profile)

# 提取结果
svm_results <- svm_profile$results
svm_variables <- svm_profile$variables

# 计算每个基因的平均重要性
gene_importance <- svm_variables %>%
  group_by(var) %>%
  summarise(Mean_Importance = mean(Overall, na.rm = TRUE)) %>%
  arrange(desc(Mean_Importance))

colnames(gene_importance) <- c("Gene", "Mean_Importance")
write.csv(gene_importance, "SVM_RFE_all_genes_importance_multiclass.csv", row.names = FALSE)

# SVM-RFE准确率图
p_svm_accuracy <- ggplot(svm_results, aes(x = Variables, y = Accuracy)) +
  geom_line(color = "blue", size = 1) +
  geom_point(color = "blue", size = 2) +
  labs(title = "SVM-RFE - Accuracy vs Feature Number (Multi-class)",
       x = "Number of Features", y = "Cross-Validation Accuracy") +
  theme_classic()
ggsave("SVM_RFE_accuracy_multiclass.pdf", p_svm_accuracy, width = 9.7, height = 4.6)

# 计算错误率 (1 - 准确率)
svm_results$ErrorRate <- 1 - svm_results$Accuracy
# 新增：Figure 5D: SVM-RFE特征数量与错误率关系
p_svm_error <- ggplot(svm_results, aes(x = Variables)) +
  geom_line(aes(y = ErrorRate, color = "Error Rate"), size = 1) +
  geom_point(aes(y = ErrorRate, color = "Error Rate"), size = 2) +
  scale_color_manual(values = c("Error Rate" = "red")) +
  labs(title = "SVM-RFE - Error Rate vs Feature Number",
       x = "Number of Features",
       y = "Cross-Validation Error Rate",
       color = "Metric") +
  theme_classic() +
  theme(legend.position = "top")
ggsave("SVM_RFE_accuracy_multiclass.pdf", p_svm_error, width = 9.7, height = 4.6)

# 3. Lasso回归 (多分类) ------------
cat("\n=== Lasso回归分析 (多分类) ===\n")
set.seed(123456)
x <- as.matrix(scale(X, center = TRUE, scale = TRUE))
y_numeric <- ifelse(Y == "treatment", 1, 0)
# 多分类Lasso (multinomial)
# 首先拟合完整的Lasso模型（不指定lambda）
lasso_model_full <- glmnet(x, y_numeric, 
                           family = "binomial", 
                           alpha = 1)

# 5倍交叉验证Lasso
cv_lasso <- cv.glmnet(x, y_numeric, 
                      family = "binomial",
                      alpha = 1,  # Lasso回归
                      nfolds = 5,
                      type.measure = "deviance")
plot(cv_lasso)

# 查看最佳lambda值
best_lambda <- cv_lasso$lambda.min
cat("最佳lambda值:", best_lambda, "\n")

lambda_1se <- cv_lasso$lambda.1se
cat("最佳lambda.1se值:", lambda_1se, "\n")

# 使用最佳lambda拟合最终模型（用于基因选择）
lasso_model <- glmnet(x, y_numeric, 
                      family = "binomial", 
                      alpha = 1, 
                      lambda = best_lambda)
# 获取系数
lasso_coef <- as.matrix(coef(lasso_model))
lasso_genes <- rownames(lasso_coef)[which(lasso_coef != 0)]
lasso_genes <- lasso_genes[-1]  # 移除截距项
cat("Lasso选择的基因数量:", length(lasso_genes), "\n")
print(lasso_genes)
# [1] "AES"      "ACTG1"    "TSC22D2"  "NFAT5"    "CCDC102B" "PPIA"     "TMSB4X"   "MT1A"    
# [9] "PHLDA1"   "HSPA1A"   "SLC26A7"

# 使用完整的Lasso模型绘制系数路径
coef_matrix <- as.matrix(coef(lasso_model_full))
# 准备系数路径数据
lambda_values <- lasso_model_full$lambda
coef_path_data <- data.frame()

for (i in 2:nrow(coef_matrix)) {  # 跳过截距项
  gene_name <- rownames(coef_matrix)[i]
  gene_coefs <- coef_matrix[i, ]
  
  temp_df <- data.frame(
    log_lambda = log(lambda_values),
    Coefficient = gene_coefs,
    Gene = gene_name
  )
  coef_path_data <- rbind(coef_path_data, temp_df)
}

# 绘制系数路径图
library(RColorBrewer)
library(ggsci)
colpalettes <-unique(c(pal_npg("nrc")(10),pal_aaas("default")(10),pal_nejm("default")(8),pal_lancet("lanonc")(9),
                      pal_jama("default")(7),pal_jco("default")(10),pal_ucscgb("default")(26),pal_d3("category10")(10),
                      pal_locuszoom("default")(7),pal_igv("default")(51),
                      pal_uchicago("default")(9),pal_startrek("uniform")(7),
                      pal_tron("legacy")(7),pal_futurama("planetexpress")(12),pal_rickandmorty("schwifty")(12),
                      pal_simpsons("springfield")(16),pal_gsea("default")(12)))

pie(1:219,col = colpalettes)
p_lasso_coef <- ggplot(coef_path_data, aes(x = log_lambda, y = Coefficient, color = Gene)) + 
  geom_vline(xintercept = log(best_lambda), size = 0.8, color = 'grey60', alpha = 0.8, linetype = 2) +
  geom_vline(xintercept = log(lambda_1se), size = 0.8, color = 'grey50', alpha = 0.8, linetype = 3) +
  geom_line(size = 1) +
  xlab("Log(Lambda)") + 
  ylab("Coefficients") + 
  theme_classic() +
  scale_x_continuous(expand = c(0.01, 0.01)) + 
  scale_y_continuous(expand = c(0.01, 0.01)) + 
  scale_color_manual(values = colpalettes)+
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size = 12, color = 'black'), 
        axis.text = element_text(size = 10, color = 'black'), 
        legend.title = element_blank(), 
        legend.text = element_text(size = 8, color = 'black'), 
        legend.position = 'right') +
  annotate('text', x = log(best_lambda), y = max(coef_path_data$Coefficient) * 0.9,
           label = paste0('lambda.min\n', round(best_lambda, 5)), 
           color = 'black', size = 3) +
  annotate('text', x = log(lambda_1se), y = max(coef_path_data$Coefficient) * 0.7,
           label = paste0('lambda.1se\n', round(lambda_1se, 5)), 
           color = 'black', size = 3)+
  guides(col = guide_legend(ncol = 4))

# 显示图形
print(p_lasso_coef)

#改进的Lasso交叉验证图
lasso_cv_data <- data.frame(
  lambda = cv_lasso$lambda,
  cvm = cv_lasso$cvm,
  cvsd = cv_lasso$cvsd,
  cvup = cv_lasso$cvup,
  cvlo = cv_lasso$cvlo,
  nzero = cv_lasso$nzero
)

lasso_cv_data$log_lambda <- log(lasso_cv_data$lambda)
lasso_cv_data$NZERO <- paste0(lasso_cv_data$nzero, ' genes')

p_lasso_cv <- ggplot(lasso_cv_data, aes(x = log_lambda, y = cvm, color = NZERO)) +
  geom_errorbar(aes(ymin = cvlo, ymax = cvup), width = 0.05, size = 0.5) +
  geom_vline(xintercept = log(best_lambda), size = 0.8, color = 'grey60', alpha = 0.8, linetype = 2) +
  geom_vline(xintercept = log(lambda_1se), size = 0.8, color = 'grey50', alpha = 0.8, linetype = 3) +
  geom_point(size = 2) +
  xlab("Log(Lambda)") +
  ylab('Binomial Deviance') +
  theme_classic() +
  scale_x_continuous(expand = c(0.02, 0.02)) +
  scale_y_continuous(expand = c(0.02, 0.02)) +
  scale_color_manual(values = colpalettes)+
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(size = 12, color = 'black'),
    axis.text = element_text(size = 10, color = 'black'),
    legend.title = element_blank(),
    legend.text = element_text(size = 8, color = 'black'),
    legend.position = 'bottom'
  ) +
  annotate('text', x = log(best_lambda), 
           y = max(lasso_cv_data$cvm) * 0.95, 
           label = paste0('lambda.min\n', round(best_lambda, 5)), 
           color = 'black', size = 3) +
  annotate('text', x = log(lambda_1se), 
           y = max(lasso_cv_data$cvm) * 0.75, 
           label = paste0('lambda.1se\n', round(lambda_1se, 5)), 
           color = 'black', size = 3)

# 保存Lasso图形
ggsave("Lasso_coef_path.pdf", p_lasso_coef, width = 10, height = 4)
ggsave("Lasso_CV.pdf", p_lasso_cv, width = 7, height = 4)

# 保存选择的基因
cat("选择的基因及其系数:\n")
for (gene in lasso_genes) {
  coef_value <- lasso_coef[gene, 1]
  cat(sprintf("%-10s: %8.4f\n", gene, coef_value))
}
# AES       :  -1.0878
# ACTG1     :  -0.4088
# TSC22D2   :  -0.7919
# NFAT5     :   0.2921
# CCDC102B  :   0.1259
# PPIA      :  -0.2119
# TMSB4X    :  -0.2433
# MT1A      :   0.0616
# PHLDA1    :   0.3161
# HSPA1A    :   0.7977
# SLC26A7   :  -0.1211

# 保存选择的基因
lasso_results <- data.frame(
  Gene = lasso_genes,
  Coefficient = as.numeric(lasso_coef[lasso_genes, ]))
write.csv(lasso_results, "lasso_selected_genes.csv", row.names = FALSE)

# cv_lasso <- cv.glmnet(x, y_numeric, 
#                       family = "binomial",
#                       alpha = 1,
#                       nfolds = 5,
#                       type.measure = "deviance")
# best_lambda <- cv_lasso$lambda.min
# lambda_1se <- cv_lasso$lambda.1se
# 
# cat("最佳lambda值:", best_lambda, "\n")
# cat("lambda.1se值:", lambda_1se, "\n")
# 
# # 使用最佳lambda拟合模型
# lasso_model <- glmnet(x, y_numeric, 
#                       family = "binomial", 
#                       alpha = 1, 
#                       lambda = best_lambda)
# # 提取非零系数基因
# lasso_coef <- coef(lasso_model)
# selected_genes <- unique(unlist(lapply(lasso_coef, function(coef_mat) {
#   genes <- rownames(coef_mat)[which(coef_mat != 0)]
#   genes <- genes[genes != "(Intercept)"]
#   return(genes)
# })))
# 
# cat("Lasso选择的基因数量:", length(selected_genes), "\n")
# print(selected_genes)
# # [1] "ACTB"     "RALGAPA1" "AES"      "CDKN1A"   "RLF"      "FOSB"     "ACTG1"    "TSC22D2" 
# # [9] "HSPA1A"   "SLC26A7"  "NFATC2"   "PPIA"     "FAM133B"  "TMSB4X"   "INTS6" 
# 
# # 系数路径图 (多分类)
# lasso_model_full <- glmnet(x, Y, family = "binomial", alpha = 1)
# 
# # 绘制系数路径图
# pdf("Lasso_coef_path_multiclass.pdf", width = 8, height = 6)
# par(mfrow = c(2, 3))
# for (i in 1:length(lasso_model_full$beta)) {
#   class_name <- names(lasso_model_full$beta)[i]
#   plot(lasso_model_full, main = paste("Coefficient Path -", class_name))
#   abline(v = log(best_lambda), col = "red", lty = 2)
#   abline(v = log(lambda_1se), col = "blue", lty = 2)
# }
# dev.off()
# 
# # 交叉验证图
# pdf("Lasso_CV_multiclass.pdf", width = 7, height = 5)
# plot(cv_lasso)
# title("Lasso Regression - Cross Validation (Multi-class)")
# dev.off()

# 4. XGBoost (多分类) ------------
cat("\n=== XGBoost分析 (多分类) ===\n")
set.seed(123456)
# 将因子标签转换为数值 (从0开始)
y_numeric <- as.numeric(Y) - 1
num_class <- length(unique(y_numeric))

# 创建XGBoost数据
x_matrix <- as.matrix(X)
dtrain <- xgb.DMatrix(data = x_matrix, label = y_numeric)

# XGBoost参数 (多分类)
params <- list(
  objective = "multi:softprob",
  eval_metric = "mlogloss",
  num_class = num_class,
  max_depth = 6,
  eta = 0.1,
  subsample = 0.8,
  colsample_bytree = 0.8)

# 训练模型
xgb_model <- xgb.train(
  params = params,
  data = dtrain,
  nrounds = 100,
  verbose = 1)

# 特征重要性
xgb_importance <- xgb.importance(feature_names = colnames(X), model = xgb_model)
xgb_genes <- data.frame(gene = xgb_importance$Feature)
# 特征重要性图
p_xgb_importance <- ggplot(xgb_importance, 
                           aes(x = reorder(Feature, Gain), y = Gain, fill = Gain)) +
  geom_bar(stat = "identity", alpha = 0.8) +
  scale_fill_gradient(low = "lightblue", high = "darkblue") +
  coord_flip() +
  labs(title = "XGBoost Feature Importance (Top 15, Multi-class)",
       x = "Gene", y = "Gain") +
  theme_classic()

ggsave("XGBoost_importance_multiclass.pdf", p_xgb_importance, width = 3.7, height = 7)

# 保存重要基因
write.csv(xgb_importance, "XGBoost_important_genes_multiclass.csv", row.names = FALSE)


# 5. 韦恩图分析 ------------
cat("\n=== 韦恩图分析 ===\n")
# 从各方法中选择重要基因
##RF选择重要性基因
# 绘制重要性下降曲线
importance_scores <- rf_importance$MeanDecreaseAccuracy
genes_ordered <- rf_importance$Gene
# 计算相对重要性（相对于最大值）
relative_importance <- importance_scores / max(importance_scores)
# 绘制重要性下降曲线
plot_data <- data.frame(
  Gene = factor(genes_ordered, levels = genes_ordered),
  Importance = relative_importance)

p_elbow <- ggplot(plot_data, aes(x = seq_along(Gene), y = Importance)) +
  geom_line(color = "blue", size = 1) +
  geom_point(color = "blue", size = 2) +
  geom_vline(xintercept = c(1, 5, 10, 15, 20,25,30,35), linetype = "dashed", alpha = 0.5) +
  labs(title = "Feature Importance - Elbow Method",
       x = "Rank (by Importance)",
       y = "Relative Importance") +
  theme_classic()
print(p_elbow)
##区肘部重要性下降的基因：从第6个斜率陡然上升
ggsave("RF_p_elbow.pdf", p_elbow, width = 7, height = 4)

#RF选择
rf_selected <- head(rf_importance$Gene[order(-rf_importance$MeanDecreaseGini)], 14)
# SVM-RFE: 选择最优特征集中的基因
svm_selected <- svm_profile$optVariables
# Lasso: 已选择的基因
lasso_selected <- lasso_genes
# XGBoost: 选择重要性前10的基因
xgb_selected <- head(xgb_importance$Feature, 10)
xgb_selected <- xgb_importance$Feature
# 创建韦恩图数据
venn_list <- list(
  RandomForest = rf_selected,
  SVM_RFE = svm_selected,
  Lasso = lasso_selected,
  XGBoost = xgb_selected)

# 绘制韦恩图
venn_plot <- ggvenn(venn_list, 
                    fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
                    stroke_size = 0.5, set_name_size = 4) +
  labs(title = "Feature Selection Venn Diagram (Multi-class)")

ggsave("Venn_diagram_multiclass.pdf", venn_plot, width = 8, height = 6)

# 找出交集基因
common_genes <- Reduce(intersect, venn_list)
# [1] "AES"     "ACTG1"   "TSC22D2" "HSPA1A"  "PHLDA1"  "SLC26A7" "MT1A"    "NFAT5"

# 保存所有选择的基因
all_selected_genes <- data.frame(
  Method = c(rep("RandomForest", length(rf_selected)),
             rep("SVM_RFE", length(svm_selected)),
             rep("Lasso", length(lasso_selected)),
             rep("XGBoost", length(xgb_selected))),
  Gene = c(rf_selected, svm_selected, lasso_selected, xgb_selected))

write.csv(all_selected_genes, "ML_all_selected_feature_genes_multiclass.csv", row.names = FALSE)
write.csv(data.frame(geneNames = common_genes),file = "common_genes.csv")
save.image(file = "/data/nas1/xiangxiaoyu_OD/1.BPD+MSCs/1.scRNA/05_Machine_learning/Machine.RData")
