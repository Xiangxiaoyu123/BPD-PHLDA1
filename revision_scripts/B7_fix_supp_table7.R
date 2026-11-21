###############################################################################
# B7  修正 Supplementary Table 7（相邻两列基因名重复）   —— 回应 R2-minor 1
#   同时检查其他补充表是否也存在“内容完全相同的重复列”
###############################################################################
source(here::here("revision_scripts", "00_setup.R"))
library(openxlsx)
od <- file.path(OUT, "B7_supptable")

SUPP_DIR <- file.path(BASE, "Supplementary_Tables")      # <- 改成补充表所在目录
files <- list.files(SUPP_DIR, pattern = "\\.xlsx$", full.names = TRUE)

report <- list()
for (f in files) {
  wb_out <- createWorkbook()
  for (sh in getSheetNames(f)) {
    df <- read.xlsx(f, sheet = sh, check.names = FALSE)
    dup <- duplicated(lapply(df, function(x) as.character(x)))        # 内容完全相同的列
    if (any(dup)) report[[length(report) + 1]] <-
      data.frame(file = basename(f), sheet = sh, removed_columns = paste(colnames(df)[dup], collapse = "; "))
    df <- df[, !dup, drop = FALSE]
    addWorksheet(wb_out, sh); writeData(wb_out, sh, df)
  }
  saveWorkbook(wb_out, file.path(od, basename(f)), overwrite = TRUE)
}
report <- bind_rows(report)
write.csv(report, file.path(od, "duplicated_columns_removed.csv"), row.names = FALSE)
print(report)
# 注意：若 Table 7 的两列只是“列名不同但内容都是基因名”，会被识别并删除；
#       若一列是 gene symbol、另一列是 Ensembl ID，则不会删除——这种情况请手工确认是否需要保留。
