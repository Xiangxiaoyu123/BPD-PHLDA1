###############################################################################
# Read the GSE275938 CellBender .h5 files and build the Seurat objects
#
# Source      : r.00_rawdata.r (lines 1-133)
# Input       : GSE275938 CellBender-filtered .h5 files (GEO)
# Output      : seurat_list.rds
#
# Code below is unchanged from the original analysis script.
###############################################################################
# -----------------------------------------------------------------------------
# PATHS: this file keeps the absolute paths of the analysis server on which the
# study was run (/data/nas1/xiangxiaoyu_OD/1.BPD+MSCs/1.scRNA/... and
# /home/xxyloveme/study/BPD+MSCs/...). Set PROJ_DIR in config/paths.R and
# replace those prefixes with your own project directory before running.
# -----------------------------------------------------------------------------


#############1.单细胞数据及处理##############
setwd("/data/nas1/xiangxiaoyu_OD/1.BPD+MSCs/1.scRNA/01_Data_Download/")
#####单个样本读取-HDF5格式（.h5/.hdf5）
#####加载包
library(Seurat)
library(purrr)
library(Matrix)
library(hdf5r)
#加载数据#
remove(list = ls()) #清除 Global Environment

list.files()  #查看当前工作目录下的文件

h5_files <- c(
  Term_infant_2 = "GSM8488401_20d-norm-out_cellbender_filtered.h5",
  Term_infant_1 = "GSM8488400_0d-norm-out_cellbender_filtered.h5",
  BPD_PH_2 = "GSM8488399_5796-2-out_cellbender_filtered.h5",
  BPD_PH_1 = "GSM8488398_5796-1-out_cellbender_filtered.h5",
  BPD_2 = "GSM8488397_7mo-2-out_cellbender_filtered.h5",
  BPD_1 = "GSM8488396_7mo-1-out_cellbender_filtered.h5",
  Acute_preterm_injury = "GSM8488395_7466-1-filtered_feature_bc_matrix.h5"
)

# 查看第一个文件的结构
h5_file <- hdf5r::H5File$new(h5_files[["Sample1"]], mode = "r")
print(h5_file)
h5_file$close_all()

# 专门针对您的CellBender文件格式的读取函数
read_custom_cellbender_h5 <- function(filename) {
  if (!file.exists(filename)) {
    stop("File not found: ", filename)
  }
  
  # 打开H5文件
  infile <- H5File$new(filename = filename, mode = 'r')
  
  # 检查是否有matrix组
  if (!"matrix" %in% names(infile)) {
    infile$close_all()
    stop("File does not contain 'matrix' group: ", filename)
  }
  
  # 获取matrix组
  matrix_group <- infile[["matrix"]]
  
  # 读取稀疏矩阵组件
  data <- matrix_group[["data"]][]
  indices <- matrix_group[["indices"]][]
  indptr <- matrix_group[["indptr"]][]
  shape <- matrix_group[["shape"]][]
  
  # 读取barcodes和features
  barcodes <- matrix_group[["barcodes"]][]
  features <- matrix_group[["features/name"]][]
  
  # 关闭文件
  infile$close_all()
  
  # 处理重复的基因名
  if (any(duplicated(features))) {
    features <- make.unique(features, sep = "_")
  }
  
  # 创建稀疏矩阵
  sparse.mat <- sparseMatrix(
    i = indices + 1,  # R使用1-based索引
    p = indptr,
    x = as.numeric(data),
    dims = shape,
    dimnames = list(features, barcodes)
  )
  
  # 转换为dgCMatrix格式
  sparse.mat <- as(sparse.mat, "dgCMatrix")
  
  return(sparse.mat)
}

# 批量读取并创建Seurat对象
seurat_list <- map(names(h5_files), ~ {
  sample_name <- .x
  file_path <- h5_files[[.x]]
  
  message("Processing sample: ", sample_name)
  
  tryCatch({
    # 检查文件类型并选择适当的读取函数
    if(grepl("cellbender", file_path, ignore.case = TRUE)) {
      # 使用自定义函数读取CellBender格式的文件
      counts <- read_custom_cellbender_h5(file_path)
    } else {
      # 对于标准10X格式的文件
      counts <- Read10X_h5(filename = file_path)
    }
    
    # 创建Seurat对象
    seurat_object <- CreateSeuratObject(
      counts = counts,
      project = sample_name,
      min.cells = 3,
      min.features = 200
    )
    
    # 添加样本信息
    seurat_object[["sample"]] <- sample_name
    
    return(seurat_object)
  }, error = function(e) {
    message("Error processing sample ", sample_name, ": ", e$message)
    return(NULL)
  })
})

# 移除可能为NULL的元素（处理失败的样本）
seurat_list <- compact(seurat_list)

# 给列表命名
names(seurat_list) <- names(h5_files)[sapply(seurat_list, function(x) !is.null(x))]

# 查看成功处理的样本信息
if(length(seurat_list) > 0) {
  map(seurat_list, ~ {
    message("Sample: ", unique(.x$sample))
    message("Cells: ", ncol(.x))
    message("Features: ", nrow(.x))
    message("-------------------")
  })
} else {
  message("No samples were successfully processed.")
}

saveRDS(seurat_list,file = "seurat_list.rds")
