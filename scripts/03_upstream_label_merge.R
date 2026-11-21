###############################################################################
# Label and merge the samples (adapted from the public pipeline of Shirazi et al. (https://github.com/SucreLab/HumanNeonatal))
#
# Source      : r.00_rawdata.r (lines 255-302)
# Input       : load_all_data.rds
# Output      : label_merge.rds
#
# Code below is unchanged from the original analysis script.
###############################################################################
# -----------------------------------------------------------------------------
# PATHS: this file keeps the absolute paths of the analysis server on which the
# study was run (/data/nas1/xiangxiaoyu_OD/1.BPD+MSCs/1.scRNA/... and
# /home/xxyloveme/study/BPD+MSCs/...). Set PROJ_DIR in config/paths.R and
# replace those prefixes with your own project directory before running.
# -----------------------------------------------------------------------------


############2_label_merge#########################
# ---
#   title: "scRNAseq Neonatal Lung Injury 1.2"
# author: Shawyon Shirazi & Nick Negretti
# date: 3/14/24
# output: rmarkdown::github_document
# ---
  
## Load libraries and helper functions
# ```{r, results="hide", message = FALSE, load-libraries}
source("./helper_functions/globals.R")
source("./helper_functions/libraries.R")
source("./helper_functions/cluster.R")
source("./helper_functions/colors.R")

N_WORKERS <- 12
plan("multicore", workers = N_WORKERS)

# Load premerged object
# ```{r load-data}
alldata <- readRDS("./load_all_data.rds")

# Function to merge list of data & split project name into descriptive labels
# ```{r merge}
merge_data <- function(alldata) {
  alldata <- Reduce(function(x, y) {merge(x, y)}, alldata)
  alldata$dataset <- "Unassigned"
  alldata$age <- "Unassigned"
  alldata$condition <- "Unassigned"
  alldata$broad_condition <- "Unassigned"
  alldata$dataset <- sapply(str_split(alldata$orig.ident, "_"),
                            function(x) {x[1]})
  alldata$age <- sapply(str_split(alldata$orig.ident, "_"),
                        function(x) {x[2]})
  alldata$condition <- sapply(str_split(alldata$orig.ident, "_"),
                              function(x) {x[3]})
  alldata$broad_condition <- sapply(str_split(alldata$orig.ident, "_"),
                                    function(x) {x[4]})
  return(alldata)
}

# Merge the data and save for next step
# ```{r merge all}
alldata <- merge_data(alldata)

saveRDS(alldata, "./label_merge.rds", compress = "bzip2")
