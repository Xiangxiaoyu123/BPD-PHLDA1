###############################################################################
# config/paths.R  —  the only file that has to be edited before running
#
#   source(here::here("config", "paths.R"))
#
# The scripts in scripts/ are kept exactly as they were run and still contain
# the absolute paths of the original analysis server. Set PROJ_DIR below and
# replace those prefixes (see README, "Paths") before running them.
###############################################################################

PROJ_DIR    <- normalizePath(".", mustWork = FALSE)   # <- project directory
DATA_DIR    <- file.path(PROJ_DIR, "data")            # GEO downloads
RESULTS_DIR <- file.path(PROJ_DIR, "results")
FIG_DIR     <- file.path(RESULTS_DIR, "figures")

for (d in c(DATA_DIR, RESULTS_DIR, FIG_DIR))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

# Prefixes that appear in scripts/ and revision_scripts/ and have to be mapped
# onto PROJ_DIR:
#   /data/nas1/xiangxiaoyu_OD/1.BPD+MSCs/1.scRNA     (original analysis)
#   /home/xxyloveme/study/BPD+MSCs                   (revision analyses)
OLD_PREFIXES <- c("/data/nas1/xiangxiaoyu_OD/1.BPD+MSCs/1.scRNA",
                  "/home/xxyloveme/study/BPD+MSCs")

set.seed(2026)
options(stringsAsFactors = FALSE)
