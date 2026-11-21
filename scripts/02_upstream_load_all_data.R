###############################################################################
# Load and QC-filter all samples (adapted from the public pipeline of Shirazi et al. (https://github.com/SucreLab/HumanNeonatal))
#
# Source      : r.00_rawdata.r (lines 137-254)
# Input       : GSE275938 .h5 files
# Output      : load_all_data.rds
#
# Code below is unchanged from the original analysis script.
###############################################################################
# -----------------------------------------------------------------------------
# PATHS: this file keeps the absolute paths of the analysis server on which the
# study was run (/data/nas1/xiangxiaoyu_OD/1.BPD+MSCs/1.scRNA/... and
# /home/xxyloveme/study/BPD+MSCs/...). Set PROJ_DIR in config/paths.R and
# replace those prefixes with your own project directory before running.
# -----------------------------------------------------------------------------


###########################文献中的代码##############
############1_load_all_data##########################
# ---
#   title: "scRNAseq Neonatal Lung Injury 1.1"
# author: Shawyon Shirazi & Nick Negretti
# date: 3/13/24
# output: rmarkdown::github_document
# ---
#   
#   ```{r}
# getwd()
# ```

## Load libraries and helper functions
# ```{r load libraries & set up helper functions}
source("./helper_functions/globals.R")
source("./helper_functions/libraries.R")
source("./helper_functions/cluster.R")
source("./helper_functions/colors.R")

N_WORKERS <- 12
plan("multicore", workers = N_WORKERS)

## Set up parallel processing futures for loading, annotating, and cleaning data
# ```{r construct-futures}

pct_mt_high <- 10
pct_mt_low <- 0
min_nfeature <- 500
max_nfeature <- 7000

terminfant_1 %<-% {CreateSeuratObject(counts = ReadCB_h5("GSM8488400_0d-norm-out_cellbender_filtered.h5"),
                                      project = "Term infant 1_0d_Term infant_Normal") %>%
    PercentageFeatureSet(., pattern = "^MT-",
                         col.name = "percent.mt") %>%
    subset(., subset = percent.mt < pct_mt_high &
             percent.mt > pct_mt_low &
             nFeature_RNA > min_nfeature &
             nFeature_RNA < max_nfeature)
}

terminfant_2 %<-% {CreateSeuratObject(counts = ReadCB_h5("GSM8488401_20d-norm-out_cellbender_filtered.h5"),
                                      project = "Term infant 2_20d_Term infant_Normal") %>%
    PercentageFeatureSet(., pattern = "^MT-",
                         col.name = "percent.mt") %>%
    subset(., subset = percent.mt < pct_mt_high &
             percent.mt > pct_mt_low &
             nFeature_RNA > min_nfeature &
             nFeature_RNA < max_nfeature)
}

bpdph_1 %<-% {CreateSeuratObject(counts = ReadCB_h5("GSM8488398_5796-1-out_cellbender_filtered.h5"),
                                 project = "BPD+PH 1_7m_BPD+PH_Severe") %>%
    PercentageFeatureSet(., pattern = "^MT-",
                         col.name = "percent.mt") %>%
    subset(., subset = percent.mt < pct_mt_high &
             percent.mt > pct_mt_low &
             nFeature_RNA > min_nfeature &
             nFeature_RNA < max_nfeature)
}

bpdph_2 %<-% {CreateSeuratObject(counts = ReadCB_h5("GSM8488399_5796-2-out_cellbender_filtered.h5"),
                                 project = "BPD+PH 2_7m_BPD+PH_Severe") %>%
    PercentageFeatureSet(., pattern = "^MT-",
                         col.name = "percent.mt") %>%
    subset(., subset = percent.mt < pct_mt_high &
             percent.mt > pct_mt_low &
             nFeature_RNA > min_nfeature &
             nFeature_RNA < max_nfeature)
}

preterm_1 %<-% {CreateSeuratObject(counts = ReadCB_h5("GSM8488395_7466-1-filtered_feature_bc_matrix.h5"),
                                   project = "Acute preterm injury 1_26wkGA_Acute preterm injury_Preterm") %>%
    PercentageFeatureSet(., pattern = "^MT-",
                         col.name = "percent.mt") %>%
    subset(., subset = percent.mt < pct_mt_high &
             percent.mt > pct_mt_low &
             nFeature_RNA > min_nfeature &
             nFeature_RNA < max_nfeature)
}

bpd_1 %<-% {CreateSeuratObject(counts = ReadCB_h5("GSM8488396_7mo-1-out_cellbender_filtered.h5"),
                               project = "BPD 1_7m_BPD_Severe") %>%
    PercentageFeatureSet(., pattern = "^MT-",
                         col.name = "percent.mt") %>%
    subset(., subset = percent.mt < pct_mt_high &
             percent.mt > pct_mt_low &
             nFeature_RNA > min_nfeature &
             nFeature_RNA < max_nfeature)
}

bpd_2 %<-% {CreateSeuratObject(counts = ReadCB_h5("GSM8488397_7mo-2-out_cellbender_filtered.h5"),
                               project = "BPD 2_7m_BPD_Severe") %>%
    PercentageFeatureSet(., pattern = "^MT-",
                         col.name = "percent.mt") %>%
    subset(., subset = percent.mt < pct_mt_high &
             percent.mt > pct_mt_low &
             nFeature_RNA > min_nfeature &
             nFeature_RNA < max_nfeature)
}

## Activate data loading and store in a list of objects
# ```{r create-data-list}
expression_objs <- list(terminfant_1 = terminfant_1,
                        terminfant_2 = terminfant_2,
                        bpdph_1 = bpdph_1,
                        bpdph_2 = bpdph_2,
                        preterm_1 = preterm_1,
                        bpd_1 = bpd_1,
                        bpd_2 = bpd_2)

## Save data for next step
# ```{r save-listobj}
saveRDS(expression_objs,
        "./load_all_data.rds",
        compress = "bzip2")
