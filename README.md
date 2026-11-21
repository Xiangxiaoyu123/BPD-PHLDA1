# Code for: Single-cell profiling identifies PHLDA1 as a candidate marker of adventitial fibroblasts associated with altered inferred immune cell composition in bronchopulmonary dysplasia

Wu Deng, Xiaoyu Xiang, Sizhe Yang, Yi Heng Zan

This repository contains the R code used for the analyses in the manuscript.
`SCRIPT_MAP.csv` lists every script with its origin, its inputs and the figures
and tables it produces. `MISSING_CODE.md` lists the steps for which no script is
included.

## Layout

```
config/            paths.R — the only file that has to be edited
scripts/           the original analysis, in execution order (01–24)
  exploratory_not_reported/   datasets examined but not reported in the manuscript
revision_scripts/  analyses added in response to the reviewers (B1–B10)
helper_functions/  helper functions from the upstream pipeline (see below)
data/              empty; download the GEO datasets here (data/README.md)
results/           empty; output of the scripts
env/               how to record the R environment
```

## Data

All datasets are public (GEO: GSE275938, GSE242310, GSE235562, GSE106910,
GSE164358). See `data/README.md` for what to download and where to put it.
The human lung tissue sections used for histology are archival autopsy
specimens and are not public.

## Run order

The scripts in `scripts/` are numbered in the order in which they were run:

| Step | Script | What it does |
|------|--------|--------------|
| 01–05 | load, QC, merge, SCTransform, Harmony, split into lineages | preprocessing of GSE275938 |
| 06–10 | clustering and lineage annotation; per-lineage annotation | 30 clusters, four lineages, seven mesenchymal subpopulations |
| 11–13 | merge the annotated lineages, re-embed, lineage figures | Fig. 2A–C; Supplementary Figs. 1–2 |
| 14–15 | mesenchymal composition and differential expression; adventitial fibroblast subclusters and enrichment | Figs. 2D–H, 3 |
| 16–17 | CellChat; Monocle 2 pseudotime | Fig. 4 |
| 18–19 | WGCNA; LASSO, Random Forest, SVM-RFE and XGBoost | Fig. 5 |
| 20 | candidate gene expression and ROC analysis | Fig. 6A–D |
| 21 | xCell and CIBERSORT deconvolution and correlation with PHLDA1 | Fig. 7 |
| 22–23 | rat hyperoxia dataset; nasopharyngeal basal cells | Fig. 6E–I |
| 24 | scTenifoldKnk (template; see `MISSING_CODE.md`) | Fig. 8A–C |

The scripts in `revision_scripts/` are run after `scripts/` and read the objects
the latter produce. They all begin with `source(here::here("revision_scripts",
"00_setup.R"))`, which holds the paths and the shared helper functions:

| Script | Reviewer comment | Produces |
|--------|------------------|----------|
| `B1_annotation_concordance.R` | R1-1 | agreement with the annotation and DEGs of the original study |
| `B2_donor_composition.R`, `B3_pseudobulk_DESeq2.R`, `B4_trajectory_LODO.R` | R1-2 | donor-level composition, pseudo-bulk DESeq2, leave-one-donor-out pseudotime |
| `B4a_rebuild_monocle2.R` | — | rebuilds the Monocle 2 object after loss of the saved object |
| `B5_sex_check.R`, `B5b_scTenifoldKnk.R`, `B5c_sex_figures.R`, `B5_replot_CNS.R` | R2-2 | donor sex, sex-chromosome sensitivity analyses |
| `B7_fix_supp_table7.R` | R2 minor 1 | removes the duplicated gene-symbol column |
| `B9_immune_scRNA.R`, `B9b_focus_AdvFB_CD8.R`, `B9c_immune_evidence_figures.R`, `B9_replot_CNS.R` | R2-3 | donor-level immune composition and adventitial FB–CD8⁺ T cell communication |
| `B10_dataset_audit.R` | editor | audit of the public datasets and of how they are described |

## Paths

The scripts are published exactly as they were run and therefore still contain
the absolute paths of the analysis servers:

```
/data/nas1/xiangxiaoyu_OD/1.BPD+MSCs/1.scRNA     (scripts/)
/home/xxyloveme/study/BPD+MSCs                   (revision_scripts/)
```

Set `PROJ_DIR` in `config/paths.R` and replace these two prefixes with your own
project directory before running anything. Nothing else in the code was
changed, with one exception: in `revision_scripts/`, the line that sourced
`00_setup_R2.R` by absolute path now reads
`source(here::here("revision_scripts", "00_setup.R"))`, and `00_setup.R` is the
file that was called `00_setup_R2.R` during the revision.

## Environment

R 4.3.2 with Seurat 4.4.0. Main packages: Seurat, harmony, sctransform,
glmGamPoi, MAST, DESeq2, speckle, CellChat, monocle (v2), WGCNA, glmnet,
randomForest, caret, e1071, xgboost, pROC, xCell, CIBERSORT, scTenifoldKnk,
clusterProfiler, ggplot2, ComplexHeatmap. See `env/README.md`.

## Upstream code

GSE275938 was generated and first analyzed by Shirazi et al. The preprocessing
and annotation of that dataset were run with the authors' public pipeline,
https://github.com/SucreLab/HumanNeonatal. The scripts whose names contain
`upstream`, the files in `helper_functions/` and `scripts/11`–`13` come from
that repository, in some cases adapted. Please cite the original study when
using them.

## Citation

Please cite the manuscript. Zenodo DOI: [to be added after the Zenodo release].

## Licence

The original analysis code in this repository is released under the MIT License
(see `LICENSE`). Note that `helper_functions/` and the `upstream` scripts come
from https://github.com/SucreLab/HumanNeonatal and are covered by that
repository's terms.
