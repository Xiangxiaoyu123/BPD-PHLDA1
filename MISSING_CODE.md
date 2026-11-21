# Code that is not in this repository

This file lists the analysis steps of the manuscript for which no script is
included, so that it is clear what the repository does and does not cover.

## 1. Image analysis: immunohistochemistry (Fig. 8H–N) and multiplex immunofluorescence (Fig. 9, Supplementary Fig. 14)

The H-score, the percentage of DAB-positive cells and the intensity grades were
produced interactively in QuPath, and the multiplex immunofluorescence
quantification (adventitial ring, PDGFRA⁺PHLDA1⁺ classification, per-vessel
proportions, sensitivity analyses and the nearest-neighbour permutation test)
was run on the cell tables exported from the Vectra/QuPath workflow. The
scripts and the QuPath project were not kept in a form that can be published.

Before uploading the repository, add whichever of the following you still have:

- the QuPath project or the `.groovy` scripts used for cell detection and the
  positive-cell thresholds;
- the per-cell measurement tables exported from QuPath;
- the script that computed the H-score, the enrichment index and the
  distance-to-CD8 curves.

The procedure itself is described in the Methods of the manuscript
("Histology and immunohistochemistry", "Multiplex immunofluorescence").

## 2. The original scTenifoldKnk run

The script that produced the published PHLDA1 virtual knockout result
(Fig. 8A–C, Supplementary Table 18) was lost.
`scripts/24_scTenifoldKnk_template.R` is the template it was adapted from, and
it still carries the gene name and the paths of the project it came from.
`revision_scripts/B5b_scTenifoldKnk.R` re-runs the analysis for PHLDA1 and is
the script to use.

## 3. CIBERSORT

`scripts/21_immune_deconvolution_xCell_CIBERSORT.R` contains
`source("CIBERSORT.R")`. The CIBERSORT source code is distributed under its own
licence by https://cibersortx.stanford.edu and cannot be redistributed here.
Request it from the authors and place `CIBERSORT.R` in the working directory.

## 4. Supplementary tables

The supplementary tables were assembled from the outputs of the scripts listed
in `SCRIPT_MAP.csv`; the spreadsheet-formatting steps are not included.
