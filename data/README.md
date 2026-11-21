# Data

No data files are kept in this repository. All datasets analyzed in the study
are public and can be downloaded from the Gene Expression Omnibus
(https://www.ncbi.nlm.nih.gov/geo/):

| Accession  | Platform      | Content                                                                 | Used in |
|------------|---------------|-------------------------------------------------------------------------|---------|
| GSE275938  | scRNA-seq     | Human neonatal lung: term (n = 2), BPD (n = 2), BPD+PH (n = 2), acute preterm injury (n = 1) | primary dataset |
| GSE242310  | scRNA-seq     | Rat lung, 95% O2 or room air from postnatal day 1 to 10, collected on day 10 | cross-species validation |
| GSE235562  | microarray    | Tracheal aspirate-derived MSCs, 24 h cyclic mechanical stretch and/or 40%/80% O2 | WGCNA, machine learning, test set |
| GSE106910  | bulk RNA-seq  | Tracheal aspirate-derived MSC lines from 19 infants                      | validation set, severity strata |
| GSE164358  | bulk RNA-seq  | Nasopharyngeal basal cells from preterm and term infants                 | additional validation |

Download each accession into `data/<accession>/`. For GSE275938 the scripts
expect the CellBender-filtered `.h5` files and, for the comparison with the
original annotation, `GSE275938_cell_metadata.csv.gz`.

The human lung tissue sections used for histology, immunohistochemistry and
multiplex immunofluorescence are archival autopsy specimens and are not public.
