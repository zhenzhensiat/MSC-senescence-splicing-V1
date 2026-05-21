# MSC Senescence Alternative Splicing Analysis

Code repository for "Alternative splicing constitutes a transcription-independent regulatory layer in replicative senescence of human umbilical cord mesenchymal stem cells."

## Repository structure

```
.
├── 00_project_config.R              # Global configuration (source this first)
├── theme_bindlab.R                  # Nature-style ggplot2 theme
├── 01_QC_and_senescence_validation.R     # PCA, SenMayo, QC metrics
├── 02_differential_expression.R          # DESeq2 DEG analysis
├── 03_alternative_splicing.R             # rMATS + limma DSE analysis
├── 04_rna_editing.R                      # RNA editing analysis
├── 05_multilayer_integration.R           # Multi-omics convergence
├── 06_GO_enrichment.R                   # GO BP enrichment
├── 07_functional_separation.R           # Fisher test, functional separation
├── 08_KEGG_enrichment.R                 # KEGG pathway enrichment
├── 09_editing_splicing_crosstalk.R      # RNA editing-splicing crosstalk
├── Fig1-6.R                             # Main figures 1–6
├── SFigure1-8.R                         # Supplementary figures S1–S8
├── README.md
└── LICENSE
```

## Execution order

Scripts should be run in numerical order (01 → 02 → 03 → 04 → 05 → 06 → 07 → 08 → 09 → Fig1-6.R → SFigure1-8.R). Each script begins with `source("00_project_config.R")`.

## Data

Raw and processed data files are available from [GEO accession / repository URL — to be added upon publication]. Place data files in the `data/` directory and raw input files in `data_raw/` before running the scripts.

## Software requirements

- R ≥ 4.2.0
- Key packages: DESeq2, rMATS, limma, ggplot2, cowplot, pheatmap, dplyr, tidyr, clusterProfiler, GSVA, reshape2, scales, ggrepel, svglite

To install all required packages, run `load_packages()` calls at the top of each script, which automatically install missing packages from CRAN and Bioconductor.

## License

MIT License — see LICENSE file.
