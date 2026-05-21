# 00_project_config.R — global configuration for MSC senescence transcriptome analysis
# Source this file at the beginning of each analysis script: source("00_project_config.R")

# Project root: auto-detect from working directory
project_root <- getwd()
setwd(project_root)

# Directory structure
dir_raw     <- file.path(project_root, "data_raw")
dir_data    <- file.path(project_root, "data")
dir_fig     <- file.path(project_root, "figures")

for (d in c(dir_raw, dir_data, dir_fig)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# Nature-style publication theme
source(file.path(project_root, "theme_bindlab.R"))

# Passage color palette (blue=young → red=senescent)
passage_colors <- c(
  "P2"  = "#4575B4", "P8"  = "#91BFDB",
  "P10" = "#FC8D59", "P12" = "#D73027"
)

# Differential expression colors
de_colors <- c("Up" = "#D73027", "Down" = "#4575B4", "NS" = "grey70")

# Alternative splicing type colors
as_type_colors <- c(
  "SE" = "#4575B4", "A3SS" = "#91BFDB", "A5SS" = "#FEE090",
  "MXE" = "#FC8D59", "RI" = "#D73027"
)

# Heatmaps: blue-white-red for z-score expression
heatmap_colors <- colorRampPalette(c("#4575B4", "white", "#D73027"))(100)

# Generic 8-color discrete palette
discrete_colors <- c("#4575B4", "#D73027", "#FF7F00", "#984EA3",
                     "#4DAF4A", "#A65628", "#F781BF", "#999999")

# Raw data file paths
f_expression <- file.path(dir_raw, "MSC_passage_senescence_study_genes.expression.xls")
f_SE_JCEC    <- file.path(dir_raw, "SE.MATS.JCEC.xls")
f_SE_JC      <- file.path(dir_raw, "SE.MATS.JC.xls")
f_MXE_JC     <- file.path(dir_raw, "MXE.MATS.JC.xls")
f_MXE_JCEC   <- file.path(dir_raw, "MXE.MATS.JCEC.xls")
f_A3SS_JC    <- file.path(dir_raw, "A3SS.MATS.JC.xls")
f_A3SS_JCEC  <- file.path(dir_raw, "A3SS.MATS.JCEC.xls")
f_A5SS_JC    <- file.path(dir_raw, "A5SS.MATS.JC.xls")
f_A5SS_JCEC  <- file.path(dir_raw, "A5SS.MATS.JCEC.xls")
f_RI_JC      <- file.path(dir_raw, "RI.MATS.JC.xls")
f_RI_JCEC    <- file.path(dir_raw, "RI.MATS.JCEC.xls")
f_editing    <- file.path(dir_raw, "editing.xls")
f_editing_type <- file.path(dir_raw, "editing.type.xls")
f_snp_annot  <- file.path(dir_raw, "snp.annot.xls")
f_snp_stat   <- file.path(dir_raw, "snp.annot.stat.xls")

# Sample metadata
sample_table <- data.frame(
  sample_id   = c("P2-1","P2-2", "P8-1","P8-2","P8-3",
                   "P10-1","P10-2","P10-3", "P12-1","P12-2","P12-3"),
  passage     = c("P2","P2", "P8","P8","P8",
                   "P10","P10","P10", "P12","P12","P12"),
  passage_num = c(2,2, 8,8,8, 10,10,10, 12,12,12),
  replicate   = c(1,2, 1,2,3, 1,2,3, 1,2,3),
  stringsAsFactors = FALSE
)
sample_table$passage <- factor(sample_table$passage,
                                levels = c("P2","P8","P10","P12"))

# Utility: save ggplot as PDF + PNG + SVG (300 DPI)
save_fig <- function(plot, filename, w = 8, h = 6, dpi = 300, svg = TRUE) {
  base_path <- file.path(dir_fig, filename)
  ggsave(base_path, plot, width = w, height = h, dpi = dpi, device = cairo_pdf)
  png_path <- sub("\\.pdf$", ".png", base_path)
  ggsave(png_path, plot, width = w, height = h, dpi = dpi)
  if (svg && requireNamespace("svglite", quietly = TRUE)) {
    svg_path <- sub("\\.pdf$", ".svg", base_path)
    ggsave(svg_path, plot, width = w, height = h, device = function(...)
      svglite::svglite(..., fix_text_size = FALSE))
  }
  cat("  saved:", filename, "\n")
}

# Utility: save pheatmap as PDF
save_pheatmap <- function(filename, w = 8, h = 7) {
  pdf(file.path(dir_fig, filename), width = w, height = h)
}

# Utility: save intermediate data as CSV
save_data <- function(df, filename) {
  fpath <- file.path(dir_data, filename)
  write.csv(df, fpath, row.names = FALSE)
  cat("  saved:", filename, "\n")
}

# Utility: install and load R packages (CRAN + Bioconductor)
load_packages <- function(pkgs) {
  for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      tryCatch(
        install.packages(pkg, repos = "https://cloud.r-project.org", quiet = TRUE),
        error = function(e) {
          if (!requireNamespace("BiocManager", quietly = TRUE))
            install.packages("BiocManager", repos = "https://cloud.r-project.org")
          BiocManager::install(pkg, ask = FALSE, update = FALSE)
        }
      )
    }
    library(pkg, character.only = TRUE)
  }
}

# Print project summary
cat("====================================================\n")
cat("  MSC Passage Senescence Transcriptome Project\n")
cat("====================================================\n")
cat("  Root    :", project_root, "\n")
cat("  Data    :", dir_data, "\n")
cat("  Figures :", dir_fig, "\n")
cat("  Samples :", nrow(sample_table),
    "(", paste(levels(sample_table$passage), collapse = " -> "), ")\n")
cat("====================================================\n\n")
