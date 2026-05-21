# 12_KEGG_enrichment_local.R — Local KEGG pathway enrichment analysis (offline)

suppressPackageStartupMessages(library(data.table))

# ---- 0. Configuration ----
PROJECT_ROOT <- getwd()
CACHE_DIR    <- file.path(PROJECT_ROOT, "stringdb_cache")
DATA_DIR     <- file.path(PROJECT_ROOT, "data")
OUT_DIR      <- DATA_DIR

GMT_FILE   <- file.path(CACHE_DIR, "c2.all.v2026.1.Hs.symbols.gmt")
LAYER_CSV  <- file.path(DATA_DIR, "05_gene_layer_categories.csv")
EXPR_FILE  <- file.path(DATA_DIR, "02_VST_expression.csv")

PADJ_THRESH      <- 0.05
MIN_GENESET_SIZE <- 5
MAX_GENESET_SIZE <- 500
MIN_OVERLAP      <- 3

# ---- 1. Parse MSigDB GMT, extract KEGG pathways ----
cat("============================================================\n")
cat("  Local KEGG Pathway Enrichment (Fully Offline)\n")
cat("============================================================\n\n")

cat("[1] Parsing MSigDB GMT:", GMT_FILE, "\n")
gmt_lines <- readLines(GMT_FILE)
cat(sprintf("  Total gene sets: %d\n", length(gmt_lines)))

# Parse GMT lines into list
parse_gmt <- function(lines) {
  gs_list <- list()
  for (line in lines) {
    parts <- strsplit(line, "\t")[[1]]
    if (length(parts) < 3) next
    gs_id    <- parts[1]
    gs_name  <- parts[2]
    gs_genes <- unique(parts[-(1:2)])
    gs_genes <- gs_genes[nchar(gs_genes) > 0]
    gs_list[[gs_id]] <- list(name = gs_name, genes = gs_genes)
  }
  return(gs_list)
}

all_gs <- parse_gmt(gmt_lines)

# Filter to KEGG pathways only, exclude MEDICUS (environmental factors)
kegg_gs <- all_gs[grepl("^KEGG_", names(all_gs))]
kegg_gs <- kegg_gs[!grepl("KEGG_MEDICUS", names(kegg_gs))]
cat(sprintf("  Core KEGG pathways: %d\n", length(kegg_gs)))

# Size filter
orig_n  <- length(kegg_gs)
kegg_gs <- Filter(function(x) {
  n <- length(x$genes)
  n >= MIN_GENESET_SIZE && n <= MAX_GENESET_SIZE
}, kegg_gs)
cat(sprintf("  After size filter [%d, %d]: %d\n", MIN_GENESET_SIZE, MAX_GENESET_SIZE, length(kegg_gs)))
cat(sprintf("  Removed %d pathways outside size range\n", orig_n - length(kegg_gs)))

# ---- 2. Build pathway gene universe ----
cat("\n[2] Building gene universe from KEGG pathways...\n")
all_kegg_genes <- unique(unlist(lapply(kegg_gs, `[[`, "genes")))
cat(sprintf("  Unique genes across all KEGG pathways: %d\n", length(all_kegg_genes)))

# ---- 3. Load background from expression data ----
cat("\n[3] Loading background gene set from expression matrix...\n")
expr <- fread(EXPR_FILE, select = "symbol")
bg_genes_all <- expr[[1]]
bg_genes_all <- unique(bg_genes_all[bg_genes_all != "" & !is.na(bg_genes_all)])
cat(sprintf("  Total detected genes: %d\n", length(bg_genes_all)))

# Universe = detected genes that are in KEGG pathways
universe <- intersect(bg_genes_all, all_kegg_genes)
cat(sprintf("  Universe (detected + in KEGG): %d\n", length(universe)))

# ---- 4. Load gene lists ----
cat("\n[4] Loading gene lists...\n")
layer_df <- fread(LAYER_CSV)

deg_all <- layer_df[layer %in% c("DEG only", "DEG & DSE"), unique(symbol)]
dse_all <- layer_df[layer %in% c("DSE only", "DEG & DSE"), unique(symbol)]
overlap  <- intersect(deg_all, dse_all)

deg_only <- setdiff(deg_all, overlap)
dse_only <- setdiff(dse_all, overlap)

cat(sprintf("  DEG all: %d | DSE all: %d | Overlap: %d\n",
            length(deg_all), length(dse_all), length(overlap)))
cat(sprintf("  DEG-only: %d | DSE-only: %d\n",
            length(deg_only), length(dse_only)))

if (length(overlap) > 0) {
  cat(sprintf("  Overlap genes: %s\n", paste(overlap, collapse = ", ")))
}

# ---- 5. Hypergeometric enrichment ----
cat("\n[5] Running hypergeometric enrichment...\n")

do_enrichment <- function(gene_list, list_name, universe, kegg_gs) {
  cat(sprintf("\n  [%s] Enrichment:\n", list_name))

  # Restrict to genes in universe
  gene_set <- intersect(gene_list, universe)
  n_annotated <- length(gene_set)
  n_total     <- length(gene_list)
  cat(sprintf("    %d/%d genes in universe\n", n_annotated, n_total))

  if (n_annotated < MIN_OVERLAP) {
    cat("    Too few annotated genes, skipping.\n")
    return(data.frame())
  }

  N <- length(universe)   # total background
  n <- n_annotated         # size of our gene set

  results <- data.frame(
    pathway_id   = character(0),
    description  = character(0),
    set_size     = integer(0),
    overlap      = integer(0),
    overlap_genes = character(0),
    enrichment   = numeric(0),
    pval         = numeric(0),
    padj         = numeric(0),
    stringsAsFactors = FALSE
  )

  tested <- 0
  for (pid in names(kegg_gs)) {
    pw_genes_all <- kegg_gs[[pid]]$genes
    pw_genes     <- intersect(pw_genes_all, universe)
    M <- length(pw_genes)
    if (M < 2) next

    # Overlap between pathway and our gene set
    overlap_genes <- intersect(pw_genes, gene_set)
    k <- length(overlap_genes)
    if (k < MIN_OVERLAP) next

    tested <- tested + 1

    # Hypergeometric test: P(X >= k)
    pval <- phyper(k - 1, M, N - M, n, lower.tail = FALSE)

    # Enrichment fold = observed / expected
    expected <- M * n / N
    enrichment <- if (expected > 0) k / expected else 0

    results <- rbind(results, data.frame(
      pathway_id    = pid,
      description   = kegg_gs[[pid]]$name,
      set_size      = M,
      overlap       = k,
      overlap_genes = paste(sort(overlap_genes), collapse = ";"),
      enrichment    = round(enrichment, 2),
      pval          = pval,
      padj          = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  cat(sprintf("    Tested: %d pathways\n", tested))

  if (nrow(results) == 0) {
    cat("    No results with sufficient overlap.\n")
    return(results)
  }

  # BH correction
  results$padj <- p.adjust(results$pval, method = "BH")
  results <- results[order(results$pval), ]

  n_sig <- sum(results$padj < PADJ_THRESH, na.rm = TRUE)
  cat(sprintf("    Significant (padj < %.2f): %d\n", PADJ_THRESH, n_sig))

  return(results)
}

# Run for each gene set
res_deg     <- do_enrichment(deg_only, "DEG-only", universe, kegg_gs)
res_dse     <- do_enrichment(dse_only, "DSE-only", universe, kegg_gs)
res_overlap <- do_enrichment(overlap,  "Overlap",  universe, kegg_gs)

# ---- 6. Save results ----
cat("\n[6] Saving results...\n")

save_csv <- function(df, fname, list_label) {
  fpath <- file.path(OUT_DIR, fname)
  fwrite(df, fpath)
  cat(sprintf("  Saved %s: %d terms, %d significant\n",
              fname, nrow(df), sum(df$padj < PADJ_THRESH, na.rm = TRUE)))
}

if (nrow(res_deg) > 0)     save_csv(res_deg,     "12_KEGG_DEG_only.csv", "DEG-only")
if (nrow(res_dse) > 0)     save_csv(res_dse,     "12_KEGG_DSE_only.csv", "DSE-only")
if (nrow(res_overlap) > 0) save_csv(res_overlap, "12_KEGG_overlap.csv",  "Overlap")

# ---- 7. Print top results ----
cat("\n=== Top KEGG Pathways ===\n")

for (set_name in c("DEG-only", "DSE-only", "Overlap")) {
  res <- switch(set_name,
    "DEG-only" = res_deg,
    "DSE-only" = res_dse,
    "Overlap"  = res_overlap
  )
  if (nrow(res) == 0) {
    cat(sprintf("\n  [%s] No results\n", set_name))
    next
  }
  n_show <- min(10, nrow(res))
  sig <- res[res$padj < PADJ_THRESH, ]
  cat(sprintf("\n  [%s] Top %d:\n", set_name, n_show))
  for (i in seq_len(n_show)) {
    r <- res[i, ]
    sig_mark <- if (!is.na(r$padj) && r$padj < PADJ_THRESH) "*" else " "
    # Strip URL prefix from description
    desc <- sub("^https?://.*geneset/", "", r$description)
    desc <- gsub("_", " ", r$pathway_id)
    cat(sprintf("    %s %s: %d genes, fold=%.2f, padj=%.2e %s\n",
                sig_mark, desc, r$overlap, r$enrichment, r$padj,
                paste0("(", substr(r$overlap_genes, 1, 80), "...)")))
  }
}

cat("\n=== KEGG Enrichment Complete ===\n")
cat(sprintf("Output files in %s:\n", OUT_DIR))
cat("  12_KEGG_DEG_only.csv\n")
cat("  12_KEGG_DSE_only.csv\n")
cat("  12_KEGG_overlap.csv\n")
