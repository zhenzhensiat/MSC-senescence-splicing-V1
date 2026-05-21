# 11_GO_enrichment_local.R — Local GO biological process enrichment analysis (offline)

suppressPackageStartupMessages(library(data.table))

# ---- 0. Configuration ----
PROJECT_ROOT <- getwd()
CACHE_DIR    <- file.path(PROJECT_ROOT, "stringdb_cache")
DATA_DIR     <- file.path(PROJECT_ROOT, "data")
OUT_DIR      <- DATA_DIR

GO_OBO       <- file.path(CACHE_DIR, "go-basic.obo")
GO_GAF       <- file.path(CACHE_DIR, "goa_human.gaf.gz")
LAYER_CSV    <- file.path(DATA_DIR, "05_gene_layer_categories.csv")
EXPR_FILE    <- file.path(DATA_DIR, "02_VST_expression.csv")

PADJ_THRESH      <- 0.05
MIN_GENESET_SIZE <- 3
MAX_GENESET_SIZE <- 500
MIN_OVERLAP      <- 3

#' ---- 1. Parse GO OBO ----
cat("============================================================\n")
cat("  Local GO BP Enrichment (Fully Offline) — R Edition\n")
cat("============================================================\n\n")

cat("[1] Parsing GO OBO:", GO_OBO, "\n")
obo_lines <- readLines(GO_OBO)

# State machine: [Term] ... blank line => one term
parse_obo <- function(lines) {
  terms <- list()
  current <- NULL
  i <- 1
  pb <- txtProgressBar(style = 3)
  for (line in lines) {
    if ((i %% 10000) == 0) setTxtProgressBar(pb, i / length(lines))
    i <- i + 1
    line <- trimws(line)
    if (line == "[Term]") {
      current <- list(
        id = NA_character_, name = NA_character_,
        namespace = NA_character_, parents = character(0),
        alt_ids = character(0)
      )
    } else if (line == "" && !is.null(current)) {
      if (!is.na(current$id)) {
        terms[[current$id]] <- current
        for (alt in current$alt_ids) {
          terms[[alt]] <- current
        }
      }
      current <- NULL
    } else if (!is.null(current)) {
      if (startsWith(line, "id: ")) {
        current$id <- substring(line, 5)
      } else if (startsWith(line, "name: ")) {
        current$name <- substring(line, 7)
      } else if (startsWith(line, "namespace: ")) {
        current$namespace <- substring(line, 12)
      } else if (startsWith(line, "is_a: ")) {
        parent <- strsplit(substring(line, 7), " !")[[1]][1]
        current$parents <- c(current$parents, trimws(parent))
      } else if (startsWith(line, "alt_id: ")) {
        current$alt_ids <- c(current$alt_ids, substring(line, 9))
      }
    }
  }
  close(pb)
  cat(sprintf("  Parsed %d terms\n", length(terms)))
  return(terms)
}

go_terms <- parse_obo(obo_lines)

#' ---- 2. Parse GOA Human GAF ----
cat(sprintf("\n[2] Parsing GOA GAF: %s ...\n", GO_GAF))
con <- gzfile(GO_GAF, "rt")
gaf_lines <- readLines(con)
close(con)

# Build symbol -> GO IDs and GO ID -> symbols
symbol_to_go <- new.env(hash = TRUE, parent = emptyenv())
go_to_symbols <- new.env(hash = TRUE, parent = emptyenv())

total_lines <- 0L
annotated_bp  <- 0L
for (line in gaf_lines) {
  if (nchar(line) == 0 || substr(line, 1, 1) == "!") next
  total_lines <- total_lines + 1
  cols <- strsplit(line, "\t")[[1]]
  if (length(cols) < 5) next

  symbol  <- cols[3]
  go_id   <- cols[5]
  # Fix GO:GO: double prefix
  if (startsWith(go_id, "GO:GO:")) go_id <- substring(go_id, 4)
  qualifier <- if (length(cols) >= 4) cols[4] else ""

  # Skip NOT qualifier annotations
  if (grepl("NOT", qualifier)) next
  # Skip outdated terms
  if (is.null(go_terms[[go_id]])) next
  # Only biological_process
  if (go_terms[[go_id]]$namespace != "biological_process") next

  # Add to symbol->GO
  existing <- symbol_to_go[[symbol]]
  if (is.null(existing)) {
    symbol_to_go[[symbol]] <- go_id
  } else {
    symbol_to_go[[symbol]] <- unique(c(existing, go_id))
  }

  # Add to GO->symbols
  existing2 <- go_to_symbols[[go_id]]
  if (is.null(existing2)) {
    go_to_symbols[[go_id]] <- symbol
  } else {
    go_to_symbols[[go_id]] <- unique(c(existing2, symbol))
  }
  annotated_bp <- annotated_bp + 1
}

cat(sprintf("  Processed %d lines\n", total_lines))
cat(sprintf("  %d genes annotated to BP\n", length(ls(symbol_to_go))))
cat(sprintf("  %d unique GO BP terms\n", length(ls(go_to_symbols))))

#' ---- 3. Load background gene set ----
cat("\n[3] Loading background gene set...\n")
expr <- fread(EXPR_FILE, select = "symbol")
bg_genes_all <- unique(expr$symbol[expr$symbol != "" & !is.na(expr$symbol)])
cat(sprintf("  Total detected genes: %d\n", length(bg_genes_all)))

# Universe = detected genes that have any GO BP annotation
go_annotated_genes <- ls(symbol_to_go)
universe <- intersect(bg_genes_all, go_annotated_genes)
cat(sprintf("  Universe (detected + GO BP annotated): %d\n", length(universe)))

#' ---- 4. Load gene lists ----
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
  cat(sprintf("  Overlap genes: %s\n", paste(sort(overlap), collapse = ",")))
}

#' ---- 5. Hypergeometric enrichment ----
cat("\n[5] Running hypergeometric enrichment...\n")

#' Log-space hypergeometric p-value (avoids overflow with large factorials)
#' phyper uses C code internally which is numerically stable, so we use it directly

do_go_enrichment <- function(gene_list, list_name, universe, go_to_symbols, go_terms) {
  cat(sprintf("\n  [%s] Enrichment:\n", list_name))

  # Restrict to genes in universe
  gene_set <- intersect(gene_list, universe)
  n_annotated <- length(gene_set)
  n_total     <- length(gene_list)
  cat(sprintf("    %d/%d genes with GO BP annotation\n", n_annotated, n_total))

  if (n_annotated < MIN_OVERLAP) {
    cat("    Too few annotated genes, skipping.\n")
    return(data.frame())
  }

  N <- length(universe)   # total background
  n <- n_annotated         # size of gene set in universe

  go_ids <- ls(go_to_symbols)
  total_terms <- length(go_ids)

  tested <- 0L
  results_list <- list()

  for (goid in go_ids) {
    pw_genes_all <- go_to_symbols[[goid]]
    pw_genes     <- intersect(pw_genes_all, universe)
    M <- length(pw_genes)
    if (M < MIN_GENESET_SIZE || M > MAX_GENESET_SIZE) next

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

    term_info <- go_terms[[goid]]
    term_name <- if (!is.null(term_info)) term_info$name else goid

    results_list[[length(results_list) + 1]] <- data.frame(
      go_id       = goid,
      term        = term_name,
      set_size    = M,
      overlap     = k,
      overlap_genes = paste(sort(overlap_genes), collapse = ";"),
      enrichment  = round(enrichment, 2),
      pval        = pval,
      stringsAsFactors = FALSE
    )
  }

  cat(sprintf("    Terms tested: %d\n", tested))

  if (length(results_list) == 0) {
    cat("    No results with sufficient overlap.\n")
    return(data.frame())
  }

  results <- rbindlist(results_list)
  results$padj <- p.adjust(results$pval, method = "BH")
  results <- results[order(results$pval)]

  n_sig <- sum(results$padj < PADJ_THRESH, na.rm = TRUE)
  cat(sprintf("    Significant (padj < %.2f): %d\n", PADJ_THRESH, n_sig))

  return(results)
}

# Run for each gene set
res_deg     <- do_go_enrichment(deg_only, "DEG-only", universe, go_to_symbols, go_terms)
res_dse     <- do_go_enrichment(dse_only, "DSE-only", universe, go_to_symbols, go_terms)
res_overlap <- do_go_enrichment(overlap,  "Overlap",  universe, go_to_symbols, go_terms)

#' ---- 6. Category classification ----
categorize_term <- function(term_name, go_id) {
  t <- tolower(term_name)
  if (grepl("immune|inflammat|defense|virus|interferon|cytokine|t cell|b cell|chemokine|complement|toll|nfr?[-\u03ba]b|nf kappa", t))
    return("Immune/Inflammatory")
  if (grepl("development|differential|morphogen|organogen|embryo|pattern|limb|axis|fate|specification|mesenchyme|osteoblast|chondro|adipogen|neurogen|angiogen", t))
    return("Development/Differentiation")
  if (grepl("signal|kinase|phosphat|receptor|akt|pi3k|mapk|erk|jnk|wnt|notch|hedgehog|tgf|bmp|gpcr|cAMP|calcium|ras|rho|stat|smad|catenin", t))
    return("Signal Transduction")
  if (grepl("adhesion|integrin|collagen|fibri|extracellular|matrix|ecm|basement|migration|motility|chemo", t))
    return("ECM/Cell Adhesion")
  if (grepl("cell cycle|dna repair|mitosis|meiosis|proliferation|replication|apoptosis|senescen|p53|checkpoint|telomere|chromatin", t))
    return("Cell Cycle/DNA Repair")
  if (grepl("metabol|glycolysis|oxidative|lipid|fatty acid|amino acid|glucose|respiration|biosynthetic|catabolic|mitochondri", t))
    return("Metabolism")
  return("Other")
}

# ---- 7. Save results ----
cat("\n[6] Saving results...\n")

save_go_csv <- function(df, fname, list_label) {
  if (nrow(df) == 0) {
    cat(sprintf("  %s: no results\n", fname))
    return()
  }
  fpath <- file.path(OUT_DIR, fname)
  fwrite(df, fpath)
  cat(sprintf("  Saved %s: %d terms, %d significant\n",
              fname, nrow(df), sum(df$padj < PADJ_THRESH, na.rm = TRUE)))
}

save_go_csv(res_deg,     "11_GO_BP_DEG_only.csv", "DEG-only")
save_go_csv(res_dse,     "11_GO_BP_DSE_only.csv", "DSE-only")
if (nrow(res_overlap) > 0) save_go_csv(res_overlap, "11_GO_BP_overlap.csv", "Overlap")

# Category comparison for significant terms
# Combine DEG significant terms with categories
if (nrow(res_deg) > 0) {
  sig_deg <- res_deg[padj < PADJ_THRESH]
  if (nrow(sig_deg) > 0) {
    sig_deg$category <- sapply(sig_deg$term, categorize_term)
    cat_tab <- table(sig_deg$category)
    cat_df <- data.frame(
      category = names(cat_tab),
      count    = as.integer(cat_tab)
    )
    cat_df <- cat_df[order(-cat_df$count), ]
    fwrite(cat_df, file.path(OUT_DIR, "11_GO_category_comparison.csv"))
    cat(sprintf("  Saved 11_GO_category_comparison.csv: %d categories\n", nrow(cat_df)))
  }
}

# Combined top terms (top 10 per gene set)
combine_top <- function(df1, name1, df2, name2, n = 10) {
  parts <- list()
  if (nrow(df1) > 0) {
    d1 <- df1[1:min(n, nrow(df1))]
    d1$gene_set <- name1
    parts[[1]] <- d1
  }
  if (nrow(df2) > 0) {
    d2 <- df2[1:min(n, nrow(df2))]
    d2$gene_set <- name2
    parts[[2]] <- d2
  }
  if (length(parts) > 0) rbindlist(parts) else NULL
}

combined <- combine_top(res_deg, "DEG-only", res_dse, "DSE-only")
if (!is.null(combined)) {
  fwrite(combined, file.path(OUT_DIR, "11_GO_combined_top.csv"))
  cat("  Saved 11_GO_combined_top.csv\n")
}

#' ---- 8. Print top results ----
cat("\n=== Top GO BP terms ===\n")

for (set_name in c("DEG-only", "DSE-only", "Overlap")) {
  res <- switch(set_name,
    "DEG-only" = res_deg,
    "DSE-only" = res_dse,
    "Overlap"  = res_overlap
  )
  if (nrow(res) == 0 || all(is.na(res$padj) | res$padj >= PADJ_THRESH)) {
    cat(sprintf("\n  [%s] No significant terms\n", set_name))
    next
  }
  sig <- res[padj < PADJ_THRESH]
  n_show <- min(15, nrow(sig))
  cat(sprintf("\n  [%s] Top %d significant terms:\n", set_name, n_show))
  for (i in seq_len(n_show)) {
    r <- sig[i]
    cat(sprintf("    %s: %d genes, fold=%.2f, padj=%.2e\n",
                r$term, r$overlap, r$enrichment, r$padj))
  }
}

cat("\n=== GO Enrichment Complete ===\n")
cat(sprintf("Output files in %s:\n", OUT_DIR))
cat("  11_GO_BP_DEG_only.csv\n")
cat("  11_GO_BP_DSE_only.csv\n")
cat("  11_GO_category_comparison.csv\n")
cat("  11_GO_combined_top.csv\n")
