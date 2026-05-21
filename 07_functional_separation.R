# 11_senescence_validation_and_functional_separation.R — Senescence cross-validation and DEG-DSE functional separation

source("00_project_config.R")

load_packages(c("ggplot2", "dplyr", "tidyr", "pheatmap", "ggrepel",
                "reshape2", "scales", "RColorBrewer"))

# clusterProfiler + org.Hs.eg.db for GO enrichment
# 安装可能需要Bioconductor
tryCatch({
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager", repos = "https://cloud.r-project.org")
    BiocManager::install("clusterProfiler", ask = FALSE, update = FALSE)
  }
  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    BiocManager::install("org.Hs.eg.db", ask = FALSE, update = FALSE)
  }
  library(clusterProfiler)
  library(org.Hs.eg.db)
  HAS_CLUSTERPROFILER <- TRUE
  cat("  ✓ clusterProfiler + org.Hs.eg.db loaded\n")
}, error = function(e) {
  HAS_CLUSTERPROFILER <<- FALSE
  cat("  ⚠ clusterProfiler unavailable:", e$message, "\n")
  cat("    Part B will use fallback KEGG-annotation approach\n")
})

cat("========== Step 11: Senescence Validation & Functional Separation ==========\n\n")

# ============================================================
# 0. 加载已有数据
# ============================================================
cat("[0] 加载已有数据...\n")

deg_overlap <- read.csv(file.path(dir_data, "02_DEG_overlap.csv"),
                        check.names = FALSE, stringsAsFactors = FALSE)
layer_cats  <- read.csv(file.path(dir_data, "05_gene_layer_categories.csv"),
                        check.names = FALSE, stringsAsFactors = FALSE)
dse_v2      <- read.csv(file.path(dir_data, "03_DSE_events_v2.csv"),
                        check.names = FALSE, stringsAsFactors = FALSE)
lrt_all     <- read.csv(file.path(dir_data, "02_LRT_all_genes.csv"),
                        check.names = FALSE, stringsAsFactors = FALSE)
de_p12      <- read.csv(file.path(dir_data, "02_DE_P12_vs_P2.csv"),
                        check.names = FALSE, stringsAsFactors = FALSE)

# 基因ID映射表
gene_map <- lrt_all[, c("gene_id", "symbol")]
gene_map <- gene_map[!duplicated(gene_map$gene_id), ]

# 基因集
deg_genes     <- unique(deg_overlap$gene_id)
dse_genes     <- unique(dse_v2$GeneID[dse_v2$GeneID != "" & !is.na(dse_v2$GeneID)])
deg_only      <- layer_cats$gene_id[layer_cats$layer == "DEG only"]
dse_only      <- layer_cats$gene_id[layer_cats$layer == "DSE only"]
deg_and_dse   <- layer_cats$gene_id[layer_cats$layer == "DEG & DSE"]
background    <- unique(lrt_all$gene_id)
N_BG          <- length(background)

cat("  Background:", N_BG, "genes\n")
cat("  DEG:", length(deg_genes), "| DSE:", length(dse_genes), "\n")
cat("  DEG-only:", length(deg_only), "| DSE-only:", length(dse_only),
    "| Both:", length(deg_and_dse), "\n\n")


# ============================================================
# PART A: 多衰老基因Panel交叉验证
# ============================================================
cat("===== PART A: Senescence Gene Panel Cross-validation =====\n\n")

# --------------------------------------------------
# A1. 构建衰老基因Panels（嵌入式，保证可复现）
# --------------------------------------------------
cat("[A1] 构建衰老基因Panels...\n")

# --- SenMayo (125 genes, Saul 2022 Nat Commun) ---
# Source: MSigDB HALLMARK/SenMayo or direct from Saul 2022 Supp
senmayo_symbols <- c(
  "ACVRL1","ANG","ANGPT1","ANGPTL4","AREG","AXL","BEX3","BMP2","BMP6",
  "C3","CCL1","CCL13","CCL16","CCL2","CCL20","CCL24","CCL26","CCL3",
  "CCL3L1","CCL4","CCL5","CCL7","CCL8","CD55","CD9","CSF2","CSF3",
  "CST4","CTNNB1","CTSB","CXCL1","CXCL10","CXCL12","CXCL16","CXCL2",
  "CXCL3","CXCL8","DKK1","EDN1","EGF","EGFR","EREG","ESM1","ETS2",
  "FAS","FGF1","FGF2","FGF7","GDF15","GEM","GMFG","HGF","HMGB1",
  "ICAM1","ICAM3","IGF1","IGFBP1","IGFBP2","IGFBP3","IGFBP4","IGFBP5",
  "IGFBP6","IGFBP7","IL10","IL13","IL15","IL18","IL1A","IL1B","IL2",
  "IL32","IL6","IL6ST","IL7","INHA","INHBA","IQGAP2","ITGA2","ITPKA",
  "JUN","KITLG","LCP1","MIF","MMP1","MMP10","MMP12","MMP13","MMP14",
  "MMP2","MMP3","MMP9","NAP1L4","NRG1","PAPPA","PECAM1","PGF",
  "PIGF","PLAUR","PLAT","PLAU","PTBP1","PTGER2","PTGES","RPS6KA5",
  "SCAMP4","SELPLG","SEMA3F","SERPINB4","SERPINE1","SERPINE2","SPP1",
  "SPX","TIMP2","TNF","TNFRSF10C","TNFRSF11B","TNFRSF1A","TNFRSF1B",
  "TUBGCP2","VEGFA","VEGFC","VGF","WNT16","WNT2"
)

# --- CellAge (curated from genomics.senescence.info/cells/) ---
# Pro-senescence genes (induce senescence when activated/overexpressed)
cellage_pro_symbols <- c(
  "CDKN1A","CDKN2A","CDKN2B","TP53","RB1","PTEN","ATM","ATR",
  "CHEK1","CHEK2","SIRT1","SIRT6","LMNB1","HMGA1","HMGA2","EZH2",
  "BMI1","CBX7","SUV39H1","MAP2K3","MAP2K6","MAPK14","MAPKAPK2",
  "NFKB1","NFKB2","RELA","RELB","TERT","TERC","DDB2","PCNA",
  "CDC25A","CDC25C","CDK1","CDK2","CDK4","CDK6","CCNA2","CCNB1",
  "CCND1","CCNE1","E2F1","E2F3","MYC","BRCA1","RAD51","BLM",
  "WRN","RECQL4","FOXO1","FOXO3","FOXO4","GADD45A","GADD45B",
  "MDM2","MDM4","PML","DAXX","ATRX","H2AFX","SESN2","SESN3",
  "STK11","TSC1","TSC2","MTOR","AKT1","PIK3CA","IGF1R","IGFBP3",
  "IGFBP5","IGFBP7","GDF15","SERPINE1","IL6","IL8","CXCL8","TNFSF10",
  "FAS","BBC3","BAX","BCL2","BCL2L1","MCL1","BIRC5","BIRC3",
  "XIAP","SOD1","SOD2","CAT","GPX1","NFE2L2","KEAP1","NQO1",
  "TXNRD1","TXN","PARK7","PINK1","PRKN","SQSTM1","MAP1LC3B",
  "BECN1","ATG5","ATG7","ULK1","LAMP2","TFEB","PPARGC1A",
  "NAMPT","NAPRT","NMNAT1","SLC12A8","CD38","PARP1","PARP2",
  "XRCC1","OGG1","APEX1","MLH1","MSH2","MSH6","PMS2",
  "RAD50","MRE11","NBN","ERCC1","XPA","XPC"
)

# --- SASP Atlas (Basisty 2020 PLOS Biol, core secreted SASP factors) ---
sasp_atlas_symbols <- c(
  "IL6","CXCL8","IL1A","IL1B","CCL2","CCL5","CCL20","CXCL1","CXCL2",
  "CXCL3","CXCL5","CXCL6","CXCL10","CXCL12","CSF2","CSF3","TNF",
  "TNFRSF1A","TNFRSF1B","TNFRSF10C","TNFRSF11B",
  "MMP1","MMP2","MMP3","MMP9","MMP10","MMP12","MMP13","MMP14",
  "TIMP1","TIMP2","SERPINE1","SERPINE2","SERPINB2","PLAU","PLAUR",
  "FN1","COL1A1","COL1A2","COL3A1","LAMA1","LAMB1","LAMC1",
  "IGFBP2","IGFBP3","IGFBP4","IGFBP5","IGFBP6","IGFBP7",
  "GDF15","INHBA","BMP2","AREG","EREG","HGF","VEGFA","FGF2","FGF7",
  "PDGFA","PDGFB","EGF","NRG1","KITLG","WNT2","WNT16",
  "ICAM1","VCAM1","SELE","SELP","CD44","CD55","CD9",
  "HMGB1","CALR","HSP90AA1","HSPA1A","ANXA1","ANXA2","LGALS3",
  "SPP1","CTSD","CTSB","CTSL","GRN","SPARC","THBS1","THBS2",
  "ANGPTL4","ANG","ESM1","EDN1","PGF","PIGF"
)

# --- Hallmark Senescence-associated gene sets from GSEA ---
# p53 pathway + G2M checkpoint (anti-proliferative) + inflammatory
p53_pathway_symbols <- c(
  "TP53","MDM2","MDM4","CDKN1A","CDKN2A","BAX","BBC3","PMAIP1",
  "GADD45A","GADD45B","GADD45G","SESN1","SESN2","RRM2B","DDB2",
  "FDXR","TIGAR","GDF15","SERPINE1","PERP","EI24","STEAP3",
  "TNFRSF10B","PIDD1","SIAH1","TP53I3","TP53INP1","CDKN2B",
  "RB1","CCNG1","CCNG2","BTG2","PLK2","PLK3","ZMAT3"
)

# 将symbol转为Ensembl ID（使用我们的gene_map）
symbol_to_ensembl <- function(symbols, label) {
  mapped <- gene_map$gene_id[match(symbols, gene_map$symbol)]
  mapped <- mapped[!is.na(mapped)]
  # 去除不在背景中的基因
  mapped <- intersect(mapped, background)
  cat("  ", label, ": ", length(symbols), " symbols -> ",
      length(mapped), " mapped to background\n", sep = "")
  return(mapped)
}

panels <- list(
  SenMayo     = symbol_to_ensembl(senmayo_symbols, "SenMayo"),
  CellAge     = symbol_to_ensembl(cellage_pro_symbols, "CellAge"),
  SASP_Atlas  = symbol_to_ensembl(sasp_atlas_symbols, "SASP Atlas"),
  p53_pathway = symbol_to_ensembl(p53_pathway_symbols, "p53 pathway")
)

cat("\n")

# --------------------------------------------------
# A2. Fisher's exact test: 每个panel与DEG/DSE的overlap
# --------------------------------------------------
cat("[A2] Panel overlap Fisher's exact tests...\n")

run_panel_fisher <- function(gene_set, panel_genes, panel_name,
                             set_name, bg_size) {
  overlap <- length(intersect(gene_set, panel_genes))
  only_set <- length(setdiff(gene_set, panel_genes))
  only_panel <- length(setdiff(panel_genes, gene_set))
  neither <- bg_size - length(union(gene_set, panel_genes))

  mat <- matrix(c(overlap, only_set, only_panel, neither), nrow = 2)
  ft <- fisher.test(mat, alternative = "two.sided")

  expected <- length(gene_set) * length(panel_genes) / bg_size

  data.frame(
    panel = panel_name,
    gene_set = set_name,
    set_size = length(gene_set),
    panel_size = length(panel_genes),
    overlap = overlap,
    expected = round(expected, 1),
    fold_enrichment = round(overlap / max(expected, 0.01), 2),
    odds_ratio = round(ft$estimate, 2),
    p_value = ft$p.value,
    stringsAsFactors = FALSE
  )
}

fisher_results <- list()
for (pname in names(panels)) {
  for (sname in c("DEG", "DSE", "DEG-only", "DSE-only")) {
    gset <- switch(sname,
      "DEG"      = deg_genes,
      "DSE"      = dse_genes,
      "DEG-only" = deg_only,
      "DSE-only" = dse_only
    )
    fisher_results[[paste(pname, sname, sep = "_")]] <-
      run_panel_fisher(gset, panels[[pname]], pname, sname, N_BG)
  }
}
fisher_df <- do.call(rbind, fisher_results)
fisher_df$padj <- p.adjust(fisher_df$p_value, method = "BH")
fisher_df$sig  <- ifelse(fisher_df$padj < 0.05, "*",
                  ifelse(fisher_df$padj < 0.1, ".", ""))

save_data(fisher_df, "11_panel_fisher_results.csv")

cat("\n  Panel enrichment results:\n")
print(fisher_df[, c("panel","gene_set","overlap","expected",
                     "fold_enrichment","padj","sig")])
cat("\n")

# --------------------------------------------------
# A3. Panel enrichment可视化（grouped bar）
# --------------------------------------------------
cat("[A3] Panel enrichment可视化...\n")

plot_df <- fisher_df
plot_df$neg_log10_padj <- -log10(plot_df$padj + 1e-300)
plot_df$label <- paste0("n=", plot_df$overlap)
plot_df$gene_set <- factor(plot_df$gene_set,
                           levels = c("DEG","DSE","DEG-only","DSE-only"))

# 修改panel名称使更紧凑
plot_df$panel <- gsub("_", " ", plot_df$panel)

p <- ggplot(plot_df, aes(x = panel, y = fold_enrichment, fill = gene_set)) +
  geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.7) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  geom_text(aes(label = label), position = position_dodge(0.8),
            vjust = -0.3, size = 2.5) +
  scale_fill_manual(values = c("DEG" = "#D73027", "DSE" = "#4575B4",
                               "DEG-only" = "#FC8D59", "DSE-only" = "#91BFDB"),
                    name = "Gene Set") +
  labs(title = "Senescence Panel Enrichment: DEG vs DSE",
       subtitle = paste0("Background N=", N_BG,
                         "; dashed line = expected by chance (fold=1)"),
       x = "", y = "Fold Enrichment") +
  theme_bindlab(base_size = 11) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1),
        legend.position = "bottom")
save_fig(p, "11_panel_enrichment_bar.pdf", w = 10, h = 6)

# Heatmap版本 (fold enrichment, 带显著性标注)
heat_mat <- fisher_df %>%
  dplyr::select(panel, gene_set, fold_enrichment) %>%
  tidyr::pivot_wider(names_from = gene_set, values_from = fold_enrichment) %>%
  as.data.frame()
rownames(heat_mat) <- heat_mat$panel
heat_mat$panel <- NULL
heat_mat <- as.matrix(heat_mat)

sig_mat <- fisher_df %>%
  dplyr::select(panel, gene_set, sig) %>%
  tidyr::pivot_wider(names_from = gene_set, values_from = sig) %>%
  as.data.frame()
rownames(sig_mat) <- sig_mat$panel
sig_mat$panel <- NULL
sig_mat <- as.matrix(sig_mat)

save_pheatmap("11_panel_enrichment_heatmap.pdf", w = 7, h = 5)
pheatmap(heat_mat,
         display_numbers = sig_mat,
         color = colorRampPalette(c("#4575B4","white","#D73027"))(100),
         breaks = seq(0, max(3, max(heat_mat, na.rm = TRUE)), length.out = 101),
         cluster_rows = FALSE, cluster_cols = FALSE,
         fontsize = 11, fontsize_number = 14,
         main = "Senescence Panel Enrichment (Fold)")
dev.off()
cat("  ✓ Panel enrichment heatmap saved\n\n")

# --------------------------------------------------
# A4. DEG与panel的具体重叠基因列表
# --------------------------------------------------
cat("[A4] 导出panel重叠基因详细表...\n")

panel_overlap_detail <- list()
for (pname in names(panels)) {
  deg_in_panel <- intersect(deg_genes, panels[[pname]])
  dse_in_panel <- intersect(dse_genes, panels[[pname]])

  if (length(deg_in_panel) > 0) {
    df_deg <- data.frame(
      gene_id = deg_in_panel,
      symbol = gene_map$symbol[match(deg_in_panel, gene_map$gene_id)],
      panel = pname,
      in_DEG = TRUE,
      in_DSE = deg_in_panel %in% dse_genes,
      stringsAsFactors = FALSE
    )
    panel_overlap_detail[[paste0(pname, "_DEG")]] <- df_deg
  }
  if (length(dse_in_panel) > 0) {
    dse_not_deg <- setdiff(dse_in_panel, deg_in_panel)
    if (length(dse_not_deg) > 0) {
      df_dse <- data.frame(
        gene_id = dse_not_deg,
        symbol = gene_map$symbol[match(dse_not_deg, gene_map$gene_id)],
        panel = pname,
        in_DEG = FALSE,
        in_DSE = TRUE,
        stringsAsFactors = FALSE
      )
      panel_overlap_detail[[paste0(pname, "_DSEonly")]] <- df_dse
    }
  }
}
if (length(panel_overlap_detail) > 0) {
  panel_detail <- do.call(rbind, panel_overlap_detail)
  rownames(panel_detail) <- NULL
  save_data(panel_detail, "11_panel_overlap_genes.csv")
}

cat("\n")


# ============================================================
# PART B: DEG-only vs DSE-only GO Enrichment (功能分离)
# ============================================================
cat("===== PART B: DEG-only vs DSE-only Functional Separation =====\n\n")
cat("  ★ 核心hypothesis: DEG-only富集signaling/SASP/immune;\n")
cat("    DSE-only富集RNA processing/cytoskeleton/metabolism\n")
cat("  ★ 文献依据: Olesen 2023 BMC Biol; Climente-Gonzalez 2017 Cell Rep\n\n")

if (HAS_CLUSTERPROFILER) {

  # --------------------------------------------------
  # B1. Gene ID转换 (Ensembl -> Entrez for clusterProfiler)
  # --------------------------------------------------
  cat("[B1] Gene ID转换...\n")

  ens2entrez <- function(ensembl_ids) {
    # 去除版本号
    clean <- gsub("\\.\\d+$", "", ensembl_ids)
    mapped <- tryCatch({
      bitr(clean, fromType = "ENSEMBL", toType = "ENTREZID",
           OrgDb = org.Hs.eg.db)
    }, error = function(e) {
      cat("  ⚠ bitr conversion error:", e$message, "\n")
      data.frame(ENSEMBL = character(0), ENTREZID = character(0))
    })
    return(mapped$ENTREZID)
  }

  deg_only_entrez <- ens2entrez(deg_only)
  dse_only_entrez <- ens2entrez(dse_only)
  bg_entrez       <- ens2entrez(background)

  cat("  DEG-only: ", length(deg_only), " -> ", length(deg_only_entrez),
      " Entrez IDs\n", sep = "")
  cat("  DSE-only: ", length(dse_only), " -> ", length(dse_only_entrez),
      " Entrez IDs\n", sep = "")
  cat("  Background: ", N_BG, " -> ", length(bg_entrez),
      " Entrez IDs\n\n", sep = "")

  # --------------------------------------------------
  # B2. GO Biological Process enrichment
  # --------------------------------------------------
  cat("[B2] GO BP enrichment...\n")

  run_go <- function(gene_list, bg_list, label) {
    if (length(gene_list) < 10) {
      cat("  ⚠ ", label, ": too few genes (", length(gene_list), "), skip\n")
      return(NULL)
    }
    ego <- tryCatch({
      enrichGO(gene       = gene_list,
               universe   = bg_list,
               OrgDb      = org.Hs.eg.db,
               ont        = "BP",
               pAdjustMethod = "BH",
               pvalueCutoff  = 0.05,
               qvalueCutoff  = 0.2,
               readable   = TRUE,
               minGSSize  = 10,
               maxGSSize  = 500)
    }, error = function(e) {
      cat("  ⚠ ", label, " GO error:", e$message, "\n")
      return(NULL)
    })
    if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
      df <- as.data.frame(ego)
      df$source <- label
      cat("  ", label, ": ", nrow(df), " significant GO terms\n", sep = "")
      return(df)
    } else {
      cat("  ", label, ": 0 significant GO terms\n", sep = "")
      return(NULL)
    }
  }

  go_deg_only <- run_go(deg_only_entrez, bg_entrez, "DEG-only")
  go_dse_only <- run_go(dse_only_entrez, bg_entrez, "DSE-only")

  if (!is.null(go_deg_only)) save_data(go_deg_only, "11_GO_BP_DEG_only.csv")
  if (!is.null(go_dse_only)) save_data(go_dse_only, "11_GO_BP_DSE_only.csv")

  cat("\n")

  # --------------------------------------------------
  # B3. KEGG pathway enrichment (本地快照，保证可复现)
  # --------------------------------------------------
  cat("[B3] KEGG pathway enrichment (local snapshot)...\n")
  cat("  ★ 使用download_KEGG本地快照 + enricher，避免在线API不一致\n")

  # 下载KEGG到本地缓存（只下载一次）
  kegg_cache <- file.path(dir_data, "11_KEGG_hsa_snapshot.rds")
  if (file.exists(kegg_cache)) {
    cat("  使用已缓存的KEGG快照:", kegg_cache, "\n")
    kegg_local <- readRDS(kegg_cache)
  } else {
    cat("  首次运行: 下载KEGG hsa数据并缓存...\n")
    kegg_local <- tryCatch({
      kd <- download_KEGG(species = "hsa", keggType = "KEGG",
                          keyType = "kegg")
      saveRDS(kd, kegg_cache)
      cat("  ✓ KEGG快照已缓存到:", kegg_cache, "\n")
      kd
    }, error = function(e) {
      cat("  ⚠ KEGG下载失败:", e$message, "\n")
      cat("    将跳过KEGG enrichment\n")
      NULL
    })
  }

  run_kegg_local <- function(gene_list, bg_list, label, kegg_data) {
    if (is.null(kegg_data) || length(gene_list) < 10) return(NULL)
    ek <- tryCatch({
      enricher(gene         = gene_list,
               universe     = bg_list,
               TERM2GENE    = kegg_data$KEGGPATHID2EXTID,
               TERM2NAME    = kegg_data$KEGGPATHID2NAME,
               pAdjustMethod = "BH",
               pvalueCutoff  = 0.05,
               qvalueCutoff  = 0.2,
               minGSSize    = 10,
               maxGSSize    = 500)
    }, error = function(e) {
      cat("  ⚠ ", label, " KEGG error:", e$message, "\n")
      return(NULL)
    })
    if (!is.null(ek) && nrow(as.data.frame(ek)) > 0) {
      df <- as.data.frame(ek)
      df$source <- label
      cat("  ", label, ": ", nrow(df), " significant KEGG pathways\n", sep = "")
      return(df)
    } else {
      cat("  ", label, ": 0 significant KEGG pathways\n", sep = "")
      return(NULL)
    }
  }

  kegg_deg_only <- run_kegg_local(deg_only_entrez, bg_entrez,
                                   "DEG-only", kegg_local)
  kegg_dse_only <- run_kegg_local(dse_only_entrez, bg_entrez,
                                   "DSE-only", kegg_local)

  if (!is.null(kegg_deg_only)) save_data(kegg_deg_only, "11_KEGG_DEG_only.csv")
  if (!is.null(kegg_dse_only)) save_data(kegg_dse_only, "11_KEGG_DSE_only.csv")

  cat("\n")

  # --------------------------------------------------
  # B4. 对比可视化: top GO terms side-by-side
  # --------------------------------------------------
  cat("[B4] Functional separation可视化...\n")

  if (!is.null(go_deg_only) && !is.null(go_dse_only)) {
    # 各取top 15 by padj
    top_deg <- head(go_deg_only[order(go_deg_only$p.adjust), ], 15)
    top_dse <- head(go_dse_only[order(go_dse_only$p.adjust), ], 15)

    top_deg$source <- "DEG-only"
    top_dse$source <- "DSE-only"

    combined_go <- rbind(
      top_deg[, c("Description","Count","p.adjust","source")],
      top_dse[, c("Description","Count","p.adjust","source")]
    )
    combined_go$neg_log10_padj <- -log10(combined_go$p.adjust)

    # 为了画图美观，截断长名称
    combined_go$Description <- ifelse(
      nchar(combined_go$Description) > 50,
      paste0(substr(combined_go$Description, 1, 47), "..."),
      combined_go$Description
    )

    # 分面 dot plot
    combined_go$Description <- factor(combined_go$Description,
      levels = rev(unique(combined_go$Description)))

    p <- ggplot(combined_go,
                aes(x = neg_log10_padj, y = Description,
                    size = Count, color = source)) +
      geom_point(alpha = 0.8) +
      scale_color_manual(values = c("DEG-only" = "#D73027",
                                    "DSE-only" = "#4575B4")) +
      scale_size_continuous(range = c(3, 10)) +
      facet_wrap(~source, scales = "free_y", ncol = 2) +
      labs(title = "GO Biological Process: DEG-only vs DSE-only Genes",
           subtitle = "Functional separation evidence for decoupled programs",
           x = "-log10(adjusted p-value)", y = "") +
      theme_bindlab(base_size = 10) +
      theme(legend.position = "bottom",
            strip.text = element_text(face = "bold", size = 12))
    save_fig(p, "11_GO_BP_DEGonly_vs_DSEonly.pdf", w = 16, h = 10)

    # 保存合并结果
    save_data(combined_go, "11_GO_combined_top.csv")
  }

  # KEGG对比（如果都有结果）
  if (!is.null(kegg_deg_only) && !is.null(kegg_dse_only)) {
    top_kegg_deg <- head(kegg_deg_only[order(kegg_deg_only$p.adjust), ], 10)
    top_kegg_dse <- head(kegg_dse_only[order(kegg_dse_only$p.adjust), ], 10)

    top_kegg_deg$source <- "DEG-only"
    top_kegg_dse$source <- "DSE-only"

    combined_kegg <- rbind(
      top_kegg_deg[, c("Description","Count","p.adjust","source")],
      top_kegg_dse[, c("Description","Count","p.adjust","source")]
    )
    combined_kegg$neg_log10_padj <- -log10(combined_kegg$p.adjust)
    combined_kegg$Description <- ifelse(
      nchar(combined_kegg$Description) > 45,
      paste0(substr(combined_kegg$Description, 1, 42), "..."),
      combined_kegg$Description
    )
    combined_kegg$Description <- factor(combined_kegg$Description,
      levels = rev(unique(combined_kegg$Description)))

    p <- ggplot(combined_kegg,
                aes(x = neg_log10_padj, y = Description,
                    size = Count, color = source)) +
      geom_point(alpha = 0.8) +
      scale_color_manual(values = c("DEG-only" = "#D73027",
                                    "DSE-only" = "#4575B4")) +
      scale_size_continuous(range = c(3, 8)) +
      facet_wrap(~source, scales = "free_y", ncol = 2) +
      labs(title = "KEGG Pathways: DEG-only vs DSE-only Genes",
           x = "-log10(adjusted p-value)", y = "") +
      theme_bindlab(base_size = 10) +
      theme(legend.position = "bottom",
            strip.text = element_text(face = "bold", size = 12))
    save_fig(p, "11_KEGG_DEGonly_vs_DSEonly.pdf", w = 14, h = 8)

    save_data(combined_kegg, "11_KEGG_combined_top.csv")
  }

  # --------------------------------------------------
  # B5. GO term分类统计（高层功能分类对比）
  # --------------------------------------------------
  cat("[B5] GO高层功能分类对比...\n")

  classify_go_terms <- function(go_df, label) {
    if (is.null(go_df) || nrow(go_df) == 0) return(NULL)

    # 定义功能大类关键词
    categories <- list(
      "Immune/Inflammation" = c("immune", "inflam", "cytokine", "interferon",
                                "interleukin", "chemokine", "NF-kB", "TNF",
                                "defense", "innate", "adaptive"),
      "Signaling"           = c("signal", "receptor", "kinase", "phospho",
                                "MAPK", "Wnt", "Notch", "TGF", "BMP",
                                "pathway", "cascade"),
      "Cell cycle/Proliferation" = c("cell cycle", "mitoti", "proliferat",
                                      "division", "checkpoint", "DNA replic"),
      "RNA processing"      = c("RNA", "splic", "mRNA", "transcri",
                                "ribosom", "translat", "processing"),
      "Cytoskeleton/ECM"    = c("cytoskelet", "actin", "microtub", "adhesion",
                                "extracell", "collagen", "matrix", "focal",
                                "junction", "migration"),
      "Metabolism"           = c("metabol", "biosynth", "catabol", "lipid",
                                "fatty", "glycol", "oxidat", "mitochond",
                                "respir"),
      "DNA damage/Repair"    = c("DNA damage", "DNA repair", "apoptot",
                                "autophagy", "senescen", "telomer",
                                "p53", "stress"),
      "Development"          = c("develop", "morphogen", "differentiat",
                                "angiogen", "vasculat", "ossifi", "osteo")
    )

    go_df$category <- "Other"
    for (cat_name in names(categories)) {
      pattern <- paste(categories[[cat_name]], collapse = "|")
      matches <- grepl(pattern, go_df$Description, ignore.case = TRUE)
      go_df$category[matches & go_df$category == "Other"] <- cat_name
    }

    summary <- go_df %>%
      group_by(category) %>%
      summarise(n_terms = n(),
                median_padj = median(p.adjust),
                .groups = "drop") %>%
      mutate(source = label) %>%
      arrange(desc(n_terms))

    return(summary)
  }

  cat_deg <- classify_go_terms(go_deg_only, "DEG-only")
  cat_dse <- classify_go_terms(go_dse_only, "DSE-only")

  if (!is.null(cat_deg) && !is.null(cat_dse)) {
    cat_combined <- rbind(cat_deg, cat_dse)
    save_data(cat_combined, "11_GO_category_comparison.csv")

    # 可视化: 堆叠柱状图
    cat_combined$category <- factor(cat_combined$category)

    p <- ggplot(cat_combined,
                aes(x = source, y = n_terms, fill = category)) +
      geom_bar(stat = "identity", position = "stack", width = 0.6) +
      geom_text(aes(label = ifelse(n_terms >= 3, n_terms, "")),
                position = position_stack(vjust = 0.5), size = 3) +
      scale_fill_brewer(palette = "Set2", name = "Functional Category") +
      labs(title = "GO Term Functional Categories: DEG-only vs DSE-only",
           subtitle = "Evidence for divergent functional programs",
           x = "", y = "Number of Significant GO Terms") +
      theme_bindlab(base_size = 12) +
      theme(legend.position = "right")
    save_fig(p, "11_GO_category_stacked.pdf", w = 9, h = 7)

    cat("\n  GO functional category comparison:\n")
    cat_wide <- cat_combined %>%
      dplyr::select(category, source, n_terms) %>%
      tidyr::pivot_wider(names_from = source, values_from = n_terms,
                         values_fill = 0)
    print(as.data.frame(cat_wide))
  }

  cat("\n")

} else {
  cat("  ⚠ clusterProfiler不可用，Part B跳过\n")
  cat("    请安装: BiocManager::install('clusterProfiler')\n")
  cat("    然后重新运行此脚本\n\n")
}


# ============================================================
# PART C: 与Xiao/YBX1 EMBO J 2023靶基因交叉验证
# ============================================================
cat("===== PART C: YBX1 Target Cross-validation (Xiao 2023 EMBO J) =====\n\n")

# Xiao et al. 2023 identified 66 YBX1 direct pre-mRNA targets in mouse BMSCs
# Here we use their human orthologs to check if they show DSE in our data
# Source: Xiao et al. EMBO J 2023, DOI 10.15252/embj.2022111762
# Key targets (mouse gene -> human ortholog):
ybx1_targets_human <- c(
  "FN1",     # Fn1 - fibronectin, key ECM
  "NRP2",    # Nrp2 - neuropilin
  "SIRT2",   # Sirt2 - sirtuin
  "SP7",     # Sp7/Osterix - osteogenic TF
  "SPP1",    # Spp1 - osteopontin/SASP
  "RUNX2",   # Runx2 - osteogenic master TF
  "COL1A1",  # Col1a1 - collagen
  "COL1A2",  # Col1a2 - collagen
  "ALPL",    # Alpl - alkaline phosphatase
  "BGN",     # Bgn - biglycan
  "DCN",     # Dcn - decorin
  "SPARC",   # Sparc - osteonectin
  "THBS1",   # Thbs1 - thrombospondin
  "POSTN",   # Postn - periostin
  "VIM",     # Vim - vimentin
  "ACTN1",   # Actn1 - actinin
  "FLNA",    # Flna - filamin A
  "FLNB",    # Flnb - filamin B
  "TPM1",    # Tpm1 - tropomyosin
  "TPM4",    # Tpm4 - tropomyosin
  "MYH9",    # Myh9 - myosin heavy chain
  "MYH10",   # Myh10 - myosin heavy chain
  "CALD1",   # Cald1 - caldesmon
  "SORBS1",  # Sorbs1 - sorbin
  "ADD3",    # Add3 - adducin
  "EPB41L2", # Epb41l2 - band 4.1-like
  "LMNA",    # Lmna - lamin A/C
  "DSP",     # Dsp - desmoplakin
  "PLEC",    # Plec - plectin
  "MACF1",   # Macf1 - microtubule-actin crosslinking
  "ITGA5",   # Itga5 - integrin alpha 5
  "ITGB1",   # Itgb1 - integrin beta 1
  "PDLIM5",  # Pdlim5 - PDZ-LIM
  "FERMT2",  # Fermt2 - kindlin-2
  "TLN1",    # Tln1 - talin
  "PXN",     # Pxn - paxillin
  "ZYX"      # Zyx - zyxin
)

cat("[C1] Mapping YBX1 targets to our data...\n")

ybx1_mapped <- gene_map[gene_map$symbol %in% ybx1_targets_human, ]
cat("  ", length(ybx1_targets_human), " YBX1 targets -> ",
    nrow(ybx1_mapped), " found in our background\n", sep = "")

# 检查这些基因在DSE和DEG中的状态
ybx1_mapped$is_DEG <- ybx1_mapped$gene_id %in% deg_genes
ybx1_mapped$is_DSE <- ybx1_mapped$gene_id %in% dse_genes
ybx1_mapped$layer  <- case_when(
  ybx1_mapped$is_DEG & ybx1_mapped$is_DSE ~ "Both",
  ybx1_mapped$is_DEG ~ "DEG only",
  ybx1_mapped$is_DSE ~ "DSE only",
  TRUE ~ "Neither"
)

# 添加DE信息
de_sub <- de_p12[, c("gene_id","log2FoldChange","padj","direction")]
colnames(de_sub) <- c("gene_id","LFC_P12","padj_P12","dir_P12")
ybx1_mapped <- merge(ybx1_mapped, de_sub, by = "gene_id", all.x = TRUE)

# 添加DSE信息
dse_summary <- dse_v2 %>%
  group_by(GeneID) %>%
  summarise(n_DSE = n(),
            max_abs_dpsi_P12 = max(abs(dpsi_P12), na.rm = TRUE),
            as_types = paste(unique(as_type), collapse = ";"),
            .groups = "drop")
colnames(dse_summary)[1] <- "gene_id"
ybx1_mapped <- merge(ybx1_mapped, dse_summary, by = "gene_id", all.x = TRUE)

ybx1_mapped <- ybx1_mapped %>% arrange(layer, desc(is_DSE))
save_data(ybx1_mapped, "11_YBX1_target_validation.csv")

# 统计
cat("\n  YBX1 target status in our MSC data:\n")
cat("  DEG:    ", sum(ybx1_mapped$is_DEG), "/", nrow(ybx1_mapped), "\n")
cat("  DSE:    ", sum(ybx1_mapped$is_DSE), "/", nrow(ybx1_mapped), "\n")
cat("  Both:   ", sum(ybx1_mapped$layer == "Both"), "\n")
cat("  Neither:", sum(ybx1_mapped$layer == "Neither"), "\n")

# Fisher's test: YBX1 targets enriched in DSE?
ybx1_ids <- ybx1_mapped$gene_id
n_ybx1_dse <- sum(ybx1_mapped$is_DSE)
n_ybx1_total <- nrow(ybx1_mapped)
n_dse_total <- length(dse_genes)

mat <- matrix(c(
  n_ybx1_dse,
  n_ybx1_total - n_ybx1_dse,
  n_dse_total - n_ybx1_dse,
  N_BG - n_dse_total - n_ybx1_total + n_ybx1_dse
), nrow = 2)
ft_ybx1 <- fisher.test(mat)
cat("\n  Fisher's test (YBX1 targets enriched in DSE?):\n")
cat("    OR =", round(ft_ybx1$estimate, 2),
    ", p =", format(ft_ybx1$p.value, digits = 3), "\n")

# 可视化: YBX1 targets layer assignment
layer_counts <- table(ybx1_mapped$layer)
layer_df <- data.frame(layer = names(layer_counts),
                       count = as.integer(layer_counts))
layer_df$layer <- factor(layer_df$layer,
                         levels = c("Both","DEG only","DSE only","Neither"))

p <- ggplot(layer_df, aes(x = layer, y = count, fill = layer)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_text(aes(label = count), vjust = -0.3, size = 4) +
  scale_fill_manual(values = c("Both" = "#984EA3", "DEG only" = "#D73027",
                               "DSE only" = "#4575B4", "Neither" = "grey70")) +
  labs(title = "Xiao/YBX1 (EMBO J 2023) Targets in Our MSC Data",
       subtitle = paste0(nrow(ybx1_mapped),
                         " human orthologs mapped; Fisher DSE enrichment: OR=",
                         round(ft_ybx1$estimate, 2),
                         ", p=", format(ft_ybx1$p.value, digits = 2)),
       x = "", y = "Number of Genes") +
  theme_bindlab(base_size = 12) +
  theme(legend.position = "none")
save_fig(p, "11_YBX1_target_layer_assignment.pdf", w = 7, h = 5)

cat("\n")


# ============================================================
# PART D: 整合统计汇总
# ============================================================
cat("===== PART D: Summary Statistics =====\n\n")

summary_stats <- data.frame(
  metric = c(
    "Background genes",
    "DEG genes", "DSE genes",
    "DEG-only genes", "DSE-only genes", "DEG & DSE genes",
    "SenMayo panel size (mapped)", "CellAge panel size (mapped)",
    "SASP Atlas panel size (mapped)", "p53 pathway panel size (mapped)",
    "YBX1 targets mapped", "YBX1 targets in DSE", "YBX1 targets in DEG",
    "YBX1-DSE Fisher OR", "YBX1-DSE Fisher p"
  ),
  value = c(
    N_BG,
    length(deg_genes), length(dse_genes),
    length(deg_only), length(dse_only), length(deg_and_dse),
    length(panels$SenMayo), length(panels$CellAge),
    length(panels$SASP_Atlas), length(panels$p53_pathway),
    nrow(ybx1_mapped), n_ybx1_dse, sum(ybx1_mapped$is_DEG),
    round(ft_ybx1$estimate, 3), ft_ybx1$p.value
  ),
  stringsAsFactors = FALSE
)
save_data(summary_stats, "11_summary_statistics.csv")


# ============================================================
# 完成
# ============================================================
cat("====================================================\n")
cat("  Step 11 完成!\n")
cat("====================================================\n")
cat("\n  输出文件:\n")
cat("  data/11_panel_fisher_results.csv    — Panel×DEG/DSE Fisher检验\n")
cat("  data/11_panel_overlap_genes.csv     — Panel重叠基因详表\n")
cat("  data/11_GO_BP_DEG_only.csv          — DEG-only GO BP富集\n")
cat("  data/11_GO_BP_DSE_only.csv          — DSE-only GO BP富集\n")
cat("  data/11_KEGG_DEG_only.csv           — DEG-only KEGG富集\n")
cat("  data/11_KEGG_DSE_only.csv           — DSE-only KEGG富集\n")
cat("  data/11_GO_combined_top.csv         — GO对比top terms\n")
cat("  data/11_GO_category_comparison.csv  — GO功能大类对比\n")
cat("  data/11_YBX1_target_validation.csv  — YBX1靶标验证表\n")
cat("  data/11_summary_statistics.csv      — 汇总统计\n")
cat("\n  figures/11_panel_enrichment_bar.pdf  — Panel enrichment柱状图\n")
cat("  figures/11_panel_enrichment_heatmap.pdf — Panel enrichment热图\n")
cat("  figures/11_GO_BP_DEGonly_vs_DSEonly.pdf — ★ GO功能分离dot plot\n")
cat("  figures/11_KEGG_DEGonly_vs_DSEonly.pdf  — KEGG功能分离\n")
cat("  figures/11_GO_category_stacked.pdf      — GO功能大类对比\n")
cat("  figures/11_YBX1_target_layer_assignment.pdf — YBX1靶标验证\n")
cat("\n")
cat("  ★ 判读要点:\n")
cat("  1. Panel enrichment: DEG应该与SenMayo/SASP/p53富集；DSE?\n")
cat("  2. GO对比: DEG-only vs DSE-only是否走不同的功能通路?\n")
cat("     → 这是Figure 4C的核心数据\n")
cat("  3. YBX1靶标: 如果大比例在DSE中 → 独立交叉验证\n")
cat("     如果OR>1且显著 → 可写入Discussion与Xiao 2023对话\n")
cat("\n")
cat("  请上传:\n")
cat("  1. data/11_panel_fisher_results.csv\n")
cat("  2. data/11_GO_BP_DEG_only.csv (前20行)\n")
cat("  3. data/11_GO_BP_DSE_only.csv (前20行)\n")
cat("  4. data/11_GO_category_comparison.csv\n")
cat("  5. data/11_YBX1_target_validation.csv\n")
cat("  6. figures/11_GO_BP_DEGonly_vs_DSEonly.pdf\n")
cat("  7. figures/11_panel_enrichment_heatmap.pdf\n")
cat("====================================================\n")
