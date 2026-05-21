# 03_alternative_splicing.R — Differential alternative splicing analysis with logit-transformed PSI

source("00_project_config.R")

load_packages(c("data.table", "limma", "ggplot2", "pheatmap",
                "reshape2", "dplyr"))

cat("=== Step 3 v3: Alternative Splicing (DSE) ===\n")
cat("  ★ v3: logit(PSI) -> limma (Schafer 2015; SUPPA2 standard)\n\n")

# ── 3.0 配置 ──
DPSI_THRESH    <- 0.05    # modT配合的effect-size阈值(原始PSI尺度)
DPSI_STRICT    <- 0.10    # 纯effect-size筛选阈值
FDR_THRESH     <- 0.05
MIN_READS      <- 10
PSI_PSEUDO     <- 0.001   # logit变换伪计数，防止log(0)

sample_info <- sample_table

as_files <- list(
  SE   = NULL,
  A3SS = f_A3SS_JCEC,
  A5SS = f_A5SS_JCEC,
  MXE  = f_MXE_JC,
  RI   = f_RI_JCEC
)

# ── 3.1 合并SE分片 ──
cat("3.1 Merging SE split files...\n")

se_merged_path <- file.path(dir_data, "SE_MATS_JCEC_merged.tsv")

if (file.exists(se_merged_path)) {
  cat("  Already merged, skipping\n")
  as_files$SE <- se_merged_path
} else {
  if (file.exists(f_SE_JCEC)) {
    cat("  Found complete SE JCEC file\n")
    as_files$SE <- f_SE_JCEC
  } else {
    se_parts <- sort(list.files(dir_raw, pattern = "SE.MATS.JCEC.part|SE_MATS_JCEC_part",
                                full.names = TRUE))
    if (length(se_parts) > 0) {
      cat("  Found", length(se_parts), "SE JCEC parts, merging...\n")
      se_list <- lapply(se_parts, function(f) fread(f, sep = "\t"))
      se_merged <- rbindlist(se_list, use.names = TRUE)
      # ★ v3: 去重
      se_merged <- unique(se_merged, by = "AS_ID")
      fwrite(se_merged, se_merged_path, sep = "\t")
      as_files$SE <- se_merged_path
      cat("  Merged SE JCEC:", nrow(se_merged), "events (after dedup)\n")
    } else {
      cat("  WARNING: No SE JCEC file found!\n")
    }
  }
}


# ── 3.2 加载+过滤 ──

load_and_filter_as <- function(filepath, as_type) {
  cat("\n  Loading", as_type, "from", basename(filepath), "...\n")
  dt <- fread(filepath, sep = "\t")
  setnames(dt, trimws(names(dt)))

  # ★ v3: AS_ID去重
  n_before <- nrow(dt)
  dt <- unique(dt, by = "AS_ID")
  if (nrow(dt) < n_before) cat("    Removed", n_before - nrow(dt), "duplicate AS_IDs\n")

  psi_cols <- paste0("IncLevel_", sample_info$sample_id)
  ijc_cols <- paste0("IJC_", sample_info$sample_id)
  sjc_cols <- paste0("SJC_", sample_info$sample_id)

  missing <- setdiff(c(psi_cols, ijc_cols, sjc_cols), names(dt))
  if (length(missing) > 0) {
    cat("    WARNING: missing columns:", paste(head(missing, 5), collapse = ", "), "\n")
    return(NULL)
  }

  psi_mat <- as.matrix(dt[, ..psi_cols])
  colnames(psi_mat) <- sample_info$sample_id

  ijc_mat <- as.matrix(dt[, ..ijc_cols])
  sjc_mat <- as.matrix(dt[, ..sjc_cols])
  total_mat <- ijc_mat + sjc_mat
  colnames(total_mat) <- sample_info$sample_id

  psi_complete <- rowSums(is.na(psi_mat)) == 0

  read_pass <- sapply(levels(sample_info$passage), function(p) {
    samps <- sample_info$sample_id[sample_info$passage == p]
    n_req <- max(1, length(samps) - 1)
    rowSums(total_mat[, samps, drop = FALSE] >= MIN_READS) >= n_req
  })
  read_ok <- rowSums(read_pass) == ncol(read_pass)

  psi_var <- apply(psi_mat, 1, sd, na.rm = TRUE) > 0.01

  keep <- psi_complete & read_ok & psi_var
  cat("    Raw:", nrow(dt), "-> Filtered:", sum(keep), "events\n")

  result <- data.table(
    AS_ID       = dt$AS_ID,
    GeneID      = dt$GeneID,
    Gene_symbol = dt$Gene_symbol,
    Chr         = dt$Chr,
    Strand      = dt$Strand,
    as_type     = as_type
  )
  result <- cbind(result, as.data.table(psi_mat))
  result <- result[keep]

  return(result)
}


# ── 3.3 ★ v3: limma DSE — logit变换后建模 ──

run_limma_dse <- function(result_dt) {
  psi_cols <- sample_info$sample_id
  psi_mat <- as.matrix(result_dt[, ..psi_cols])

  # ★ v3核心修正: logit变换
  # PSI ∈ [0,1] → 限定到 (0,1) → logit → ℝ
  # 文献依据: Schafer et al. 2015 Bioinformatics; SUPPA2 默认做法
  # 目的: 方差稳定化，消除边界效应对limma的影响
  psi_bounded <- pmin(pmax(psi_mat, PSI_PSEUDO), 1 - PSI_PSEUDO)
  psi_logit   <- log(psi_bounded / (1 - psi_bounded))

  # limma设计矩阵
  design <- model.matrix(~ 0 + passage, data = sample_info)
  colnames(design) <- levels(sample_info$passage)

  contrasts <- makeContrasts(
    P8_vs_P2  = P8 - P2,
    P10_vs_P2 = P10 - P2,
    P12_vs_P2 = P12 - P2,
    levels = design
  )

  # ★ 在logit尺度上拟合limma
  fit <- lmFit(psi_logit, design)
  fit2 <- contrasts.fit(fit, contrasts)
  fit2 <- eBayes(fit2)

  comps <- c("P8_vs_P2", "P10_vs_P2", "P12_vs_P2")
  for (comp in comps) {
    tt <- topTable(fit2, coef = comp, number = Inf, sort.by = "none")
    comp_short <- gsub("_vs_P2", "", comp)

    # ★ dPSI在原始PSI尺度上计算（非logit空间）
    # 保证effect-size阈值和报告值与文献惯例一致
    p2_samps   <- sample_info$sample_id[sample_info$passage == "P2"]
    comp_samps <- sample_info$sample_id[sample_info$passage == comp_short]
    dpsi_raw   <- rowMeans(psi_mat[, comp_samps, drop = FALSE]) -
                  rowMeans(psi_mat[, p2_samps, drop = FALSE])

    result_dt[[paste0("dpsi_", comp_short)]]  <- dpsi_raw
    result_dt[[paste0("pval_", comp_short)]]   <- tt$P.Value   # 来自logit模型
    result_dt[[paste0("padj_", comp_short)]]   <- tt$adj.P.Val # 来自logit模型
    result_dt[[paste0("tstat_", comp_short)]]  <- tt$t

    # Strategy A: p-value(logit模型) + effect-size(原始PSI)
    result_dt[[paste0("dse_modT_", comp_short)]] <-
      (tt$adj.P.Val < FDR_THRESH) & (abs(dpsi_raw) >= DPSI_THRESH)

    # Strategy B: 纯effect-size（原始PSI尺度）
    result_dt[[paste0("dse_effect_", comp_short)]] <-
      abs(dpsi_raw) >= DPSI_STRICT
  }

  # 分组均值PSI（原始尺度）
  for (p in levels(sample_info$passage)) {
    samps <- sample_info$sample_id[sample_info$passage == p]
    result_dt[[paste0("psi_", p)]] <- rowMeans(psi_mat[, samps, drop = FALSE])
  }

  return(result_dt)
}


# ── 3.4 主分析循环 ──
cat("\n3.4 Running DSE analysis...\n")

all_results <- list()
summary_rows <- list()

for (as_type in names(as_files)) {
  fpath <- as_files[[as_type]]
  if (is.null(fpath) || !file.exists(fpath)) {
    cat("  SKIP:", as_type, "- file not found\n")
    next
  }

  dt <- load_and_filter_as(fpath, as_type)
  if (is.null(dt) || nrow(dt) == 0) next

  dt <- run_limma_dse(dt)

  row <- data.frame(
    as_type      = as_type,
    total_events = nrow(dt),
    DSE_modT_P8  = sum(dt$dse_modT_P8, na.rm = TRUE),
    DSE_modT_P10 = sum(dt$dse_modT_P10, na.rm = TRUE),
    DSE_modT_P12 = sum(dt$dse_modT_P12, na.rm = TRUE),
    DSE_effect_P8  = sum(dt$dse_effect_P8, na.rm = TRUE),
    DSE_effect_P10 = sum(dt$dse_effect_P10, na.rm = TRUE),
    DSE_effect_P12 = sum(dt$dse_effect_P12, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  row$DSE_any_modT   <- nrow(dt[dse_modT_P8 == TRUE | dse_modT_P10 == TRUE | dse_modT_P12 == TRUE])
  row$DSE_any_effect <- nrow(dt[dse_effect_P8 == TRUE | dse_effect_P10 == TRUE | dse_effect_P12 == TRUE])

  summary_rows[[as_type]] <- row
  all_results[[as_type]] <- dt

  cat("  ", as_type, "DSE (modT):", row$DSE_any_modT,
      " | DSE (effect):", row$DSE_any_effect, "\n")
}

combined <- rbindlist(all_results, fill = TRUE)
summary_df <- rbindlist(summary_rows)

totals <- data.frame(
  as_type = "TOTAL",
  total_events   = sum(summary_df$total_events),
  DSE_modT_P8    = sum(summary_df$DSE_modT_P8),
  DSE_modT_P10   = sum(summary_df$DSE_modT_P10),
  DSE_modT_P12   = sum(summary_df$DSE_modT_P12),
  DSE_effect_P8  = sum(summary_df$DSE_effect_P8),
  DSE_effect_P10 = sum(summary_df$DSE_effect_P10),
  DSE_effect_P12 = sum(summary_df$DSE_effect_P12),
  DSE_any_modT   = sum(summary_df$DSE_any_modT),
  DSE_any_effect = sum(summary_df$DSE_any_effect),
  stringsAsFactors = FALSE
)
summary_df <- rbind(summary_df, totals, fill = TRUE)

cat("\n=== DSE Summary ===\n")
print(summary_df)

# ── 3.5 保存 ──
cat("\n3.5 Saving results...\n")

save_data(summary_df, "03_DSE_summary_v2.csv")

dse_any <- combined[dse_modT_P8 == TRUE | dse_modT_P10 == TRUE | dse_modT_P12 == TRUE]

out_cols <- c("AS_ID", "as_type", "GeneID", "Gene_symbol", "Chr", "Strand",
              "psi_P2", "psi_P8", "psi_P10", "psi_P12",
              "dpsi_P8", "dpsi_P10", "dpsi_P12",
              "padj_P8", "padj_P10", "padj_P12",
              "dse_modT_P8", "dse_modT_P10", "dse_modT_P12")
out_cols <- intersect(out_cols, names(dse_any))
save_data(dse_any[, ..out_cols], "03_DSE_events_v2.csv")

core_genes <- combined[(dse_modT_P8 == TRUE) + (dse_modT_P10 == TRUE) + (dse_modT_P12 == TRUE) >= 2]
save_data(core_genes[, ..out_cols], "03_DSE_core_genes.csv")

cat("  DSE events saved:", nrow(dse_any), "\n")
cat("  Core DSE genes (>=2 comparisons):", uniqueN(core_genes$Gene_symbol), "\n")

# ── 3.6 剪接因子表达 ──
cat("\n3.6 Splicing factor expression...\n")

# 加载表达数据计算SF的log2FC
sf_genes <- c("SRSF1","SRSF2","SRSF3","SRSF4","SRSF5","SRSF6","SRSF7","SRSF9","SRSF10","SRSF11",
              "HNRNPA1","HNRNPA2B1","HNRNPC","HNRNPD","HNRNPF","HNRNPH1","HNRNPK","HNRNPL","HNRNPM","HNRNPU",
              "PTBP1","PTBP2","RBFOX1","RBFOX2","MBNL1","MBNL2","CELF1","CELF2",
              "TRA2A","TRA2B","ESRP1","ESRP2","QKI","NOVA1","NOVA2",
              "SF3B1","U2AF1","U2AF2","SF1","PRPF8","SNRNP200",
              "KHDRBS1","KHDRBS2","KHDRBS3","YBX1","FUS","TARDBP","TIA1","TIAL1","ELAVL1")

raw_expr <- read.delim(f_expression, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
tpm_cols <- grep("_tpm$", colnames(raw_expr), value = TRUE)
tpm_mat  <- as.matrix(raw_expr[, tpm_cols])
rownames(tpm_mat) <- raw_expr$id
colnames(tpm_mat) <- gsub("_tpm$", "", tpm_cols)
tpm_mat <- tpm_mat[, sample_info$sample_id]
gene_info <- raw_expr[, c("id", "Symbol")]

sf_rows <- list()
for (g in sf_genes) {
  idx <- which(gene_info$Symbol == g)
  if (length(idx) == 0) next
  eid <- gene_info$id[idx[1]]
  vals <- tpm_mat[eid, ]
  p2_mean  <- mean(vals[sample_info$passage == "P2"])
  p8_mean  <- mean(vals[sample_info$passage == "P8"])
  p10_mean <- mean(vals[sample_info$passage == "P10"])
  p12_mean <- mean(vals[sample_info$passage == "P12"])
  sf_rows[[g]] <- data.frame(
    gene = g, P2 = p2_mean, P8 = p8_mean, P10 = p10_mean, P12 = p12_mean,
    lfc_P8  = log2((p8_mean + 1)  / (p2_mean + 1)),
    lfc_P10 = log2((p10_mean + 1) / (p2_mean + 1)),
    lfc_P12 = log2((p12_mean + 1) / (p2_mean + 1)),
    stringsAsFactors = FALSE
  )
}
sf_df <- do.call(rbind, sf_rows)
save_data(sf_df, "03_splicing_factor_expression.csv")

cat("  SF expression:", nrow(sf_df), "genes\n")


# ── 3.7 可视化 ──
cat("\n3.7 Generating figures...\n")

# Convert to stable data.frame for visualization (avoids data.table version issues)
combined <- as.data.frame(combined)

# Fig 3A: DSE barplot
plot_data <- summary_df[summary_df$as_type != "TOTAL", ]
plot_long <- reshape2::melt(plot_data,
                  id.vars = "as_type",
                  measure.vars = c("DSE_modT_P8", "DSE_modT_P10", "DSE_modT_P12"),
                  variable.name = "comparison", value.name = "count")
plot_long$comparison <- gsub("DSE_modT_", "", plot_long$comparison)
plot_long$comparison <- factor(plot_long$comparison, levels = c("P8", "P10", "P12"))
plot_long$as_type <- factor(plot_long$as_type, levels = c("SE", "A3SS", "A5SS", "MXE", "RI"))

p1 <- ggplot(plot_long, aes(x = as_type, y = count, fill = comparison)) +
  geom_bar(stat = "identity", position = position_dodge(0.7), width = 0.6) +
  geom_text(aes(label = count), position = position_dodge(0.7), vjust = -0.3, size = 2.5) +
  scale_fill_manual(values = passage_colors[c("P8","P10","P12")],
                    labels = c("P8 vs P2", "P10 vs P2", "P12 vs P2")) +
  labs(x = "AS type", y = "Number of DSE",
       title = "Differential splicing events (FDR < 0.05 & |dPSI| >= 0.05)") +
  theme_bindlab() + theme(legend.title = element_blank())
save_fig(p1, "03_DSE_barplot_v2.pdf", w = 7, h = 4)

# Fig 3B: Volcano
volcano_list <- list()
for (comp in c("P8", "P10", "P12")) {
  vdf <- data.frame(
    dpsi    = combined[[paste0("dpsi_", comp)]],
    log10p  = -log10(pmax(combined[[paste0("padj_", comp)]], 1e-50)),
    sig     = ifelse(combined[[paste0("dse_modT_", comp)]], "DSE (FDR & dPSI)", "NS"),
    gene    = combined$Gene_symbol,
    comparison = paste0(comp, " vs P2"),
    stringsAsFactors = FALSE
  )
  volcano_list[[comp]] <- vdf
}
vdf_all <- do.call(rbind, volcano_list)
vdf_all$sig <- factor(vdf_all$sig, levels = c("NS", "DSE (FDR & dPSI)"))
vdf_all$comparison <- factor(vdf_all$comparison, levels = c("P8 vs P2", "P10 vs P2", "P12 vs P2"))

p2 <- ggplot(vdf_all, aes(x = dpsi, y = log10p, color = sig)) +
  geom_point(size = 0.3, alpha = 0.4) +
  scale_color_manual(values = c("NS" = "grey70", "DSE (FDR & dPSI)" = "#D73027")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50", linewidth = 0.3) +
  geom_vline(xintercept = c(-DPSI_THRESH, DPSI_THRESH), linetype = "dashed", color = "grey50", linewidth = 0.3) +
  facet_wrap(~ comparison, nrow = 1) +
  labs(x = "dPSI (original scale)", y = "-log10(FDR from logit model)", color = NULL) +
  theme_bindlab() + theme(legend.position = "bottom", legend.key.size = unit(0.3, "cm"))
save_fig(p2, "03_DSE_volcano_v2.pdf", w = 10, h = 4)

# Fig 3C: Top DSE heatmap
combined$max_abs_dpsi <- pmax(abs(combined$dpsi_P8), abs(combined$dpsi_P10), abs(combined$dpsi_P12), na.rm = TRUE)
top30 <- head(combined[order(-combined$max_abs_dpsi), ], 30)

hm_mat <- as.matrix(top30[, c("dpsi_P8", "dpsi_P10", "dpsi_P12")])
rownames(hm_mat) <- paste0(top30$Gene_symbol, " (", top30$as_type, ")")
colnames(hm_mat) <- c("dPSI P8", "dPSI P10", "dPSI P12")

save_pheatmap("03_DSE_heatmap_top30_v2.pdf", w = 5, h = 8)
pheatmap(hm_mat, color = heatmap_colors, breaks = seq(-0.5, 0.5, length.out = 101),
         cluster_cols = FALSE, cluster_rows = TRUE, fontsize_row = 7,
         main = "Top 30 DSE by |dPSI|")
dev.off()

png(file.path(dir_fig, "03_DSE_heatmap_top30_v2.png"), width = 5, height = 8, units = "in", res = 300)
pheatmap(hm_mat, color = heatmap_colors, breaks = seq(-0.5, 0.5, length.out = 101),
         cluster_cols = FALSE, cluster_rows = TRUE, fontsize_row = 7,
         main = "Top 30 DSE by |dPSI|")
dev.off()

# Fig 3D: Direction
dir_data_plot <- data.frame()
for (comp in c("P8", "P10", "P12")) {
  subset_dt <- combined[combined[[paste0("dse_modT_", comp)]] == TRUE, ]
  n_inc <- sum(subset_dt[[paste0("dpsi_", comp)]] > 0, na.rm = TRUE)
  n_exc <- sum(subset_dt[[paste0("dpsi_", comp)]] < 0, na.rm = TRUE)
  dir_data_plot <- rbind(dir_data_plot, data.frame(
    comparison = paste0(comp, " vs P2"),
    direction = c("Inclusion", "Exclusion"), count = c(n_inc, n_exc)))
}

p3 <- ggplot(dir_data_plot, aes(x = comparison, y = count, fill = direction)) +
  geom_bar(stat = "identity", position = position_dodge(0.7), width = 0.6) +
  geom_text(aes(label = count), position = position_dodge(0.7), vjust = -0.3, size = 3) +
  scale_fill_manual(values = c("Inclusion" = "#D73027", "Exclusion" = "#4575B4")) +
  labs(x = NULL, y = "Number of DSE", title = "DSE direction: inclusion vs exclusion") +
  theme_bindlab() + theme(legend.title = element_blank())
save_fig(p3, "03_DSE_direction_v2.pdf", w = 5, h = 4)

# Fig 3E: SF heatmap
sf_top <- sf_df[order(-abs(sf_df$lfc_P8)), ][1:min(20, nrow(sf_df)), ]
sf_hm <- as.matrix(sf_top[, c("lfc_P8", "lfc_P10", "lfc_P12")])
rownames(sf_hm) <- sf_top$gene
colnames(sf_hm) <- c("LFC P8", "LFC P10", "LFC P12")

save_pheatmap("03_splicing_factor_heatmap_v2.pdf", w = 5, h = 6)
pheatmap(sf_hm, color = heatmap_colors, breaks = seq(-1, 1, length.out = 101),
         cluster_cols = FALSE, fontsize_row = 8,
         main = "Splicing factor expression changes (LFC vs P2)")
dev.off()

png(file.path(dir_fig, "03_splicing_factor_heatmap_v2.png"), width = 5, height = 6, units = "in", res = 300)
pheatmap(sf_hm, color = heatmap_colors, breaks = seq(-1, 1, length.out = 101),
         cluster_cols = FALSE, fontsize_row = 8,
         main = "Splicing factor expression changes (LFC vs P2)")
dev.off()

cat("\n=== Step 3 v3 complete! ===\n")
cat("Key change: logit(PSI) transformation before limma modeling\n")
cat("  p-values from logit model; dPSI on original PSI scale\n")
cat("Output files: same as v2 (overwritten in place)\n")
