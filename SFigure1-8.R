# SFigure1-8_v5.R — Supplementary figures S1–S8

source("00_project_config.R")
load_packages(c("ggplot2", "cowplot", "reshape2", "scales", "dplyr", "tidyr", "pheatmap", "grid", "ggrepel"))

cat("\n========== Generating Supplementary Figures S1–S8 ==========\n")

# ---- SFigure1_v5.R ----

cat("\n========== SFig1_v5: QC ==========\n")

PW_FULL <- 6.7; PW_HALF <- 3.35
FIG_W <- PW_FULL; FIG_H <- 4.8
TXT <- 7

# ============ 1a: Library size boxplot ============
cat("  1a: Library size boxplot...\n")
qc_stats <- read.csv(file.path(dir_data, "01_QC_statistics.csv"))
qc_stats$passage <- factor(qc_stats$passage, levels = c("P2", "P8", "P10", "P12"))

figS1a <- ggplot(qc_stats, aes(x = passage, y = total_count / 1e6, fill = passage)) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.4, linewidth = 0.4, color = "grey30") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, linewidth = 0.4, color = "grey30") +
  geom_jitter(width = 0.12, size = 1.6, alpha = 0.7, color = "grey30") +
  scale_fill_manual(values = passage_colors, guide = "none") +
  labs(x = "Passage", y = "Library size (M reads)", subtitle = "RNA-seq library sizes") +
  theme_bindlab(base_size = TXT)
cat("  1a done\n")

# ============ 1b: Detected genes ============
cat("  1b: Detected genes...\n")
figS1b <- ggplot(qc_stats, aes(x = passage, y = genes_detected_count, fill = passage)) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.4, linewidth = 0.4, color = "grey30") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, linewidth = 0.4, color = "grey30") +
  geom_jitter(width = 0.12, size = 1.6, alpha = 0.7, color = "grey30") +
  scale_fill_manual(values = passage_colors, guide = "none") +
  labs(x = "Passage", y = "Detected genes", subtitle = "Number of detected genes") +
  theme_bindlab(base_size = TXT)
cat("  1b done\n")

# ============ 1c: Sample correlation heatmap ============
cat("  1c: Sample correlation heatmap...\n")
vst <- read.csv(file.path(dir_data, "02_VST_expression.csv"), stringsAsFactors = FALSE, check.names = FALSE)
sample_cols <- grep("^P\\d+-\\d+$", colnames(vst), value = TRUE)
if (length(sample_cols) == 0) sample_cols <- grep("^P\\d+\\.\\d+$", colnames(vst), value = TRUE)
vst_mat <- as.matrix(vst[, sample_cols])
cor_mat <- cor(vst_mat, method = "pearson")

sample_anno <- data.frame(
  Passage = c("P2","P2","P8","P8","P8","P10","P10","P10","P12","P12","P12"),
  row.names = colnames(cor_mat))
anno_colors <- list(Passage = passage_colors)

cor_hm <- pheatmap(cor_mat,
  color = colorRampPalette(c("#4575B4", "white", "#D73027"))(100),
  annotation_col = sample_anno, annotation_row = sample_anno,
  annotation_colors = anno_colors,
  cluster_rows = TRUE, cluster_cols = TRUE,
  show_colnames = TRUE, show_rownames = TRUE,
  fontsize = TXT - 1, border_color = "grey90",
  legend = TRUE, main = "Sample Pearson correlation",
  silent = TRUE)
figS1c <- cor_hm$gtable
figS1c <- grid::grobTree(grid::rectGrob(gp = grid::gpar(fill = "white", col = NA)), figS1c)
cat("  1c done\n")

# ============ Assembly ============
cat("  Assembling SFigure1_v5...\n")
SFig1_row1 <- plot_grid(figS1a, figS1b, labels = c("a", "b"), label_size = 9,
                         rel_widths = c(1, 1), ncol = 2)
SFig1 <- plot_grid(SFig1_row1, figS1c, labels = c("", "c"), label_size = 9,
                    ncol = 1, rel_heights = c(0.7, 1))
save_fig(SFig1, "SFigure1_v5.pdf", w = FIG_W, h = FIG_H, dpi = 300)
cat("  SFigure1_v5 complete.\n\n")

# ---- SFigure2_v5.R ----

cat("\n========== SFig2_v5: DEG火山图 ==========\n")

PW_FULL <- 6.7; PW_HALF <- 3.35
FIG_W <- PW_FULL; FIG_H <- 3.5
TXT <- 7

make_volcano <- function(df, passage_label, fc_cut = 1, fdr_cut = 0.05) {
  # Cap FDR floor at 1e-20 for sensible y-axis
  df$neg_log10_padj <- pmin(-log10(pmax(df$padj, 1e-20)), 20)
  df$diffexpressed <- "NS"
  df$diffexpressed[df$log2FoldChange > fc_cut & df$padj < fdr_cut] <- "Up"
  df$diffexpressed[df$log2FoldChange < -fc_cut & df$padj < fdr_cut] <- "Down"
  df$diffexpressed <- factor(df$diffexpressed, levels = c("Up", "Down", "NS"))

  n_up <- sum(df$diffexpressed == "Up")
  n_down <- sum(df$diffexpressed == "Down")

  p <- ggplot(df, aes(x = log2FoldChange, y = neg_log10_padj, color = diffexpressed)) +
    geom_point(size = 0.8, alpha = 0.5, stroke = 0) +
    scale_color_manual(values = c("Up" = "#D73027", "Down" = "#4575B4", "NS" = "grey70"), guide = "none") +
    geom_hline(yintercept = -log10(fdr_cut), linetype = "dashed", color = "grey60", linewidth = 0.3) +
    geom_vline(xintercept = c(-fc_cut, fc_cut), linetype = "dashed", color = "grey60", linewidth = 0.3) +
    scale_x_continuous(limits = c(-4, 4)) +
    coord_cartesian(ylim = c(0, 20)) +
    annotate("text", x = -3.5, y = 19,
             label = paste0("Down: ", n_down), size = 2.4, color = "#4575B4", fontface = "bold", hjust = 0) +
    annotate("text", x = 1.5, y = 19,
             label = paste0("Up: ", n_up), size = 2.4, color = "#D73027", fontface = "bold", hjust = 0) +
    labs(x = expression(log[2](fold~change)), y = expression(-log[10](FDR)),
         title = passage_label) +
    theme_bindlab(base_size = TXT) +
    theme(plot.title = element_text(size = TXT, face = "bold", hjust = 0.5))
  p
}

de_p8  <- read.csv(file.path(dir_data, "02_DE_P8_vs_P2.csv"), stringsAsFactors = FALSE)
de_p10 <- read.csv(file.path(dir_data, "02_DE_P10_vs_P2.csv"), stringsAsFactors = FALSE)
de_p12 <- read.csv(file.path(dir_data, "02_DE_P12_vs_P2.csv"), stringsAsFactors = FALSE)

figS2a <- make_volcano(de_p8,  "P8 vs P2")
figS2b <- make_volcano(de_p10, "P10 vs P2")
figS2c <- make_volcano(de_p12, "P12 vs P2")

SFig2 <- plot_grid(figS2a, figS2b, figS2c, labels = c("a", "b", "c"),
                    label_size = 9, ncol = 3)
save_fig(SFig2, "SFigure2_v5.pdf", w = FIG_W, h = FIG_H, dpi = 300)
cat("  SFigure2_v5 complete.\n\n")

# ---- SFigure3_v5.R ----

cat("\n========== SFig3_v5: GSEA Hallmark ==========\n")

PW_FULL <- 6.7; PW_HALF <- 3.35
FIG_W <- PW_FULL; FIG_H <- 5.0
TXT <- 7

# ============ Data ============
gsea <- read.csv(file.path(dir_data, "02_fgsea_hallmark_all.csv"), stringsAsFactors = FALSE)

# ============ 3a: Top pathways dot plot ============
cat("  3a: Top GSEA dot plot...\n")
gsea_top <- gsea %>%
  group_by(comparison) %>%
  slice_min(padj, n = 8) %>%
  ungroup()
gsea_top$Description <- factor(gsea_top$pathway, levels = rev(unique(gsea_top$pathway)))
gsea_top$comparison <- factor(gsea_top$comparison,
  levels = c("P8_vs_P2", "P10_vs_P2", "P12_vs_P2"),
  labels = c("P8 vs P2", "P10 vs P2", "P12 vs P2"))

figS3a <- ggplot(gsea_top, aes(x = NES, y = Description, color = comparison, size = -log10(padj))) +
  geom_point(alpha = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.3) +
  scale_color_manual(values = c("P8 vs P2" = "#91BFDB", "P10 vs P2" = "#FC8D59", "P12 vs P2" = "#D73027"),
                     name = "") +
  scale_size_continuous(range = c(1.5, 4), name = expression(-log[10](p[adj]))) +
  labs(x = "Normalized Enrichment Score (NES)", y = "",
       subtitle = "Hallmark GSEA: top enriched pathways") +
  theme_bindlab(base_size = TXT) +
  theme(plot.subtitle = element_text(size = TXT, face = "bold", hjust = 0.5),
        legend.position = "bottom", legend.key.size = unit(0.18, "cm"),
        legend.text = element_text(size = TXT - 2),
        axis.text.y = element_text(size = TXT - 1))
cat("  3a done\n")

# ============ 3b: Key pathway NES trends ============
cat("  3b: Key pathway NES trends...\n")
key_pathways <- c("HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_INFLAMMATORY_RESPONSE",
                  "HALLMARK_P53_PATHWAY", "HALLMARK_APOPTOSIS", "HALLMARK_INTERFERON_GAMMA_RESPONSE",
                  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")
gsea_key <- gsea %>%
  filter(pathway %in% key_pathways) %>%
  mutate(
    passage_num = case_when(
      comparison == "P8_vs_P2" ~ 8,
      comparison == "P10_vs_P2" ~ 10,
      comparison == "P12_vs_P2" ~ 12
    ),
    pathway_label = gsub("_", " ", pathway),
    pathway_label = tools::toTitleCase(tolower(pathway_label))
  )

figS3b <- ggplot(gsea_key, aes(x = passage_num, y = NES, color = pathway_label, group = pathway_label)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2.0) +
  scale_x_continuous(breaks = c(8, 10, 12), labels = c("P8", "P10", "P12")) +
  scale_color_manual(values = c("#D73027", "#4575B4", "#4DAF4A", "#FF7F00", "#984EA3", "#A65628"),
                     name = "") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.3) +
  labs(x = "Passage", y = "NES", subtitle = "Key senescence-related pathway trends") +
  theme_bindlab(base_size = TXT) +
  theme(plot.subtitle = element_text(size = TXT, face = "bold", hjust = 0.5),
        legend.position = "bottom", legend.key.size = unit(0.18, "cm"),
        legend.text = element_text(size = TXT - 2))
cat("  3b done\n")

# ============ Assembly ============
cat("  Assembling SFigure3_v5...\n")
SFig3 <- plot_grid(figS3a, figS3b, labels = c("a", "b"), label_size = 9,
                    rel_widths = c(1.3, 1), ncol = 2)
save_fig(SFig3, "SFigure3_v5.pdf", w = FIG_W, h = FIG_H, dpi = 300)
cat("  SFigure3_v5 complete.\n\n")

# ---- SFigure4_v5.R ----

cat("\n========== SFig4_v5: DSE filtering comparison ==========\n")

PW_FULL <- 6.7; PW_HALF <- 3.35
FIG_W <- PW_FULL; FIG_H <- 5.0
TXT <- 7

# ============ Data ============
dse_summary <- read.csv(file.path(dir_data, "03_DSE_summary_v2.csv"), stringsAsFactors = FALSE)

# ============ 4a: AND vs effect-size comparison barplot ============
cat("  4a: AND vs effect-size comparison...\n")
dse_type_sub <- subset(dse_summary, as_type != "TOTAL")
dse_type_sub$as_type <- factor(dse_type_sub$as_type, levels = c("SE", "MXE", "A3SS", "A5SS", "RI"))

cmp_long <- melt(dse_type_sub,
  id.vars = "as_type",
  measure.vars = c("DSE_any_modT", "DSE_any_effect"),
  variable.name = "filter", value.name = "n_events")
cmp_long$filter <- factor(cmp_long$filter,
  levels = c("DSE_any_effect", "DSE_any_modT"),
  labels = c("Effect-size only\n(|dPSI| >= 0.05)", "AND\n(FDR < 0.05 AND |dPSI| >= 0.05)"))

figS4a <- ggplot(cmp_long, aes(x = as_type, y = n_events, fill = filter)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, color = "white", linewidth = 0.15) +
  scale_fill_manual(values = c("Effect-size only\n(|dPSI| >= 0.05)" = "grey70",
                                "AND\n(FDR < 0.05 AND |dPSI| >= 0.05)" = "#4575B4"),
                    name = "Filtering strategy") +
  geom_text(aes(label = n_events), position = position_dodge(width = 0.7),
            size = 2.2, vjust = -0.3, color = "grey30", fontface = "bold") +
  scale_x_discrete(labels = c("SE" = "SE", "MXE" = "MXE", "A3SS" = "A3SS", "A5SS" = "A5SS", "RI" = "RI")) +
  labs(x = "AS type", y = "Number of DSE events",
       subtitle = paste0("AND: ", sum(dse_type_sub$DSE_any_modT),
                         " | Effect-size: ", sum(dse_type_sub$DSE_any_effect),
                         " | Reduction: ", round((1 - sum(dse_type_sub$DSE_any_modT) /
                         sum(dse_type_sub$DSE_any_effect)) * 100, 1), "%")) +
  theme_bindlab(base_size = TXT) +
  theme(plot.subtitle = element_text(size = TXT, face = "bold", hjust = 0.5),
        legend.position = "bottom", legend.key.size = unit(0.22, "cm"),
        legend.text = element_text(size = TXT - 2))
cat("  4a done\n")

# ============ 4b: DSE per passage (AND-filtered only) ============
cat("  4b: DSE per passage...\n")
dse_totals <- dse_summary[dse_summary$as_type == "TOTAL", ]
pass_data <- data.frame(
  passage = factor(c("P8", "P10", "P12"), levels = c("P8", "P10", "P12")),
  n_and = c(dse_totals$DSE_modT_P8, dse_totals$DSE_modT_P10, dse_totals$DSE_modT_P12),
  n_eff = c(dse_totals$DSE_effect_P8, dse_totals$DSE_effect_P10, dse_totals$DSE_effect_P12)
)

figS4b <- ggplot(pass_data, aes(x = passage)) +
  geom_col(aes(y = n_eff, fill = "Effect-size"), width = 0.4, alpha = 0.35, color = "grey70", linewidth = 0.2) +
  geom_col(aes(y = n_and, fill = "AND"), width = 0.35, color = "white", linewidth = 0.15) +
  geom_text(aes(y = n_and, label = n_and), size = 2.5, fontface = "bold", vjust = -0.4, color = "#2F5496") +
  geom_text(aes(y = n_eff, label = n_eff), size = 2.2, vjust = -1.8, color = "grey50") +
  scale_fill_manual(values = c("AND" = "#4575B4", "Effect-size" = "grey70"),
                    name = "Filter", labels = c("AND" = "AND (FDR+PSI)", "Effect-size" = "Effect-size only")) +
  labs(x = "Passage", y = "Number of DSE events", subtitle = "High-confidence DSE per passage") +
  ylim(0, max(pass_data$n_eff) * 1.25) +
  theme_bindlab(base_size = TXT) +
  theme(plot.subtitle = element_text(size = TXT, face = "bold", hjust = 0.5),
        legend.position = c(0.15, 0.88),
        legend.background = element_rect(fill = alpha("white", 0.9), color = "grey85", linewidth = 0.25),
        legend.key.size = unit(0.18, "cm"),
        legend.text = element_text(size = TXT - 2))
cat("  4b done\n")

# ============ Assembly ============
SFig4 <- plot_grid(figS4a, figS4b, labels = c("a", "b"), label_size = 9,
                    rel_widths = c(1.4, 1), ncol = 2)
save_fig(SFig4, "SFigure4_v5.pdf", w = FIG_W, h = FIG_H, dpi = 300)
cat("  SFigure4_v5 complete.\n\n")

# ---- SFigure5_v5.R ----

cat("\n========== SFig5_v5: 三阶段衰老程序 ==========\n")

PW_FULL <- 6.7; PW_HALF <- 3.35
FIG_W <- PW_FULL; FIG_H <- 3.5
TXT <- 7

# ============ Data ============
dse_summary <- read.csv(file.path(dir_data, "03_DSE_summary_v2.csv"), stringsAsFactors = FALSE)
sm <- read.csv(file.path(dir_data, "01_SenMayo_scores.csv"), stringsAsFactors = FALSE)
sf_expr <- read.csv(file.path(dir_data, "03_splicing_factor_expression.csv"), stringsAsFactors = FALSE)
sm$passage <- factor(sm$passage, levels = c("P2", "P8", "P10", "P12"))

dse_total <- dse_summary[dse_summary$as_type == "TOTAL", ]
dse_all <- dse_summary[dse_summary$as_type != "TOTAL", ]

# Track 1: DSE events per passage
dse_phase <- data.frame(
  passage = factor(c("P8", "P10", "P12"), levels = c("P8", "P10", "P12")),
  n_modT = c(dse_total$DSE_modT_P8, dse_total$DSE_modT_P10, dse_total$DSE_modT_P12),
  n_effect = c(dse_total$DSE_effect_P8, dse_total$DSE_effect_P10, dse_total$DSE_effect_P12)
)

track1 <- ggplot(dse_phase, aes(x = passage, y = n_modT)) +
  geom_col(fill = "#4575B4", width = 0.5, alpha = 0.85) +
  geom_text(aes(label = n_modT, y = n_modT / 2), size = 3.0, fontface = "bold", color = "white") +
  geom_text(aes(y = n_modT + 5, label = paste0("(", n_effect, ")")), size = 2.0, color = "grey50") +
  labs(x = "", y = "High-confidence DSE", title = "DSE events (AND)") +
  theme_bindlab(base_size = TXT) +
  theme(plot.title = element_text(size = TXT, face = "bold", hjust = 0.5))

# Track 2: SenMayo score trend
track2 <- ggplot(sm, aes(x = passage, y = SenMayo)) +
  stat_summary(fun = mean, geom = "line", aes(group = 1), linewidth = 0.9, color = "#D73027") +
  stat_summary(fun = mean, geom = "point", size = 2.8, color = "#D73027") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.12, linewidth = 0.5, color = "#D73027") +
  geom_point(aes(fill = passage), shape = 21, size = 1.8, color = "white",
             stroke = 0.3, alpha = 0.5, position = position_jitter(width = 0.08)) +
  scale_fill_manual(values = passage_colors, guide = "none") +
  coord_cartesian(ylim = c(6.5, 7.8)) +
  labs(x = "", y = "SenMayo score", title = "SenMayo") +
  theme_bindlab(base_size = TXT) +
  theme(plot.title = element_text(size = TXT, face = "bold", hjust = 0.5))

# Track 3: Key SF expression
sf_key <- sf_expr[sf_expr$gene %in% c("PTBP1", "SRSF2", "HNRNPL", "MBNL1"), ]
sf_long <- melt(sf_key, id.vars = "gene",
                measure.vars = c("P2", "P8", "P10", "P12"),
                variable.name = "passage", value.name = "tpm")
sf_long$passage <- factor(sf_long$passage, levels = c("P2", "P8", "P10", "P12"))

track3 <- ggplot(sf_long, aes(x = passage, y = tpm, color = gene, group = gene)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.8) +
  scale_color_manual(values = c("PTBP1" = "#D73027", "SRSF2" = "#4575B4",
                                 "HNRNPL" = "#4DAF4A", "MBNL1" = "#FF7F00"), name = "") +
  labs(x = "Passage", y = "TPM", title = "Key SF expression") +
  theme_bindlab(base_size = TXT) +
  theme(plot.title = element_text(size = TXT, face = "bold", hjust = 0.5),
        legend.position = "right", legend.key.size = unit(0.18, "cm"),
        legend.text = element_text(size = TXT - 2))

# Assembly
SFig5 <- plot_grid(track1, track2, track3, ncol = 3,
                    rel_widths = c(1, 1, 1.25),
                    labels = c("a", "b", "c"), label_size = 9)
save_fig(SFig5, "SFigure5_v5.pdf", w = FIG_W, h = FIG_H, dpi = 300)
cat("  SFigure5_v5 complete.\n\n")

# ---- SFigure6_v5.R ----

cat("\n========== SFig6_v5: Multi-layer convergence ==========\n")

PW_FULL <- 6.7; PW_HALF <- 3.35
FIG_W <- PW_FULL; FIG_H <- 5.5
TXT <- 7

# ============ Data ============
gene_cat <- read.csv(file.path(dir_data, "05_gene_layer_categories.csv"), stringsAsFactors = FALSE)
trends <- read.csv(file.path(dir_data, "05_multilayer_trends.csv"), stringsAsFactors = FALSE)

# ============ 6a: Layer category distribution ============
cat("  6a: Layer category distribution...\n")
layer_dist <- gene_cat %>%
  count(layer) %>%
  mutate(pct = round(n / sum(n) * 100, 1)) %>%
  arrange(desc(n))
layer_dist$layer <- factor(layer_dist$layer, levels = rev(layer_dist$layer))

layer_colors <- c("DEG only" = "#D73027", "DSE only" = "#4575B4",
                  "DEG & DSE" = "#4DAF4A", "editing_only" = "#FF7F00")

figS6a <- ggplot(layer_dist, aes(x = layer, y = n, fill = layer)) +
  geom_col(width = 0.55, color = "white", linewidth = 0.2) +
  geom_text(aes(label = paste0(n, " (", pct, "%)")), hjust = -0.05, size = 2.8, color = "grey30") +
  coord_flip() +
  scale_fill_manual(values = layer_colors, guide = "none") +
  labs(x = "", y = "Number of genes",
       subtitle = "Gene distribution across omics layers") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  theme_bindlab(base_size = TXT) +
  theme(plot.subtitle = element_text(size = TXT, face = "bold", hjust = 0.5))
cat("  6a done\n")

# ============ 6b: Top multi-layer genes with layer annotation ============
cat("  6b: Top multi-layer genes...\n")
# Use actual layer names: "DEG & DSE"
multi_layer <- gene_cat %>%
  filter(layer == "DEG & DSE") %>%
  distinct(symbol, .keep_all = TRUE) %>%
  arrange(layer) %>%
  head(20)

if (nrow(multi_layer) > 0) {
  multi_layer$symbol <- factor(multi_layer$symbol, levels = rev(multi_layer$symbol))

  figS6b <- ggplot(multi_layer, aes(x = layer, y = symbol, fill = layer)) +
    geom_tile(color = "white", linewidth = 0.5, width = 0.7, height = 0.8) +
    scale_fill_manual(values = c("DEG & DSE" = "#4DAF4A"), guide = "none") +
    labs(x = "", y = "",
         subtitle = sprintf("Genes perturbed in both DEG \u0026 DSE (n=%d)", nrow(multi_layer))) +
    theme_bindlab(base_size = TXT) +
    theme(plot.subtitle = element_text(size = TXT, face = "bold", hjust = 0.5),
          axis.text.y = element_text(size = TXT - 1),
          axis.text.x = element_text(size = TXT - 1, face = "bold", color = "#4DAF4A"))
} else {
  figS6b <- ggplot() + annotate("text", x = 0.5, y = 0.5,
    label = "Multi-layer genes:\nconvergent across transcriptome + splicing",
    size = 3, color = "grey40") + theme_void()
}
cat("  6b done\n")

# ============ Assembly ============
SFig6 <- plot_grid(figS6a, figS6b, labels = c("a", "b"), label_size = 9,
                    rel_widths = c(1, 1.1), ncol = 2)
save_fig(SFig6, "SFigure6_v5.pdf", w = FIG_W, h = FIG_H, dpi = 300)
cat("  SFigure6_v5 complete.\n\n")

# ---- SFigure7_v5.R ----

cat("\n========== SFig7_v5: RNA editing ==========\n")

PW_FULL <- 6.7; PW_HALF <- 3.35
FIG_W <- PW_FULL; FIG_H <- 4.0
TXT <- 7

# ============ Data ============
editing_summary <- read.csv(file.path(dir_data, "04_editing_summary.csv"), stringsAsFactors = FALSE)
editing_type_long <- read.csv(file.path(dir_data, "04_editing_type_long.csv"), stringsAsFactors = FALSE)

# ============ 7a: Editing type distribution ============
cat("  7a: Editing type distribution...\n")
if (nrow(editing_type_long) > 0) {
  # Aggregate across samples
  etype_agg <- editing_type_long %>%
    group_by(editing_type) %>%
    summarise(count = sum(count), .groups = "drop") %>%
    arrange(desc(count))
  n_types <- nrow(etype_agg)
  etype_agg$editing_type <- factor(etype_agg$editing_type, levels = rev(etype_agg$editing_type))

  # Generate enough colors (12 editing types A->G, A->C, etc.)
  etype_colors <- c("#4575B4","#D73027","#4DAF4A","#FF7F00","#984EA3",
                    "#A65628","#F781BF","#999999","#E41A1C","#377EB8",
                    "#FFFF33","#66C2A5")

  figS7a <- ggplot(etype_agg, aes(x = editing_type, y = count, fill = editing_type)) +
    geom_col(width = 0.55, color = "white", linewidth = 0.2) +
    geom_text(aes(label = scales::comma(count), y = count + max(count) * 0.03),
              size = 2.6, fontface = "bold", color = "grey30") +
    scale_fill_manual(values = etype_colors[1:n_types], guide = "none") +
    labs(x = "Editing type", y = "Number of sites",
         subtitle = "RNA editing site distribution") +
    theme_bindlab(base_size = TXT) +
    theme(plot.subtitle = element_text(size = TXT, face = "bold", hjust = 0.5),
          axis.text.x = element_text(angle = 20, hjust = 1))
} else {
  figS7a <- ggplot() + annotate("text", x = 0.5, y = 0.5,
    label = "RNA editing type distribution", size = 3, color = "grey50") + theme_void()
}
cat("  7a done\n")

# ============ 7b: Editing functional consequence distribution ============
cat("  7b: Editing functional consequence...\n")
edit_func <- read.csv(file.path(dir_data, "04_editing_function_dist.csv"), stringsAsFactors = FALSE)

if (nrow(edit_func) > 0) {
  # Aggregate across passages, use function_type and N
  func_agg <- edit_func %>%
    group_by(function_type) %>%
    summarise(count = sum(N), .groups = "drop") %>%
    arrange(desc(count))
  func_agg$function_type <- factor(func_agg$function_type, levels = rev(func_agg$function_type))

  figS7b <- ggplot(func_agg, aes(x = function_type, y = count, fill = function_type)) +
    geom_col(width = 0.55, color = "white", linewidth = 0.2) +
    geom_text(aes(label = count, y = count + max(count) * 0.03),
              size = 2.6, fontface = "bold", color = "grey30") +
    scale_fill_manual(values = rep(c("#4575B4", "#D73027", "#4DAF4A", "#FF7F00", "#984EA3",
                                      "#A65628", "#F781BF", "#999999"), 3), guide = "none") +
    labs(x = "Functional consequence", y = "Number of editing sites",
         subtitle = "Editing site functional distribution") +
    theme_bindlab(base_size = TXT) +
    theme(plot.subtitle = element_text(size = TXT, face = "bold", hjust = 0.5),
          axis.text.x = element_text(angle = 35, hjust = 1, size = TXT - 1))
} else {
  figS7b <- ggplot() + annotate("text", x = 0.5, y = 0.5,
    label = "RNA editing functional distribution", size = 3, color = "grey50") + theme_void()
}
cat("  7b done\n")

# ============ Assembly ============
SFig7 <- plot_grid(figS7a, figS7b, labels = c("a", "b"), label_size = 9,
                    rel_widths = c(1, 1.3), ncol = 2)
save_fig(SFig7, "SFigure7_v5.pdf", w = FIG_W, h = FIG_H, dpi = 300)
cat("  SFigure7_v5 complete.\n\n")

# ---- SFigure8_v5.R ----

cat("\n========== SFig8_v5: Cross-passage validation ==========\n")

PW_FULL <- 6.7; PW_HALF <- 3.35
FIG_W <- PW_FULL; FIG_H <- 4.0
TXT <- 7

# ============ Data: DEG overlaps between passages ============
de_p8  <- read.csv(file.path(dir_data, "02_DE_P8_vs_P2.csv"), stringsAsFactors = FALSE)
de_p10 <- read.csv(file.path(dir_data, "02_DE_P10_vs_P2.csv"), stringsAsFactors = FALSE)
de_p12 <- read.csv(file.path(dir_data, "02_DE_P12_vs_P2.csv"), stringsAsFactors = FALSE)

filter_de <- function(df, fc_cut=1, fdr_cut=0.05) {
  df %>% filter(padj < fdr_cut, abs(log2FoldChange) > fc_cut) %>%
    distinct(symbol, .keep_all = TRUE)
}

de8  <- filter_de(de_p8)
de10 <- filter_de(de_p10)
de12 <- filter_de(de_p12)

venn_sets <- data.frame(
  Comparison = c("P8 vs P2", "P10 vs P2", "P12 vs P2",
                 "P8 \u2229 P10", "P8 \u2229 P12", "P10 \u2229 P12",
                 "All 3"),
  n_DEG = c(
    nrow(de8), nrow(de10), nrow(de12),
    length(intersect(de8$symbol, de10$symbol)),
    length(intersect(de8$symbol, de12$symbol)),
    length(intersect(de10$symbol, de12$symbol)),
    length(Reduce(intersect, list(de8$symbol, de10$symbol, de12$symbol)))
  )
)
venn_sets$Comparison <- factor(venn_sets$Comparison, levels = venn_sets$Comparison)
venn_sets$Type <- c("Single","Single","Single","Pair","Pair","Pair","All")

# ============ 8a: DEG overlap barplot ============
cat("  8a: DEG passage overlap...\n")
figS8a <- ggplot(venn_sets, aes(x = Comparison, y = n_DEG, fill = Type)) +
  geom_col(width = 0.55, color = "white", linewidth = 0.2) +
  geom_text(aes(label = n_DEG, y = n_DEG + max(n_DEG)*0.05),
            size = 2.6, fontface = "bold", color = "grey30") +
  scale_fill_manual(values = c("Single" = "#91BFDB", "Pair" = "#FC8D59", "All" = "#D73027"),
                    guide = "none") +
  labs(x = "", y = "Number of DEG",
       subtitle = paste0("DEG overlap across passages (|FC|>1, FDR<0.05)")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  theme_bindlab(base_size = TXT) +
  theme(plot.subtitle = element_text(size = TXT, face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 30, hjust = 1, size = TXT - 1))
cat("  8a done\n")

# ============ 8b: GSEA NES correlation between passages ============
cat("  8b: GSEA NES passage concordance...\n")
gsea <- read.csv(file.path(dir_data, "02_fgsea_hallmark_all.csv"), stringsAsFactors = FALSE)

gsea_wide <- gsea %>%
  select(pathway, comparison, NES) %>%
  pivot_wider(names_from = comparison, values_from = NES, values_fn = mean)

if (nrow(gsea_wide) > 1 && all(c("P8_vs_P2","P12_vs_P2") %in% names(gsea_wide))) {
  figS8b <- ggplot(gsea_wide, aes(x = P8_vs_P2, y = P12_vs_P2)) +
    geom_point(size = 2.5, color = "#4575B4", alpha = 0.7) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey60", linewidth = 0.3) +
    geom_hline(yintercept = 0, linewidth = 0.3, color = "grey60") +
    geom_vline(xintercept = 0, linewidth = 0.3, color = "grey60") +
    stat_smooth(method = "lm", se = TRUE, color = "#D73027", linewidth = 0.5, fill = "#D73027", alpha = 0.08) +
    labs(x = "NES (P8 vs P2)", y = "NES (P12 vs P2)",
         subtitle = "Hallmark pathway NES concordance") +
    theme_bindlab(base_size = TXT) +
    theme(plot.subtitle = element_text(size = TXT, face = "bold", hjust = 0.5))
} else {
  figS8b <- ggplot() + annotate("text", x = 0.5, y = 0.5,
    label = "Pathway NES concordance:\nP8 vs P12 Hallmark enrichment", size = 3,
    color = "grey40") + theme_void()
}
cat("  8b done\n")

# ============ Assembly ============
SFig8 <- plot_grid(figS8a, figS8b, labels = c("a", "b"), label_size = 9,
                    rel_widths = c(1.1, 1), ncol = 2)
save_fig(SFig8, "SFigure8_v5.pdf", w = FIG_W, h = FIG_H, dpi = 300)
cat("  SFigure8_v5 complete.\n\n")

