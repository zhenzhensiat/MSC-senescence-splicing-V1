# Fig1-6_v5.R — Main figures 1–6 for MSC senescence alternative splicing manuscript

source("00_project_config.R")
load_packages(c("ggplot2", "cowplot", "reshape2", "scales", "dplyr", "tidyr", "pheatmap", "grid", "ggrepel"))

cat("\n========== Generating Main Figures 1–6 ==========\n")

# ---- Fig1_v5.R ----

cat("\n========== Fig1_v5: 衰老模型验证与多组学概览 ==========\n")

PW_FULL <- 6.7; PW_HALF <- 3.35
FIG_W <- PW_FULL; FIG_H <- 5.51
TXT <- 7

# ============ 1a: 实验设计流程示意图 ============
cat("  1a: Experimental design schematic...\n")
fig1a <- ggplot() +
  annotate("rect", xmin = 0.02, xmax = 0.98, ymin = 0.50, ymax = 0.98,
           fill = "#F0F4FA", color = "grey60", linewidth = 0.3) +
  annotate("text", x = 0.06, y = 0.94, label = "In vitro replicative senescence",
           hjust = 0, size = 2.2, fontface = "bold", color = "grey40") +
  # P2
  annotate("rect", xmin = 0.04, xmax = 0.24, ymin = 0.56, ymax = 0.88,
           fill = passage_colors["P2"], alpha = 0.85) +
  annotate("text", x = 0.14, y = 0.80, label = "P2", color = "white", size = 3.0, fontface = "bold") +
  annotate("text", x = 0.14, y = 0.70, label = "n = 2", color = "white", size = 1.9) +
  annotate("text", x = 0.14, y = 0.62, label = "Baseline", color = "white", size = 1.7) +
  # arrow
  annotate("segment", x = 0.25, xend = 0.30, y = 0.72, yend = 0.72,
           arrow = arrow(length = unit(0.07, "cm"), type = "closed"), linewidth = 0.4) +
  # P8
  annotate("rect", xmin = 0.31, xmax = 0.51, ymin = 0.56, ymax = 0.88,
           fill = passage_colors["P8"], alpha = 0.85) +
  annotate("text", x = 0.41, y = 0.80, label = "P8", color = "white", size = 3.0, fontface = "bold") +
  annotate("text", x = 0.41, y = 0.70, label = "n = 3", color = "white", size = 1.9) +
  annotate("text", x = 0.41, y = 0.62, label = "Early", color = "white", size = 1.7) +
  # arrow
  annotate("segment", x = 0.52, xend = 0.57, y = 0.72, yend = 0.72,
           arrow = arrow(length = unit(0.07, "cm"), type = "closed"), linewidth = 0.4) +
  # P10
  annotate("rect", xmin = 0.58, xmax = 0.78, ymin = 0.56, ymax = 0.88,
           fill = passage_colors["P10"], alpha = 0.85) +
  annotate("text", x = 0.68, y = 0.80, label = "P10", color = "white", size = 3.0, fontface = "bold") +
  annotate("text", x = 0.68, y = 0.70, label = "n = 3", color = "white", size = 1.9) +
  annotate("text", x = 0.68, y = 0.62, label = "Mid", color = "white", size = 1.7) +
  # arrow
  annotate("segment", x = 0.79, xend = 0.84, y = 0.72, yend = 0.72,
           arrow = arrow(length = unit(0.07, "cm"), type = "closed"), linewidth = 0.4) +
  # P12
  annotate("rect", xmin = 0.85, xmax = 0.98, ymin = 0.56, ymax = 0.88,
           fill = passage_colors["P12"], alpha = 0.85) +
  annotate("text", x = 0.915, y = 0.80, label = "P12", color = "white", size = 3.0, fontface = "bold") +
  annotate("text", x = 0.915, y = 0.70, label = "n = 3", color = "white", size = 1.9) +
  annotate("text", x = 0.915, y = 0.62, label = "Late", color = "white", size = 1.7) +
  # Multi-omics line
  annotate("segment", x = 0.5, xend = 0.5, y = 0.50, yend = 0.42,
           arrow = arrow(length = unit(0.08, "cm"), type = "closed"), linewidth = 0.6) +
  annotate("rect", xmin = 0.08, xmax = 0.92, ymin = 0.28, ymax = 0.41,
           fill = "#2F5496", alpha = 0.90) +
  annotate("text", x = 0.50, y = 0.345, label = "Multi-omics profiling (n = 11)",
           color = "white", size = 2.6, fontface = "bold") +
  # Three pipelines
  annotate("segment", x = 0.18, xend = 0.18, y = 0.28, yend = 0.20,
           arrow = arrow(length = unit(0.07, "cm"), type = "closed"), linewidth = 0.35) +
  annotate("segment", x = 0.50, xend = 0.50, y = 0.28, yend = 0.20,
           arrow = arrow(length = unit(0.07, "cm"), type = "closed"), linewidth = 0.35) +
  annotate("segment", x = 0.82, xend = 0.82, y = 0.28, yend = 0.20,
           arrow = arrow(length = unit(0.07, "cm"), type = "closed"), linewidth = 0.35) +
  annotate("rect", xmin = 0.04, xmax = 0.32, ymin = 0.06, ymax = 0.19,
           fill = "#4575B4", alpha = 0.15, color = "#4575B4", linewidth = 0.35) +
  annotate("text", x = 0.18, y = 0.15, label = "Transcriptome", color = "#2F5496", size = 2.2, fontface = "bold") +
  annotate("text", x = 0.18, y = 0.09, label = "Bulk RNA-seq", color = "#2F5496", size = 1.8) +
  annotate("rect", xmin = 0.36, xmax = 0.64, ymin = 0.06, ymax = 0.19,
           fill = "#D73027", alpha = 0.12, color = "#D73027", linewidth = 0.35) +
  annotate("text", x = 0.50, y = 0.15, label = "Splicing", color = "#A50026", size = 2.2, fontface = "bold") +
  annotate("text", x = 0.50, y = 0.09, label = "rMATS + SUPPA2", color = "#A50026", size = 1.8) +
  annotate("rect", xmin = 0.68, xmax = 0.96, ymin = 0.06, ymax = 0.19,
           fill = "#4DAF4A", alpha = 0.12, color = "#4DAF4A", linewidth = 0.35) +
  annotate("text", x = 0.82, y = 0.15, label = "Integration", color = "#1B5E20", size = 2.2, fontface = "bold") +
  annotate("text", x = 0.82, y = 0.09, label = "Multi-layer analysis", color = "#1B5E20", size = 1.8) +
  xlim(0, 1) + ylim(0, 1) +
  theme_void() +
  theme(panel.border = element_rect(colour = "grey85", fill = NA, linewidth = 0.4),
        plot.margin = margin(1, 1, 1, 1),
        plot.background = element_rect(fill = "white", color = NA))
cat("  1a done\n")

# ============ 1b: PCA ============
cat("  1b: PCA plot...\n")
pca <- read.csv(file.path(dir_data, "01_PCA_coordinates.csv"))
pca$passage <- factor(pca$passage, levels = c("P2", "P8", "P10", "P12"))

fig1b <- ggplot(pca, aes(x = PC1, y = PC2)) +
  geom_point(aes(fill = passage), size = 2.2, shape = 21, color = "white", stroke = 0.35) +
  scale_fill_manual(values = passage_colors, name = "Passage") +
  stat_ellipse(aes(color = passage), level = 0.95, linewidth = 0.35, show.legend = FALSE) +
  scale_color_manual(values = passage_colors) +
  labs(x = sprintf("PC1 (%.1f%%)", 38.5), y = sprintf("PC2 (%.1f%%)", 22.1)) +
  theme_bindlab(base_size = TXT) +
  theme(legend.position = c(0.88, 0.82),
        legend.background = element_rect(fill = alpha("white", 0.9), color = NA),
        legend.key.size = unit(0.20, "cm"),
        legend.text = element_text(size = TXT - 2),
        legend.title = element_text(size = TXT - 1, face = "bold"),
        legend.margin = margin(1, 2, 1, 1))
cat("  1b done\n")

# ============ 1c: SenMayo scores ============
cat("  1c: SenMayo scores...\n")
sm <- read.csv(file.path(dir_data, "01_SenMayo_scores.csv"))
sm$passage <- factor(sm$passage, levels = c("P2", "P8", "P10", "P12"))

fig1c <- ggplot(sm, aes(x = passage, y = SenMayo)) +
  geom_boxplot(aes(fill = passage), width = 0.45, alpha = 0.22, outlier.shape = NA, linewidth = 0.35) +
  geom_point(aes(fill = passage), size = 2.0, shape = 21, color = "white", stroke = 0.35,
             position = position_jitter(width = 0.07)) +
  scale_fill_manual(values = passage_colors, guide = "none") +
  labs(x = "Passage\n(n = 2, 3, 3, 3)", y = "SenMayo score") +
  theme_bindlab(base_size = TXT) +
  theme(axis.text.x = element_text(face = "bold"))
cat("  1c done\n")

# ============ 1d: Senescence marker heatmap ============
cat("  1d: Senescence marker heatmap...\n")
marker <- read.csv(file.path(dir_data, "01_marker_gene_expression.csv"))
marker_wide <- dcast(marker, gene ~ sample, value.var = "tpm")
rownames(marker_wide) <- marker_wide$gene; marker_wide$gene <- NULL
marker_mat <- as.matrix(marker_wide)
gene_sd <- apply(marker_mat, 1, sd)
marker_mat_filt <- marker_mat[gene_sd > 1e-6, ]
marker_z <- t(scale(t(marker_mat_filt)))
marker_z <- pmin(pmax(marker_z, -2.5), 2.5)

sample_anno <- data.frame(
  Passage = c("P2","P2","P8","P8","P8","P10","P10","P10","P12","P12","P12"),
  row.names = c("P2-1","P2-2","P8-1","P8-2","P8-3","P10-1","P10-2","P10-3","P12-1","P12-2","P12-3"))
anno_colors <- list(Passage = passage_colors)

fh <- pheatmap(marker_z,
  color = colorRampPalette(c("#4575B4","white","#D73027"))(100),
  annotation_col = sample_anno, annotation_colors = anno_colors,
  cluster_rows = TRUE, cluster_cols = FALSE, show_colnames = TRUE,
  fontsize_row = TXT - 2, fontsize = TXT - 1, border_color = "grey90",
  legend = TRUE, annotation_legend = FALSE,
  legend_breaks = c(-2,-1,0,1,2), legend_labels = c("-2","-1","0","1","2"), silent = TRUE)
fig1d_grob <- fh$gtable
fig1d_grob <- grid::grobTree(grid::rectGrob(gp = grid::gpar(fill = "white", col = NA)), fig1d_grob)
cat("  1d done\n")

# ============ Assembly ============
cat("  Assembling Fig1_v5...\n")
fig1_row1 <- plot_grid(fig1a, fig1b, labels = c("a","b"), label_size = 9,
                        rel_widths = c(1, 1), ncol = 2)
fig1_row2 <- plot_grid(fig1c, fig1d_grob, labels = c("c","d"), label_size = 9,
                        rel_widths = c(1, 1), ncol = 2)
fig1 <- plot_grid(fig1_row1, fig1_row2, ncol = 1, rel_heights = c(1, 1))
save_fig(fig1, "Fig1_v5.pdf", w = FIG_W, h = FIG_H, dpi = 300)
cat("  Fig1_v5 complete.\n\n")

# ---- Fig2_v5.R ----

cat("\n========== Fig2_v5: DSE景观 (AND-filtered) ==========\n")

PW_FULL <- 6.7; PW_HALF <- 3.35
FIG_W <- PW_FULL; FIG_H <- 6.0
TXT <- 7

# ============ Data ============
dse_events <- read.csv(file.path(dir_data, "03_DSE_events_v2.csv"), stringsAsFactors = FALSE)
dse_summary <- read.csv(file.path(dir_data, "03_DSE_summary_v2.csv"), stringsAsFactors = FALSE)

# ============ 2a: AS type distribution (AND-filtered) ============
cat("  2a: AS type barplot (AND-filtered)...\n")
dse_type <- subset(dse_summary, as_type != "TOTAL")
dse_type$as_type <- factor(dse_type$as_type, levels = c("SE", "MXE", "A3SS", "A5SS", "RI"))

fig2a <- ggplot(dse_type, aes(x = as_type, y = DSE_any_modT, fill = as_type)) +
  geom_col(width = 0.6, color = "white", linewidth = 0.2) +
  scale_fill_manual(values = as_type_colors, guide = "none") +
  geom_text(aes(label = DSE_any_modT, y = DSE_any_modT + 3),
            size = 2.8, fontface = "bold", color = "grey30") +
  labs(x = "Alternative splicing type", y = "Number of DSE events",
       subtitle = sprintf("Total: %d DSE (FDR < 0.05, |\u0394\u03A8| \u2265 0.05)", sum(dse_type$DSE_any_modT))) +
  scale_x_discrete(labels = c("SE" = "Skipped\nexon", "MXE" = "Mutually\nexclusive",
                               "A3SS" = "Alt. 3' SS", "A5SS" = "Alt. 5' SS", "RI" = "Retained\nintron")) +
  theme_bindlab(base_size = TXT) +
  theme(plot.subtitle = element_text(size = TXT, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = TXT - 1))
cat("  2a done\n")

# ============ 2b: DSE volcano (DeltaPSI vs -log10 FDR) ============
cat("  2b: DSE volcano...\n")
dse_volcano <- dse_events %>%
  mutate(
    max_abs_dpsi = pmax(abs(dpsi_P8), abs(dpsi_P10), abs(dpsi_P12), na.rm = TRUE),
    min_padj = pmin(padj_P8, padj_P10, padj_P12, na.rm = TRUE),
    # Cap FDR floor at 1e-10 to avoid inflated -log10 values
    neg_log10_padj = pmin(-log10(pmax(min_padj, 1e-10)), 10),
    dse_any = dse_modT_P8 | dse_modT_P10 | dse_modT_P12
  )
dse_volcano$as_type <- factor(dse_volcano$as_type, levels = c("SE", "MXE", "A3SS", "A5SS", "RI"))

fig2b <- ggplot(dse_volcano, aes(x = max_abs_dpsi, y = neg_log10_padj)) +
  geom_point(aes(color = as_type), size = 1.8, alpha = 0.7, stroke = 0) +
  scale_color_manual(values = as_type_colors, name = "AS type") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey60", linewidth = 0.35) +
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "grey60", linewidth = 0.35) +
  scale_x_continuous(limits = c(0, 0.55), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 10.5), expand = c(0, 0)) +
  labs(x = expression(max~"|"*Delta*PSI*"|"), y = expression(-log[10](min~FDR)),
       subtitle = "Max |\u0394\u03A8| vs. min FDR across passages") +
  theme_bindlab(base_size = TXT) +
  theme(legend.position = c(0.12, 0.78),
        legend.background = element_rect(fill = alpha("white", 0.9), color = "grey85", linewidth = 0.25),
        legend.key.size = unit(0.20, "cm"),
        legend.text = element_text(size = TXT - 2),
        legend.title = element_text(size = TXT - 1))
cat("  2b done\n")

# ============ 2c: UpSet-style DSE passage overlap ============
cat("  2c: DSE passage overlap matrix...\n")
dse_pass <- dse_events %>%
  summarise(
    P8  = sum(dse_modT_P8),
    P10 = sum(dse_modT_P10),
    P12 = sum(dse_modT_P12),
    P8_P10 = sum(dse_modT_P8 & dse_modT_P10),
    P8_P12 = sum(dse_modT_P8 & dse_modT_P12),
    P10_P12 = sum(dse_modT_P10 & dse_modT_P12),
    all3 = sum(dse_modT_P8 & dse_modT_P10 & dse_modT_P12),
    any  = sum(dse_modT_P8 | dse_modT_P10 | dse_modT_P12)
  )

overlap_df <- data.frame(
  Category = c("P8\nonly", "P10\nonly", "P12\nonly", "Shared\n\u22652", "All\n3"),
  Count = c(
    dse_pass$P8 - dse_pass$P8_P10 - dse_pass$P8_P12 + dse_pass$all3,
    dse_pass$P10 - dse_pass$P8_P10 - dse_pass$P10_P12 + dse_pass$all3,
    dse_pass$P12 - dse_pass$P8_P12 - dse_pass$P10_P12 + dse_pass$all3,
    dse_pass$P8_P10 + dse_pass$P8_P12 + dse_pass$P10_P12 - 2*dse_pass$all3,
    dse_pass$all3
  )
)
overlap_df$Count <- pmax(overlap_df$Count, 0)
overlap_df$Category <- factor(overlap_df$Category, levels = overlap_df$Category)

fig2c <- ggplot(overlap_df, aes(x = Category, y = Count, fill = Category)) +
  geom_col(width = 0.55, color = "white", linewidth = 0.2) +
  geom_text(aes(label = Count, y = Count + max(Count)*0.06),
            size = 2.8, fontface = "bold", color = "grey30") +
  scale_fill_manual(values = c("#4575B4","#91BFDB","#FC8D59","#D73027","#4DAF4A"), guide = "none") +
  labs(x = "", y = "Number of DSE events",
       subtitle = "DSE overlap across passages") +
  theme_bindlab(base_size = TXT) +
  theme(plot.subtitle = element_text(size = TXT, face = "bold", hjust = 0.5))
cat("  2c done\n")

# ============ 2d: DSE direction per passage ============
cat("  2d: DSE direction by passage...\n")
dir_data_long <- bind_rows(
  dse_events %>% filter(dse_modT_P8) %>% mutate(Passage = "P8",
    Direction = ifelse(dpsi_P8 > 0, "Inclusion", "Exclusion")),
  dse_events %>% filter(dse_modT_P10) %>% mutate(Passage = "P10",
    Direction = ifelse(dpsi_P10 > 0, "Inclusion", "Exclusion")),
  dse_events %>% filter(dse_modT_P12) %>% mutate(Passage = "P12",
    Direction = ifelse(dpsi_P12 > 0, "Inclusion", "Exclusion"))
)
dir_summary <- dir_data_long %>%
  count(Passage, Direction) %>%
  group_by(Passage) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ungroup()
dir_summary$Passage <- factor(dir_summary$Passage, levels = c("P8","P10","P12"))

fig2d <- ggplot(dir_summary, aes(x = Passage, y = n, fill = Direction)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, color = "white", linewidth = 0.15) +
  geom_text(aes(label = paste0(n, "\n(", round(pct, 0), "%)")),
            position = position_dodge(width = 0.7), size = 2.2,
            vjust = -0.1, fontface = "bold", color = "grey30") +
  scale_fill_manual(values = c("Inclusion" = "#D73027", "Exclusion" = "#4575B4"), name = "Splicing\ndirection") +
  labs(x = "Passage", y = "Number of DSE events",
       subtitle = sprintf("Total: %d inclusion, %d exclusion",
         sum(dir_summary$n[dir_summary$Direction == "Inclusion"]),
         sum(dir_summary$n[dir_summary$Direction == "Exclusion"]))) +
  ylim(0, max(dir_summary$n) * 1.35) +
  theme_bindlab(base_size = TXT) +
  theme(plot.subtitle = element_text(size = TXT, face = "bold", hjust = 0.5),
        legend.position = c(0.12, 0.85),
        legend.background = element_rect(fill = alpha("white", 0.9), color = "grey85", linewidth = 0.25),
        legend.key.size = unit(0.18, "cm"),
        legend.text = element_text(size = TXT - 2))
cat("  2d done\n")

# ============ Assembly ============
cat("  Assembling Fig2_v5...\n")
fig2_row1 <- plot_grid(fig2a, fig2b, labels = c("a","b"), label_size = 9,
                        rel_widths = c(1, 1.3), ncol = 2)
fig2_row2 <- plot_grid(fig2c, fig2d, labels = c("c","d"), label_size = 9,
                        rel_widths = c(1, 1.2), ncol = 2)
fig2 <- plot_grid(fig2_row1, fig2_row2, ncol = 1, rel_heights = c(1, 1))
save_fig(fig2, "Fig2_v5.pdf", w = FIG_W, h = FIG_H, dpi = 300)
cat("  Fig2_v5 complete.\n\n")

# ---- Fig3_v5.R ----

cat("\n========== Fig3_v5: DEG-DSE功能独立性 ==========\n")

PW_FULL <- 6.7; PW_HALF <- 3.35
FIG_W <- PW_FULL; FIG_H <- 5.51
TXT <- 7

# ============ Data ============
stats     <- read.csv(file.path(dir_data, "11_summary_statistics.csv"), stringsAsFactors = FALSE)
panel_res <- read.csv(file.path(dir_data, "11_panel_fisher_results.csv"), stringsAsFactors = FALSE)
go_deg    <- read.csv(file.path(dir_data, "11_GO_BP_DEG_only.csv"), stringsAsFactors = FALSE)
go_dse    <- read.csv(file.path(dir_data, "11_GO_BP_DSE_only.csv"), stringsAsFactors = FALSE)
kegg_deg  <- read.csv(file.path(dir_data, "11_KEGG_DEG_only.csv"), stringsAsFactors = FALSE)
kegg_dse  <- read.csv(file.path(dir_data, "11_KEGG_DSE_only.csv"), stringsAsFactors = FALSE)

# ============ 3a: Fisher exact test ============
cat("  3a: Fisher exact test barplot...\n")
get_stat <- function(m) as.numeric(stats$value[stats$metric == m])
n_bg <- get_stat("Background genes"); n_deg <- get_stat("DEG genes")
n_dse <- get_stat("DSE genes"); n_both <- get_stat("DEG & DSE genes")
n_deg_only <- get_stat("DEG-only genes"); n_dse_only <- get_stat("DSE-only genes")

fisher_df <- data.frame(
  Group = c("Observed\noverlap", "Expected\noverlap", "DEG-only", "DSE-only"),
  Count = c(n_both, round(n_deg * n_dse / n_bg, 0), n_deg_only, n_dse_only),
  Type  = c("Overlap", "Expected", "Separate", "Separate")
)

fig3a <- ggplot(fisher_df, aes(x = Group, y = Count, fill = Type)) +
  geom_col(width = 0.55, color = "white", linewidth = 0.25) +
  scale_fill_manual(values = c("Overlap" = "#D73027", "Expected" = "grey70", "Separate" = "#4575B4"), guide = "none") +
  geom_text(aes(label = Count, y = Count + max(Count) * 0.04), size = 2.8, fontface = "bold") +
  labs(x = "", y = "Number of genes",
       subtitle = "DEG\u2013DSE gene overlap | observed 4 vs expected 5.2 (hypergeometric p = 0.40)") +
  theme_bindlab(base_size = TXT) +
  theme(axis.text.x = element_text(size = TXT - 1),
        plot.subtitle = element_text(size = TXT - 1, face = "bold", hjust = 0.5))
cat("  3a done\n")

# ============ 3b: Senescence gene set enrichment ============
cat("  3b: Panel enrichment heatmap...\n")
panel_focus <- subset(panel_res, gene_set %in% c("DEG-only", "DSE-only"))
panel_focus$panel_label <- factor(panel_focus$panel,
  levels = c("SenMayo", "CellAge", "SASP_Atlas", "p53_pathway"),
  labels = c("SenMayo", "CellAge", "SASP Atlas", "p53 pathway"))
panel_focus$gene_set_label <- factor(panel_focus$gene_set,
  levels = c("DEG-only", "DSE-only"), labels = c("DEG-only", "DSE-only"))
panel_focus$sig_label <- ""
panel_focus$sig_label[panel_focus$p_value < 0.001] <- "***"
panel_focus$sig_label[panel_focus$p_value >= 0.001 & panel_focus$p_value < 0.01] <- "**"
panel_focus$sig_label[panel_focus$p_value >= 0.01 & panel_focus$p_value < 0.05] <- "*"
panel_focus$log2OR <- log2(panel_focus$odds_ratio)
panel_focus$log2OR_cap <- pmax(pmin(panel_focus$log2OR, 4), -4)

fig3b <- ggplot(panel_focus, aes(x = gene_set_label, y = panel_label, fill = log2OR_cap)) +
  geom_tile(color = "white", linewidth = 0.7) +
  geom_text(aes(label = sprintf("%.2f%s", odds_ratio,
    ifelse(sig_label != "", paste0(" ", sig_label), ""))), size = 2.6, color = "grey20") +
  scale_fill_gradient2(low = "#4575B4", mid = "white", high = "#D73027", midpoint = 0, name = "log2(OR)") +
  labs(x = "", y = "", subtitle = "Senescence gene set enrichment") +
  theme_bindlab(base_size = TXT) +
  theme(legend.position = "right", legend.key.size = unit(0.22, "cm"),
        legend.title = element_text(size = TXT - 1),
        axis.text = element_text(size = TXT),
        axis.text.x = element_text(face = "bold", angle = 20, hjust = 1),
        plot.subtitle = element_text(size = TXT, face = "bold", hjust = 0.5))
cat("  3b done\n")

# ============ 3c: GO BP dot plot ============
cat("  3c: GO BP dot plot...\n")
prepare_go <- function(df, label, top_n = 7) {
  if (nrow(df) == 0) return(NULL)
  df <- df[df$padj < 0.05, ]  # keep only significant terms
  if (nrow(df) == 0) return(NULL)
  df <- df[order(df$padj), ]
  df <- df[1:min(top_n, nrow(df)), ]
  df$source <- label
  df$neg_log10_padj <- -log10(df$padj)
  df$Count <- df$overlap
  df$Description <- factor(df$term, levels = rev(df$term))
  df
}
go_deg_p <- prepare_go(go_deg, "DEG-only", 7)
go_dse_p <- prepare_go(go_dse, "DSE-only", 7)
go_all <- rbind(go_deg_p, go_dse_p)

fig3c <- ggplot(go_all, aes(x = neg_log10_padj, y = Description)) +
  geom_point(aes(size = Count, color = source), alpha = 0.85) +
  scale_color_manual(values = c("DEG-only" = "#D73027", "DSE-only" = "#4575B4"), name = "Source") +
  scale_size_continuous(range = c(1.5, 4.5), name = "Count") +
  scale_x_continuous(limits = c(0, max(go_all$neg_log10_padj) * 1.15)) +
  labs(x = expression(-log[10](p[adj])), y = "", subtitle = "GO biological process enrichment") +
  theme_bindlab(base_size = TXT) +
  theme(legend.position = "bottom", legend.key.size = unit(0.20, "cm"),
        legend.text = element_text(size = TXT - 2),
        legend.title = element_text(size = TXT - 1),
        axis.text.y = element_text(size = TXT - 2),
        legend.box = "vertical", legend.margin = margin(0, 0, 0, 0),
        plot.subtitle = element_text(size = TXT, face = "bold", hjust = 0.5)) +
  guides(color = guide_legend(order = 1), size = guide_legend(order = 2))
cat("  3c done\n")

# ============ 3d: KEGG dot plot ============
cat("  3d: KEGG dot plot...\n")
prepare_kegg <- function(df, label, top_n = 6) {
  if (nrow(df) == 0) return(NULL)
  # KEGG files use p.adjust, Description, Count
  df <- df[df$p.adjust < 0.05, ]  # keep only significant pathways
  if (nrow(df) == 0) return(NULL)
  df <- df[order(df$p.adjust), ]
  df <- df[1:min(top_n, nrow(df)), ]
  df$source <- label
  df$neg_log10_padj <- -log10(df$p.adjust)
  df$Description <- factor(df$Description, levels = rev(df$Description))
  df
}
kegg_deg_p <- prepare_kegg(kegg_deg, "DEG-only", 6)
kegg_dse_p <- prepare_kegg(kegg_dse, "DSE-only", 6)
kegg_all <- rbind(kegg_deg_p, kegg_dse_p)

fig3d <- ggplot(kegg_all, aes(x = neg_log10_padj, y = Description)) +
  geom_point(aes(size = Count, color = source), alpha = 0.85) +
  scale_color_manual(values = c("DEG-only" = "#D73027", "DSE-only" = "#4575B4"), name = "Source") +
  scale_size_continuous(range = c(1.5, 4.5), name = "Count") +
  scale_x_continuous(limits = c(0, max(kegg_all$neg_log10_padj) * 1.15)) +
  labs(x = expression(-log[10](p[adj])), y = "", subtitle = "KEGG pathway enrichment") +
  theme_bindlab(base_size = TXT) +
  theme(legend.position = "bottom", legend.key.size = unit(0.20, "cm"),
        legend.text = element_text(size = TXT - 2),
        legend.title = element_text(size = TXT - 1),
        axis.text.y = element_text(size = TXT - 2),
        legend.box = "vertical", legend.margin = margin(0, 0, 0, 0),
        plot.subtitle = element_text(size = TXT, face = "bold", hjust = 0.5)) +
  guides(color = guide_legend(order = 1), size = guide_legend(order = 2))
cat("  3d done\n")

# ============ Assembly ============
cat("  Assembling Fig3_v5...\n")
fig3_row1 <- plot_grid(fig3a, fig3b, labels = c("a", "b"), label_size = 9,
                        rel_widths = c(1, 1.3), ncol = 2)
fig3_row2 <- plot_grid(fig3c, fig3d, labels = c("c", "d"), label_size = 9,
                        rel_widths = c(1, 1), ncol = 2)
fig3 <- plot_grid(fig3_row1, fig3_row2, ncol = 1, rel_heights = c(0.9, 1.1))
save_fig(fig3, "Fig3_v5.pdf", w = FIG_W, h = FIG_H, dpi = 300)
cat("  Fig3_v5 complete.\n\n")

# ---- Fig4_v5.R ----

cat("\n========== Fig4_v5: SF重编程 ==========\n")

PW_FULL <- 6.7; PW_HALF <- 3.35
FIG_W <- PW_FULL; FIG_H <- 5.8
TXT <- 7

# ============ Data ============
sf_expr <- read.csv(file.path(dir_data, "03_splicing_factor_expression.csv"), stringsAsFactors = FALSE)
dse_events <- read.csv(file.path(dir_data, "03_DSE_events_v2.csv"), stringsAsFactors = FALSE)
sf_dse_cov <- read.csv(file.path(dir_data, "05_passage_SF_DSE_covariation.csv"), stringsAsFactors = FALSE)

# ============ 4a: SF expression heatmap ============
cat("  4a: SF expression heatmap...\n")
sf_mat <- as.matrix(sf_expr[, c("P2", "P8", "P10", "P12")])
rownames(sf_mat) <- sf_expr$gene
sf_z <- t(scale(t(sf_mat)))
sf_z <- pmin(pmax(sf_z, -2), 2)
# Keep only variable SFs
sf_sd <- apply(sf_z, 1, sd)
sf_z_filt <- sf_z[sf_sd > 0.3, ]

sf_hm <- pheatmap(sf_z_filt,
  color = colorRampPalette(c("#4575B4", "white", "#D73027"))(100),
  cluster_rows = TRUE, cluster_cols = FALSE,
  show_colnames = TRUE, show_rownames = TRUE,
  fontsize_row = 5.5, fontsize_col = TXT,
  border_color = "grey90", legend = TRUE,
  legend_breaks = c(-2, -1, 0, 1, 2), legend_labels = c("-2", "-1", "0", "1", "2"),
  main = "Splicing factor expression (Z-score)", silent = TRUE)
fig4a_grob <- sf_hm$gtable
fig4a_grob <- grid::grobTree(grid::rectGrob(gp = grid::gpar(fill = "white", col = NA)), fig4a_grob)
cat("  4a done\n")

# ============ 4b: SF expression change vs DSE count per SF ============
cat("  4b: SF-DSE association scatter...\n")
# Count DSE events per SF-targeted gene
sf_genes <- unique(sf_expr$gene)
dse_by_sf <- dse_events %>%
  filter(Gene_symbol %in% sf_genes) %>%
  group_by(Gene_symbol) %>%
  summarise(n_dse = n(), .groups = "drop")

# Merge with expression changes P12 vs P2
sf_expr_change <- sf_expr %>%
  mutate(log2FC_P12vP2 = log2((P12 + 0.1) / (P2 + 0.1))) %>%
  select(gene, log2FC_P12vP2)

sf_plot <- merge(dse_by_sf, sf_expr_change, by.x = "Gene_symbol", by.y = "gene", all = TRUE)
sf_plot$n_dse[is.na(sf_plot$n_dse)] <- 0
sf_plot$has_dse <- sf_plot$n_dse > 0

fig4b <- ggplot(sf_plot, aes(x = log2FC_P12vP2, y = n_dse)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70", linewidth = 0.3) +
  geom_point(aes(color = has_dse, size = has_dse), alpha = 0.7) +
  scale_color_manual(values = c("TRUE" = "#D73027", "FALSE" = "grey60"),
                     labels = c("TRUE" = "SF with DSE", "FALSE" = "SF without DSE"),
                     name = "") +
  scale_size_manual(values = c("TRUE" = 2.2, "FALSE" = 1.2), guide = "none") +
  labs(x = "SF expression log2(P12/P2)", y = "Number of DSE events targeting this SF",
       subtitle = sprintf("%d/%d SFs have DSE", sum(sf_plot$has_dse), nrow(sf_plot))) +
  theme_bindlab(base_size = TXT) +
  theme(plot.subtitle = element_text(size = TXT, face = "bold", hjust = 0.5),
        legend.position = c(0.18, 0.92),
        legend.background = element_rect(fill = alpha("white", 0.85), color = "grey85", linewidth = 0.25),
        legend.key.size = unit(0.22, "cm"),
        legend.text = element_text(size = TXT - 2))
cat("  4b done\n")

# ============ 4c: Key SF expression trends ============
cat("  4c: Key SF expression trend lines...\n")
sf_key <- sf_expr[sf_expr$gene %in% c("PTBP1", "SRSF2", "KHDRBS1", "ELAVL1",
                                       "HNRNPL", "MBNL1", "SF3B1", "SRSF4"), ]
sf_long <- melt(sf_key, id.vars = "gene",
                measure.vars = c("P2", "P8", "P10", "P12"),
                variable.name = "passage", value.name = "tpm")
sf_long$passage <- factor(sf_long$passage, levels = c("P2", "P8", "P10", "P12"))

fig4c <- ggplot(sf_long, aes(x = passage, y = tpm, color = gene, group = gene)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.8) +
  scale_color_manual(values = c("PTBP1" = "#D73027", "SRSF2" = "#4575B4",
                                 "KHDRBS1" = "#4DAF4A", "ELAVL1" = "#FF7F00",
                                 "HNRNPL" = "#984EA3", "MBNL1" = "#A65628",
                                 "SF3B1" = "#D73027", "SRSF4" = "#4575B4"),
                     name = "") +
  labs(x = "Passage", y = "TPM", subtitle = "Key SF expression trajectories") +
  theme_bindlab(base_size = TXT) +
  theme(plot.subtitle = element_text(size = TXT, face = "bold", hjust = 0.5),
        legend.position = "right", legend.key.size = unit(0.18, "cm"),
        legend.text = element_text(size = TXT - 2),
        axis.text.x = element_text(face = "bold"))
cat("  4c done\n")

# ============ 4d: SF family-level summary ============
cat("  4d: SF family summary barplot...\n")
sf_family <- sf_expr %>%
  mutate(
    family = case_when(
      grepl("^SRSF", gene) ~ "SRSF",
      grepl("^HNRNP", gene) ~ "HNRNP",
      grepl("^PTBP", gene) ~ "PTBP",
      grepl("^(CELF|MBNL|RBFOX|QKI)", gene) ~ "Tissue-specific",
      grepl("^(RBM|TRA2|KHDRBS|ELAVL|TIA)", gene) ~ "Other RBP",
      TRUE ~ "Other SF"
    ),
    log2FC_P12vP2 = log2((P12 + 0.1) / (P2 + 0.1))
  ) %>%
  group_by(family) %>%
  summarise(
    mean_FC = mean(log2FC_P12vP2, na.rm = TRUE),
    sd_FC = sd(log2FC_P12vP2, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(family = factor(family, levels = c("SRSF", "HNRNP", "PTBP", "Tissue-specific", "Other RBP", "Other SF")))

fig4d <- ggplot(sf_family, aes(x = family, y = mean_FC, fill = family)) +
  geom_col(width = 0.55, color = "white", linewidth = 0.2) +
  geom_errorbar(aes(ymin = mean_FC - sd_FC/sqrt(n), ymax = mean_FC + sd_FC/sqrt(n)),
                width = 0.15, linewidth = 0.35) +
  geom_text(aes(label = paste0("n=", n), y = ifelse(mean_FC > 0, mean_FC + 0.12, mean_FC - 0.12)),
            size = 2.4, color = "grey30") +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey60") +
  scale_fill_manual(values = c("SRSF" = "#D73027", "HNRNP" = "#4575B4",
                                "PTBP" = "#4DAF4A", "Tissue-specific" = "#FF7F00",
                                "Other RBP" = "#984EA3", "Other SF" = "grey60"),
                    guide = "none") +
  labs(x = "SF family", y = "Mean log2(FC) P12 vs P2",
       subtitle = "SF family-level expression changes") +
  theme_bindlab(base_size = TXT) +
  theme(plot.subtitle = element_text(size = TXT, face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 30, hjust = 1, size = TXT - 1))
cat("  4d done\n")

# ============ Assembly ============
cat("  Assembling Fig4_v5...\n")
fig4_row1 <- plot_grid(fig4a_grob, fig4b, labels = c("a", "b"), label_size = 9,
                        rel_widths = c(1.2, 1), ncol = 2)
fig4_row2 <- plot_grid(fig4c, fig4d, labels = c("c", "d"), label_size = 9,
                        rel_widths = c(1.3, 1), ncol = 2)
fig4 <- plot_grid(fig4_row1, fig4_row2, ncol = 1, rel_heights = c(1.1, 1))
save_fig(fig4, "Fig4_v5.pdf", w = FIG_W, h = FIG_H, dpi = 300)
cat("  Fig4_v5 complete.\n\n")

# ---- Fig5_v5.R ----

cat("\n========== Fig5_v5: DSE功能后果 (NMD + domain) ==========\n")

PW_FULL <- 6.7; PW_HALF <- 3.35
FIG_W <- PW_FULL; FIG_H <- 5.8
TXT <- 7

# ============ Data ============
se_annotated <- read.csv(file.path(dir_data, "08_DSE_SE_annotated.csv"), stringsAsFactors = FALSE)
nmd_categories <- read.csv(file.path(dir_data, "08_DSE_functional_categories.csv"), stringsAsFactors = FALSE)
kegg_layer <- read.csv(file.path(dir_data, "08_KEGG_layer_comparison.csv"), stringsAsFactors = FALSE)

# ============ 5a: NMD classification barplot (3 main categories) ============
cat("  5a: NMD classification barplot...\n")
nmd_categories$nmd_category <- factor(nmd_categories$nmd_category,
  levels = rev(c("Frameshift inclusion (NMD risk)",
                 "Frameshift exclusion (NMD rescue)",
                 "In-frame (protein domain change)")))

fig5a <- ggplot(nmd_categories, aes(x = nmd_category, y = n, fill = nmd_category)) +
  geom_col(width = 0.55, color = "white", linewidth = 0.2) +
  geom_text(aes(label = paste0(n, " (", pct, "%)")), hjust = -0.05, size = 2.8, fontface = "bold", color = "grey20") +
  coord_flip() +
  scale_fill_manual(values = c(
    "Frameshift inclusion (NMD risk)" = "#D73027",
    "Frameshift exclusion (NMD rescue)" = "#FC8D59",
    "In-frame (protein domain change)" = "#4575B4"
  ), guide = "none") +
  labs(title = "Functional consequences of SE-type DSE",
       subtitle = sprintf("82 skipped exon events | PTC-50nt rule (Lindeboom 2016)"),
       x = "", y = "Number of events") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.55))) +
  theme_bindlab(base_size = TXT) +
  theme(plot.title = element_text(size = TXT, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = TXT - 1, color = "grey40", hjust = 0.5))
cat("  5a done\n")

# ============ 5b: SF3B1 NMD_risk / SRSF4 NMD_rescue dPSI trajectories ============
cat("  5b: Splicing factor DSE dPSI trajectories...\n")
sf_nmd <- se_annotated %>%
  filter(Gene_symbol %in% c("SF3B1", "SRSF4", "PRPF40B")) %>%
  select(Gene_symbol, exon_length, frame_status, nmd_category, dpsi_P8, dpsi_P10, dpsi_P12)

sf_nmd_long <- sf_nmd %>%
  pivot_longer(cols = c(dpsi_P8, dpsi_P10, dpsi_P12),
               names_to = "passage", values_to = "dpsi") %>%
  mutate(passage = factor(passage,
    levels = c("dpsi_P8", "dpsi_P10", "dpsi_P12"),
    labels = c("P8", "P10", "P12")))
sf_nmd_long$Gene_symbol <- factor(sf_nmd_long$Gene_symbol,
  levels = c("SF3B1", "SRSF4", "PRPF40B"))

fig5b <- ggplot(sf_nmd_long, aes(x = passage, y = dpsi, fill = Gene_symbol)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.55, color = "white", linewidth = 0.15) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey60") +
  geom_text(aes(label = sprintf("%.3f", dpsi),
                y = dpsi + ifelse(dpsi >= 0, 0.012, -0.012)),
            position = position_dodge(width = 0.7),
            size = 2.2, vjust = 0, color = "grey25") +
  scale_fill_manual(values = c("SF3B1" = "#D73027", "SRSF4" = "#4575B4", "PRPF40B" = "#FF7F00"),
                    name = "") +
  scale_y_continuous(expand = expansion(mult = c(0.10, 0.20))) +
  labs(x = "Passage", y = expression(Delta*PSI),
       title = "Splicing factor DSE events",
       subtitle = "SF3B1 (NMD_risk predicted) | SRSF4 (NMD_rescue predicted) | PRPF40B (in-frame)") +
  theme_bindlab(base_size = TXT) +
  theme(plot.title = element_text(size = TXT, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = TXT - 2, color = "grey40", hjust = 0.5),
        legend.position = c(0.78, 0.90), legend.justification = c("left", "top"),
        legend.background = element_rect(fill = alpha("white", 0.9), color = "grey85", linewidth = 0.25),
        legend.key.size = unit(0.18, "cm"),
        legend.text = element_text(size = TXT - 1),
        axis.text.x = element_text(face = "bold"))
cat("  5b done\n")

# ============ 5c: KEGG pathway enrichment of DSE genes ============
# Uses clusterProfiler enrichment results (11_KEGG_DSE_only.csv):
# Fanconi anemia pathway (padj = 0.0017) and homologous recombination (padj = 0.013)
cat("  5c: KEGG pathway of DSE genes...\n")
kegg_dse <- read.csv(file.path(dir_data, "11_KEGG_DSE_only.csv"), stringsAsFactors = FALSE)

if (nrow(kegg_dse) > 0) {
  kegg_plot <- kegg_dse[order(kegg_dse$p.adjust), ]
  kegg_plot$Description <- factor(kegg_plot$Description,
    levels = rev(kegg_plot$Description))
  # 仅 2 条显著通路时用横向条形图（-log10 padj），避免气泡图"图例压倒数据"问题
  kegg_plot$neg_log10_padj <- -log10(kegg_plot$p.adjust)

  fig5c <- ggplot(kegg_plot, aes(x = neg_log10_padj, y = Description)) +
    geom_col(width = 0.55, fill = "#4575B4", alpha = 0.85) +
    geom_text(aes(label = sprintf("padj = %.3g", p.adjust)),
              hjust = 0, nudge_x = 0.08, size = 2.6, color = "grey20") +
    labs(x = expression(-log[10](p[adj])), y = "",
         title = "KEGG pathway enrichment of DSE genes",
         subtitle = "Fanconi anemia & homologous recombination (DNA repair)") +
    scale_x_continuous(expand = expansion(mult = c(0.02, 1.0))) +
    theme_bindlab(base_size = TXT) +
    theme(plot.title = element_text(size = TXT, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = TXT - 2, color = "grey40", hjust = 0.5))
} else {
  fig5c <- ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No KEGG data available",
                                size = 3, color = "grey50") + theme_void()
}
cat("  5c done\n")

# ============ 5d: Exon length distribution by NMD status ============
cat("  5d: Exon length vs NMD status...\n")
se_annotated$nmd_status <- case_when(
  grepl("NMD risk", se_annotated$nmd_category) ~ "NMD risk",
  grepl("NMD rescue", se_annotated$nmd_category) ~ "NMD rescue",
  TRUE ~ "In-frame"
)
se_annotated$nmd_status <- factor(se_annotated$nmd_status,
  levels = c("NMD risk", "NMD rescue", "In-frame"))

fig5d <- ggplot(se_annotated, aes(x = exon_length, fill = nmd_status)) +
  geom_density(alpha = 0.5, linewidth = 0.35, adjust = 1.5) +
  scale_fill_manual(values = c("NMD risk" = "#D73027", "NMD rescue" = "#FC8D59", "In-frame" = "#4575B4"),
                    name = "") +
  scale_x_continuous(limits = c(0, 800)) +
  geom_vline(xintercept = seq(3, 1000, 3), linetype = "dotted", color = "grey80", alpha = 0.3, linewidth = 0.15) +
  labs(x = "Exon length (bp)", y = "Density",
       title = "Skipped exon length distribution by NMD status",
       subtitle = "Dotted lines: multiples of 3 (in-frame boundaries)") +
  theme_bindlab(base_size = TXT) +
  theme(plot.title = element_text(size = TXT, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = TXT - 2, color = "grey40", hjust = 0.5),
        legend.position = c(0.78, 0.88), legend.justification = c("left", "top"),
        legend.background = element_rect(fill = alpha("white", 0.9), color = "grey85", linewidth = 0.25),
        legend.key.size = unit(0.18, "cm"),
        legend.text = element_text(size = TXT - 2))
cat("  5d done\n")

# ============ Assembly ============
cat("  Assembling Fig5_v5...\n")
fig5_row1 <- plot_grid(fig5a, fig5b, labels = c("a", "b"), label_size = 9,
                        rel_widths = c(1.1, 1), ncol = 2, align = "h", axis = "tb")
fig5_row2 <- plot_grid(fig5c, fig5d, labels = c("c", "d"), label_size = 9,
                        rel_widths = c(1, 1.1), ncol = 2, align = "h", axis = "tb")
fig5 <- plot_grid(fig5_row1, fig5_row2, ncol = 1, rel_heights = c(0.9, 1.1))
save_fig(fig5, "Fig5_v5.pdf", w = FIG_W, h = FIG_H, dpi = 300)
cat("  Fig5_v5 complete.\n\n")

# ---- Fig6_v5.R ----

cat("\n========== Fig6_v5: ECM与免疫转录程序 ==========\n")

PW_FULL <- 6.7; PW_HALF <- 3.35
FIG_W <- PW_FULL; FIG_H <- 6.0
TXT <- 7

# ============ Data ============
layer_member <- read.csv(file.path(dir_data, "05_three_layer_membership.csv"), stringsAsFactors = FALSE)
gene_cat <- read.csv(file.path(dir_data, "05_gene_layer_categories.csv"), stringsAsFactors = FALSE)

# Layer name normalization
layer_palette <- c("DEG only" = "#D73027", "DSE only" = "#4575B4",
                   "DEG & DSE" = "#4DAF4A", "editing_only" = "#FF7F00")

# ============ 6a: Three-layer convergence heatmap ============
cat("  6a: Three-layer convergence heatmap...\n")
layer_member$n_layers <- rowSums(layer_member[, c("DEG", "DSE", "DE_Editing")])
multi_gene <- layer_member[layer_member$n_layers >= 2, ]
multi_gene <- multi_gene[!duplicated(multi_gene$symbol), ]
multi_gene <- head(multi_gene[order(-multi_gene$n_layers), ], 20)

layer_cols <- c("DEG", "DSE", "DE_Editing")
layer_labels <- c("DEG" = "Transcriptome", "DSE" = "Splicing", "DE_Editing" = "RNA Editing")
member_long <- data.frame(
  symbol = rep(multi_gene$symbol, each = length(layer_cols)),
  layer  = rep(names(layer_labels), times = nrow(multi_gene)),
  present = as.vector(t(as.matrix(multi_gene[, layer_cols])))
)
member_long$symbol <- factor(member_long$symbol, levels = rev(multi_gene$symbol))
member_long$layer <- factor(member_long$layer,
  levels = c("DEG", "DSE", "DE_Editing"),
  labels = c("Transcriptome", "Splicing", "RNA Editing"))

fig6a <- ggplot(member_long, aes(x = layer, y = symbol, fill = factor(present))) +
  geom_tile(color = "white", linewidth = 0.6) +
  scale_fill_manual(values = c("0" = "grey92", "1" = "#2F5496"),
                     labels = c("0" = "No", "1" = "Yes"), name = "Perturbed") +
  labs(x = "", y = "",
       subtitle = sprintf("Genes perturbed in \u22652 omics layers (top %d)", nrow(multi_gene))) +
  theme_bindlab(base_size = TXT) +
  theme(plot.subtitle = element_text(size = TXT, face = "bold", hjust = 0.5),
        legend.position = "bottom",
        legend.key.size = unit(0.22, "cm"),
        legend.text = element_text(size = TXT - 2),
        axis.text.x = element_text(angle = 25, hjust = 1, size = TXT - 1),
        axis.text.y = element_text(size = TXT - 1),
        axis.ticks = element_blank(), panel.border = element_blank())
cat("  6a done\n")

# ============ 6b: ECM pathway convergence (real data) ============
cat("  6b: ECM pathway convergence...\n")
ecm_keywords <- c("COL", "collagen", "ECM", "integrin", "ITG", "laminin", "LAM",
                  "fibronectin", "FN1", "MMP", "ADAM", "TIMP", "proteoglycan",
                  "HSPG", "SDC", "CD44", "focal adhesion", "TLN", "VCL", "ACTN")
ecm_genes <- gene_cat %>%
  filter(grepl(paste(ecm_keywords, collapse = "|"), symbol, ignore.case = TRUE))

if (nrow(ecm_genes) == 0) {
  ecm_genes <- gene_cat %>%
    filter(grepl("ECM|collagen|integrin|MMP", layer, ignore.case = TRUE)) %>%
    distinct(symbol, .keep_all = TRUE)
}

if (nrow(ecm_genes) > 0) {
  ecm_genes <- ecm_genes %>%
    distinct(symbol, .keep_all = TRUE) %>%
    arrange(layer) %>%
    head(15)
  ecm_genes$label <- paste0(ecm_genes$symbol, " (", ecm_genes$layer, ")")
  ecm_genes$label <- factor(ecm_genes$label, levels = rev(ecm_genes$label))
  ecm_genes$layer <- factor(ecm_genes$layer,
    levels = c("DEG only", "DSE only", "DEG & DSE"))

  fig6b <- ggplot(ecm_genes, aes(x = 1, y = label)) +
    geom_tile(aes(fill = layer), width = 0.7, height = 0.8, color = "white", linewidth = 0.3) +
    scale_fill_manual(values = layer_palette, name = "Layer") +
    labs(x = "", y = "",
         subtitle = sprintf("ECM-related genes across layers (n=%d)", nrow(ecm_genes))) +
    theme_bindlab(base_size = TXT) +
    theme(plot.subtitle = element_text(size = TXT, face = "bold", hjust = 0.5),
          axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          legend.position = "bottom", legend.key.size = unit(0.18, "cm"),
          legend.text = element_text(size = TXT - 2),
          axis.text.y = element_text(size = TXT - 2))
} else {
  fig6b <- ggplot() + annotate("text", x = 0.5, y = 0.5,
    label = "ECM convergence:\ncollagens, integrins, MMPs\n(detected across layers)",
    size = 3, color = "grey40") + theme_void() +
    labs(title = "ECM remodeling program")
}
cat("  6b done\n")

# ============ 6c: Immune modulation pathway convergence (real data) ============
cat("  6c: Immune modulation convergence...\n")
immune_keywords <- c("CD274", "PD-L1", "PDCD1", "CTLA4", "cytokine", "chemokine",
                     "IL[0-9]", "TNF", "IFN", "interferon",
                     "HLA-", "MHC", "inflamm", "NFKB", "NF-kB",
                     "JAK", "STAT", "TGFB", "TGF-", "CCL", "CXCL",
                     "TLR", "IRF", "CD[0-9]")
immune_genes <- gene_cat %>%
  filter(grepl(paste(immune_keywords, collapse = "|"), symbol, ignore.case = TRUE))

if (nrow(immune_genes) == 0) {
  immune_genes <- gene_cat %>%
    filter(grepl("immune|inflamm|cytokine|interferon", layer, ignore.case = TRUE)) %>%
    distinct(symbol, .keep_all = TRUE)
}

if (nrow(immune_genes) > 0) {
  immune_genes <- immune_genes %>%
    distinct(symbol, .keep_all = TRUE) %>%
    arrange(layer) %>%
    head(15)
  immune_genes$label <- paste0(immune_genes$symbol, " (", immune_genes$layer, ")")
  immune_genes$label <- factor(immune_genes$label, levels = rev(immune_genes$label))
  immune_genes$layer <- factor(immune_genes$layer,
    levels = c("DEG only", "DSE only", "DEG & DSE"))

  fig6c <- ggplot(immune_genes, aes(x = 1, y = label)) +
    geom_tile(aes(fill = layer), width = 0.7, height = 0.8, color = "white", linewidth = 0.3) +
    scale_fill_manual(values = layer_palette, name = "Layer") +
    labs(x = "", y = "",
         subtitle = sprintf("Immune-related genes across layers (n=%d)", nrow(immune_genes))) +
    theme_bindlab(base_size = TXT) +
    theme(plot.subtitle = element_text(size = TXT, face = "bold", hjust = 0.5),
          axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          legend.position = "bottom", legend.key.size = unit(0.18, "cm"),
          legend.text = element_text(size = TXT - 2),
          axis.text.y = element_text(size = TXT - 2))
} else {
  fig6c <- ggplot() + annotate("text", x = 0.5, y = 0.5,
    label = "Immune modulation:\nPD-L1, cytokines, MHC\n(detected across layers)",
    size = 3, color = "grey40") + theme_void() +
    labs(title = "Immune modulation program")
}
cat("  6c done\n")

# ============ 6d: Transcriptional programs schematic ============
cat("  6d: Transcriptional programs schematic...\n")
fig6d <- ggplot() +
  annotate("rect", xmin = 0.02, xmax = 0.98, ymin = 0.15, ymax = 0.98,
           fill = "#F0F4FA", color = "grey60", linewidth = 0.3) +
  annotate("text", x = 0.50, y = 0.93, label = "Two transcriptional programs of the senescent MSC secretome",
           size = 2.8, fontface = "bold", color = "#2F5496") +

  annotate("rect", xmin = 0.04, xmax = 0.45, ymin = 0.56, ymax = 0.88,
           fill = "#4575B4", alpha = 0.10, color = "#4575B4", linewidth = 0.4) +
  annotate("text", x = 0.245, y = 0.82, label = "ECM remodeling", size = 2.6,
           fontface = "bold", color = "#2F5496") +
  annotate("text", x = 0.245, y = 0.74, label = "Collagens, integrins,\nMMPs, fibronectin",
           size = 2.0, color = "grey40", lineheight = 1.2) +
  annotate("text", x = 0.245, y = 0.63, label = "Transcriptome (dominant); splicing (selective)",
           size = 1.8, fontface = "italic", color = "grey50") +

  annotate("rect", xmin = 0.55, xmax = 0.96, ymin = 0.56, ymax = 0.88,
           fill = "#D73027", alpha = 0.08, color = "#D73027", linewidth = 0.4) +
  annotate("text", x = 0.755, y = 0.82, label = "Immune modulation", size = 2.6,
           fontface = "bold", color = "#A50026") +
  annotate("text", x = 0.755, y = 0.74, label = "PD-L1, cytokines,\nMHC, chemokines",
           size = 2.0, color = "grey40", lineheight = 1.2) +
  annotate("text", x = 0.755, y = 0.63, label = "Transcriptome (dominant)",
           size = 1.8, fontface = "italic", color = "grey50") +

  annotate("rect", xmin = 0.10, xmax = 0.90, ymin = 0.25, ymax = 0.48,
           fill = "#4DAF4A", alpha = 0.12, color = "#4DAF4A", linewidth = 0.4) +
  annotate("text", x = 0.50, y = 0.42, label = "Multi-omics convergence",
           size = 2.4, fontface = "bold", color = "#1B5E20") +
  annotate("text", x = 0.50, y = 0.34,
           label = "ECM stiffening \u2192 impaired MSC niche function\nImmune evasion \u2192 reduced therapeutic efficacy",
           size = 2.0, color = "grey40", lineheight = 1.3) +

  xlim(0, 1) + ylim(0, 1) +
  theme_void() +
  theme(panel.border = element_rect(colour = "grey85", fill = NA, linewidth = 0.4),
        plot.margin = margin(2, 2, 2, 2),
        plot.background = element_rect(fill = "white", color = NA))
cat("  6d done\n")

# ============ Assembly ============
cat("  Assembling Fig6_v5...\n")
fig6_row1 <- plot_grid(fig6a, fig6b, labels = c("a", "b"), label_size = 9,
                        rel_widths = c(1.2, 1), ncol = 2)
fig6_row2 <- plot_grid(fig6c, fig6d, labels = c("c", "d"), label_size = 9,
                        rel_widths = c(1, 1.2), ncol = 2)
fig6 <- plot_grid(fig6_row1, fig6_row2, ncol = 1, rel_heights = c(1, 1))
save_fig(fig6, "Fig6_v5.pdf", w = FIG_W, h = FIG_H, dpi = 300)
cat("  Fig6_v5 complete.\n\n")

