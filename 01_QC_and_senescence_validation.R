# 01_QC_and_senescence_validation.R — Data quality assessment and senescence trajectory validation

source("00_project_config.R")

load_packages(c("ggplot2", "pheatmap", "RColorBrewer", "ggrepel",
                "reshape2", "dplyr", "tidyr", "GSVA", "matrixStats"))

cat("========== Step 1: QC & Senescence Validation ==========\n\n")

# --------------------------------------------------
# 1.1 读取数据
# --------------------------------------------------
cat("[1/7] 读取表达矩阵...\n")

raw <- read.delim(f_expression, sep = "\t", header = TRUE,
                  stringsAsFactors = FALSE, check.names = FALSE)

tpm_cols   <- grep("_tpm$", colnames(raw), value = TRUE)
count_cols <- grep("_count$", colnames(raw), value = TRUE)

tpm_mat   <- as.matrix(raw[, tpm_cols])
count_mat <- as.matrix(raw[, count_cols])
rownames(tpm_mat)   <- raw$id
rownames(count_mat) <- raw$id

gene_info <- raw[, c("id", "Symbol")]

# 统一列名为样品ID
colnames(tpm_mat)   <- gsub("_tpm$", "", tpm_cols)
colnames(count_mat) <- gsub("_count$", "", count_cols)

# 确保列顺序与sample_table一致
tpm_mat   <- tpm_mat[, sample_table$sample_id]
count_mat <- count_mat[, sample_table$sample_id]

cat("  矩阵维度:", nrow(tpm_mat), "genes ×", ncol(tpm_mat), "samples\n\n")

# --------------------------------------------------
# 1.2 基本QC统计
# --------------------------------------------------
cat("[2/7] 基本QC统计...\n")

qc_stats <- data.frame(
  sample      = sample_table$sample_id,
  passage     = sample_table$passage,
  total_count = colSums(count_mat),
  genes_detected_count = colSums(count_mat > 0),
  genes_tpm_ge1  = colSums(tpm_mat >= 1),
  genes_tpm_ge5  = colSums(tpm_mat >= 5),
  median_tpm_expressed = apply(tpm_mat, 2, function(x) median(x[x > 0])),
  stringsAsFactors = FALSE
)
save_data(qc_stats, "01_QC_statistics.csv")

# 图: 文库大小
p <- ggplot(qc_stats, aes(x = sample, y = total_count / 1e6, fill = passage)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = passage_colors) +
  labs(title = "Library Size", x = "", y = "Total Counts (million)") +
  theme_bindlab(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_fig(p, "01_library_size.pdf", w = 8, h = 5)

# 图: 检测基因数
p <- ggplot(qc_stats, aes(x = sample, y = genes_tpm_ge1, fill = passage)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = passage_colors) +
  labs(title = "Genes Detected (TPM >= 1)", x = "", y = "Number of Genes") +
  theme_bindlab(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_fig(p, "01_genes_detected.pdf", w = 8, h = 5)

cat("\n")

# --------------------------------------------------
# 1.3 PCA分析
# --------------------------------------------------
cat("[3/7] PCA分析...\n")

expressed <- rowSums(tpm_mat >= 1) >= 2
log_tpm   <- log2(tpm_mat[expressed, ] + 1)
cat("  用于PCA的基因数:", nrow(log_tpm), "\n")

pca_res <- prcomp(t(log_tpm), center = TRUE, scale. = TRUE)
var_pct <- round(summary(pca_res)$importance[2, 1:5] * 100, 1)

pca_df <- data.frame(
  sample  = sample_table$sample_id,
  passage = sample_table$passage,
  PC1 = pca_res$x[, 1],
  PC2 = pca_res$x[, 2],
  PC3 = pca_res$x[, 3]
)
save_data(pca_df, "01_PCA_coordinates.csv")

# 图: PC1 vs PC2
p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = passage)) +
  geom_point(size = 5) +
  geom_text_repel(aes(label = sample), size = 3.5, max.overlaps = 20) +
  scale_color_manual(values = passage_colors) +
  labs(title = "PCA: PC1 vs PC2",
       x = paste0("PC1 (", var_pct[1], "%)"),
       y = paste0("PC2 (", var_pct[2], "%)")) +
  theme_bindlab(base_size = 12)
save_fig(p, "01_PCA_PC1_PC2.pdf", w = 8, h = 6)

# 图: PC1 vs PC3
p <- ggplot(pca_df, aes(x = PC1, y = PC3, color = passage)) +
  geom_point(size = 5) +
  geom_text_repel(aes(label = sample), size = 3.5, max.overlaps = 20) +
  scale_color_manual(values = passage_colors) +
  labs(title = "PCA: PC1 vs PC3",
       x = paste0("PC1 (", var_pct[1], "%)"),
       y = paste0("PC3 (", var_pct[3], "%)")) +
  theme_bindlab(base_size = 12)
save_fig(p, "01_PCA_PC1_PC3.pdf", w = 8, h = 6)

cat("  方差解释: PC1=", var_pct[1], "% PC2=", var_pct[2],
    "% PC3=", var_pct[3], "%\n\n")

# --------------------------------------------------
# 1.4 样品相关性热图
# --------------------------------------------------
cat("[4/7] 样品相关性热图...\n")

cor_mat <- cor(log_tpm, method = "pearson")

anno_col <- data.frame(Passage = sample_table$passage,
                       row.names = sample_table$sample_id)
anno_colors <- list(Passage = passage_colors)

save_pheatmap("01_correlation_heatmap.pdf", w = 8, h = 7)
pheatmap(cor_mat,
         clustering_method = "ward.D2",
         color = cor_heatmap_colors,
         annotation_col = anno_col,
         annotation_colors = anno_colors,
         display_numbers = TRUE, number_format = "%.3f",
         fontsize_number = 8, fontsize = 10,
         main = "Sample Pearson Correlation (log2 TPM)")
dev.off()
cat("  ✓ 相关性热图已保存\n\n")

# --------------------------------------------------
# 1.5 经典衰老标志基因表达
# --------------------------------------------------
cat("[5/7] 经典衰老标志基因检查...\n")

marker_genes <- list(
  `Cell Cycle Arrest (expect UP)` = c("CDKN1A","CDKN2A","CDKN2B","TP53","RB1"),
  `SASP Factors (expect UP)`      = c("IL6","CXCL8","CCL2","CXCL1","MMP1",
                                       "MMP3","SERPINE1","TIMP1","IGFBP3","IGFBP7"),
  `Proliferation (expect DOWN)`   = c("MCM3","MYBL2","RRM2","CDC20","STMN1",
                                       "E2F1","MKI67","TOP2A"),
  `Other Senescence Markers`      = c("GPNMB","STAT1","PRNP","SOD2",
                                       "LMNB1","HMGA2","GATA6"),
  `Splicing Factors`              = c("SRSF3","HNRNPA1","YBX1","SRSF1",
                                       "HNRNPD","PTBP1","RBFOX2"),
  `RNA Editing Enzymes`           = c("ADAR","ADARB1","ADARB2")
)

# 提取所有标志基因表达
marker_data_list <- list()
for (cat_name in names(marker_genes)) {
  symbols <- marker_genes[[cat_name]]
  mapped  <- find_ensembl(symbols, gene_info)
  
  for (i in seq_len(nrow(mapped))) {
    if (is.na(mapped$ensembl[i])) {
      cat("  ⚠ 未找到:", mapped$symbol[i], "\n")
      next
    }
    vals <- tpm_mat[mapped$ensembl[i], ]
    marker_data_list[[length(marker_data_list) + 1]] <- data.frame(
      gene      = mapped$symbol[i],
      category  = cat_name,
      sample    = sample_table$sample_id,
      passage   = sample_table$passage,
      passage_num = sample_table$passage_num,
      tpm       = as.numeric(vals),
      stringsAsFactors = FALSE
    )
  }
}
marker_expr <- do.call(rbind, marker_data_list)
save_data(marker_expr, "01_marker_gene_expression.csv")

# 分类别画图
for (cat_name in names(marker_genes)) {
  sub <- marker_expr[marker_expr$category == cat_name, ]
  if (nrow(sub) == 0) next
  
  means <- sub %>%
    group_by(gene, passage, passage_num) %>%
    summarise(mean_tpm = mean(tpm), .groups = "drop")
  
  n_genes <- length(unique(sub$gene))
  gene_colors <- colorRampPalette(
    c("#4575B4","#D73027","#FF7F00","#984EA3",
      "#4DAF4A","#A65628","#F781BF","#66C2A5"))(n_genes)
  
  p <- ggplot(sub, aes(x = passage, y = tpm, color = gene)) +
    geom_jitter(width = 0.1, alpha = 0.5, size = 2) +
    geom_line(data = means, aes(x = passage, y = mean_tpm, group = gene),
              linewidth = 1) +
    geom_point(data = means, aes(x = passage, y = mean_tpm), size = 3) +
    facet_wrap(~gene, scales = "free_y", ncol = 4) +
    scale_color_manual(values = gene_colors) +
    labs(title = cat_name, x = "Passage", y = "TPM") +
    theme_bindlab(base_size = 12) +
    theme(legend.position = "none")
  
  safe_name <- gsub("[^A-Za-z0-9]", "_", cat_name)
  save_fig(p, paste0("01_markers_", safe_name, ".pdf"),
           w = 12, h = max(4, ceiling(n_genes / 4) * 3))
}
cat("\n")

# --------------------------------------------------
# 1.6 SenMayo衰老评分
# --------------------------------------------------
cat("[6/7] SenMayo衰老评分...\n")

senmayo_genes <- c(
  "ACVR1B","ANG","ANGPT1","ANGPTL4","AREG","AXL","BEX3",
  "BMP2","BMP6","C3","CCL1","CCL13","CCL16","CCL2","CCL20",
  "CCL24","CCL26","CCL3","CCL3L1","CCL4","CCL5","CCL7","CCL8",
  "CD55","CD9","CSF1","CSF2","CSF2RB","CST10","CTNNB1","CTSB",
  "CXCL1","CXCL10","CXCL12","CXCL13","CXCL14","CXCL16","CXCL2",
  "CXCL3","CXCL8","CXCR2","DKK1","EDN1","EGF","EGFR","EREG",
  "ESM1","ETS2","FAS","FGF1","FGF2","FGF7","GDF15","GEM","GMFG",
  "HGF","HMGB1","ICAM1","ICAM3","IGF1","IGFBP1","IGFBP2",
  "IGFBP3","IGFBP4","IGFBP5","IGFBP6","IGFBP7","IL10","IL13",
  "IL15","IL18","IL1A","IL1B","IL2","IL32","IL6","IL6ST","IL7",
  "INHA","IQGAP2","ITGA2","ITPKA","JUN","KITLG","LCP1","MIF",
  "MMP1","MMP10","MMP12","MMP13","MMP14","MMP2","MMP3","MMP9",
  "NAP1L4","NRG1","PAPPA","PECAM1","PGF","PIGF","PLAT","PLAU",
  "PLAUR","PTBP1","PTGER2","PTGES","RPS6KA5","SCAMP4","SELPLG",
  "SEMA3F","SERPINB4","SERPINE1","SERPINE2","SPP1","SPX","TIMP2",
  "TNF","TNFRSF10C","TNFRSF11B","TNFRSF1A","TNFRSF1B","TUBGCP2",
  "VEGFA","VEGFC","VGF","WNT16","WNT2"
)

senmayo_map   <- find_ensembl(senmayo_genes, gene_info)
senmayo_found <- senmayo_map$ensembl[!is.na(senmayo_map$ensembl)]
senmayo_in_mat <- intersect(senmayo_found, rownames(log_tpm))
cat("  SenMayo基因: 总", length(senmayo_genes),
    " | 匹配", length(senmayo_in_mat), "\n")

senmayo_list <- list(SenMayo = senmayo_in_mat)

# ssGSEA评分
tryCatch({
  param <- ssgseaParam(log_tpm, senmayo_list)
  scores <- gsva(param)
}, error = function(e) {
  scores <<- gsva(log_tpm, senmayo_list, method = "ssgsea", verbose = FALSE)
})

score_df <- data.frame(
  sample      = colnames(scores),
  passage     = sample_table$passage,
  passage_num = sample_table$passage_num,
  SenMayo     = as.numeric(scores["SenMayo", ]),
  stringsAsFactors = FALSE
)
save_data(score_df, "01_SenMayo_scores.csv")

# 统计
cor_test <- cor.test(score_df$passage_num, score_df$SenMayo, method = "spearman")
cat("  Spearman rho =", round(cor_test$estimate, 4),
    ", p =", format(cor_test$p.value, digits = 4), "\n")

score_summary <- score_df %>%
  group_by(passage) %>%
  summarise(mean = mean(SenMayo), sd = sd(SenMayo), .groups = "drop")
cat("  代次均值:\n")
print(as.data.frame(score_summary))

# 图: 箱线图
p <- ggplot(score_df, aes(x = passage, y = SenMayo, fill = passage)) +
  geom_boxplot(width = 0.5, alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 3, alpha = 0.8) +
  scale_fill_manual(values = passage_colors) +
  labs(title = paste0("SenMayo Score (rho=", round(cor_test$estimate, 3),
                      ", p=", format(cor_test$p.value, digits = 3), ")"),
       x = "Passage", y = "SenMayo ssGSEA Score") +
  theme_bindlab(base_size = 12) +
  theme(legend.position = "none")
save_fig(p, "01_SenMayo_boxplot.pdf", w = 6, h = 5)

# 图: 趋势线
p <- ggplot(score_df, aes(x = passage_num, y = SenMayo)) +
  geom_point(aes(color = passage), size = 4) +
  geom_smooth(method = "lm", se = TRUE, color = "grey30", linetype = "dashed") +
  scale_color_manual(values = passage_colors) +
  labs(title = "SenMayo Score Trajectory",
       x = "Passage Number", y = "ssGSEA Score") +
  theme_bindlab(base_size = 12)
save_fig(p, "01_SenMayo_trajectory.pdf", w = 7, h = 5)

cat("\n")

# --------------------------------------------------
# 1.7 表达分布总览
# --------------------------------------------------
cat("[7/7] 表达分布总览...\n")

tpm_long <- reshape2::melt(log2(tpm_mat + 1))
colnames(tpm_long) <- c("gene", "sample", "log2TPM")
tpm_long$passage <- sample_table$passage[match(tpm_long$sample,
                                                sample_table$sample_id)]

# 图: 箱线图
p <- ggplot(tpm_long[tpm_long$log2TPM > 0, ],
            aes(x = sample, y = log2TPM, fill = passage)) +
  geom_boxplot(width = 0.6, outlier.size = 0.3) +
  scale_fill_manual(values = passage_colors) +
  labs(title = "Expression Distribution", x = "", y = "log2(TPM+1)") +
  theme_bindlab(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_fig(p, "01_expression_boxplot.pdf", w = 8, h = 5)

# 图: 密度图
p <- ggplot(tpm_long[tpm_long$log2TPM > 0, ],
            aes(x = log2TPM, color = sample)) +
  geom_density(linewidth = 0.7) +
  labs(title = "Expression Density (expressed genes)",
       x = "log2(TPM+1)", y = "Density") +
  theme_bindlab(base_size = 12)
save_fig(p, "01_expression_density.pdf", w = 10, h = 5)

# --------------------------------------------------
# 完成
# --------------------------------------------------
cat("\n")
cat("====================================================\n")
cat("  Step 1 完成!\n")
cat("====================================================\n")
cat("  输出图表:", dir_fig, "/01_*.pdf\n")
cat("  输出数据:", dir_data, "/01_*.csv\n")
cat("\n")
cat("  ★ 判读指南:\n")
cat("  1. PCA: 同代次聚在一起? PC1沿P2→P12排列?\n")
cat("  2. 相关性: 同代次内 > 不同代次间?\n")
cat("  3. 标志基因: CDKN1A/SERPINE1↑, MKI67/TOP2A↓?\n")
cat("  4. SenMayo: rho>0.7 且 p<0.05 = 衰老轨迹成立\n")
cat("\n")
cat("  请把figures目录下的01_开头PDF发给我判读!\n")
cat("====================================================\n")
