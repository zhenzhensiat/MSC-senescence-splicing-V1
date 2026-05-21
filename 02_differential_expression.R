# 02_differential_expression.R — DESeq2 differential expression and enrichment analysis

source("00_project_config.R")

load_packages(c("DESeq2", "ggplot2", "pheatmap", "RColorBrewer",
                "ggrepel", "dplyr", "reshape2", "fgsea"))

cat("========== Step 2: Differential Expression ==========\n\n")

# --------------------------------------------------
# 2.1 构建DESeq2对象
# --------------------------------------------------
cat("[1/6] 构建DESeq2对象...\n")

raw <- read.delim(f_expression, sep = "\t", header = TRUE,
                  stringsAsFactors = FALSE, check.names = FALSE)

count_cols <- grep("_count$", colnames(raw), value = TRUE)
count_mat  <- as.matrix(raw[, count_cols])
rownames(count_mat) <- raw$id
colnames(count_mat) <- gsub("_count$", "", count_cols)

# 取整（DESeq2要求整数）
count_mat <- round(count_mat)
count_mat <- count_mat[, sample_table$sample_id]

gene_info <- raw[, c("id", "Symbol")]

# 构建colData
col_data <- data.frame(
  row.names   = sample_table$sample_id,
  passage     = sample_table$passage,
  passage_num = sample_table$passage_num
)

# 过滤低表达基因：至少3个样品中count >= 10
keep <- rowSums(count_mat >= 10) >= 3
count_filtered <- count_mat[keep, ]
cat("  过滤前:", nrow(count_mat), "基因\n")
cat("  过滤后:", nrow(count_filtered), "基因\n")

# 构建DESeq2 Dataset
dds <- DESeqDataSetFromMatrix(
  countData = count_filtered,
  colData   = col_data,
  design    = ~ passage
)
dds$passage <- relevel(dds$passage, ref = "P2")

cat("✓ DESeq2对象构建完成\n\n")

# --------------------------------------------------
# 2.2 LRT时间序列检验（全局检验）
# --------------------------------------------------
cat("[2/6] LRT时间序列检验...\n")

dds_lrt <- DESeq(dds, test = "LRT", reduced = ~ 1)
res_lrt <- results(dds_lrt, alpha = 0.05)

lrt_df <- as.data.frame(res_lrt)
lrt_df$gene_id <- rownames(lrt_df)
lrt_df$symbol  <- gene_info$Symbol[match(lrt_df$gene_id, gene_info$id)]
lrt_df <- lrt_df[order(lrt_df$padj), ]

n_sig_lrt <- sum(lrt_df$padj < 0.05, na.rm = TRUE)
cat("  LRT显著基因（padj < 0.05）:", n_sig_lrt, "\n")

save_data(lrt_df, "02_LRT_all_genes.csv")

cat("✓ LRT检验完成\n\n")

# --------------------------------------------------
# 2.3 配对比较（Wald检验）
# --------------------------------------------------
cat("[3/6] 配对比较（P8/P10/P12 vs P2）...\n")

dds_wald <- DESeq(dds, test = "Wald")

contrasts <- list(
  "P8_vs_P2"  = c("passage", "P8",  "P2"),
  "P10_vs_P2" = c("passage", "P10", "P2"),
  "P12_vs_P2" = c("passage", "P12", "P2")
)

de_results <- list()
de_summary <- data.frame()

for (name in names(contrasts)) {
  res <- results(dds_wald, contrast = contrasts[[name]], alpha = 0.05)
  df  <- as.data.frame(res)
  df$gene_id <- rownames(df)
  df$symbol  <- gene_info$Symbol[match(df$gene_id, gene_info$id)]
  df$comparison <- name
  
  # 分类
  df$direction <- "NS"
  df$direction[df$padj < 0.05 & df$log2FoldChange > 1]  <- "Up"
  df$direction[df$padj < 0.05 & df$log2FoldChange < -1] <- "Down"
  
  n_up   <- sum(df$direction == "Up", na.rm = TRUE)
  n_down <- sum(df$direction == "Down", na.rm = TRUE)
  
  de_results[[name]] <- df
  de_summary <- rbind(de_summary, data.frame(
    comparison = name, up = n_up, down = n_down,
    total = n_up + n_down, stringsAsFactors = FALSE
  ))
  
  save_data(df, paste0("02_DE_", name, ".csv"))
  cat("  ", name, ": Up=", n_up, " Down=", n_down, "\n")
}

save_data(de_summary, "02_DE_summary.csv")
cat("✓ 配对比较完成\n\n")

# --------------------------------------------------
# 2.4 火山图
# --------------------------------------------------
cat("[4/6] 绑制火山图...\n")

for (name in names(de_results)) {
  df <- de_results[[name]]
  df <- df[!is.na(df$padj), ]
  
  # 标注top基因
  top_genes <- df %>%
    filter(direction != "NS") %>%
    arrange(padj) %>%
    head(20)
  
  p <- ggplot(df, aes(x = log2FoldChange, y = -log10(padj), color = direction)) +
    geom_point(size = 0.8, alpha = 0.6) +
    geom_point(data = df[df$direction != "NS", ], size = 1.2, alpha = 0.8) +
    scale_color_manual(values = de_colors) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey50", linewidth = 0.3) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50", linewidth = 0.3) +
    geom_text_repel(data = top_genes, aes(label = symbol),
                    size = 3, max.overlaps = 15, color = "black",
                    segment.color = "grey60", segment.size = 0.3) +
    labs(title = gsub("_", " ", name),
         subtitle = paste0("Up: ", sum(df$direction == "Up"),
                           "  Down: ", sum(df$direction == "Down")),
         x = "log2 Fold Change", y = "-log10(adjusted p-value)") +
    theme_bindlab(base_size = 12) +
    theme(legend.position = "none")
  
  save_fig(p, paste0("02_volcano_", name, ".pdf"), w = 7, h = 6)
}

cat("✓ 火山图完成\n\n")

# --------------------------------------------------
# 2.5 热图: Top DEGs
# --------------------------------------------------
cat("[5/6] DEG热图...\n")

# 获取VST标准化数据
vsd <- vst(dds_wald, blind = FALSE)
vsd_mat <- assay(vsd)

# 从LRT中取top100显著基因
top100_lrt <- lrt_df %>%
  filter(padj < 0.05) %>%
  head(100)

top_genes_ids <- top100_lrt$gene_id
top_mat <- vsd_mat[top_genes_ids, ]

# Z-score标准化（按行）
top_mat_z <- t(scale(t(top_mat)))

# 行名用gene symbol
row_labels <- gene_info$Symbol[match(rownames(top_mat_z), gene_info$id)]
row_labels[row_labels == "-" | is.na(row_labels)] <- rownames(top_mat_z)[row_labels == "-" | is.na(row_labels)]
rownames(top_mat_z) <- make.unique(row_labels)

# 注释
anno_col <- data.frame(Passage = col_data$passage, row.names = rownames(col_data))
anno_colors <- list(Passage = passage_colors)

save_pheatmap("02_heatmap_top100_LRT.pdf", w = 10, h = 14)
pheatmap(top_mat_z,
         color = heatmap_colors,
         clustering_method = "ward.D2",
         cluster_cols = FALSE,
         annotation_col = anno_col,
         annotation_colors = anno_colors,
         show_rownames = TRUE,
         fontsize_row = 6,
         fontsize = 10,
         breaks = seq(-2.5, 2.5, length.out = 101),
         main = "Top 100 LRT DEGs (Z-score)")
dev.off()
cat("  ✓ saved: 02_heatmap_top100_LRT.pdf\n")

# 也存一个PNG版热图
png(file.path(dir_fig, "02_heatmap_top100_LRT.png"), width = 10, height = 14,
    units = "in", res = 300)
pheatmap(top_mat_z,
         color = heatmap_colors,
         clustering_method = "ward.D2",
         cluster_cols = FALSE,
         annotation_col = anno_col,
         annotation_colors = anno_colors,
         show_rownames = TRUE,
         fontsize_row = 6,
         fontsize = 10,
         breaks = seq(-2.5, 2.5, length.out = 101),
         main = "Top 100 LRT DEGs (Z-score)")
dev.off()

# DEG overlap统计（UpSet风格用Venn数据）
deg_lists <- list()
for (name in names(de_results)) {
  df <- de_results[[name]]
  deg_lists[[name]] <- df$gene_id[df$direction != "NS"]
}

# 交集统计
all_degs <- unique(unlist(deg_lists))
overlap_mat <- sapply(deg_lists, function(x) all_degs %in% x)
rownames(overlap_mat) <- all_degs

overlap_summary <- data.frame(
  gene_id = all_degs,
  symbol = gene_info$Symbol[match(all_degs, gene_info$id)],
  P8_vs_P2  = overlap_mat[, "P8_vs_P2"],
  P10_vs_P2 = overlap_mat[, "P10_vs_P2"],
  P12_vs_P2 = overlap_mat[, "P12_vs_P2"],
  n_comparisons = rowSums(overlap_mat),
  stringsAsFactors = FALSE
)
overlap_summary <- overlap_summary[order(-overlap_summary$n_comparisons), ]
save_data(overlap_summary, "02_DEG_overlap.csv")

# 交集统计图
n_only_p8  <- sum(overlap_summary$n_comparisons == 1 & overlap_summary$P8_vs_P2)
n_only_p10 <- sum(overlap_summary$n_comparisons == 1 & overlap_summary$P10_vs_P2)
n_only_p12 <- sum(overlap_summary$n_comparisons == 1 & overlap_summary$P12_vs_P2)
n_all3     <- sum(overlap_summary$n_comparisons == 3)
n_p8_p10   <- sum(overlap_summary$P8_vs_P2 & overlap_summary$P10_vs_P2 & !overlap_summary$P12_vs_P2)
n_p8_p12   <- sum(overlap_summary$P8_vs_P2 & !overlap_summary$P10_vs_P2 & overlap_summary$P12_vs_P2)
n_p10_p12  <- sum(!overlap_summary$P8_vs_P2 & overlap_summary$P10_vs_P2 & overlap_summary$P12_vs_P2)

overlap_bar <- data.frame(
  category = c("P8 only", "P10 only", "P12 only",
               "P8∩P10", "P8∩P12", "P10∩P12", "All three"),
  count = c(n_only_p8, n_only_p10, n_only_p12,
            n_p8_p10, n_p8_p12, n_p10_p12, n_all3),
  stringsAsFactors = FALSE
)
overlap_bar$category <- factor(overlap_bar$category,
                               levels = rev(overlap_bar$category))

p <- ggplot(overlap_bar, aes(x = category, y = count)) +
  geom_bar(stat = "identity", fill = "#4575B4", width = 0.6) +
  geom_text(aes(label = count), hjust = -0.2, size = 3.5) +
  coord_flip() +
  labs(title = "DEG overlap across comparisons",
       subtitle = paste0("Total unique DEGs: ", nrow(overlap_summary)),
       x = "", y = "Number of DEGs") +
  theme_bindlab(base_size = 12)

save_fig(p, "02_DEG_overlap_bar.pdf", w = 7, h = 5)

cat("✓ 热图和overlap分析完成\n\n")

# --------------------------------------------------
# 2.6 fgsea通路富集分析
# --------------------------------------------------
cat("[6/6] fgsea通路富集...\n")

# 从本地缓存加载Hallmark基因集（需先运行 00_download_databases.R）
hallmark_rds <- file.path(dir_data, "msigdb_hallmark.rds")
if (!file.exists(hallmark_rds)) {
  cat("  ✗ 未找到本地基因集缓存!\n")
  cat("  请先运行: source('scripts/00_download_databases.R')\n")
  stop("缺少 msigdb_hallmark.rds，请先运行数据库下载脚本")
}

hallmark_db  <- readRDS(hallmark_rds)
hallmark_list <- hallmark_db$ensembl
# 如果ensembl版不可用（手动GMT方案），用symbol版
use_symbol <- is.null(hallmark_list)
if (use_symbol) {
  hallmark_list <- hallmark_db$symbol
  cat("  ℹ 使用Gene Symbol版基因集\n")
} else {
  cat("  ℹ 使用Ensembl ID版基因集\n")
}
cat("  基因集数量:", length(hallmark_list), "\n")

# 对每个配对比较做GSEA
gsea_all <- list()

for (name in names(de_results)) {
  df <- de_results[[name]]
  df <- df[!is.na(df$stat), ]
  
  # 构建排序统计量
  if (use_symbol) {
    # Symbol模式：用gene symbol作为key
    df_ranked <- df[!is.na(df$symbol) & df$symbol != "-", ]
    ranks <- setNames(df_ranked$stat, df_ranked$symbol)
    # 去重（同symbol取绝对值最大的）
    ranks <- tapply(ranks, names(ranks), function(x) x[which.max(abs(x))])
  } else {
    ranks <- setNames(df$stat, df$gene_id)
  }
  ranks <- sort(ranks, decreasing = TRUE)
  
  # 运行fgsea
  set.seed(42)
  fgsea_res <- fgsea(pathways = hallmark_list,
                     stats = ranks,
                     minSize = 15,
                     maxSize = 500)
  
  fgsea_df <- as.data.frame(fgsea_res)
  fgsea_df$comparison <- name
  fgsea_df$pathway_short <- gsub("^HALLMARK_", "", fgsea_df$pathway)
  
  # 去掉leadingEdge列（list列不方便存CSV）
  fgsea_df$leadingEdge <- sapply(fgsea_df$leadingEdge, function(x) paste(x, collapse = ";"))
  
  gsea_all[[name]] <- fgsea_df
  
  n_up_pw   <- sum(fgsea_df$padj < 0.05 & fgsea_df$NES > 0, na.rm = TRUE)
  n_down_pw <- sum(fgsea_df$padj < 0.05 & fgsea_df$NES < 0, na.rm = TRUE)
  cat("  ", name, ": 上调通路=", n_up_pw, " 下调通路=", n_down_pw, "\n")
}

gsea_combined <- do.call(rbind, gsea_all)
save_data(gsea_combined, "02_fgsea_hallmark_all.csv")

# GSEA汇总热图：NES值
sig_pathways <- gsea_combined %>%
  filter(padj < 0.05) %>%
  pull(pathway) %>%
  unique()

if (length(sig_pathways) > 0) {
  nes_matrix <- gsea_combined %>%
    filter(pathway %in% sig_pathways) %>%
    mutate(pathway_short = gsub("^HALLMARK_", "", pathway)) %>%
    dplyr::select(pathway_short, comparison, NES) %>%
    tidyr::pivot_wider(names_from = comparison, values_from = NES) %>%
    as.data.frame()
  
  rownames(nes_matrix) <- nes_matrix$pathway_short
  nes_matrix$pathway_short <- NULL
  nes_mat <- as.matrix(nes_matrix)
  
  # padj标记矩阵
  padj_matrix <- gsea_combined %>%
    filter(pathway %in% sig_pathways) %>%
    mutate(pathway_short = gsub("^HALLMARK_", "", pathway)) %>%
    dplyr::select(pathway_short, comparison, padj) %>%
    tidyr::pivot_wider(names_from = comparison, values_from = padj) %>%
    as.data.frame()
  rownames(padj_matrix) <- padj_matrix$pathway_short
  padj_matrix$pathway_short <- NULL
  padj_mat <- as.matrix(padj_matrix)
  
  # 星号标注
  star_mat <- matrix("", nrow = nrow(padj_mat), ncol = ncol(padj_mat))
  star_mat[padj_mat < 0.05] <- "*"
  star_mat[padj_mat < 0.01] <- "**"
  star_mat[padj_mat < 0.001] <- "***"
  
  # 显示数字 = NES + 星号
  display_mat <- matrix(paste0(round(nes_mat, 2), star_mat),
                        nrow = nrow(nes_mat), ncol = ncol(nes_mat))
  
  h_fig <- max(8, nrow(nes_mat) * 0.35 + 3)
  
  save_pheatmap("02_GSEA_hallmark_heatmap.pdf", w = 9, h = h_fig)
  pheatmap(nes_mat,
           color = heatmap_colors,
           breaks = seq(-3, 3, length.out = 101),
           cluster_cols = FALSE,
           clustering_method = "ward.D2",
           display_numbers = display_mat,
           number_color = "black",
           fontsize_number = 8,
           fontsize = 10,
           fontsize_row = 9,
           main = "GSEA Hallmark Pathways (NES)")
  dev.off()
  cat("  ✓ saved: 02_GSEA_hallmark_heatmap.pdf\n")
  
  # PNG版
  png(file.path(dir_fig, "02_GSEA_hallmark_heatmap.png"),
      width = 9, height = h_fig, units = "in", res = 300)
  pheatmap(nes_mat,
           color = heatmap_colors,
           breaks = seq(-3, 3, length.out = 101),
           cluster_cols = FALSE,
           clustering_method = "ward.D2",
           display_numbers = display_mat,
           number_color = "black",
           fontsize_number = 8,
           fontsize = 10,
           fontsize_row = 9,
           main = "GSEA Hallmark Pathways (NES)")
  dev.off()
}

# 每个比较的top通路条形图
for (name in names(gsea_all)) {
  df <- gsea_all[[name]] %>%
    filter(padj < 0.05) %>%
    arrange(NES) %>%
    mutate(pathway_short = factor(pathway_short, levels = pathway_short))
  
  if (nrow(df) == 0) next
  
  # 取top10上调 + top10下调
  top_up   <- df %>% filter(NES > 0) %>% tail(10)
  top_down <- df %>% filter(NES < 0) %>% head(10)
  plot_df  <- rbind(top_down, top_up)
  plot_df$pathway_short <- factor(plot_df$pathway_short,
                                  levels = plot_df$pathway_short)
  
  p <- ggplot(plot_df, aes(x = pathway_short, y = NES,
                            fill = ifelse(NES > 0, "Up", "Down"))) +
    geom_bar(stat = "identity", width = 0.7) +
    coord_flip() +
    scale_fill_manual(values = de_colors[c("Up", "Down")], guide = "none") +
    labs(title = paste0("Hallmark GSEA: ", gsub("_", " ", name)),
         x = "", y = "Normalized Enrichment Score (NES)") +
    theme_bindlab(base_size = 11)
  
  save_fig(p, paste0("02_GSEA_barplot_", name, ".pdf"), w = 9, h = 7)
}

cat("✓ fgsea富集分析完成\n\n")

# --------------------------------------------------
# 完成
# --------------------------------------------------

# 保存VST矩阵供后续步骤使用
vsd_df <- as.data.frame(vsd_mat)
vsd_df$gene_id <- rownames(vsd_df)
vsd_df$symbol  <- gene_info$Symbol[match(vsd_df$gene_id, gene_info$id)]
save_data(vsd_df, "02_VST_expression.csv")

cat("====================================================\n")
cat("  Step 2 完成!\n")
cat("====================================================\n")
cat("\n  输出文件:\n")
cat("  data/02_LRT_all_genes.csv        — LRT全局检验结果\n")
cat("  data/02_DE_P8_vs_P2.csv          — P8 vs P2差异基因\n")
cat("  data/02_DE_P10_vs_P2.csv         — P10 vs P2差异基因\n")
cat("  data/02_DE_P12_vs_P2.csv         — P12 vs P2差异基因\n")
cat("  data/02_DE_summary.csv           — DEG数量汇总\n")
cat("  data/02_DEG_overlap.csv          — DEG交集\n")
cat("  data/02_fgsea_hallmark_all.csv   — Hallmark GSEA结果\n")
cat("  data/02_VST_expression.csv       — VST标准化表达（供后续用）\n")
cat("  figures/02_volcano_*.pdf         — 火山图\n")
cat("  figures/02_heatmap_top100_LRT.pdf— Top DEG热图\n")
cat("  figures/02_DEG_overlap_bar.pdf   — DEG交集柱形图\n")
cat("  figures/02_GSEA_*.pdf            — GSEA结果图\n")
cat("\n  请把 data/02_DE_summary.csv 和\n")
cat("  data/02_fgsea_hallmark_all.csv 上传给我判读!\n")
cat("====================================================\n")
