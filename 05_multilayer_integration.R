# 05_multilayer_integration.R — Multi-layer integration: DEG × DSE × RNA editing cross-network

source("00_project_config.R")

load_packages(c("ggplot2", "dplyr", "tidyr", "pheatmap", "RColorBrewer",
  "ggrepel", "reshape2", "UpSetR",
  "corrplot", "scales"))

cat("========== Step 5: Multi-layer Integration (v3.1) ==========\n\n")

# ============================================================
# 0. 加载所有已有数据
# ============================================================
cat("[0] 加载Step2/3/4结果数据...\n")

deg_overlap  <- read.csv(file.path(dir_data, "02_DEG_overlap.csv"), stringsAsFactors = FALSE)
de_p8  <- read.csv(file.path(dir_data, "02_DE_P8_vs_P2.csv"), stringsAsFactors = FALSE)
de_p10  <- read.csv(file.path(dir_data, "02_DE_P10_vs_P2.csv"), stringsAsFactors = FALSE)
de_p12  <- read.csv(file.path(dir_data, "02_DE_P12_vs_P2.csv"), stringsAsFactors = FALSE)
de_summary  <- read.csv(file.path(dir_data, "02_DE_summary.csv"), stringsAsFactors = FALSE)
lrt_all  <- read.csv(file.path(dir_data, "02_LRT_all_genes.csv"), stringsAsFactors = FALSE)
vst_expr  <- read.csv(file.path(dir_data, "02_VST_expression.csv"), stringsAsFactors = FALSE,
  check.names = FALSE)  # ★ v3.1修复: 保留原始列名P2-1而非P2.1
gsea_all  <- read.csv(file.path(dir_data, "02_fgsea_hallmark_all.csv"), stringsAsFactors = FALSE)
dse_v2  <- read.csv(file.path(dir_data, "03_DSE_events_v2.csv"), stringsAsFactors = FALSE)
dse_core  <- read.csv(file.path(dir_data, "03_DSE_core_genes.csv"), stringsAsFactors = FALSE)
sf_expr  <- read.csv(file.path(dir_data, "03_splicing_factor_expression.csv"), stringsAsFactors = FALSE)
edit_all  <- read.csv(file.path(dir_data, "04_DE_editing_all.csv"), stringsAsFactors = FALSE)
edit_summary <- read.csv(file.path(dir_data, "04_editing_summary.csv"), stringsAsFactors = FALSE)
edit_gained_lost <- read.csv(file.path(dir_data, "04_editing_gained_lost.csv"), stringsAsFactors = FALSE)
edit_type  <- read.csv(file.path(dir_data, "04_editing_type_long.csv"), stringsAsFactors = FALSE)
edit_func  <- read.csv(file.path(dir_data, "04_editing_function_dist.csv"), stringsAsFactors = FALSE)
edit_struct  <- read.csv(file.path(dir_data, "04_editing_structure_dist.csv"), stringsAsFactors = FALSE)
senmayo_scores <- read.csv(file.path(dir_data, "01_SenMayo_scores.csv"), stringsAsFactors = FALSE)

# 基因注释
gene_info_full <- read.delim(f_expression, sep = "\t", header = TRUE,
  stringsAsFactors = FALSE, check.names = FALSE)
gene_info_map <- gene_info_full[, c("id", "Symbol")]

# 样品列名
# ★ read.csv(check.names=TRUE) 会把 "P2-1" 转成 "P2.1"
# 需要同时匹配连字符和点号格式
sample_cols <- grep("^P\\d+-\\d+$", colnames(vst_expr), value = TRUE)
if (length(sample_cols) == 0) {
  sample_cols <- grep("^P\\d+\\.\\d+$", colnames(vst_expr), value = TRUE)
  if (length(sample_cols) > 0)
  cat("  ★ 注意: VST列名含点号(read.csv转换), 已自动匹配\n")
}
if (length(sample_cols) == 0)
  stop("无法在VST矩阵中匹配样品列! 检查02_VST_expression.csv的列名格式")

# 背景基因数（DESeq2过滤后的expressed genes）
N_BACKGROUND <- nrow(lrt_all)

cat("  Background genes (expressed): ", N_BACKGROUND, "\n")
cat("  DEGs: ", nrow(deg_overlap), "\n")
cat("  DSE events (v2): ", nrow(dse_v2), "\n")
cat("  Editing records: ", nrow(edit_all), "\n")
cat("  ✓ 数据加载完成\n\n")

# ============================================================
# 通用工具: Fisher's exact test for gene set overlap
# ★ v3.1修正: alternative="two.sided" 替代 "greater"
#   依据: ?fisher.test — "greater"对反富集(observed<expected)返回p≈1.0(假阴性)
#   "two.sided"对富集和反富集方向均给出正确p值
#   文献: ENCODE重叠分析惯例 (PMID:22955616); R stats::fisher.test 文档
# ============================================================
#' 对两个基因集做Fisher's exact test
#' @param set_a, set_b 两个基因ID向量
#' @param background 背景基因总数
#' @return data.frame with OR, p-value, counts
fisher_overlap_test <- function(set_a, set_b, background,
  label_a = "A", label_b = "B") {
  a_only  <- length(setdiff(set_a, set_b))
  b_only  <- length(setdiff(set_b, set_a))
  both  <- length(intersect(set_a, set_b))
  neither <- background - a_only - b_only - both
  # 防止负值（基因不全在background中）
  neither <- max(0, neither)

  mat <- matrix(c(both, b_only, a_only, neither), nrow = 2,
  dimnames = list(c(paste0(label_a, "+"), paste0(label_a, "-")),
  c(paste0(label_b, "+"), paste0(label_b, "-"))))

  ft <- fisher.test(mat, alternative = "two.sided")  # ★ v3.1修正: greater → two.sided
  data.frame(
  set_A = label_a, set_B = label_b,
  n_A = length(set_a), n_B = length(set_b),
  overlap = both,
  expected = round(length(set_a) * length(set_b) / background, 1),
  fold_enrichment = round(both / max(length(set_a) * length(set_b) / background, 0.1), 2),
  odds_ratio = round(as.numeric(ft$estimate), 2),
  p_value = ft$p.value,
  stringsAsFactors = FALSE
  )
}

# ============================================================
# PART A: Step3 DSE 深入分析
# ============================================================
cat("===== PART A: 可变剪接深入分析 =====\n\n")

# --------------------------------------------------
# A1. DSE代次分布统计
# --------------------------------------------------
cat("[A1] DSE分布统计...\n")

dse_stats <- dse_v2 %>%
  group_by(as_type) %>%
  summarise(
  total_events = n(),
  modT_P8  = sum(dse_modT_P8 == TRUE, na.rm = TRUE),
  modT_P10 = sum(dse_modT_P10 == TRUE, na.rm = TRUE),
  modT_P12 = sum(dse_modT_P12 == TRUE, na.rm = TRUE),
  modT_any = sum(dse_modT_P8 | dse_modT_P10 | dse_modT_P12, na.rm = TRUE),
  .groups = "drop"
  )

total_row <- dse_stats %>% summarise(across(where(is.numeric), sum)) %>%
  mutate(as_type = "TOTAL")
dse_stats <- rbind(dse_stats, total_row)
save_data(dse_stats, "05_DSE_summary_corrected.csv")

cat("  DSE统计 (modT: FDR<0.05 & |dPSI|>=0.05):\n")
print(as.data.frame(dse_stats[, c("as_type","total_events","modT_P8","modT_P10","modT_P12","modT_any")]))

# 图: DSE类型×代次分布
plot_dse_types <- dse_stats %>%
  filter(as_type != "TOTAL") %>%
  dplyr::select(as_type, modT_P8, modT_P10, modT_P12) %>%
  reshape2::melt(id.vars = "as_type", variable.name = "comparison", value.name = "count") %>%
  mutate(comparison = factor(gsub("modT_", "", comparison), levels = c("P8","P10","P12")))

p <- ggplot(plot_dse_types, aes(x = as_type, y = count, fill = comparison)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  scale_fill_manual(values = passage_colors[c("P8","P10","P12")], name = "vs P2") +
  labs(title = "Differential Splicing Events by Type and Passage",
  subtitle = "And criterion: FDR < 0.05 & |dPSI| >= 0.05",
  x = "AS Type", y = "Number of DSE Events") +
  theme_bindlab(base_size = 12)
save_fig(p, "05_DSE_type_distribution.pdf", w = 8, h = 5)
cat("\n")

# --------------------------------------------------
# A2. DSE方向性: inclusion vs exclusion
# --------------------------------------------------
cat("[A2] DSE方向性...\n")

dse_direction <- dse_v2 %>%
  mutate(
  dir_P8  = case_when(dse_modT_P8  & dpsi_P8  > 0 ~ "Inclusion",
  dse_modT_P8  & dpsi_P8  < 0 ~ "Exclusion", TRUE ~ NA_character_),
  dir_P10 = case_when(dse_modT_P10 & dpsi_P10 > 0 ~ "Inclusion",
  dse_modT_P10 & dpsi_P10 < 0 ~ "Exclusion", TRUE ~ NA_character_),
  dir_P12 = case_when(dse_modT_P12 & dpsi_P12 > 0 ~ "Inclusion",
  dse_modT_P12 & dpsi_P12 < 0 ~ "Exclusion", TRUE ~ NA_character_)
  )

dir_stats <- data.frame(
  comparison = rep(c("P8_vs_P2","P10_vs_P2","P12_vs_P2"), each = 2),
  direction  = rep(c("Inclusion","Exclusion"), 3),
  count = c(sum(dse_direction$dir_P8  == "Inclusion", na.rm = TRUE),
  sum(dse_direction$dir_P8  == "Exclusion", na.rm = TRUE),
  sum(dse_direction$dir_P10 == "Inclusion", na.rm = TRUE),
  sum(dse_direction$dir_P10 == "Exclusion", na.rm = TRUE),
  sum(dse_direction$dir_P12 == "Inclusion", na.rm = TRUE),
  sum(dse_direction$dir_P12 == "Exclusion", na.rm = TRUE)),
  stringsAsFactors = FALSE
)
save_data(dir_stats, "05_DSE_direction_stats.csv")

dir_stats$passage <- factor(gsub("_vs_P2", "", dir_stats$comparison), levels = c("P8","P10","P12"))
p <- ggplot(dir_stats, aes(x = passage, y = count, fill = direction)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.6) +
  scale_fill_manual(values = c("Inclusion" = "#D73027", "Exclusion" = "#4575B4")) +
  labs(title = "DSE Direction: Exon Inclusion vs Exclusion",
  x = "Passage (vs P2)", y = "Number of DSE Events") +
  theme_bindlab(base_size = 12)
save_fig(p, "05_DSE_direction.pdf", w = 7, h = 5)
cat("\n")

# --------------------------------------------------
# A3. dPSI分布密度图
# --------------------------------------------------
cat("[A3] dPSI分布...\n")

dpsi_long <- dse_v2 %>%
  dplyr::select(AS_ID, as_type, dpsi_P8, dpsi_P10, dpsi_P12) %>%
  reshape2::melt(id.vars = c("AS_ID","as_type"), variable.name = "comparison", value.name = "dpsi") %>%
  mutate(comparison = factor(gsub("dpsi_", "", comparison), levels = c("P8","P10","P12"))) %>%
  filter(is.finite(dpsi))

p <- ggplot(dpsi_long, aes(x = dpsi, color = comparison)) +
  geom_density(linewidth = 0.8) +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", color = "grey50") +
  scale_color_manual(values = passage_colors[c("P8","P10","P12")], name = "vs P2") +
  xlim(-0.5, 0.5) +
  labs(title = "Distribution of dPSI Across Passages", x = "dPSI (vs P2)", y = "Density") +
  theme_bindlab(base_size = 12)
save_fig(p, "05_DSE_dpsi_density.pdf", w = 7, h = 5)

p <- p + facet_wrap(~as_type, scales = "free_y", ncol = 3) +
  labs(title = "dPSI Distribution by AS Type")
save_fig(p, "05_DSE_dpsi_by_type.pdf", w = 12, h = 6)
cat("\n")

# --------------------------------------------------
# A4. 剪接因子表达变化热图
# --------------------------------------------------
cat("[A4] 剪接因子表达热图...\n")

sf_mat <- as.matrix(sf_expr[, c("lfc_P8","lfc_P10","lfc_P12")])
rownames(sf_mat) <- sf_expr$gene
colnames(sf_mat) <- c("P8 vs P2", "P10 vs P2", "P12 vs P2")

sf_changed_mask <- apply(abs(sf_mat), 1, max) > 0.2
sf_mat_plot <- sf_mat[sf_changed_mask, , drop = FALSE]
sf_mat_plot[!is.finite(sf_mat_plot)] <- 0

if (nrow(sf_mat_plot) > 3) {
  save_pheatmap("05_splicing_factor_heatmap.pdf", w = 6, h = max(6, nrow(sf_mat_plot) * 0.3 + 2))
  pheatmap(sf_mat_plot,
  color = colorRampPalette(c("#4575B4","white","#D73027"))(100),
  breaks = seq(-1.5, 1.5, length.out = 101),
  cluster_cols = FALSE, clustering_method = "ward.D2",
  display_numbers = TRUE, number_format = "%.2f",
  fontsize_number = 8, fontsize = 10,
  main = "Splicing Factor Expression Changes (log2FC)")
  dev.off()
  cat("  ✓ saved: 05_splicing_factor_heatmap.pdf (", nrow(sf_mat_plot), " SFs)\n")
}
cat("\n")

# ============================================================
# PART B: Step4 Editing 深入分析
# ============================================================
cat("===== PART B: RNA Editing 深入分析 =====\n\n")

# B1. 编辑类型代次热图
cat("[B1] 编辑类型代次热图...\n")

edit_type_avg <- edit_type %>%
  group_by(passage, editing_type) %>%
  summarise(mean_prop = mean(proportion), .groups = "drop")

edit_type_mat <- edit_type_avg %>%
  tidyr::pivot_wider(names_from = passage, values_from = mean_prop) %>%
  as.data.frame()
rownames(edit_type_mat) <- edit_type_mat$editing_type
edit_type_mat$editing_type <- NULL
edit_type_mat <- as.matrix(edit_type_mat[, c("P2","P8","P10","P12")])

save_pheatmap("05_editing_type_heatmap.pdf", w = 6, h = 7)
pheatmap(edit_type_mat, color = colorRampPalette(c("white","#FDD49E","#D7301F"))(100),
  cluster_cols = FALSE, clustering_method = "ward.D2",
  display_numbers = TRUE, number_format = "%.3f",
  fontsize_number = 9, fontsize = 11,
  main = "RNA Editing Type Proportions by Passage")
dev.off()
cat("  ✓ saved\n")

# A->G趋势
ag_trend <- edit_type %>%
  filter(editing_type == "A->G") %>%
  group_by(passage, passage_num) %>%
  summarise(mean_count = mean(count), sd_count = sd(count), .groups = "drop") %>%
  mutate(passage = factor(passage, levels = c("P2","P8","P10","P12")))

p <- ggplot(ag_trend, aes(x = passage, y = mean_count)) +
  geom_bar(stat = "identity", fill = "#D73027", width = 0.6, alpha = 0.8) +
  geom_errorbar(aes(ymin = mean_count - sd_count, ymax = mean_count + sd_count), width = 0.2) +
  labs(title = "A-to-G (ADAR-mediated) Editing Sites by Passage",
  x = "Passage", y = "Number of A->G Sites") +
  theme_bindlab(base_size = 12)
save_fig(p, "05_AG_editing_trend.pdf", w = 6, h = 5)

# B2. gained/lost
cat("[B2] gained/lost...\n")
edit_gl <- edit_gained_lost
edit_gl$comparison <- factor(edit_gl$comparison, levels = c("P8 vs P2","P10 vs P2","P12 vs P2"))
gl_long <- reshape2::melt(edit_gl, id.vars = "comparison", variable.name = "type", value.name = "count")
p <- ggplot(gl_long, aes(x = comparison, y = count, fill = type)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.6) +
  scale_fill_manual(values = c("gained" = "#D73027", "lost" = "#4575B4")) +
  labs(title = "Editing Sites: Gained vs Lost (vs P2)", x = "", y = "Number of Sites") +
  theme_bindlab(base_size = 12)
save_fig(p, "05_editing_gained_lost.pdf", w = 7, h = 5)

# B3. 功能分类
cat("[B3] 功能分类...\n")
edit_func$passage <- factor(edit_func$passage, levels = c("P2","P8","P10","P12"))
p <- ggplot(edit_func, aes(x = passage, y = proportion, fill = function_type)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  scale_fill_manual(values = c("synonymous SNV" = "#91BFDB", "nonsynonymous SNV" = "#FC8D59",
  "stopgain" = "#D73027", "unknown" = "grey70")) +
  labs(title = "Functional Distribution of Exonic Editing Sites", x = "Passage", y = "Proportion") +
  theme_bindlab(base_size = 12)
save_fig(p, "05_editing_function_stacked.pdf", w = 7, h = 5)
cat("\n")

# ============================================================
# PART C: 三层交叉整合 (★ v3.1核心修正)
# ============================================================
cat("===== PART C: 三层交叉整合 (DEG x DSE x DE-Editing) =====\n\n")

# --------------------------------------------------
# C1. 构建基因集 — 每层使用"差异变化"基因集
# --------------------------------------------------
cat("[C1] 构建三层基因集...\n")
cat("  ★ v3原则: 每层使用差异变化基因，而非所有检测到的基因\n")
cat("  ★ 文献依据: Ritchie 2015 Nat Rev Genet — 多组学整合要求差异层\n\n")

# Layer 1: DEGs (|LFC|>1 & padj<0.05 in any comparison)
deg_genes <- unique(deg_overlap$gene_id)

# Layer 2: DSE genes (effect-size |dPSI|>=0.1 in any comparison)
dse_genes <- unique(dse_v2$GeneID[dse_v2$GeneID != "" & !is.na(dse_v2$GeneID)])

# Layer 2b: Core DSE genes (effect in ALL 3 comparisons)
dse_core_genes <- unique(dse_core$GeneID[dse_core$GeneID != "" & !is.na(dse_core$GeneID)])

# Layer 3: ★ 差异编辑基因 (|delta_freq|>=0.2 in any comparison)
# 不用"所有编辑基因"(5,222=40.6%背景) — 那会产生虚假富集
edit_de_effect <- edit_all[edit_all$de_effect == TRUE, ]
edit_de_gene_ids <- unique(unlist(regmatches(edit_de_effect$structure_gene,
  gregexpr("ENSG\\d+", edit_de_effect$structure_gene))))

# 参考: 所有编辑基因（仅用于描述统计，不参与交叉检验）
edit_all_gene_ids <- unique(unlist(regmatches(edit_all$structure_gene,
  gregexpr("ENSG\\d+", edit_all$structure_gene))))

cat("  Layer 1 — DEG genes:  ", length(deg_genes),
  " (", round(length(deg_genes)/N_BACKGROUND*100,1), "% of background)\n")
cat("  Layer 2 — DSE genes:  ", length(dse_genes),
  " (", round(length(dse_genes)/N_BACKGROUND*100,1), "%)\n")
cat("  Layer 2b— DSE core genes:  ", length(dse_core_genes),
  " (", round(length(dse_core_genes)/N_BACKGROUND*100,1), "%)\n")
cat("  Layer 3 — DE-Editing genes:  ", length(edit_de_gene_ids),
  " (", round(length(edit_de_gene_ids)/N_BACKGROUND*100,1), "%)\n")
cat("  (ref)  — All editing genes:  ", length(edit_all_gene_ids),
  " (", round(length(edit_all_gene_ids)/N_BACKGROUND*100,1), "% — NOT used for overlap)\n")

cat("\n")

# --------------------------------------------------
# C2. 交叉计数 + Fisher's exact test (★ v3.1: two.sided)
# --------------------------------------------------
cat("[C2] 基因交叉 + Fisher's exact test...\n")
cat("  ★ v3.1修正: Fisher检验使用双侧alternative (two.sided)\n")
cat("  ★ 文献依据: ?fisher.test推荐; ENCODE重叠分析惯例 (PMID:22955616)\n")
cat("  Background N = ", N_BACKGROUND, "\n\n")

# 逐对检验
ft_results <- rbind(
  fisher_overlap_test(deg_genes, dse_genes,  N_BACKGROUND, "DEG", "DSE"),
  fisher_overlap_test(deg_genes, dse_core_genes,  N_BACKGROUND, "DEG", "DSE_core"),
  fisher_overlap_test(deg_genes, edit_de_gene_ids, N_BACKGROUND, "DEG", "DE_Editing"),
  fisher_overlap_test(dse_genes, edit_de_gene_ids, N_BACKGROUND, "DSE", "DE_Editing"),
  fisher_overlap_test(dse_core_genes, edit_de_gene_ids, N_BACKGROUND, "DSE_core", "DE_Editing")
)
ft_results$padj <- p.adjust(ft_results$p_value, method = "BH")
ft_results$significant <- ft_results$padj < 0.05

save_data(ft_results, "05_overlap_fisher_tests.csv")

cat("  Fisher's exact test results:\n")
print(as.data.frame(ft_results))
cat("\n")

# 三层交叉
triple <- Reduce(intersect, list(deg_genes, dse_genes, edit_de_gene_ids))
cat("  DEG ∩ DSE ∩ DE-Editing (三层): ", length(triple), "\n")

# 三层Fisher: 先DEG∩DSE，再测(DEG∩DSE) ∩ DE-Editing
deg_dse_set <- intersect(deg_genes, dse_genes)
deg_dse_core <- intersect(deg_genes, dse_core_genes)  # ★ v3.1修复: 此变量在统计汇总中引用
cat("  DEG n DSE: ", length(deg_dse_set), "\n")
cat("  DEG n DSE (core): ", length(deg_dse_core), "\n")
ft_triple <- fisher_overlap_test(deg_dse_set, edit_de_gene_ids, N_BACKGROUND,
  "DEG_n_DSE", "DE_Editing")
cat("  三层交叉Fisher: overlap=", ft_triple$overlap,
  " expected=", ft_triple$expected,
  " fold=", ft_triple$fold_enrichment,
  " p=", format(ft_triple$p_value, digits = 3), "\n\n")

# --------------------------------------------------
# C3. UpSet图: 三层交叉（使用差异编辑层）
# --------------------------------------------------
cat("[C3] UpSet图...\n")

all_genes_union <- unique(c(deg_genes, dse_genes, edit_de_gene_ids))
upset_df <- data.frame(
  gene_id  = all_genes_union,
  DEG  = as.integer(all_genes_union %in% deg_genes),
  DSE  = as.integer(all_genes_union %in% dse_genes),
  DE_Editing = as.integer(all_genes_union %in% edit_de_gene_ids),
  stringsAsFactors = FALSE
)
upset_df$symbol <- gene_info_map$Symbol[match(upset_df$gene_id, gene_info_map$id)]
save_data(upset_df, "05_three_layer_membership.csv")

pdf(file.path(dir_fig, "05_UpSet_three_layers.pdf"), width = 8, height = 5)
upset(upset_df[, c("DEG","DSE","DE_Editing")],
  sets = c("DE_Editing", "DSE", "DEG"), keep.order = TRUE, order.by = "freq",
  main.bar.color = "#4575B4",
  sets.bar.color = c("#FF7F00", "#984EA3", "#D73027"),
  point.size = 3, line.size = 1,
  text.scale = c(1.5, 1.2, 1.2, 1, 1.5, 1.2),
  mainbar.y.label = "Number of Genes", sets.x.label = "Genes per Layer")
dev.off()
cat("  ✓ saved: 05_UpSet_three_layers.pdf\n")

png(file.path(dir_fig, "05_UpSet_three_layers.png"), width = 8, height = 5, units = "in", res = 300)
upset(upset_df[, c("DEG","DSE","DE_Editing")],
  sets = c("DE_Editing", "DSE", "DEG"), keep.order = TRUE, order.by = "freq",
  main.bar.color = "#4575B4",
  sets.bar.color = c("#FF7F00", "#984EA3", "#D73027"),
  point.size = 3, line.size = 1,
  text.scale = c(1.5, 1.2, 1.2, 1, 1.5, 1.2),
  mainbar.y.label = "Number of Genes", sets.x.label = "Genes per Layer")
dev.off()
cat("\n")

# --------------------------------------------------
# C4. 七区域统计 + Fisher enrichment bar
# --------------------------------------------------
cat("[C4] Venn区域统计...\n")

only_deg  <- sum(upset_df$DEG == 1 & upset_df$DSE == 0 & upset_df$DE_Editing == 0)
only_dse  <- sum(upset_df$DEG == 0 & upset_df$DSE == 1 & upset_df$DE_Editing == 0)
only_edit <- sum(upset_df$DEG == 0 & upset_df$DSE == 0 & upset_df$DE_Editing == 1)
deg_dse_only  <- sum(upset_df$DEG == 1 & upset_df$DSE == 1 & upset_df$DE_Editing == 0)
deg_edit_only <- sum(upset_df$DEG == 1 & upset_df$DSE == 0 & upset_df$DE_Editing == 1)
dse_edit_only <- sum(upset_df$DEG == 0 & upset_df$DSE == 1 & upset_df$DE_Editing == 1)
all_three  <- sum(upset_df$DEG == 1 & upset_df$DSE == 1 & upset_df$DE_Editing == 1)

venn_stats <- data.frame(
  region = c("DEG only","DSE only","DE-Editing only",
  "DEG n DSE","DEG n DE-Edit","DSE n DE-Edit","All three"),
  count = c(only_deg, only_dse, only_edit,
  deg_dse_only, deg_edit_only, dse_edit_only, all_three)
)
save_data(venn_stats, "05_venn_statistics.csv")
cat("  Venn区域统计:\n")
print(venn_stats)

venn_stats$region <- factor(venn_stats$region, levels = rev(venn_stats$region))
p <- ggplot(venn_stats, aes(x = region, y = count)) +
  geom_bar(stat = "identity", fill = "#4575B4", width = 0.65) +
  geom_text(aes(label = count), hjust = -0.15, size = 3.5) +
  coord_flip() +
  labs(title = "Three-layer Overlap: DEG x DSE x DE-Editing",
  subtitle = paste0("Background: ", N_BACKGROUND, " expressed genes; Editing = differential editing only"),
  x = "", y = "Number of Genes") +
  theme_bindlab(base_size = 12)
save_fig(p, "05_venn_barplot.pdf", w = 8, h = 5)

# Fisher enrichment可视化
ft_results$label <- paste0(ft_results$set_A, " n ", ft_results$set_B)
ft_results$neg_log10_padj <- -log10(ft_results$padj + 1e-300)
ft_results$label <- factor(ft_results$label, levels = rev(ft_results$label))

p <- ggplot(ft_results, aes(x = label, y = fold_enrichment, fill = significant)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_text(aes(label = paste0("n=", overlap, " p=", format(padj, digits = 2))),
  hjust = -0.05, size = 3) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  coord_flip() +
  scale_fill_manual(values = c("TRUE" = "#D73027", "FALSE" = "grey60")) +
  labs(title = "Overlap Significance (Fisher's Exact Test)",
  subtitle = paste0("Background N=", N_BACKGROUND, "; dashed line = expected by chance"),
  x = "", y = "Fold Enrichment (observed / expected)") +
  theme_bindlab(base_size = 12) +
  theme(legend.position = "none")
save_fig(p, "05_fisher_enrichment.pdf", w = 9, h = 5)

cat("\n")

# --------------------------------------------------
# C5. DEG∩DSE详细表
# --------------------------------------------------
cat("[C5] DEG∩DSE详细表...\n")

deg_dse_detail <- data.frame(gene_id = deg_dse_set, stringsAsFactors = FALSE)
deg_dse_detail$symbol <- gene_info_map$Symbol[match(deg_dse_detail$gene_id, gene_info_map$id)]

de_p12_sub <- de_p12[, c("gene_id","log2FoldChange","padj","direction")]
colnames(de_p12_sub) <- c("gene_id","DE_lfc_P12","DE_padj_P12","DE_direction_P12")
deg_dse_detail <- merge(deg_dse_detail, de_p12_sub, by = "gene_id", all.x = TRUE)

dse_summary_per_gene <- dse_v2 %>%
  mutate(max_abs_dpsi = pmax(abs(dpsi_P8), abs(dpsi_P10), abs(dpsi_P12), na.rm = TRUE)) %>%
  group_by(GeneID) %>%
  summarise(n_DSE_events = n(), as_types = paste(unique(as_type), collapse = ";"),
  max_dpsi_P8 = dpsi_P8[which.max(max_abs_dpsi)],
  max_dpsi_P10 = dpsi_P10[which.max(max_abs_dpsi)],
  max_dpsi_P12 = dpsi_P12[which.max(max_abs_dpsi)], .groups = "drop")
colnames(dse_summary_per_gene)[1] <- "gene_id"
deg_dse_detail <- merge(deg_dse_detail, dse_summary_per_gene, by = "gene_id", all.x = TRUE)

deg_dse_detail$has_DE_editing <- deg_dse_detail$gene_id %in% edit_de_gene_ids
deg_dse_detail <- deg_dse_detail %>% arrange(DE_padj_P12)
save_data(deg_dse_detail, "05_DEG_DSE_overlap_detail.csv")

triple_detail <- deg_dse_detail %>%
  filter(has_DE_editing == TRUE) %>%
  dplyr::select(gene_id, symbol, DE_lfc_P12, DE_direction_P12,
  n_DSE_events, as_types, max_dpsi_P12, has_DE_editing)
save_data(triple_detail, "05_triple_overlap_genes.csv")

cat("  DEG n DSE: ", nrow(deg_dse_detail), "\n")
cat("  DEG n DSE n DE-Editing: ", nrow(triple_detail), "\n\n")

# ============================================================
# PART D: 多层热图 + 散点图
# ============================================================
cat("===== PART D: 多层可视化 =====\n\n")

# --------------------------------------------------
# D1. DEG∩DSE核心基因: 表达+剪接双层热图
# --------------------------------------------------
cat("[D1] DEG∩DSE(core) 双层热图...\n")

deg_dse_core_detail <- deg_dse_detail %>%
  filter(gene_id %in% dse_core_genes) %>%
  mutate(importance = abs(DE_lfc_P12) + abs(max_dpsi_P12)) %>%
  arrange(desc(importance))

top_n <- min(40, nrow(deg_dse_core_detail))
if (top_n >= 5) {
  top_genes <- deg_dse_core_detail$gene_id[1:top_n]

  vst_sub <- vst_expr[vst_expr$gene_id %in% top_genes, ]
  vst_mat <- as.matrix(vst_sub[, sample_cols])
  rownames(vst_mat) <- vst_sub$symbol
  # 零方差过滤
  vst_mat <- vst_mat[apply(vst_mat, 1, var, na.rm = TRUE) > 1e-10, , drop = FALSE]
  vst_mat_z <- t(scale(t(vst_mat)))
  vst_mat_z[!is.finite(vst_mat_z)] <- 0

  dpsi_per_gene <- dse_v2 %>%
  filter(GeneID %in% top_genes) %>%
  mutate(max_abs = pmax(abs(dpsi_P8), abs(dpsi_P10), abs(dpsi_P12), na.rm = TRUE)) %>%
  group_by(GeneID) %>% slice_max(max_abs, n = 1, with_ties = FALSE) %>% ungroup()
  dpsi_mat <- as.matrix(dpsi_per_gene[, c("dpsi_P8","dpsi_P10","dpsi_P12")])
  rownames(dpsi_mat) <- dpsi_per_gene$Gene_symbol
  colnames(dpsi_mat) <- c("dPSI P8","dPSI P10","dPSI P12")
  dpsi_mat[!is.finite(dpsi_mat)] <- 0

  common_genes <- intersect(rownames(vst_mat_z), rownames(dpsi_mat))
  if (length(common_genes) >= 10) {
  gene_order <- deg_dse_core_detail$symbol[deg_dse_core_detail$symbol %in% common_genes]
  gene_order <- gene_order[1:min(35, length(gene_order))]

  anno_col <- data.frame(Passage = rep(c("P2","P8","P10","P12"), c(2,3,3,3)),
  row.names = sample_cols)
  anno_colors <- list(Passage = passage_colors)

  save_pheatmap("05_DEG_DSE_expression_heatmap.pdf", w = 8,
  h = max(7, length(gene_order) * 0.3 + 3))
  pheatmap(vst_mat_z[gene_order, ], color = heatmap_colors,
  breaks = seq(-2.5, 2.5, length.out = 101),
  cluster_cols = FALSE, clustering_method = "ward.D2",
  annotation_col = anno_col, annotation_colors = anno_colors,
  fontsize_row = 8, fontsize = 10,
  main = "DEG n DSE(core): Expression (Z-score)")
  dev.off()

  save_pheatmap("05_DEG_DSE_dpsi_heatmap.pdf", w = 5,
  h = max(7, length(gene_order) * 0.3 + 3))
  pheatmap(dpsi_mat[gene_order, ], color = colorRampPalette(c("#4575B4","white","#D73027"))(100),
  breaks = seq(-0.4, 0.4, length.out = 101),
  cluster_cols = FALSE, clustering_method = "ward.D2",
  display_numbers = TRUE, number_format = "%.2f",
  fontsize_number = 7, fontsize_row = 8, fontsize = 10,
  main = "DEG n DSE(core): Splicing Changes (dPSI)")
  dev.off()
  cat("  ✓ saved 双层热图\n")
  }
}
cat("\n")

# --------------------------------------------------
# D2. ★ 三组散点图: 表达变化 vs 剪接变化 (v3修正: 展示三个比较)
# --------------------------------------------------
cat("[D2] DEG vs DSE散点图 (三组比较)...\n")

# 合并所有三组DE结果与DSE
de_list <- list(
  P8  = de_p8[,  c("gene_id","symbol","log2FoldChange","padj","direction")],
  P10 = de_p10[, c("gene_id","symbol","log2FoldChange","padj","direction")],
  P12 = de_p12[, c("gene_id","symbol","log2FoldChange","padj","direction")]
)

scatter_all <- lapply(names(de_list), function(pass) {
  df <- de_list[[pass]]
  dpsi_col <- paste0("max_dpsi_", pass)
  merged <- merge(df, dse_summary_per_gene[, c("gene_id", dpsi_col)],
  by = "gene_id", all.x = FALSE)
  colnames(merged)[colnames(merged) == dpsi_col] <- "max_dpsi"
  merged$is_DEG <- merged$direction != "NS"
  merged$is_DSE <- abs(merged$max_dpsi) >= 0.1
  merged$category <- case_when(
  merged$is_DEG & merged$is_DSE ~ "Both",
  merged$is_DEG ~ "DEG only",
  merged$is_DSE ~ "DSE only",
  TRUE ~ "Neither")
  merged$comparison <- paste0(pass, " vs P2")
  merged
})
scatter_combined <- do.call(rbind, scatter_all)
scatter_combined$comparison <- factor(scatter_combined$comparison,
  levels = c("P8 vs P2","P10 vs P2","P12 vs P2"))

# 每面板top标注
top_labels <- scatter_combined %>%
  filter(category == "Both") %>%
  group_by(comparison) %>%
  mutate(score = abs(log2FoldChange) + abs(max_dpsi)) %>%
  slice_max(score, n = 8, with_ties = FALSE) %>%
  ungroup()

p <- ggplot(scatter_combined, aes(x = log2FoldChange, y = max_dpsi, color = category)) +
  geom_point(data = scatter_combined[scatter_combined$category == "Neither", ],
  size = 0.3, alpha = 0.15) +
  geom_point(data = scatter_combined[scatter_combined$category != "Neither", ],
  size = 1, alpha = 0.6) +
  geom_text_repel(data = top_labels, aes(label = symbol),
  size = 2.5, max.overlaps = 12, color = "black") +
  scale_color_manual(values = c("Both" = "#D73027", "DEG only" = "#FC8D59",
  "DSE only" = "#91BFDB", "Neither" = "grey80")) +
  geom_hline(yintercept = c(-0.1, 0.1), linetype = "dashed", color = "grey50", linewidth = 0.3) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey50", linewidth = 0.3) +
  facet_wrap(~comparison, ncol = 3) +
  labs(title = "Expression vs Splicing Changes Across Passages",
  x = "log2 Fold Change (expression)", y = "max dPSI (splicing)") +
  theme_bindlab(base_size = 11) +
  theme(legend.position = "bottom")
save_fig(p, "05_DEG_vs_DSE_scatter_3panel.pdf", w = 14, h = 6)

cat("\n")

# ============================================================
# PART E: 剪接因子-靶标关联 (★ v3重构: SF表达 × 靶事件PSI)
# ============================================================
cat("===== PART E: 剪接因子-DSE关联 (v3修正) =====\n\n")
cat("  ★ v3原则: 关联SF表达与靶事件PSI值，而非靶基因表达\n")
cat("  ★ 文献依据: Hartmann 2021 NAR; Shen 2014 PNAS\n")
cat("  ★ 原方法(SF_expr x DSE_gene_expr)混淆passage共变效应\n\n")

# --------------------------------------------------
# E1. SF表达变化与DSE数量/方向的代次协变（描述性）
# --------------------------------------------------
cat("[E1] SF变化与DSE模式代次协变...\n")

# SF变化汇总 (取显著变化的SF)
sf_changed_genes <- sf_expr$gene[apply(abs(sf_mat), 1, max) > 0.2]

# 构建代次级SF-DSE协变表
passage_sf_dse <- data.frame(
  passage = c("P8","P10","P12"),
  n_DSE_inclusion = c(dir_stats$count[dir_stats$comparison == "P8_vs_P2" & dir_stats$direction == "Inclusion"],
  dir_stats$count[dir_stats$comparison == "P10_vs_P2" & dir_stats$direction == "Inclusion"],
  dir_stats$count[dir_stats$comparison == "P12_vs_P2" & dir_stats$direction == "Inclusion"]),
  n_DSE_exclusion = c(dir_stats$count[dir_stats$comparison == "P8_vs_P2" & dir_stats$direction == "Exclusion"],
  dir_stats$count[dir_stats$comparison == "P10_vs_P2" & dir_stats$direction == "Exclusion"],
  dir_stats$count[dir_stats$comparison == "P12_vs_P2" & dir_stats$direction == "Exclusion"]),
  stringsAsFactors = FALSE
)
passage_sf_dse$inclusion_ratio <- with(passage_sf_dse,
  n_DSE_inclusion / (n_DSE_inclusion + n_DSE_exclusion))
save_data(passage_sf_dse, "05_passage_SF_DSE_covariation.csv")

cat("  代次级DSE方向变化:\n")
print(passage_sf_dse)

# --------------------------------------------------
# E2. SF表达 × DSE core事件PSI 相关性 (代次均值)
# --------------------------------------------------
cat("\n[E2] SF表达 × DSE core PSI 相关性...\n")

# 用代次均值PSI做SF-PSI相关（避免样品层面的passage混杂）
# 这是4个数据点，不做正式检验，仅做描述性热图

# SF代次均值表达（从sf_expr的TPM）
sf_passage_expr <- sf_expr[sf_expr$gene %in% sf_changed_genes, c("gene","P2","P8","P10","P12")]
rownames(sf_passage_expr) <- sf_passage_expr$gene
sf_passage_mat <- as.matrix(sf_passage_expr[, c("P2","P8","P10","P12")])

# DSE core事件PSI代次均值
dse_core_psi <- dse_core[, c("Gene_symbol","psi_P2","psi_P8","psi_P10","psi_P12")]
# 每基因取变化最大的事件
dse_core_ranked <- dse_core %>%
  mutate(max_dpsi = pmax(abs(dpsi_P8), abs(dpsi_P10), abs(dpsi_P12), na.rm = TRUE)) %>%
  group_by(Gene_symbol) %>% slice_max(max_dpsi, n = 1, with_ties = FALSE) %>% ungroup()

top_dse_n <- min(40, nrow(dse_core_ranked))
top_dse <- head(dse_core_ranked, top_dse_n)
top_dse_psi <- as.matrix(top_dse[, c("psi_P2","psi_P8","psi_P10","psi_P12")])
rownames(top_dse_psi) <- top_dse$Gene_symbol
colnames(top_dse_psi) <- c("P2","P8","P10","P12")

# Spearman相关: SF表达(4点) × DSE PSI(4点)
# 注意: 4点的相关系数只有描述意义，不做p-value声明
if (nrow(sf_passage_mat) >= 2 && nrow(top_dse_psi) >= 5) {
  cor_sf_psi <- cor(t(sf_passage_mat), t(top_dse_psi),
  method = "spearman", use = "pairwise.complete.obs")
  # 清理
  cor_sf_psi[!is.finite(cor_sf_psi)] <- 0
  good_r <- rowSums(cor_sf_psi != 0) > 0
  good_c <- colSums(cor_sf_psi != 0) > 0
  cor_sf_psi <- cor_sf_psi[good_r, good_c, drop = FALSE]

  if (nrow(cor_sf_psi) >= 2 && ncol(cor_sf_psi) >= 5) {
  cat("  SF-PSI相关性矩阵:", nrow(cor_sf_psi), "SF x", ncol(cor_sf_psi), "DSE events\n")

  save_pheatmap("05_SF_PSI_correlation.pdf", w = 14,
  h = max(5, nrow(cor_sf_psi) * 0.4 + 2))
  pheatmap(cor_sf_psi,
  color = colorRampPalette(c("#4575B4","white","#D73027"))(100),
  breaks = seq(-1, 1, length.out = 101),
  cluster_rows = TRUE, cluster_cols = TRUE,
  clustering_method = "ward.D2",
  fontsize_row = 9, fontsize_col = 7, fontsize = 10,
  main = "SF Expression x DSE PSI Correlation (passage-level Spearman)\nn=4 passages; descriptive only, no formal p-value inference")
  dev.off()
  cat("  ✓ saved: 05_SF_PSI_correlation.pdf\n")
  cat("  注: n=4数据点的Spearman仅具描述意义，论文中需注明\n")
  }
}

# SF互相关 — tryCatch保护，防止空矩阵崩溃
tryCatch({
  sf_in_vst <- vst_expr[vst_expr$symbol %in% sf_changed_genes, ]
  cat("  SF互相关: sf_in_vst =", nrow(sf_in_vst), "rows\n")
  if (nrow(sf_in_vst) >= 4) {
  sf_expr_mat_vst <- as.matrix(sf_in_vst[, sample_cols])
  rownames(sf_expr_mat_vst) <- sf_in_vst$symbol
  sf_row_var <- apply(sf_expr_mat_vst, 1, var, na.rm = TRUE)
  sf_expr_mat_vst <- sf_expr_mat_vst[sf_row_var > 1e-10, , drop = FALSE]
  cat("  过滤后:", nrow(sf_expr_mat_vst), "rows\n")

  if (nrow(sf_expr_mat_vst) >= 4) {
  sf_cor <- cor(t(sf_expr_mat_vst), method = "pearson", use = "pairwise.complete.obs")
  sf_cor[!is.finite(sf_cor)] <- 0

  save_pheatmap("05_SF_intercorrelation.pdf", w = 8, h = 7)
  pheatmap(sf_cor,
  color = colorRampPalette(c("#4575B4","white","#D73027"))(100),
  breaks = seq(-1, 1, length.out = 101),
  clustering_method = "ward.D2",
  display_numbers = TRUE, number_format = "%.2f",
  fontsize_number = 7, fontsize = 10,
  main = "Splicing Factor Inter-correlation (VST, n=11)\nCaution: passage is a confounder")
  dev.off()
  cat("  ✓ saved: 05_SF_intercorrelation.pdf\n")
  } else {
  cat("  ⚠ 过滤后SF不足4个，跳过互相关\n")
  }
  } else {
  cat("  ⚠ VST中SF不足4个，跳过互相关\n")
  }
}, error = function(e) {
  cat("  ⚠ SF互相关出错，跳过:", e$message, "\n")
  tryCatch(dev.off(), error = function(e2) {})
})
cat("\n")

# ============================================================
# PART F: 三阶段模型多层证据 (★ v3: 区分绝对/差异指标)
# ============================================================
cat("===== PART F: 三阶段模型多层证据 =====\n\n")

cat("[F1] 多层代次趋势图...\n")

trend_data <- data.frame(
  passage_num = c(2, 8, 10, 12),
  passage = factor(c("P2","P8","P10","P12"), levels = c("P2","P8","P10","P12")),
  stringsAsFactors = FALSE
)

# 绝对指标
senmayo_avg <- senmayo_scores %>%
  group_by(passage_num) %>% summarise(SenMayo = mean(SenMayo), .groups = "drop")
trend_data <- merge(trend_data, senmayo_avg, by = "passage_num", all.x = TRUE)
trend_data$editing_sites <- edit_summary$mean_total_sites
trend_data$AG_sites  <- edit_summary$mean_AG_sites

# 差异指标（vs P2）
de_sum <- de_summary
trend_data$DEG_total <- c(0, de_sum$total[de_sum$comparison == "P8_vs_P2"],
  de_sum$total[de_sum$comparison == "P10_vs_P2"],
  de_sum$total[de_sum$comparison == "P12_vs_P2"])
trend_data$DSE_total <- c(0, sum(dse_v2$dse_modT_P8 == TRUE),
  sum(dse_v2$dse_modT_P10 == TRUE),
  sum(dse_v2$dse_modT_P12 == TRUE))

save_data(trend_data, "05_multilayer_trends.csv")

# ★ v3: 分面区分绝对量和差异量
normalize_01 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (rng[2] == rng[1]) return(rep(0.5, length(x)))
  (x - rng[1]) / (rng[2] - rng[1])
}

trend_abs <- data.frame(
  passage = trend_data$passage, passage_num = trend_data$passage_num,
  `SenMayo Score` = normalize_01(trend_data$SenMayo),
  `Editing Sites (total)` = normalize_01(trend_data$editing_sites),
  `A-to-G Sites` = normalize_01(trend_data$AG_sites),
  check.names = FALSE
)
trend_abs_long <- reshape2::melt(trend_abs, id.vars = c("passage","passage_num"),
  variable.name = "Metric", value.name = "normalized")
trend_abs_long$panel <- "Absolute Metrics"

trend_diff <- data.frame(
  passage = trend_data$passage, passage_num = trend_data$passage_num,
  `DEGs (vs P2)` = normalize_01(trend_data$DEG_total),
  `DSEs (vs P2)` = normalize_01(trend_data$DSE_total),
  check.names = FALSE
)
trend_diff_long <- reshape2::melt(trend_diff, id.vars = c("passage","passage_num"),
  variable.name = "Metric", value.name = "normalized")
trend_diff_long$panel <- "Differential Metrics (vs P2)"

trend_all_long <- rbind(trend_abs_long, trend_diff_long)
trend_all_long$panel <- factor(trend_all_long$panel,
  levels = c("Absolute Metrics","Differential Metrics (vs P2)"))

metric_colors <- c("SenMayo Score" = "#D73027",
  "Editing Sites (total)" = "#984EA3",
  "A-to-G Sites" = "#4DAF4A",
  "DEGs (vs P2)" = "#FC8D59",
  "DSEs (vs P2)" = "#4575B4")

p <- ggplot(trend_all_long, aes(x = passage, y = normalized, color = Metric, group = Metric)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3.5) +
  scale_color_manual(values = metric_colors) +
  facet_wrap(~panel, ncol = 2) +
  labs(title = "Multi-layer Passage Dynamics (Normalized 0-1)",
  x = "Passage", y = "Normalized Value") +
  theme_bindlab(base_size = 12) +
  theme(legend.position = "bottom", strip.text = element_text(face = "bold"))
save_fig(p, "05_multilayer_trends.pdf", w = 12, h = 5.5)

cat("\n")

# --------------------------------------------------
# F2. 三阶段特征汇总表
# --------------------------------------------------
cat("[F2] 三阶段汇总...\n")

phase_summary <- data.frame(
  Phase = c("Phase 1 (P8)", "Phase 2 (P10)", "Phase 3 (P12)"),
  Transcriptome = c(
  paste0("IFN burst; 7:1 up:down; ", de_sum$total[1], " DEGs"),
  paste0("SASP peak; metabolic decline; ", de_sum$total[2], " DEGs"),
  paste0("p53 activation; SASP exhaustion; ", de_sum$total[3], " DEGs")),
  SenMayo = round(senmayo_avg$SenMayo[match(c(8,10,12), senmayo_avg$passage_num)], 2),
  Splicing = c(
  paste0(sum(dse_v2$dse_modT_P8 == TRUE), " DSE; early splicing response"),
  paste0(sum(dse_v2$dse_modT_P10 == TRUE), " DSE; sustained"),
  paste0(sum(dse_v2$dse_modT_P12 == TRUE), " DSE; late splicing response")),
  Editing = c(
  paste0(round(edit_summary$mean_total_sites[2]), " sites (+",
  round((edit_summary$mean_total_sites[2]/edit_summary$mean_total_sites[1]-1)*100,1), "%)"),
  paste0(round(edit_summary$mean_total_sites[3]), " sites"),
  paste0(round(edit_summary$mean_total_sites[4]), " sites")),
  stringsAsFactors = FALSE
)
save_data(phase_summary, "05_three_phase_summary.csv")

# --------------------------------------------------
# F3. 多层Heatmap
# --------------------------------------------------
cat("[F3] 多层联合热图...\n")

if (length(triple) >= 15) {
  focus_genes <- triple; focus_label <- "DEG n DSE n DE-Editing"
} else if (length(deg_dse_set) >= 15) {
  focus_genes <- deg_dse_set; focus_label <- "DEG n DSE"
} else {
  focus_genes <- head(lrt_all$gene_id[lrt_all$padj < 0.05], 50)
  focus_label <- "Top LRT genes"
}

focus_lrt <- lrt_all %>% filter(gene_id %in% focus_genes) %>% arrange(padj)
top_ids <- head(focus_lrt$gene_id, 50)

focus_vst <- vst_expr[vst_expr$gene_id %in% top_ids, ]
# ★ 过滤掉无symbol的基因（避免NA行名崩溃）
focus_vst <- focus_vst[!is.na(focus_vst$symbol) & focus_vst$symbol != "" & focus_vst$symbol != "-", ]
focus_mat <- as.matrix(focus_vst[, sample_cols])
rownames(focus_mat) <- make.unique(focus_vst$symbol)
focus_mat <- focus_mat[apply(focus_mat, 1, var, na.rm = TRUE) > 1e-10, , drop = FALSE]
focus_mat_z <- t(scale(t(focus_mat)))
focus_mat_z[!is.finite(focus_mat_z)] <- 0

if (nrow(focus_mat_z) >= 5) {
  # 构建行注释: 用gene_id直接匹配（更稳健）
  focus_vst_filt <- focus_vst[match(rownames(focus_mat_z), make.unique(focus_vst$symbol)), ]
  filt_ids <- focus_vst_filt$gene_id

  anno_row <- data.frame(
  row.names = rownames(focus_mat_z),
  DE_P12 = de_p12$direction[match(filt_ids, de_p12$gene_id)],
  has_DSE = ifelse(filt_ids %in% dse_genes, "Yes", "No"),
  has_DE_Edit = ifelse(filt_ids %in% edit_de_gene_ids, "Yes", "No")
  )
  anno_row$DE_P12[is.na(anno_row$DE_P12)] <- "NS"

  anno_col2 <- data.frame(Passage = rep(c("P2","P8","P10","P12"), c(2,3,3,3)),
  row.names = sample_cols)
  anno_colors2 <- list(Passage = passage_colors,
  DE_P12 = c("Up" = "#D73027", "Down" = "#4575B4", "NS" = "grey70"),
  has_DSE = c("Yes" = "#984EA3", "No" = "grey90"),
  has_DE_Edit = c("Yes" = "#FF7F00", "No" = "grey90"))

  h_fig <- max(8, nrow(focus_mat_z) * 0.25 + 3)
  save_pheatmap("05_multilayer_heatmap.pdf", w = 10, h = h_fig)
  pheatmap(focus_mat_z, color = heatmap_colors, breaks = seq(-2.5, 2.5, length.out = 101),
  cluster_cols = FALSE, clustering_method = "ward.D2",
  annotation_col = anno_col2, annotation_row = anno_row,
  annotation_colors = anno_colors2, fontsize_row = 7, fontsize = 10,
  main = paste0("Multi-layer Heatmap: Top ", nrow(focus_mat_z), " ", focus_label, " Genes"))
  dev.off()

  png(file.path(dir_fig, "05_multilayer_heatmap.png"), width = 10, height = h_fig, units = "in", res = 300)
  pheatmap(focus_mat_z, color = heatmap_colors, breaks = seq(-2.5, 2.5, length.out = 101),
  cluster_cols = FALSE, clustering_method = "ward.D2",
  annotation_col = anno_col2, annotation_row = anno_row,
  annotation_colors = anno_colors2, fontsize_row = 7, fontsize = 10,
  main = paste0("Multi-layer Heatmap: Top ", nrow(focus_mat_z), " ", focus_label, " Genes"))
  dev.off()
  cat("  ✓ saved: 05_multilayer_heatmap.pdf\n")
}
cat("\n")

# ============================================================
# PART G: 通路层面整合
# ============================================================
cat("===== PART G: 通路层面整合 =====\n\n")

dse_in_lrt <- lrt_all %>% filter(gene_id %in% dse_genes, padj < 0.05)
cat("  DSE基因 ∩ LRT显著:", nrow(dse_in_lrt), "/", length(dse_genes), "\n")

deg_only_set <- setdiff(deg_genes, dse_genes)
dse_only_set <- setdiff(dse_genes, deg_genes)
both_set  <- intersect(deg_genes, dse_genes)

layer_category_df <- data.frame(
  gene_id = c(deg_only_set, dse_only_set, both_set),
  layer = c(rep("DEG only", length(deg_only_set)),
  rep("DSE only", length(dse_only_set)),
  rep("DEG & DSE", length(both_set))),
  stringsAsFactors = FALSE
)
layer_category_df$symbol <- gene_info_map$Symbol[match(layer_category_df$gene_id, gene_info_map$id)]
save_data(layer_category_df, "05_gene_layer_categories.csv")

cat("  DEG only:", length(deg_only_set), " | DSE only:", length(dse_only_set),
  " | Both:", length(both_set), "\n")

layer_bar <- data.frame(
  category = factor(c("DEG only", "DEG n DSE", "DSE only"),
  levels = c("DEG only", "DEG n DSE", "DSE only")),
  count = c(length(deg_only_set), length(both_set), length(dse_only_set))
)
p <- ggplot(layer_bar, aes(x = category, y = count, fill = category)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_text(aes(label = count), vjust = -0.3, size = 4) +
  scale_fill_manual(values = c("DEG only" = "#D73027", "DEG n DSE" = "#984EA3", "DSE only" = "#4575B4")) +
  labs(title = "Gene Regulation Layers: Expression vs Splicing", x = "", y = "Number of Genes") +
  theme_bindlab(base_size = 12) + theme(legend.position = "none")
save_fig(p, "05_DEG_DSE_layer_comparison.pdf", w = 6, h = 5)

cat("\n")

# ============================================================
# 保存统计汇总
# ============================================================

integration_stats <- data.frame(
  metric = c("Background (expressed genes)",
  "DEG genes", "DSE events (|dPSI|>=0.1)", "DSE genes",
  "DSE core events (all 3)", "DSE core genes",
  "DE-Editing genes (|delta|>=0.2)",
  "All editing genes (ref, NOT used for overlap)",
  "DEG n DSE", "DEG n DSE (core)",
  "DEG n DE-Editing", "DSE n DE-Editing",
  "DEG n DSE n DE-Editing",
  "Fisher DEG-DSE fold", "Fisher DEG-DSE padj",
  "Fisher DEG-DEedit fold", "Fisher DEG-DEedit padj"),
  value = c(N_BACKGROUND,
  length(deg_genes), nrow(dse_v2), length(dse_genes),
  nrow(dse_core), length(dse_core_genes),
  length(edit_de_gene_ids),
  length(edit_all_gene_ids),
  length(deg_dse_set), length(deg_dse_core),
  length(intersect(deg_genes, edit_de_gene_ids)),
  length(intersect(dse_genes, edit_de_gene_ids)),
  length(triple),
  ft_results$fold_enrichment[1], ft_results$padj[1],
  ft_results$fold_enrichment[3], ft_results$padj[3]),
  stringsAsFactors = FALSE
)
save_data(integration_stats, "05_integration_statistics.csv")

cat("====================================================\n")
cat("  Step 5 完成! (v3.1 — Fisher检验方向修正)\n")
cat("====================================================\n")
cat("\n  v3.1修正:\n")
cat("  1. Fisher's exact test: alternative='greater' → 'two.sided'\n")
cat("  2. 依据: ?fisher.test推荐; 反富集场景(observed<expected)下'greater'返回p≈1.0\n")
cat("\n  v3关键修正:\n")
cat("  1. Editing层: 差异编辑基因(917) 替代 全部编辑基因(5,222)\n")
cat("  2. Fisher's exact test: 所有交叉都有统计检验\n")
cat("  3. SF关联: SF表达×靶事件PSI (非SF表达×靶基因表达)\n")
cat("  4. 趋势图: 绝对指标/差异指标分面展示\n")
cat("  5. 散点图: 三组比较面板(P8/P10/P12)\n")
cat("\n")
cat("  请上传以下文件给我判读:\n")
cat("  1. data/05_integration_statistics.csv\n")
cat("  2. data/05_overlap_fisher_tests.csv  ★ 已修正 (two.sided)\n")
cat("  3. data/05_venn_statistics.csv\n")
cat("  4. figures/05_fisher_enrichment.pdf  ★ 已修正\n")
cat("  5. figures/05_multilayer_trends.pdf\n")
cat("  6. figures/05_UpSet_three_layers.pdf\n")
cat("  7. figures/05_DEG_vs_DSE_scatter_3panel.pdf\n")
cat("====================================================\n")
