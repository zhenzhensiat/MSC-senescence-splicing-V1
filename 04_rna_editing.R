# 04_rna_editing.R — RNA editing analysis with limma-based differential editing detection

source("00_project_config.R")

load_packages(c("data.table", "limma", "ggplot2", "pheatmap", "reshape2", "dplyr"))

cat("=== Step 4 v2: RNA Editing Analysis ===\n")
cat("  ★ v2: limma on logit(freq) + Kruskal-Wallis global test\n\n")

# ── 4.0 配置 ──
PSI_PSEUDO <- 0.001   # logit伪计数
sample_to_passage <- setNames(as.character(sample_table$passage), sample_table$sample_id)
passage_order <- c("P2", "P8", "P10", "P12")
passage_num_map <- c(P2 = 2, P8 = 8, P10 = 10, P12 = 12)
samples_list <- split(sample_table$sample_id, sample_table$passage)

# ── 4.1 加载数据 ──
cat("\n4.1 Loading data...\n")

editing <- fread(f_editing, sep = "\t")
setnames(editing, trimws(names(editing)))
cat("  Total editing events:", nrow(editing), "\n")

editing[, passage := sample_to_passage[sample_id]]
editing[, passage_num := passage_num_map[passage]]

etype <- fread(f_editing_type, sep = "\t")
setnames(etype, trimws(names(etype)))
etype[, passage := sample_to_passage[sample]]
etype[, passage_num := passage_num_map[passage]]

cat("  Editing type table:", nrow(etype), "rows\n")

# ── 4.2 全局编辑景观 ──
cat("\n4.2 Global editing landscape...\n")

type_cols <- c("A->T","A->C","A->G","T->A","T->C","T->G",
               "C->A","C->T","C->G","G->A","G->T","G->C")

for (p in passage_order) {
  p_data <- etype[passage == p]
  total_m <- mean(p_data[["all"]])
  total_s <- sd(p_data[["all"]])
  ag_m  <- mean(p_data[["A->G"]])
  ag_pct <- ag_m / total_m * 100
  cat(sprintf("    %s: total=%.0f+/-%.0f, A->G=%.0f (%.1f%%)\n",
              p, total_m, total_s, ag_m, ag_pct))
}

# ── 4.3 编辑类型长格式 ──
etype_long <- melt(etype, id.vars = c("sample", "passage", "passage_num", "all"),
                   measure.vars = type_cols,
                   variable.name = "editing_type", value.name = "count")
etype_long[, proportion := count / all]
save_data(etype_long, "04_editing_type_long.csv")


# ── 4.4 位点x样品频率矩阵 ──
cat("\n4.4 Building site x sample frequency matrix...\n")

editing[, site_id := paste(chr, loci, editing_type, sep = ":")]

freq_wide <- dcast(editing, site_id ~ sample_id, value.var = "frequency",
                   fun.aggregate = function(x) x[1])
cat("  Unique sites:", nrow(freq_wide), "\n")

site_meta <- unique(editing[, .(site_id, chr, loci, editing_type,
                                 structure_type, structure_gene,
                                 function_type, function_gene)], by = "site_id")


# ── 4.5 位点动态 ──
cat("\n4.5 Passage-specific site dynamics...\n")

# ★ v2注释: presence/absence基于频率非NA，但NA可能来自低覆盖度
# 而非真正的"无编辑"。论文中需讨论此局限（Ramaswami & Li 2014 Methods）。
presence <- data.table(site_id = freq_wide$site_id)
for (p in passage_order) {
  samps <- intersect(samples_list[[p]], names(freq_wide))
  mat <- as.matrix(freq_wide[, ..samps])
  presence[[p]] <- rowSums(!is.na(mat)) > 0
}

for (p in passage_order) cat(sprintf("  %s: %d sites\n", p, sum(presence[[p]])))

shared_all <- sum(presence$P2 & presence$P8 & presence$P10 & presence$P12)
cat("  Shared across all:", shared_all, "\n")

gained_lost <- data.frame()
for (comp in c("P8", "P10", "P12")) {
  gained <- sum(presence[[comp]] & !presence$P2)
  lost   <- sum(!presence[[comp]] & presence$P2)
  gained_lost <- rbind(gained_lost, data.frame(
    comparison = paste0(comp, " vs P2"), gained = gained, lost = lost))
  cat(sprintf("  %s vs P2: gained=%d, lost=%d\n", comp, gained, lost))
}
save_data(gained_lost, "04_editing_gained_lost.csv")


# ── 4.6 ★ v2: 差异编辑 — 三种方法 ──
cat("\n4.6 Differential editing analysis...\n")

# ----------------------------------------------------------------
# ★ v2注释: Wilcoxon统计不可能性
# P2仅2个生物学重复。Wilcoxon秩和检验的最小双侧p-value:
#   n1=2, n2=3: min p = 0.2 (C(5,2)=10种排列，最极端占1/10×2)
#   n1=2, n2=2: min p = 0.333
# 对~15,000位点做BH校正后，FDR<0.05数学上不可达。
# 因此Wilcoxon de_strict结果恒为0，这不是bug而是数学必然。
#
# 解决方案:
# (a) limma on logit(freq): empirical Bayes跨位点收缩方差，效力更高
# (b) Kruskal-Wallis 4组全局检验: n=11，效力好于n1=2 vs n2=3
# (c) Wilcoxon保留作为参考，但不作为主要结果
# ----------------------------------------------------------------

# === 方法1: Wilcoxon (保留，仅作参考) ===
cat("  [Method 1] Wilcoxon rank-sum (reference only)...\n")

de_results_wilcox <- list()

for (comp in c("P8", "P10", "P12")) {
  p2_samps   <- intersect(samples_list[["P2"]], names(freq_wide))
  comp_samps <- intersect(samples_list[[comp]], names(freq_wide))

  p2_mat   <- as.matrix(freq_wide[, ..p2_samps])
  comp_mat <- as.matrix(freq_wide[, ..comp_samps])

  testable <- (rowSums(!is.na(p2_mat)) >= 1) & (rowSums(!is.na(comp_mat)) >= 2)

  p2_mean   <- rowMeans(p2_mat, na.rm = TRUE)
  comp_mean <- rowMeans(comp_mat, na.rm = TRUE)
  delta_freq <- comp_mean - p2_mean

  pvals <- rep(NA_real_, nrow(freq_wide))
  test_idx <- which(testable)
  for (i in test_idx) {
    g1 <- p2_mat[i, ][!is.na(p2_mat[i, ])]
    g2 <- comp_mat[i, ][!is.na(comp_mat[i, ])]
    if (length(g1) >= 2 & length(g2) >= 2) {
      tryCatch({ pvals[i] <- wilcox.test(g1, g2, exact = FALSE)$p.value },
               error = function(e) {})
    }
  }
  padj <- p.adjust(pvals, method = "BH")

  comp_df <- data.table(
    site_id = freq_wide$site_id, freq_P2 = p2_mean,
    delta_freq = delta_freq, pval_wilcox = pvals, padj_wilcox = padj,
    comparison = paste0(comp, "_vs_P2"))
  set(comp_df, j = paste0("freq_", comp), value = comp_mean)

  n_strict_w <- sum(padj < 0.05 & abs(delta_freq) >= 0.1, na.rm = TRUE)
  cat("    ", comp, "vs P2: Wilcoxon FDR<0.05 =", n_strict_w,
      "(expected: 0, see comment above)\n")

  de_results_wilcox[[comp]] <- comp_df
}

de_wilcox <- rbindlist(de_results_wilcox, fill = TRUE)


# === 方法2: ★ limma on logit(freq) (主要方法) ===
cat("\n  [Method 2] limma on logit(freq) (PRIMARY method)...\n")

# 构建频率矩阵（所有11样品）
all_samps <- sample_table$sample_id
all_samps_avail <- intersect(all_samps, names(freq_wide))
freq_mat <- as.matrix(freq_wide[, ..all_samps_avail])

# 过滤: 至少6个样品有非NA值（过半）
testable_limma <- rowSums(!is.na(freq_mat)) >= 6
cat("    Testable sites (>=6 non-NA):", sum(testable_limma), "\n")

# logit变换
freq_bounded <- pmin(pmax(freq_mat, PSI_PSEUDO), 1 - PSI_PSEUDO)
freq_logit   <- log(freq_bounded / (1 - freq_bounded))
# NA保持为NA（limma可以通过na处理）
freq_logit[is.na(freq_mat)] <- NA

# 构建设计矩阵
passage_factor <- factor(sample_table$passage[match(all_samps_avail, sample_table$sample_id)],
                         levels = passage_order)
design <- model.matrix(~ 0 + passage_factor)
colnames(design) <- levels(passage_factor)

contrasts <- limma::makeContrasts(
  P8_vs_P2  = P8  - P2,
  P10_vs_P2 = P10 - P2,
  P12_vs_P2 = P12 - P2,
  levels = design
)

# 只拟合testable位点
freq_logit_test <- freq_logit[testable_limma, ]

# limma拟合（处理NA: 对每个基因用available样品）
cat("    Running limma on", nrow(freq_logit_test), "sites...\n")
fit <- lmFit(freq_logit_test, design)
fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2)

# 提取结果
de_results_limma <- list()
for (comp in c("P8", "P10", "P12")) {
  comp_name <- paste0(comp, "_vs_P2")
  tt <- topTable(fit2, coef = comp_name, number = Inf, sort.by = "none")

  # delta_freq在原始尺度
  p2_samps   <- intersect(samples_list[["P2"]], all_samps_avail)
  comp_samps <- intersect(samples_list[[comp]], all_samps_avail)
  p2_mean   <- rowMeans(freq_mat[testable_limma, p2_samps, drop = FALSE], na.rm = TRUE)
  comp_mean <- rowMeans(freq_mat[testable_limma, comp_samps, drop = FALSE], na.rm = TRUE)

  limma_df <- data.table(
    site_id    = freq_wide$site_id[testable_limma],
    freq_P2    = p2_mean,
    delta_freq = comp_mean - p2_mean,
    pval_limma = tt$P.Value,
    padj_limma = tt$adj.P.Val,
    comparison = comp_name
  )
  set(limma_df, j = paste0("freq_", comp), value = comp_mean)

  # DE标准（使用limma p-value + 原始尺度effect-size）
  limma_df[, de_strict := (padj_limma < 0.05) & (abs(delta_freq) >= 0.1)]
  limma_df[, de_effect := abs(delta_freq) >= 0.2]

  n_strict <- sum(limma_df$de_strict, na.rm = TRUE)
  n_effect <- sum(limma_df$de_effect, na.rm = TRUE)
  cat("    ", comp, "vs P2: limma FDR<0.05+|df|>=0.1 =", n_strict,
      " | |df|>=0.2 =", n_effect, "\n")

  de_results_limma[[comp]] <- limma_df
}

de_limma <- rbindlist(de_results_limma, fill = TRUE)


# === 方法3: Kruskal-Wallis全局检验 ===
cat("\n  [Method 3] Kruskal-Wallis global test (4 groups)...\n")

kw_pvals <- rep(NA_real_, nrow(freq_mat))
kw_testable <- rowSums(!is.na(freq_mat)) >= 6

test_idx_kw <- which(kw_testable)
passage_vec <- sample_table$passage[match(all_samps_avail, sample_table$sample_id)]

for (i in test_idx_kw) {
  vals <- freq_mat[i, ]
  non_na <- !is.na(vals)
  if (sum(non_na) >= 6 && length(unique(passage_vec[non_na])) >= 3) {
    tryCatch({
      kw_pvals[i] <- kruskal.test(vals[non_na] ~ passage_vec[non_na])$p.value
    }, error = function(e) {})
  }
}

kw_padj <- p.adjust(kw_pvals, method = "BH")
n_kw_sig <- sum(kw_padj < 0.05, na.rm = TRUE)
cat("    KW significant (FDR<0.05):", n_kw_sig, "\n")

kw_df <- data.table(
  site_id  = freq_wide$site_id,
  kw_pval  = kw_pvals,
  kw_padj  = kw_padj,
  kw_sig   = kw_padj < 0.05
)


# ── 4.6b 合并输出 ──
cat("\n  Merging results...\n")

# 主输出使用limma结果（最有统计效力）
de_annotated <- merge(de_limma, site_meta, by = "site_id", all.x = TRUE)

# 合并KW结果
de_annotated <- merge(de_annotated, kw_df[, .(site_id, kw_pval, kw_padj, kw_sig)],
                      by = "site_id", all.x = TRUE)

save_data(de_annotated, "04_DE_editing_all.csv")

# 统计汇总
de_method_summary <- data.frame(
  method = c("Wilcoxon (n1=2,n2=3)", "limma logit (n=11)", "Kruskal-Wallis (n=11)"),
  DE_strict_P8 = c(
    sum(de_wilcox[comparison == "P8_vs_P2"]$padj_wilcox < 0.05 &
        abs(de_wilcox[comparison == "P8_vs_P2"]$delta_freq) >= 0.1, na.rm = TRUE),
    sum(de_limma[comparison == "P8_vs_P2"]$de_strict, na.rm = TRUE),
    n_kw_sig
  ),
  note = c("min p=0.2, FDR impossible", "empirical Bayes, primary",
           "global 4-group test"),
  stringsAsFactors = FALSE
)
save_data(de_method_summary, "04_DE_method_comparison.csv")
cat("\n  Method comparison:\n")
print(de_method_summary)


# ── 4.7 功能分布 ──
cat("\n4.7 Functional distribution...\n")

struct_dist <- editing[, .N, by = .(passage, structure_type)]
struct_dist[, total := sum(N), by = passage]
struct_dist[, proportion := N / total]
save_data(struct_dist, "04_editing_structure_dist.csv")

exonic <- editing[structure_type == "exonic"]
func_dist <- exonic[, .N, by = .(passage, function_type)]
func_dist[, total := sum(N), by = passage]
func_dist[, proportion := N / total]
save_data(func_dist, "04_editing_function_dist.csv")


# ── 4.8 ADAR分析 ──
cat("\n4.8 ADAR analysis...\n")

adar <- editing[editing_type == "A->G"]
adar_stats <- adar[, .(n_events = .N, mean_freq = mean(frequency),
                         median_freq = median(frequency)), by = passage]
print(adar_stats[match(passage_order, passage)])


# ── 4.9 汇总表 ──
cat("\n4.9 Summary...\n")

summary_rows <- list()
for (p in passage_order) {
  p_etype <- etype[passage == p]; p_data <- editing[passage == p]
  row <- data.frame(
    passage = p, passage_num = passage_num_map[p],
    n_replicates = nrow(p_etype),
    mean_total_sites = mean(p_etype[["all"]]),
    sd_total_sites = sd(p_etype[["all"]]),
    mean_AG_sites = mean(p_etype[["A->G"]]),
    mean_AG_proportion = mean(p_etype[["A->G"]] / p_etype[["all"]]),
    mean_TC_sites = mean(p_etype[["T->C"]]),
    mean_CT_sites = mean(p_etype[["C->T"]]),
    mean_GA_sites = mean(p_etype[["G->A"]]),
    unique_sites = sum(presence[[p]]),
    mean_freq_all = mean(p_data$frequency),
    mean_freq_AG = mean(p_data[editing_type == "A->G"]$frequency),
    stringsAsFactors = FALSE)
  summary_rows[[p]] <- row
}
save_data(do.call(rbind, summary_rows), "04_editing_summary.csv")


# ── 4.10 可视化 ──
cat("\n4.10 Figures...\n")

etype[, passage_f := factor(passage, levels = passage_order)]

# Fig 4A: Total sites
p1 <- ggplot(etype, aes(x = passage_f, y = all, fill = passage_f)) +
  geom_bar(stat = "summary", fun = "mean", width = 0.6, color = "black", linewidth = 0.3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
  geom_jitter(width = 0.1, size = 1.5, alpha = 0.6) +
  scale_fill_manual(values = passage_colors) +
  labs(x = "Passage", y = "Total editing sites", title = "Total RNA editing sites") +
  theme_bindlab() + theme(legend.position = "none")
save_fig(p1, "04_total_editing_sites.pdf", w = 4, h = 4)

# Fig 4B: A->G
ag_per_sample <- etype[, .(sample, passage_f = factor(passage, levels = passage_order), AG = `A->G`)]
p2 <- ggplot(ag_per_sample, aes(x = passage_f, y = AG)) +
  geom_bar(stat = "summary", fun = "mean", aes(fill = passage_f),
           width = 0.6, color = "black", linewidth = 0.3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
  geom_jitter(width = 0.1, size = 1.5, alpha = 0.6) +
  scale_fill_manual(values = passage_colors) +
  labs(x = "Passage", y = "A-to-G sites", title = "ADAR-mediated A-to-G editing") +
  theme_bindlab() + theme(legend.position = "none")
save_fig(p2, "04_ADAR_AG_barplot.pdf", w = 4, h = 4)

# Fig 4C: Type composition
major_types <- c("A->G", "T->C", "C->T", "G->A")
type_colors_map <- c("A->G" = "#D73027", "T->C" = "#FC8D59",
                      "C->T" = "#4575B4", "G->A" = "#91BFDB", "Other" = "#CCCCCC")
etype_plot <- copy(etype)
etype_plot[, Other := all - `A->G` - `T->C` - `C->T` - `G->A`]
etype_comp <- melt(etype_plot, id.vars = c("sample", "passage", "all"),
                   measure.vars = c("A->G", "T->C", "C->T", "G->A", "Other"),
                   variable.name = "type", value.name = "count")
etype_comp[, proportion := count / all]
etype_comp[, passage_f := factor(passage, levels = passage_order)]
etype_mean <- etype_comp[, .(proportion = mean(proportion)), by = .(passage_f, type)]

p3 <- ggplot(etype_mean, aes(x = passage_f, y = proportion, fill = type)) +
  geom_bar(stat = "identity", width = 0.6, color = "white", linewidth = 0.2) +
  scale_fill_manual(values = type_colors_map) +
  labs(x = "Passage", y = "Proportion", title = "Editing type composition") +
  theme_bindlab() + theme(legend.title = element_blank())
save_fig(p3, "04_editing_type_composition.pdf", w = 5, h = 4)

# Fig 4D: Structure
struct_plot <- struct_dist[structure_type %in% c("exonic", "UTR3", "UTR5")]
struct_plot[, passage_f := factor(passage, levels = passage_order)]
p4 <- ggplot(struct_plot, aes(x = passage_f, y = proportion, fill = structure_type)) +
  geom_bar(stat = "identity", position = position_dodge(0.7), width = 0.6) +
  scale_fill_manual(values = c("exonic" = "#4575B4", "UTR3" = "#FC8D59", "UTR5" = "#91BFDB")) +
  labs(x = "Passage", y = "Proportion", title = "Genomic region distribution") +
  theme_bindlab() + theme(legend.title = element_blank())
save_fig(p4, "04_editing_structure_dist.pdf", w = 5, h = 4)

# Fig 4E: Gained/lost
gl_long <- reshape2::melt(gained_lost, id.vars = "comparison",
                          variable.name = "direction", value.name = "count")
gl_long[gl_long$direction == "lost", "count"] <- -gl_long[gl_long$direction == "lost", "count"]
p5 <- ggplot(gl_long, aes(x = comparison, y = count, fill = direction)) +
  geom_bar(stat = "identity", position = "identity", width = 0.5, color = "black", linewidth = 0.3) +
  geom_hline(yintercept = 0, linewidth = 0.5) +
  geom_text(aes(label = abs(count)), vjust = ifelse(gl_long$count >= 0, -0.3, 1.3), size = 3) +
  scale_fill_manual(values = c("gained" = "#D73027", "lost" = "#4575B4")) +
  labs(x = NULL, y = "Number of sites", title = "Editing site dynamics (vs P2)") +
  theme_bindlab() + theme(legend.title = element_blank())
save_fig(p5, "04_editing_gain_loss.pdf", w = 5, h = 4)

# Fig 4F: A->G frequency
adar_plot <- editing[editing_type == "A->G"]
adar_plot[, passage_f := factor(passage, levels = passage_order)]
p6 <- ggplot(adar_plot, aes(x = frequency, fill = passage_f)) +
  geom_histogram(binwidth = 0.05, alpha = 0.5, position = "identity", color = "white", linewidth = 0.1) +
  scale_fill_manual(values = passage_colors) +
  labs(x = "Editing frequency", y = "Count", title = "A-to-G frequency distribution") +
  theme_bindlab() + theme(legend.title = element_blank())
save_fig(p6, "04_ADAR_freq_distribution.pdf", w = 5, h = 4)

# Fig 4G: Exonic function
func_exonic <- func_dist[function_type %in% c("synonymous SNV", "nonsynonymous SNV", "stopgain")]
func_exonic[, passage_f := factor(passage, levels = passage_order)]
p7 <- ggplot(func_exonic, aes(x = passage_f, y = proportion, fill = function_type)) +
  geom_bar(stat = "identity", width = 0.6, color = "white", linewidth = 0.2) +
  scale_fill_manual(values = c("synonymous SNV" = "#91BFDB",
                                "nonsynonymous SNV" = "#FC8D59", "stopgain" = "#D73027")) +
  labs(x = "Passage", y = "Proportion", title = "Exonic editing: functional impact") +
  theme_bindlab() + theme(legend.title = element_blank())
save_fig(p7, "04_editing_exonic_function.pdf", w = 5, h = 4)


cat("\n=== Step 4 v2 complete! ===\n")
cat("Key changes:\n")
cat("  1. limma on logit(freq) as PRIMARY DE method\n")
cat("  2. Kruskal-Wallis global test as supplementary\n")
cat("  3. Wilcoxon retained as reference (FDR=0 is expected)\n")
cat("\nOutput: same filenames, updated content\n")
cat("  ★ NEW: data/04_DE_method_comparison.csv\n")
cat("\n  请上传 04_DE_method_comparison.csv 和 04_editing_summary.csv 判读!\n")
