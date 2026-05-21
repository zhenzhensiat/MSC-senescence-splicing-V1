# 12_editing_splicing_crosstalk.R — RNA editing and splicing site crosstalk analysis

source("00_project_config.R")

load_packages(c("ggplot2", "dplyr", "tidyr", "data.table",
                "reshape2", "scales", "ggrepel"))

cat("========== Step 12: Editing-Splicing Crosstalk ==========\n\n")


# ============================================================
# PART A: Splice site坐标提取
# ============================================================
cat("===== PART A: Splice site坐标 =====\n\n")

cat("[A1] 读取rMATS外显子坐标...\n")

# --- SE ---
se_parts <- list.files(dir_raw, pattern = "SE.MATS.JC", full.names = TRUE)
if (length(se_parts) == 0)
  se_parts <- list.files(".", pattern = "SE_MATS_JC_part.*\\.xls$", full.names = TRUE)

se_list <- lapply(se_parts, function(f) {
  tryCatch(fread(f, sep = "\t", select = c("AS_ID","Chr","Strand",
    "exonStart_0base","exonEnd","upstreamES","upstreamEE",
    "downstreamES","downstreamEE")), error = function(e) NULL)
})
se_coords <- rbindlist(se_list[!sapply(se_list, is.null)], fill = TRUE)
se_coords <- unique(se_coords, by = "AS_ID")

se_ss <- rbind(
  data.table(chr=se_coords$Chr, pos=se_coords$exonStart_0base+1, AS_ID=se_coords$AS_ID, site_type="exon_start", as_type="SE"),
  data.table(chr=se_coords$Chr, pos=se_coords$exonEnd, AS_ID=se_coords$AS_ID, site_type="exon_end", as_type="SE"),
  data.table(chr=se_coords$Chr, pos=se_coords$upstreamEE, AS_ID=se_coords$AS_ID, site_type="upstream_end", as_type="SE"),
  data.table(chr=se_coords$Chr, pos=se_coords$downstreamES+1, AS_ID=se_coords$AS_ID, site_type="downstream_start", as_type="SE")
)
cat("  SE:", nrow(se_ss), "splice sites from", nrow(se_coords), "events\n")

# --- A3SS ---
a3ss <- tryCatch(fread(f_A3SS_JC, sep="\t", select=c("AS_ID","Chr",
  "longExonStart_0base","longExonEnd","shortES","shortEE","flankingES","flankingEE")), error=function(e) NULL)
a3ss_ss <- NULL
if (!is.null(a3ss) && nrow(a3ss) > 0) {
  a3ss_ss <- rbind(
    data.table(chr=a3ss$Chr, pos=a3ss$longExonStart_0base+1, AS_ID=a3ss$AS_ID, site_type="long_start", as_type="A3SS"),
    data.table(chr=a3ss$Chr, pos=a3ss$longExonEnd, AS_ID=a3ss$AS_ID, site_type="long_end", as_type="A3SS"),
    data.table(chr=a3ss$Chr, pos=a3ss$shortES+1, AS_ID=a3ss$AS_ID, site_type="short_start", as_type="A3SS"),
    data.table(chr=a3ss$Chr, pos=a3ss$shortEE, AS_ID=a3ss$AS_ID, site_type="short_end", as_type="A3SS"))
  cat("  A3SS:", nrow(a3ss_ss), "\n")
}

# --- A5SS ---
a5ss <- tryCatch(fread(f_A5SS_JC, sep="\t", select=c("AS_ID","Chr",
  "longExonStart_0base","longExonEnd","shortES","shortEE","flankingES","flankingEE")), error=function(e) NULL)
a5ss_ss <- NULL
if (!is.null(a5ss) && nrow(a5ss) > 0) {
  a5ss_ss <- rbind(
    data.table(chr=a5ss$Chr, pos=a5ss$longExonStart_0base+1, AS_ID=a5ss$AS_ID, site_type="long_start", as_type="A5SS"),
    data.table(chr=a5ss$Chr, pos=a5ss$longExonEnd, AS_ID=a5ss$AS_ID, site_type="long_end", as_type="A5SS"),
    data.table(chr=a5ss$Chr, pos=a5ss$shortES+1, AS_ID=a5ss$AS_ID, site_type="short_start", as_type="A5SS"),
    data.table(chr=a5ss$Chr, pos=a5ss$shortEE, AS_ID=a5ss$AS_ID, site_type="short_end", as_type="A5SS"))
  cat("  A5SS:", nrow(a5ss_ss), "\n")
}

# --- RI ---
ri <- tryCatch(fread(f_RI_JC, sep="\t", select=c("AS_ID","Chr",
  "riExonStart_0base","riExonEnd","upstreamES","upstreamEE",
  "downstreamES","downstreamEE")), error=function(e) NULL)
ri_ss <- NULL
if (!is.null(ri) && nrow(ri) > 0) {
  ri_ss <- rbind(
    data.table(chr=ri$Chr, pos=ri$upstreamEE, AS_ID=ri$AS_ID, site_type="upstream_end", as_type="RI"),
    data.table(chr=ri$Chr, pos=ri$downstreamES+1, AS_ID=ri$AS_ID, site_type="downstream_start", as_type="RI"))
  cat("  RI:", nrow(ri_ss), "\n")
}

# --- MXE (列名: 1stExonStart_0base等，R读入时可能加反引号) ---
mxe_ss <- NULL
tryCatch({
  mxe_raw <- fread(f_MXE_JC, sep = "\t")
  cn <- colnames(mxe_raw)
  e1s <- cn[grep("1stExonStart", cn)][1]
  e1e <- cn[grep("1stExonEnd", cn)][1]
  e2s <- cn[grep("2ndExonStart", cn)][1]
  e2e <- cn[grep("2ndExonEnd", cn)][1]
  if (!is.na(e1s) && !is.na(e2s)) {
    mxe_ss <- rbind(
      data.table(chr=mxe_raw$Chr, pos=mxe_raw[[e1s]]+1, AS_ID=mxe_raw$AS_ID, site_type="exon1_start", as_type="MXE"),
      data.table(chr=mxe_raw$Chr, pos=mxe_raw[[e1e]],   AS_ID=mxe_raw$AS_ID, site_type="exon1_end", as_type="MXE"),
      data.table(chr=mxe_raw$Chr, pos=mxe_raw[[e2s]]+1, AS_ID=mxe_raw$AS_ID, site_type="exon2_start", as_type="MXE"),
      data.table(chr=mxe_raw$Chr, pos=mxe_raw[[e2e]],   AS_ID=mxe_raw$AS_ID, site_type="exon2_end", as_type="MXE"))
    cat("  MXE:", nrow(mxe_ss), "\n")
  }
}, error = function(e) cat("  MXE读取跳过:", conditionMessage(e), "\n"))

# 合并
all_ss <- rbindlist(list(se_ss, a3ss_ss, a5ss_ss, ri_ss, mxe_ss), fill = TRUE)
all_ss <- all_ss[!is.na(pos) & !is.na(chr)]
all_ss[, chr := sub("^chr", "", chr)]

cat("\n  总splice sites:", nrow(all_ss), "\n")

# 标记DSE
cat("\n[A2] 标记DSE...\n")
dse_v2 <- read.csv(file.path(dir_data, "03_DSE_events_v2.csv"), stringsAsFactors = FALSE)
dse_ids <- unique(dse_v2$AS_ID)
all_ss[, is_DSE := AS_ID %in% dse_ids]
cat("  DSE splice sites:", sum(all_ss$is_DSE), "/ 总:", nrow(all_ss), "\n\n")


# ============================================================
# PART B: Editing site到最近splice site的距离
# ============================================================
cat("===== PART B: 距离计算 =====\n\n")

edit_all <- fread(file.path(dir_data, "04_DE_editing_all.csv"))
edit_unique <- edit_all[!duplicated(site_id)]
edit_unique[, chr := sub("^chr", "", chr)]

# 标记DE-editing (跨3个比较取OR)
de_sites <- edit_all[de_effect == TRUE | de_effect == "TRUE", unique(site_id)]
edit_unique[, is_DE := site_id %in% de_sites]

# KW和AG标记
edit_unique[, is_KW := !is.na(kw_padj) & kw_padj < 0.05]
edit_unique[, is_AG := editing_type == "A->G"]

cat("  Editing位点:", nrow(edit_unique), "\n")
cat("  DE-editing:", sum(edit_unique$is_DE), "\n")
cat("  A->G:", sum(edit_unique$is_AG), "\n")

# 距离计算
cat("\n[B2] 计算距离 (可能1-2分钟)...\n")

chromosomes <- intersect(unique(edit_unique$chr), unique(all_ss$chr))

edit_unique[, dist_any_ss := NA_real_]
edit_unique[, dist_DSE_ss := NA_real_]
edit_unique[, nearest_AS_ID := NA_character_]
edit_unique[, nearest_is_DSE := NA]

for (ch in chromosomes) {
  idx <- which(edit_unique$chr == ch)
  ss_ch <- all_ss[chr == ch]
  dse_ss_ch <- ss_ch[is_DSE == TRUE]
  if (nrow(ss_ch) == 0 || length(idx) == 0) next

  ep <- edit_unique$loci[idx]
  sp <- ss_ch$pos

  for (i in seq_along(idx)) {
    dists <- abs(ep[i] - sp)
    mi <- which.min(dists)
    set(edit_unique, idx[i], "dist_any_ss", dists[mi])
    set(edit_unique, idx[i], "nearest_AS_ID", ss_ch$AS_ID[mi])
    set(edit_unique, idx[i], "nearest_is_DSE", ss_ch$is_DSE[mi])
    if (nrow(dse_ss_ch) > 0)
      set(edit_unique, idx[i], "dist_DSE_ss", min(abs(ep[i] - dse_ss_ch$pos)))
  }
  if (which(chromosomes == ch) %% 5 == 0)
    cat("  chr", ch, "done\n")
}

cat("  完成。有距离:", sum(!is.na(edit_unique$dist_any_ss)), "/", nrow(edit_unique), "\n\n")


# ============================================================
# PART C: 距离分布 + 富集检验
# ============================================================
cat("===== PART C: 距离分布 + 富集 =====\n\n")

edit_wd <- edit_unique[!is.na(dist_any_ss)]

breaks_dist <- c(0, 50, 200, 500, 2000, 10000, Inf)
labels_dist <- c("0-50bp", "50-200bp", "200-500bp", "0.5-2kb", "2-10kb", ">10kb")
edit_wd[, dist_bin := cut(dist_any_ss, breaks = breaks_dist, labels = labels_dist, right = FALSE)]

# --- C1. 距离分布表 ---
cat("[C1] 距离分布...\n")
dist_by_de <- edit_wd[, .N, by = .(is_DE, dist_bin)]
dist_total <- edit_wd[, .N, by = is_DE]
dist_by_de <- merge(dist_by_de, dist_total, by = "is_DE", suffixes = c("", "_total"))
dist_by_de[, pct := round(N / N_total * 100, 2)]
save_data(as.data.frame(dist_by_de), "12_editing_splice_distance_distribution.csv")

# 安全打印: 不用dcast重命名
cat("  DE-editing距离分布:\n")
print(dist_by_de[is_DE == TRUE, .(dist_bin, N, pct)])
cat("\n  Non-DE距离分布:\n")
print(dist_by_de[is_DE == FALSE, .(dist_bin, N, pct)])

# --- C2. 密度图 ---
cat("\n[C2] 密度图...\n")
plot_data <- edit_wd[dist_any_ss <= 10000]
plot_data[, de_label := ifelse(is_DE, "DE-editing", "Non-DE editing")]

p1 <- ggplot(plot_data, aes(x = dist_any_ss, fill = de_label, color = de_label)) +
  geom_density(alpha = 0.3, linewidth = 0.8) +
  scale_x_log10(labels = scales::comma, breaks = c(10,50,200,500,2000,10000)) +
  geom_vline(xintercept = c(50, 200), linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = c("DE-editing"="#D73027", "Non-DE editing"="grey60")) +
  scale_color_manual(values = c("DE-editing"="#D73027", "Non-DE editing"="grey60")) +
  labs(title = "Distance from Editing Sites to Nearest Splice Site",
       subtitle = "Dashed lines: 50bp and 200bp thresholds",
       x = "Distance to Nearest Splice Site (bp, log10)", y = "Density", fill = "", color = "") +
  theme_bindlab(base_size = 12) +
  theme(legend.position = c(0.8, 0.8))
save_fig(p1, "12_editing_splice_distance_density.pdf", w = 9, h = 6)

# A->G专门
plot_ag <- edit_wd[dist_any_ss <= 10000 & is_AG == TRUE]
plot_ag[, de_label := ifelse(is_DE, "DE A-to-G", "Non-DE A-to-G")]
p1b <- ggplot(plot_ag, aes(x = dist_any_ss, fill = de_label, color = de_label)) +
  geom_density(alpha = 0.3, linewidth = 0.8) +
  scale_x_log10(labels = scales::comma) +
  geom_vline(xintercept = c(50, 200), linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = c("DE A-to-G"="#D73027", "Non-DE A-to-G"="#91BFDB")) +
  scale_color_manual(values = c("DE A-to-G"="#D73027", "Non-DE A-to-G"="#91BFDB")) +
  labs(title = "A-to-G (ADAR) Editing Distance to Splice Sites",
       x = "Distance (bp, log10)", y = "Density", fill = "", color = "") +
  theme_bindlab(base_size = 12) +
  theme(legend.position = c(0.8, 0.8))
save_fig(p1b, "12_AG_editing_splice_distance.pdf", w = 9, h = 6)

# --- C3. 富集检验 ---
cat("\n[C3] 富集检验...\n")

thresholds <- c(50, 100, 200, 500)
enrichment_results <- list()
for (thr in thresholds) {
  pDE  <- sum(edit_wd$is_DE == TRUE  & edit_wd$dist_any_ss <= thr, na.rm=TRUE)
  dDE  <- sum(edit_wd$is_DE == TRUE  & edit_wd$dist_any_ss > thr, na.rm=TRUE)
  pNDE <- sum(edit_wd$is_DE == FALSE & edit_wd$dist_any_ss <= thr, na.rm=TRUE)
  dNDE <- sum(edit_wd$is_DE == FALSE & edit_wd$dist_any_ss > thr, na.rm=TRUE)
  ft <- fisher.test(matrix(c(pDE, dDE, pNDE, dNDE), nrow=2))
  enrichment_results[[as.character(thr)]] <- data.frame(
    threshold_bp=thr, proximal_DE=pDE, distal_DE=dDE,
    proximal_nonDE=pNDE, distal_nonDE=dNDE,
    pct_DE_prox=round(100*pDE/(pDE+dDE),2),
    pct_nonDE_prox=round(100*pNDE/(pNDE+dNDE),2),
    OR=round(ft$estimate,3), p_value=ft$p.value, stringsAsFactors=FALSE)
}
enrich_df <- do.call(rbind, enrichment_results)
enrich_df$padj <- p.adjust(enrich_df$p_value, method = "BH")
enrich_df$sig <- ifelse(enrich_df$padj < 0.05, "*", "")

cat("  Splice proximity enrichment (DE vs non-DE):\n")
print(enrich_df)
save_data(enrich_df, "12_splice_proximity_enrichment.csv")

# --- C4. DSE-specific ---
cat("\n[C4] DSE-specific检验...\n")
for (thr in c(200, 500)) {
  nDE  <- sum(edit_wd$is_DE==TRUE  & !is.na(edit_wd$dist_DSE_ss) & edit_wd$dist_DSE_ss<=thr)
  fDE  <- sum(edit_wd$is_DE==TRUE  & (is.na(edit_wd$dist_DSE_ss) | edit_wd$dist_DSE_ss>thr))
  nNDE <- sum(edit_wd$is_DE==FALSE & !is.na(edit_wd$dist_DSE_ss) & edit_wd$dist_DSE_ss<=thr)
  fNDE <- sum(edit_wd$is_DE==FALSE & (is.na(edit_wd$dist_DSE_ss) | edit_wd$dist_DSE_ss>thr))
  if (nDE + nNDE > 0) {
    ft2 <- fisher.test(matrix(c(nDE, fDE, nNDE, fNDE), nrow=2))
    cat(sprintf("  DSE-proximal (%dbp): DE=%d nonDE=%d OR=%.2f p=%.4f\n",
                thr, nDE, nNDE, ft2$estimate, ft2$p.value))
  }
}

cat("\n")


# ============================================================
# PART D: 共定位基因 + 距离bin柱状图
# ============================================================
cat("===== PART D: 共定位基因 =====\n\n")

proximal_edits <- edit_wd[is_DE == TRUE & !is.na(dist_DSE_ss) & dist_DSE_ss <= 200]
cat("  DSE 200bp内的DE-editing:", nrow(proximal_edits), "\n")

if (nrow(proximal_edits) > 0) {
  proximal_as_ids <- unique(proximal_edits$nearest_AS_ID)
  dse_gene_map <- dse_v2[, c("AS_ID", "GeneID", "Gene_symbol")]
  proximal_genes <- unique(dse_gene_map$Gene_symbol[dse_gene_map$AS_ID %in% proximal_as_ids])
  proximal_genes <- proximal_genes[!is.na(proximal_genes) & proximal_genes != ""]

  cat("  共定位基因:", length(proximal_genes), "\n")
  if (length(proximal_genes) <= 50)
    cat("  列表:", paste(head(proximal_genes, 30), collapse=", "), "\n")

  proximal_detail <- merge(
    proximal_edits[, .(site_id, chr, loci, editing_type, dist_DSE_ss, nearest_AS_ID,
                       structure_type, function_type, function_gene)],
    dse_gene_map, by.x="nearest_AS_ID", by.y="AS_ID", all.x=TRUE)
  save_data(as.data.frame(proximal_detail), "12_proximal_editing_DSE_detail.csv")

  deg_overlap <- read.csv(file.path(dir_data, "02_DEG_overlap.csv"), stringsAsFactors=FALSE)
  n_also_deg <- sum(proximal_genes %in% deg_overlap$symbol)
  cat("  其中也是DEG:", n_also_deg, "/", length(proximal_genes), "\n")
} else {
  cat("  无DSE-proximal DE-editing位点\n")
}

# 距离bin柱状图
cat("\n[D2] 距离bin柱状图...\n")
dist_summary <- edit_wd[!is.na(dist_bin), .(
  total = .N, n_DE = sum(is_DE), n_AG = sum(is_AG), n_DE_AG = sum(is_DE & is_AG)
), by = dist_bin]
dist_summary[, pct_DE := n_DE / total * 100]
dist_summary$dist_bin <- factor(dist_summary$dist_bin, levels = labels_dist)

genome_avg <- 100 * sum(edit_wd$is_DE) / nrow(edit_wd)

p2 <- ggplot(dist_summary, aes(x = dist_bin, y = pct_DE)) +
  geom_bar(stat = "identity", fill = "#D73027", alpha = 0.8, width = 0.6) +
  geom_hline(yintercept = genome_avg, linetype = "dashed", color = "grey40") +
  annotate("text", x = 5.5, y = genome_avg + 1, label = "Genome-wide avg",
           size = 3, color = "grey40") +
  labs(title = "DE-editing Fraction by Distance to Nearest Splice Site",
       subtitle = "If editing affects splicing, expect enrichment at short distances",
       x = "Distance to Nearest Splice Site", y = "% DE-editing Sites") +
  theme_bindlab(base_size = 12)
save_fig(p2, "12_DE_editing_by_distance_bin.pdf", w = 8, h = 5)

save_data(as.data.frame(dist_summary), "12_distance_bin_summary.csv")

cat("\n")


# ============================================================
# PART E: Editing类型 x splice proximity
# ============================================================
cat("===== PART E: Editing类型 x proximity =====\n\n")

type_proximal <- edit_wd[, .(
  n_total = .N,
  n_prox_200 = sum(dist_any_ss <= 200, na.rm=TRUE),
  pct_prox = round(100 * sum(dist_any_ss <= 200, na.rm=TRUE) / .N, 2)
), by = editing_type]
type_proximal <- type_proximal[order(-pct_prox)]
cat("  按editing类型的splice-proximal比例 (N>=50):\n")
print(type_proximal[n_total >= 50])
save_data(as.data.frame(type_proximal), "12_editing_type_proximity.csv")

# A->G vs others Fisher
ag_n  <- sum(edit_wd$is_AG & edit_wd$dist_any_ss <= 200, na.rm=TRUE)
ag_f  <- sum(edit_wd$is_AG & edit_wd$dist_any_ss > 200, na.rm=TRUE)
ot_n  <- sum(!edit_wd$is_AG & edit_wd$dist_any_ss <= 200, na.rm=TRUE)
ot_f  <- sum(!edit_wd$is_AG & edit_wd$dist_any_ss > 200, na.rm=TRUE)
ft_ag <- fisher.test(matrix(c(ag_n, ag_f, ot_n, ot_f), nrow=2))
cat("\n  A->G vs other (200bp): OR=", round(ft_ag$estimate,3),
    " p=", format(ft_ag$p.value, digits=3), "\n")

cat("\n")


# ============================================================
# 完成
# ============================================================
cat("====================================================\n")
cat("  Step 12 完成!\n")
cat("====================================================\n")
cat("\n  请上传:\n")
cat("  1. data/12_splice_proximity_enrichment.csv\n")
cat("  2. data/12_distance_bin_summary.csv\n")
cat("  3. data/12_editing_type_proximity.csv\n")
cat("  4. data/12_proximal_editing_DSE_detail.csv (如有)\n")
cat("  5. figures/12_editing_splice_distance_density.pdf\n")
cat("  6. figures/12_DE_editing_by_distance_bin.pdf\n")
cat("====================================================\n")
