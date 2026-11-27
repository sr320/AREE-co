#!/usr/bin/env Rscript
# Tank-aware mixed-model differential expression analysis for PRJNA735889.
# Uses length-scaled TPM counts from Salmon/tximport, DESeq2 median-of-ratios
# effective library sizes, limma-voom precision weights, and variancePartition::dream
# with a random tank intercept and Satterthwaite degrees of freedom:
#   ~ gametogenesis + condition + (1 | tank)

usage <- paste(
  "usage: Rscript run_prjna735889_tank_mixed_model.R QUANT_DIR DESIGN.csv TX2GENE.tsv OUTPUT_DIR",
  "[--threads=4] [--min-samples=3] [--orig-dds=PATH] [--orig-results=PATH] [--run-unadjusted=no]"
)

args <- commandArgs(trailingOnly = TRUE)
is_option <- grepl("^--", args)
positional <- args[!is_option]
options <- list(threads = "4", min_samples = "3", orig_dds = NA, orig_results = NA, run_unadjusted = "no")
for (arg in args[is_option]) {
  parts <- regmatches(arg, regexec("^--([a-z-]+)=(.*)$", arg))[[1]]
  key <- if (length(parts) == 3) gsub("-", "_", parts[[2]]) else ""
  if (!key %in% names(options)) {
    stop(paste("unknown or malformed option:", arg, "\n", usage))
  }
  options[[key]] <- parts[[3]]
}
if (length(positional) != 4) {
  stop(usage)
}

quant_dir <- positional[[1]]
design_path <- positional[[2]]
tx2gene_path <- positional[[3]]
output_dir <- positional[[4]]
threads <- as.integer(options$threads)
min_samples <- as.integer(options$min_samples)
run_unadjusted <- identical(tolower(options$run_unadjusted), "yes")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(tximport)
  library(DESeq2)
  library(limma)
  library(statmod)
  library(variancePartition)
  library(BiocParallel)
})

samples <- read.csv(design_path, stringsAsFactors = FALSE, check.names = FALSE)
required_cols <- c("sample", "condition", "replicate", "run_accession", "tank", "gametogenesis")
if (!all(required_cols %in% colnames(samples))) {
  stop("design sheet is missing required columns")
}
if (anyDuplicated(samples$sample)) {
  stop("sample names must be unique")
}
samples$condition <- factor(samples$condition, levels = c("control", "acidified"))
samples$gametogenesis <- factor(samples$gametogenesis, levels = c("inactive", "active"))
samples$tank <- factor(samples$tank)
rownames(samples) <- samples$sample

quant_files <- file.path(quant_dir, samples$sample, "quant.sf")
names(quant_files) <- samples$sample
if (!all(file.exists(quant_files))) {
  stop(paste("missing Salmon quantification:", paste(quant_files[!file.exists(quant_files)], collapse = ", ")))
}

tx2gene_raw <- read.delim(tx2gene_path, stringsAsFactors = FALSE, header = FALSE, check.names = FALSE)
if (tx2gene_raw[1, 1] %in% c("transcript_id", "TXNAME")) {
  tx2gene <- read.delim(tx2gene_path, stringsAsFactors = FALSE, header = TRUE, check.names = FALSE)
  tx2gene <- tx2gene[, c("transcript_id", "gene_id")]
} else {
  tx2gene <- tx2gene_raw[, 1:2]
  colnames(tx2gene) <- c("transcript_id", "gene_id")
}
if (!all(grepl("^NCBI:GeneID:", tx2gene$gene_id))) {
  tx2gene$gene_id <- ifelse(
    grepl("^NCBI:GeneID:", tx2gene$gene_id),
    tx2gene$gene_id,
    paste0("NCBI:GeneID:", tx2gene$gene_id)
  )
}

txi_raw <- tximport(
  quant_files,
  type = "salmon",
  tx2gene = tx2gene,
  ignoreTxVersion = FALSE,
  dropInfReps = TRUE,
  countsFromAbundance = "no"
)
txi_ls <- tximport(
  quant_files,
  type = "salmon",
  tx2gene = tx2gene,
  ignoreTxVersion = FALSE,
  dropInfReps = TRUE,
  countsFromAbundance = "lengthScaledTPM"
)

keep_genes <- rownames(txi_raw$counts)[rowSums(round(txi_raw$counts) >= 10) >= min_samples]
if (!is.na(options$orig_dds) && nzchar(options$orig_dds)) {
  dds_orig <- readRDS(options$orig_dds)
  orig_genes <- rownames(dds_orig)
  if (!setequal(keep_genes, orig_genes)) {
    stop(sprintf(
      "Filtered gene set (%d) differs from original DESeq2 object (%d); overlap=%d",
      length(keep_genes), length(orig_genes), length(intersect(keep_genes, orig_genes))
    ))
  }
  keep_genes <- orig_genes
  rm(dds_orig)
}

raw_mat <- txi_raw$counts[keep_genes, samples$sample, drop = FALSE]
counts_mat <- txi_ls$counts[keep_genes, samples$sample, drop = FALSE]
rm(txi_raw, txi_ls)
invisible(gc())

sf <- DESeq2::estimateSizeFactorsForMatrix(round(counts_mat))
eff_lib <- colSums(counts_mat) * (sf / exp(mean(log(sf))))

param <- if (threads > 1) SnowParam(threads, "FORK", progressbar = FALSE) else SerialParam(progressbar = FALSE)

# 1. Primary model: dream ~ gametogenesis + condition + (1 | tank)
form_primary <- ~ gametogenesis + condition + (1 | tank)
design_fixed <- model.matrix(~ gametogenesis + condition, data = samples)
cached_fit_path <- file.path(output_dir, "dream_fit_primary.rds")
if (file.exists(cached_fit_path)) {
  cat("Loading cached dream fit:", cached_fit_path, "\n")
  fit_primary <- readRDS(cached_fit_path)
} else {
  t0 <- proc.time()[3]
  vobj_primary <- limma::voom(counts_mat, design_fixed, lib.size = eff_lib)
  fit_primary <- dream(vobj_primary, form_primary, samples, BPPARAM = param, quiet = TRUE)
  fit_primary <- variancePartition::eBayes(fit_primary, robust = TRUE)
  saveRDS(fit_primary, cached_fit_path)
  cat(sprintf("Completed primary dream fit in %.1f s\n", proc.time()[3] - t0))
}

tt_primary <- variancePartition::topTable(fit_primary, coef = "conditionacidified", number = Inf, sort.by = "none")
if (!is.null(fit_primary$s2.post) && !is.null(fit_primary$stdev.unscaled) &&
    "conditionacidified" %in% colnames(fit_primary$stdev.unscaled)) {
  se_primary <- sqrt(fit_primary$s2.post) * fit_primary$stdev.unscaled[, "conditionacidified"]
  nz <- abs(tt_primary$t) > 1e-8 & is.finite(tt_primary$t)
  max_t_diff <- max(abs(tt_primary$logFC[nz] / se_primary[nz] - tt_primary$t[nz]), na.rm = TRUE)
  cat(sprintf("Primary dream SE vs t-stat max_t_diff: %.3e\n", max_t_diff))
  if (!is.finite(max_t_diff) || max_t_diff > 1e-4) {
    stop(sprintf("SE extraction mismatch with dream t-statistic: max_t_diff=%g", max_t_diff))
  }
} else {
  se_primary <- abs(tt_primary$logFC / tt_primary$t)
}

df_primary <- if (is.matrix(fit_primary$df.total) && "conditionacidified" %in% colnames(fit_primary$df.total)) {
  fit_primary$df.total[, "conditionacidified"]
} else {
  as.numeric(fit_primary$df.total)
}

out_primary <- data.frame(
  feature_id_standardized = rownames(tt_primary),
  log2FoldChange = tt_primary$logFC,
  lfcSE = se_primary,
  stat = tt_primary$t,
  pvalue = tt_primary$P.Value,
  padj = tt_primary$adj.P.Val,
  satterthwaite_df = df_primary,
  ave_expr = tt_primary$AveExpr,
  stringsAsFactors = FALSE
)
out_primary <- out_primary[order(out_primary$pvalue, na.last = TRUE), ]
if (any(!is.finite(out_primary$log2FoldChange)) || any(!is.finite(out_primary$lfcSE)) || any(out_primary$lfcSE <= 0)) {
  stop("Non-finite or non-positive effect/SE in dream output")
}

write.table(
  out_primary,
  file.path(output_dir, "acidified_vs_control_dream_all_genes.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = ""
)
writeLines(paste(deparse(form_primary), collapse = ""), file.path(output_dir, "design_formula.txt"))

# 2. Sensitivity 1: limma-voom + duplicateCorrelation (2-pass)
v1_dc <- limma::voom(counts_mat, design_fixed, lib.size = eff_lib)
cor1 <- limma::duplicateCorrelation(v1_dc, design_fixed, block = samples$tank)
v2_dc <- limma::voom(counts_mat, design_fixed, lib.size = eff_lib, block = samples$tank, correlation = cor1$consensus.correlation)
cor2 <- limma::duplicateCorrelation(v2_dc, design_fixed, block = samples$tank)
fit_dc <- limma::eBayes(limma::lmFit(v2_dc, design_fixed, block = samples$tank, correlation = cor2$consensus.correlation), robust = TRUE)
tt_dc <- limma::topTable(fit_dc, coef = "conditionacidified", number = Inf, sort.by = "none")
se_dc <- sqrt(fit_dc$s2.post) * fit_dc$stdev.unscaled[, "conditionacidified"]

# 3. Sensitivity 2 & 3: tank-level pseudobulk DESeq2 (10 tanks)
tanks <- levels(samples$tank)
pb_counts <- sapply(tanks, function(tk) round(rowSums(raw_mat[, samples$tank == tk, drop = FALSE])))
colnames(pb_counts) <- tanks
pb_meta <- data.frame(
  tank = factor(tanks),
  condition = factor(sapply(tanks, function(tk) as.character(samples$condition[samples$tank == tk][1])), levels = c("control", "acidified")),
  gametogenesis_frac = sapply(tanks, function(tk) mean(samples$gametogenesis[samples$tank == tk] == "active")),
  row.names = tanks
)
dds_pb_adj <- DESeqDataSetFromMatrix(pb_counts, colData = pb_meta, design = ~ gametogenesis_frac + condition)
dds_pb_adj <- DESeq(dds_pb_adj, quiet = TRUE)
res_pb_adj <- results(dds_pb_adj, contrast = c("condition", "acidified", "control"), alpha = 0.05)

dds_pb_unadj <- DESeqDataSetFromMatrix(pb_counts, colData = pb_meta, design = ~ condition)
dds_pb_unadj <- DESeq(dds_pb_unadj, quiet = TRUE)
res_pb_unadj <- results(dds_pb_unadj, contrast = c("condition", "acidified", "control"), alpha = 0.05)

rownames(out_primary) <- out_primary$feature_id_standardized
genes <- out_primary$feature_id_standardized
compare_row <- function(name, lfc, se, padj, df_med) {
  sig <- !is.na(padj) & padj < 0.05
  prim_sig <- out_primary[genes, "padj"] < 0.05
  both_sig <- sig & prim_sig
  data.frame(
    model = name,
    genes_tested = sum(!is.na(padj)),
    significant_fdr_lt_0_05 = sum(sig),
    higher_in_acidified = sum(sig & lfc > 0, na.rm = TRUE),
    lower_in_acidified = sum(sig & lfc < 0, na.rm = TRUE),
    median_se = round(median(se, na.rm = TRUE), 4),
    median_df = round(df_med, 2),
    pearson_r_lfc_vs_primary = round(cor(lfc, out_primary[genes, "log2FoldChange"], use = "complete.obs"), 4),
    overlap_with_primary_fdr_0_05 = sum(both_sig, na.rm = TRUE),
    direction_agreement_overlap = if (sum(both_sig, na.rm = TRUE) > 0) {
      round(mean(sign(lfc[both_sig]) == sign(out_primary[genes, "log2FoldChange"][both_sig])), 4)
    } else NA_real_,
    stringsAsFactors = FALSE
  )
}

sens_rows <- list(
  compare_row("dream_gametogenesis_plus_condition_tank_re", out_primary[genes, "log2FoldChange"], out_primary[genes, "lfcSE"], out_primary[genes, "padj"], median(out_primary$satterthwaite_df)),
  compare_row("limma_voom_duplicateCorrelation", tt_dc[genes, "logFC"], se_dc[match(genes, rownames(tt_dc))], tt_dc[genes, "adj.P.Val"], median(fit_dc$df.total)),
  compare_row("tank_pseudobulk_deseq2_gametogenesis_frac", res_pb_adj[genes, "log2FoldChange"], res_pb_adj[genes, "lfcSE"], res_pb_adj[genes, "padj"], 7.0),
  compare_row("tank_pseudobulk_deseq2_unadjusted", res_pb_unadj[genes, "log2FoldChange"], res_pb_unadj[genes, "lfcSE"], res_pb_unadj[genes, "padj"], 8.0)
)

if (run_unadjusted) {
  form_unadj <- ~ condition + (1 | tank)
  design_unadj <- model.matrix(~ condition, data = samples)
  vobj_unadj <- limma::voom(counts_mat, design_unadj, lib.size = eff_lib)
  fit_unadj <- variancePartition::eBayes(dream(vobj_unadj, form_unadj, samples, BPPARAM = param, quiet = TRUE), robust = TRUE)
  tt_unadj <- variancePartition::topTable(fit_unadj, coef = "conditionacidified", number = Inf, sort.by = "none")
  se_unadj <- sqrt(fit_unadj$s2.post) * fit_unadj$stdev.unscaled[, "conditionacidified"]
  df_unadj <- if (is.matrix(fit_unadj$df.total) && "conditionacidified" %in% colnames(fit_unadj$df.total)) {
    fit_unadj$df.total[, "conditionacidified"]
  } else {
    as.numeric(fit_unadj$df.total)
  }
  sens_rows[[length(sens_rows) + 1]] <- compare_row(
    "dream_unadjusted_condition_tank_re",
    tt_unadj[genes, "logFC"], se_unadj[match(genes, rownames(tt_unadj))], tt_unadj[genes, "adj.P.Val"], median(df_unadj)
  )
}

if (!is.na(options$orig_results) && nzchar(options$orig_results)) {
  orig <- read.delim(options$orig_results, stringsAsFactors = FALSE)
  rownames(orig) <- orig$feature_id_standardized
  sens_rows[[length(sens_rows) + 1]] <- compare_row(
    "oyster_deseq2_unmodeled_tank_prior",
    orig[genes, "log2FoldChange"], orig[genes, "lfcSE"], orig[genes, "padj"], 47.0
  )
}

sens_df <- do.call(rbind, sens_rows)
write.table(sens_df, file.path(output_dir, "sensitivity_comparison.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

summary_df <- data.frame(
  metric = c(
    "genes_tested", "significant_padj_lt_0.05", "acidified_higher", "acidified_lower",
    "median_lfc_se", "median_satterthwaite_df", "min_satterthwaite_df", "max_satterthwaite_df",
    "dupcor_consensus_rho"
  ),
  value = c(
    nrow(out_primary),
    sum(out_primary$padj < 0.05, na.rm = TRUE),
    sum(out_primary$padj < 0.05 & out_primary$log2FoldChange > 0, na.rm = TRUE),
    sum(out_primary$padj < 0.05 & out_primary$log2FoldChange < 0, na.rm = TRUE),
    median(out_primary$lfcSE),
    median(out_primary$satterthwaite_df),
    min(out_primary$satterthwaite_df),
    max(out_primary$satterthwaite_df),
    cor2$consensus.correlation
  )
)
write.table(summary_df, file.path(output_dir, "analysis_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
cat("Wrote tank-aware dream results to", output_dir, "\n")
print(sens_df)
