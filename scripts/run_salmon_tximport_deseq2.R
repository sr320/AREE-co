usage <- paste(
  "usage: Rscript run_salmon_tximport_deseq2.R QUANT_DIR DESIGN.csv TX2GENE.tsv OUTPUT_DIR",
  "[--reference=control] [--test=selected] [--replicates=3] [--min-samples=3] [--ma-title=TITLE]"
)
args <- commandArgs(trailingOnly = TRUE)
is_option <- grepl("^--", args)
positional <- args[!is_option]
# Defaults reproduce the original PRJNA694496 selected-versus-control analysis.
options <- list(reference = "control", test = "selected", replicates = "3", min_samples = "3", ma_title = NA)
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
reference <- options$reference
test <- options$test
# Exact per-condition replicate count to require; 0 accepts any count of at least two.
replicates <- as.integer(options$replicates)
min_samples <- as.integer(options$min_samples)
contrast <- paste0(test, "_vs_", reference)
ma_title <- if (is.na(options$ma_title)) paste(tools::toTitleCase(test), "versus", reference) else options$ma_title

suppressPackageStartupMessages({
  library(DESeq2)
  library(tximport)
})

quant_dir <- positional[[1]]
design_path <- positional[[2]]
tx2gene_path <- positional[[3]]
output_dir <- positional[[4]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

samples <- read.csv(design_path, stringsAsFactors = FALSE, check.names = FALSE)
required_sample_columns <- c("sample", "condition", "replicate", "run_accession")
if (!all(required_sample_columns %in% colnames(samples))) {
  stop("design sheet is missing required columns")
}
if (anyDuplicated(samples$sample)) {
  stop("sample names must be unique")
}
condition_counts <- table(samples$condition)
levels_present <- all(c(reference, test) %in% names(condition_counts))
group_sizes <- if (levels_present) as.integer(condition_counts[c(reference, test)]) else integer()
if (!levels_present ||
    (replicates > 0 && !identical(group_sizes, c(replicates, replicates))) ||
    (replicates == 0 && any(group_sizes < 2))) {
  stop(sprintf("expected %s %s and %s %s samples",
               if (replicates > 0) replicates else "at least two", reference,
               if (replicates > 0) replicates else "at least two", test))
}

quant_files <- file.path(quant_dir, samples$sample, "quant.sf")
names(quant_files) <- samples$sample
if (!all(file.exists(quant_files))) {
  stop(paste("missing Salmon quantification:", paste(quant_files[!file.exists(quant_files)], collapse = ", ")))
}

tx2gene <- read.delim(tx2gene_path, stringsAsFactors = FALSE, check.names = FALSE)
if (!all(c("transcript_id", "gene_id") %in% colnames(tx2gene))) {
  stop("tx2gene table must contain transcript_id and gene_id")
}
txi <- tximport(
  quant_files,
  type = "salmon",
  tx2gene = tx2gene[, c("transcript_id", "gene_id")],
  ignoreTxVersion = FALSE,
  # Salmon was run without bootstraps; avoid requiring jsonlite solely to
  # inspect inferential-replicate metadata that do not exist for this run.
  dropInfReps = TRUE
)

samples$condition <- factor(samples$condition, levels = c(reference, test))
rownames(samples) <- samples$sample
dds <- DESeqDataSetFromTximport(txi, colData = samples, design = ~ condition)
genes_imported <- nrow(dds)
keep <- rowSums(counts(dds) >= 10) >= min_samples
dds <- dds[keep, ]
dds <- DESeq(dds)

result <- results(dds, contrast = c("condition", test, reference), alpha = 0.05)
result_table <- data.frame(
  feature_id_standardized = rownames(result),
  as.data.frame(result),
  row.names = NULL,
  check.names = FALSE
)
result_table <- result_table[order(result_table$pvalue, na.last = TRUE), ]
write.table(
  result_table,
  file.path(output_dir, paste0(contrast, "_deseq2_all_genes.tsv")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = ""
)

normalized <- as.data.frame(counts(dds, normalized = TRUE), check.names = FALSE)
normalized$feature_id_standardized <- rownames(normalized)
normalized <- normalized[, c("feature_id_standardized", samples$sample)]
write.table(
  normalized,
  file.path(output_dir, paste0(contrast, "_normalized_counts.tsv")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# vst() subsamples 1000 genes to fit the dispersion trend; smaller sets (e.g. pipeline tests) need the full transform.
vst_data <- if (nrow(dds) >= 1000) vst(dds, blind = FALSE) else varianceStabilizingTransformation(dds, blind = FALSE)
pca_data <- plotPCA(vst_data, intgroup = "condition", returnData = TRUE)
write.csv(pca_data, file.path(output_dir, "sample_pca.csv"), row.names = FALSE)
write.csv(data.frame(component = c("PC1", "PC2"),
                     variance_fraction = attr(pca_data, "percentVar")),
          file.path(output_dir, "sample_pca_variance.csv"), row.names = FALSE)
correlations <- cor(assay(vst_data))
write.csv(correlations, file.path(output_dir, "sample_vst_correlations.csv"), quote = FALSE)

# Retain diagnostics for every library; flags require review, not automatic removal.
cooks <- assays(dds)[["cooks"]]
cooks_cutoff <- qf(0.99, 2, ncol(dds) - 2)
sample_qc <- data.frame(sample = samples$sample, condition = samples$condition,
                        normalization_factor_median = apply(
                          normalizationFactors(dds), 2, median, na.rm = TRUE
                        ),
                        genes_above_cooks_cutoff = colSums(cooks > cooks_cutoff, na.rm = TRUE))
write.table(sample_qc, file.path(output_dir, "sample_deseq2_qc.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
gene_qc <- data.frame(feature_id_standardized = rownames(dds),
                      max_cooks = apply(cooks, 1, max, na.rm = TRUE))
write.table(gene_qc, file.path(output_dir, "gene_cooks_distances.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
summary_table <- data.frame(
  metric = c("genes_imported", "genes_after_count_filter", "genes_with_pvalue",
             "genes_with_padj", "significant_padj_lt_0.05", paste0(test, "_higher"), paste0(test, "_lower"),
             "cooks_cutoff"),
  value = c(genes_imported, nrow(dds), sum(!is.na(result$pvalue)), sum(!is.na(result$padj)),
            sum(result$padj < 0.05, na.rm = TRUE),
            sum(result$padj < 0.05 & result$log2FoldChange > 0, na.rm = TRUE),
            sum(result$padj < 0.05 & result$log2FoldChange < 0, na.rm = TRUE), cooks_cutoff))
write.table(summary_table, file.path(output_dir, "analysis_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

pdf(file.path(output_dir, "sample_pca.pdf"), width = 7, height = 5)
print(plotPCA(vst_data, intgroup = "condition"))
dev.off()

png(file.path(output_dir, "sample_pca.png"), width = 1400, height = 1000, res = 180)
print(plotPCA(vst_data, intgroup = "condition"))
dev.off()
png(file.path(output_dir, paste0(contrast, "_MA.png")), width = 1400, height = 1000, res = 180)
plotMA(result, main = ma_title, alpha = 0.05)
dev.off()
png(file.path(output_dir, "sample_correlations.png"), width = 1200, height = 1200, res = 180)
heatmap(correlations, symm = TRUE, margins = c(10, 10),
        main = "VST sample correlations", scale = "none")
dev.off()

if (requireNamespace("apeglm", quietly = TRUE)) {
  coefficient <- grep(paste0("condition_", contrast), resultsNames(dds), value = TRUE, fixed = TRUE)
  if (length(coefficient) == 1) {
    shrunk <- lfcShrink(dds, coef = coefficient, type = "apeglm")
    shrunk_table <- data.frame(
      feature_id_standardized = rownames(shrunk),
      as.data.frame(shrunk),
      row.names = NULL,
      check.names = FALSE
    )
    write.table(
      shrunk_table,
      file.path(output_dir, paste0(contrast, "_deseq2_apeglm_ranking.tsv")),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE,
      na = ""
    )
  }
}

saveRDS(dds, file.path(output_dir, paste0(contrast, "_dds.rds")))
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
