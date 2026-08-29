args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("usage: Rscript run_salmon_tximport_deseq2.R QUANT_DIR DESIGN.csv TX2GENE.tsv OUTPUT_DIR")
}

suppressPackageStartupMessages({
  library(DESeq2)
  library(tximport)
})

quant_dir <- args[[1]]
design_path <- args[[2]]
tx2gene_path <- args[[3]]
output_dir <- args[[4]]
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
if (!all(c("control", "selected") %in% names(condition_counts)) ||
    !identical(as.integer(condition_counts[c("control", "selected")]), c(3L, 3L))) {
  stop("expected three control and three selected samples")
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
  ignoreTxVersion = FALSE
)

samples$condition <- factor(samples$condition, levels = c("control", "selected"))
rownames(samples) <- samples$sample
dds <- DESeqDataSetFromTximport(txi, colData = samples, design = ~ condition)
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]
dds <- DESeq(dds)

result <- results(dds, contrast = c("condition", "selected", "control"), alpha = 0.05)
result_table <- data.frame(
  feature_id_standardized = rownames(result),
  as.data.frame(result),
  row.names = NULL,
  check.names = FALSE
)
result_table <- result_table[order(result_table$pvalue, na.last = TRUE), ]
write.table(
  result_table,
  file.path(output_dir, "selected_vs_control_deseq2_all_genes.tsv"),
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
  file.path(output_dir, "selected_vs_control_normalized_counts.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

vst_data <- vst(dds, blind = FALSE)
pca_data <- plotPCA(vst_data, intgroup = "condition", returnData = TRUE)
write.csv(pca_data, file.path(output_dir, "sample_pca.csv"), row.names = FALSE)
correlations <- cor(assay(vst_data))
write.csv(correlations, file.path(output_dir, "sample_vst_correlations.csv"), quote = FALSE)

pdf(file.path(output_dir, "sample_pca.pdf"), width = 7, height = 5)
print(plotPCA(vst_data, intgroup = "condition"))
dev.off()

if (requireNamespace("apeglm", quietly = TRUE)) {
  coefficient <- grep("condition_selected_vs_control", resultsNames(dds), value = TRUE)
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
      file.path(output_dir, "selected_vs_control_deseq2_apeglm_ranking.tsv"),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE,
      na = ""
    )
  }
}

saveRDS(dds, file.path(output_dir, "selected_vs_control_dds.rds"))
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
