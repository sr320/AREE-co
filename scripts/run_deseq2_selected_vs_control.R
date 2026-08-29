args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("usage: Rscript run_deseq2_selected_vs_control.R GENE_COUNTS.tsv SAMPLE_SHEET.csv OUTPUT_DIR")
}

suppressPackageStartupMessages(library(DESeq2))

counts_path <- args[[1]]
samples_path <- args[[2]]
output_dir <- args[[3]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

counts_table <- read.delim(counts_path, check.names = FALSE)
samples <- read.csv(samples_path, check.names = FALSE, stringsAsFactors = FALSE)
if (!all(c("sample", "condition") %in% colnames(samples))) {
  stop("sample sheet must contain sample and condition columns")
}

if (!all(samples$sample %in% colnames(counts_table))) {
  stop("not every sample-sheet sample is present in the count matrix")
}
gene_id_column <- if ("gene_id" %in% colnames(counts_table)) "gene_id" else colnames(counts_table)[[1]]
gene_ids <- counts_table[[gene_id_column]]
if (anyDuplicated(gene_ids)) {
  stop("gene identifiers in the count table must be unique")
}
count_matrix <- as.matrix(counts_table[, samples$sample, drop = FALSE])
storage.mode(count_matrix) <- "numeric"
rownames(count_matrix) <- gene_ids
count_matrix <- round(count_matrix)

samples$condition <- factor(samples$condition, levels = c("control", "selected"))
rownames(samples) <- samples$sample
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = samples,
  design = ~ condition
)
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]
dds <- DESeq(dds)

result <- results(dds, contrast = c("condition", "selected", "control"), alpha = 0.05)
result_table <- data.frame(
  feature_id_original = rownames(result),
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

if (requireNamespace("apeglm", quietly = TRUE)) {
  coefficient <- grep("condition_selected_vs_control", resultsNames(dds), value = TRUE)
  if (length(coefficient) == 1) {
    shrunk <- lfcShrink(dds, coef = coefficient, type = "apeglm")
    shrunk_table <- data.frame(
      feature_id_original = rownames(shrunk),
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
