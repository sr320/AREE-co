args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("usage: Rscript run_deseq2_leave_one_out.R DDS.rds EXCLUDED_SAMPLE FULL_RESULTS.tsv OUTPUT_DIR")
}

suppressPackageStartupMessages(library(DESeq2))

dds_path <- args[[1]]
excluded_sample <- args[[2]]
full_results_path <- args[[3]]
output_dir <- args[[4]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

dds <- readRDS(dds_path)
if (!(excluded_sample %in% colnames(dds))) {
  stop(paste("excluded sample is not present in DDS:", excluded_sample))
}
dds_sensitivity <- dds[, colnames(dds) != excluded_sample]
condition_counts <- table(colData(dds_sensitivity)$condition)
if (length(condition_counts) != 2 || any(condition_counts < 2)) {
  stop("leave-one-out analysis requires at least two samples in each condition")
}

# Re-estimate dispersions and fit the model after sample exclusion while
# retaining the primary analysis gene set and tximport normalization factors.
dds_sensitivity <- DESeq(dds_sensitivity)
result <- results(
  dds_sensitivity,
  contrast = c("condition", "selected", "control"),
  alpha = 0.05
)
result_table <- data.frame(
  feature_id_standardized = rownames(result),
  as.data.frame(result),
  row.names = NULL,
  check.names = FALSE
)
result_table <- result_table[order(result_table$pvalue, na.last = TRUE), ]
write.table(
  result_table,
  file.path(output_dir, paste0("exclude_", excluded_sample, "_deseq2_all_genes.tsv")),
  sep = "\t", quote = FALSE, row.names = FALSE, na = ""
)

full <- read.delim(full_results_path, stringsAsFactors = FALSE, check.names = FALSE)
required <- c("feature_id_standardized", "log2FoldChange", "pvalue", "padj")
if (!all(required %in% colnames(full))) {
  stop("full-results table is missing required columns")
}
comparison <- merge(
  full[, required], result_table[, required],
  by = "feature_id_standardized", suffixes = c("_full", "_leave_one_out")
)
comparison$same_direction <- with(
  comparison,
  sign(log2FoldChange_full) == sign(log2FoldChange_leave_one_out)
)
comparison$significant_full <- comparison$padj_full < 0.05
comparison$significant_leave_one_out <- comparison$padj_leave_one_out < 0.05
write.table(
  comparison,
  file.path(output_dir, paste0("exclude_", excluded_sample, "_comparison.tsv")),
  sep = "\t", quote = FALSE, row.names = FALSE, na = ""
)

finite_effects <- is.finite(comparison$log2FoldChange_full) &
  is.finite(comparison$log2FoldChange_leave_one_out)
full_significant <- comparison$significant_full %in% TRUE
summary_table <- data.frame(
  metric = c(
    "excluded_sample", "remaining_control_samples", "remaining_selected_samples",
    "genes_compared", "effect_pearson_correlation", "same_direction_all",
    "significant_full", "significant_leave_one_out",
    "full_significant_same_direction", "full_significant_retained_at_fdr_0.05"
  ),
  value = c(
    excluded_sample, condition_counts[["control"]], condition_counts[["selected"]],
    nrow(comparison),
    cor(comparison$log2FoldChange_full[finite_effects],
        comparison$log2FoldChange_leave_one_out[finite_effects]),
    sum(comparison$same_direction %in% TRUE, na.rm = TRUE),
    sum(full_significant, na.rm = TRUE),
    sum(comparison$significant_leave_one_out %in% TRUE, na.rm = TRUE),
    sum(full_significant & comparison$same_direction %in% TRUE, na.rm = TRUE),
    sum(full_significant & comparison$significant_leave_one_out %in% TRUE, na.rm = TRUE)
  )
)
write.table(
  summary_table,
  file.path(output_dir, paste0("exclude_", excluded_sample, "_summary.tsv")),
  sep = "\t", quote = FALSE, row.names = FALSE
)

saveRDS(
  dds_sensitivity,
  file.path(output_dir, paste0("exclude_", excluded_sample, "_dds.rds"))
)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
