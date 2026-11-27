usage <- paste(
  "usage: Rscript run_methylation_limma.R COVERAGE.tsv DESIGN.csv OUTPUT_DIR",
  "[--reference=control] [--test=heat] [--feature-type=cpg_site]",
  "[--min-coverage=5] [--min-samples=3] [--covariate=COL[+COL...]]"
)
args <- commandArgs(trailingOnly = TRUE)
is_opt <- grepl("^--", args)
positional <- args[!is_opt]
options <- list(
  reference = "control",
  test = "heat",
  feature_type = "cpg_site",
  min_coverage = "5",
  min_samples = "3",
  covariate = NA
)
for (arg in args[is_opt]) {
  parts <- regmatches(arg, regexec("^--([a-z-]+)=(.*)$", arg))[[1]]
  key <- if (length(parts) == 3) gsub("-", "_", parts[[2]]) else ""
  if (!key %in% names(options)) {
    stop(paste("unknown option:", arg, "\n", usage))
  }
  options[[key]] <- parts[[3]]
}
if (length(positional) != 3) {
  stop(usage)
}

suppressPackageStartupMessages({
  library(limma)
})

coverage_path <- positional[[1]]
design_path <- positional[[2]]
output_dir <- positional[[3]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

reference <- options$reference
test <- options$test
min_cov <- as.numeric(options$min_coverage)
min_samples <- as.integer(options$min_samples)
covariates <- if (is.na(options$covariate) || !nzchar(options$covariate)) character() else
  strsplit(options$covariate, "+", fixed = TRUE)[[1]]

samples <- read.csv(design_path, stringsAsFactors = FALSE, check.names = FALSE)
if (!all(c("sample", "condition") %in% colnames(samples))) {
  stop("design sheet must contain sample and condition columns")
}
if (!all(c(reference, test) %in% samples$condition)) {
  stop("design sheet missing test or reference condition")
}
samples$condition <- factor(samples$condition, levels = c(reference, test))
rownames(samples) <- samples$sample

cov_df <- read.delim(coverage_path, stringsAsFactors = FALSE, check.names = FALSE)
if (!"feature_id" %in% colnames(cov_df)) {
  stop("coverage table must contain feature_id column")
}

m_cols <- paste0(samples$sample, "_M")
u_cols <- paste0(samples$sample, "_U")
beta_cols <- paste0(samples$sample, "_beta")

if (all(m_cols %in% colnames(cov_df)) && all(u_cols %in% colnames(cov_df))) {
  M_mat <- as.matrix(cov_df[, m_cols, drop = FALSE])
  U_mat <- as.matrix(cov_df[, u_cols, drop = FALSE])
  tot_mat <- M_mat + U_mat
  keep <- rowSums(tot_mat >= min_cov, na.rm = TRUE) >= min_samples
  cov_df <- cov_df[keep, , drop = FALSE]
  M_mat <- M_mat[keep, , drop = FALSE]
  U_mat <- U_mat[keep, , drop = FALSE]
  beta_mat <- (M_mat + 0.5) / (M_mat + U_mat + 1.0)
  colnames(beta_mat) <- samples$sample
} else {
  target_cols <- if (all(beta_cols %in% colnames(cov_df))) beta_cols else samples$sample
  if (!all(target_cols %in% colnames(cov_df))) {
    stop("coverage table missing sample methylation columns")
  }
  beta_mat <- as.matrix(cov_df[, target_cols, drop = FALSE])
  colnames(beta_mat) <- samples$sample
  keep <- rowSums(!is.na(beta_mat)) >= min_samples
  cov_df <- cov_df[keep, , drop = FALSE]
  beta_mat <- beta_mat[keep, , drop = FALSE]
}

beta_clamped <- pmin(pmax(beta_mat, 1e-4), 1 - 1e-4)
mval_mat <- log2(beta_clamped / (1 - beta_clamped))

model_terms <- c(covariates, "condition")
for (cv in covariates) {
  if (!cv %in% colnames(samples)) stop(paste("missing covariate:", cv))
  samples[[cv]] <- factor(samples[[cv]])
}
design_formula <- as.formula(paste("~", paste(model_terms, collapse = " + ")))
design_mat <- model.matrix(design_formula, data = samples)
if (qr(design_mat)$rank < ncol(design_mat)) {
  stop("design matrix is not full rank")
}
coef_name <- paste0("condition", test)

fit <- lmFit(mval_mat, design_mat)
fit <- eBayes(fit, robust = TRUE)
tt <- topTable(fit, coef = coef_name, number = Inf, sort.by = "none")
se_vec <- sqrt(fit$s2.post) * fit$stdev.unscaled[, coef_name]

ctrl_idx <- samples$condition == reference
test_idx <- samples$condition == test
delta_beta <- rowMeans(beta_mat[, test_idx, drop = FALSE], na.rm = TRUE) -
  rowMeans(beta_mat[, ctrl_idx, drop = FALSE], na.rm = TRUE)

out <- data.frame(
  feature_id = cov_df$feature_id,
  effect_size = tt$logFC,
  standard_error = se_vec,
  p_value = tt$P.Value,
  adjusted_p_value = tt$adj.P.Val,
  delta_beta = delta_beta,
  gene_symbol = if ("gene_symbol" %in% colnames(cov_df)) cov_df$gene_symbol else cov_df$feature_id,
  description = if ("description" %in% colnames(cov_df)) cov_df$description else paste("Methylation locus", cov_df$feature_id),
  feature_id_standardized = if ("feature_id_standardized" %in% colnames(cov_df)) cov_df$feature_id_standardized else cov_df$feature_id,
  stringsAsFactors = FALSE
)
out <- out[order(out$p_value), ]

res_file <- file.path(output_dir, paste0(test, "_vs_", reference, "_differential_methylation.tsv"))
write.table(out, res_file, sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(
  sprintf("limma_%s_eBayes_M_value_%s", as.character(packageVersion("limma")), options$feature_type),
  file.path(output_dir, "analysis_method.txt")
)
