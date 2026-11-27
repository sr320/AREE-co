usage <- paste(
  "usage: Rscript run_proteomics_limma.R ABUNDANCE.tsv DESIGN.csv OUTPUT_DIR",
  "[--reference=control] [--test=heat] [--min-samples=3] [--log2-transformed=no] [--covariate=COL[+COL...]]"
)
args <- commandArgs(trailingOnly = TRUE)
is_opt <- grepl("^--", args)
positional <- args[!is_opt]
options <- list(
  reference = "control",
  test = "heat",
  min_samples = "3",
  log2_transformed = "no",
  covariate = NA
)
for (arg in args[is_opt]) {
  parts <- regmatches(arg, regexec("^--([a-z0-9-]+)=(.*)$", arg))[[1]]
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

abund_path <- positional[[1]]
design_path <- positional[[2]]
output_dir <- positional[[3]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

reference <- options$reference
test <- options$test
min_samples <- as.integer(options$min_samples)
already_log2 <- identical(tolower(options$log2_transformed), "yes")
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

df <- read.delim(abund_path, stringsAsFactors = FALSE, check.names = FALSE)
if (!"feature_id" %in% colnames(df)) {
  stop("abundance table must contain feature_id column")
}
if (!all(samples$sample %in% colnames(df))) {
  stop("abundance table missing sample columns from design")
}

mat <- as.matrix(df[, samples$sample, drop = FALSE])
mode(mat) <- "numeric"
keep <- rowSums(!is.na(mat) & mat > 0) >= min_samples
df <- df[keep, , drop = FALSE]
mat <- mat[keep, , drop = FALSE]

if (!already_log2) {
  mat <- log2(pmax(mat, 0) + 1)
}
# Sample-wise median centering preserves fold-change direction while removing loading bias
sample_medians <- apply(mat, 2, median, na.rm = TRUE)
global_median <- median(sample_medians, na.rm = TRUE)
mat_norm <- sweep(mat, 2, sample_medians - global_median, "-")

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

fit <- lmFit(mat_norm, design_mat)
fit <- eBayes(fit, robust = TRUE)
tt <- topTable(fit, coef = coef_name, number = Inf, sort.by = "none")
se_vec <- sqrt(fit$s2.post) * fit$stdev.unscaled[, coef_name]

out <- data.frame(
  feature_id = df$feature_id,
  effect_size = tt$logFC,
  standard_error = se_vec,
  p_value = tt$P.Value,
  adjusted_p_value = tt$adj.P.Val,
  gene_symbol = if ("gene_symbol" %in% colnames(df)) df$gene_symbol else df$feature_id,
  description = if ("description" %in% colnames(df)) df$description else paste("Protein", df$feature_id),
  feature_id_standardized = if ("feature_id_standardized" %in% colnames(df)) df$feature_id_standardized else df$feature_id,
  stringsAsFactors = FALSE
)
out <- out[order(out$p_value), ]

res_file <- file.path(output_dir, paste0(test, "_vs_", reference, "_differential_proteins.tsv"))
write.table(out, res_file, sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(
  sprintf("limma_%s_median_centered_eBayes_robust", as.character(packageVersion("limma"))),
  file.path(output_dir, "analysis_method.txt")
)
