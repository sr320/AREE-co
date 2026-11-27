usage <- paste(
  "usage: Rscript run_metabolomics_diff.R ABUNDANCE.tsv DESIGN.csv OUTPUT_DIR",
  "[--reference=control] [--test=heat] [--min-samples=3] [--covariate=COL[+COL...]]"
)
args <- commandArgs(trailingOnly = TRUE)
is_opt <- grepl("^--", args)
positional <- args[!is_opt]
options <- list(
  reference = "control",
  test = "heat",
  min_samples = "3",
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
mat <- pmax(mat[keep, , drop = FALSE], 0)

# Probabilistic Quotient Normalization (PQN) against median spectrum
ref_spec <- apply(mat, 1, median, na.rm = TRUE)
ref_spec[ref_spec <= 0] <- 1
quotients <- sweep(mat, 1, ref_spec, "/")
pqn_factors <- apply(quotients, 2, median, na.rm = TRUE)
pqn_factors[is.na(pqn_factors) | pqn_factors <= 0] <- 1
mat_pqn <- sweep(mat, 2, pqn_factors, "/")
mat_log2 <- log2(mat_pqn + 1)

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

fit <- lmFit(mat_log2, design_mat)
fit <- eBayes(fit, robust = TRUE)
tt <- topTable(fit, coef = coef_name, number = Inf, sort.by = "none")
se_vec <- sqrt(fit$s2.post) * fit$stdev.unscaled[, coef_name]

ctrl_idx <- samples$condition == reference
test_idx <- samples$condition == test
wilcox_p <- apply(mat_log2, 1, function(vals) {
  suppressWarnings(wilcox.test(vals[test_idx], vals[ctrl_idx], exact = FALSE)$p.value)
})

out <- data.frame(
  feature_id = df$feature_id,
  effect_size = tt$logFC,
  standard_error = se_vec,
  p_value = tt$P.Value,
  adjusted_p_value = tt$adj.P.Val,
  wilcox_p_value = wilcox_p,
  gene_symbol = if ("gene_symbol" %in% colnames(df)) df$gene_symbol else df$feature_id,
  description = if ("description" %in% colnames(df)) df$description else paste("Metabolite", df$feature_id),
  feature_id_standardized = if ("feature_id_standardized" %in% colnames(df)) df$feature_id_standardized else df$feature_id,
  stringsAsFactors = FALSE
)
out <- out[order(out$p_value), ]

res_file <- file.path(output_dir, paste0(test, "_vs_", reference, "_differential_metabolites.tsv"))
write.table(out, res_file, sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(
  sprintf("PQN_log2_limma_%s_eBayes_robust", as.character(packageVersion("limma"))),
  file.path(output_dir, "analysis_method.txt")
)
