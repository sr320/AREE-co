#!/usr/bin/env Rscript
# Cross-check AREE Python random-effects pooling against R metafor::rma.uni and metafor::rma.mv.
# Usage: Rscript crosscheck_meta_analysis.R INPUT_TSV OUTPUT_TSV

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("usage: Rscript crosscheck_meta_analysis.R INPUT_TSV OUTPUT_TSV")
}

suppressPackageStartupMessages({
  library(metafor)
})

input_path <- args[[1]]
output_path <- args[[2]]

df <- read.delim(input_path, stringsAsFactors = FALSE)
case_ids <- unique(df$case_id)
rows <- vector("list", length(case_ids))

for (idx in seq_along(case_ids)) {
  cid <- case_ids[[idx]]
  sub <- df[df$case_id == cid, ]
  method <- sub$method[[1]]
  rho <- sub$rho[[1]]
  yi <- as.numeric(sub$effect_size)
  sei <- as.numeric(sub$standard_error)
  grp <- as.character(sub$dependence_group)
  has_shared <- any(duplicated(grp)) && !is.na(rho)

  if (has_shared) {
    # Exact sampling covariance matrix V(rho) + between-group block reduction
    V <- diag(sei^2, nrow = length(sei))
    for (g in unique(grp)) {
      g_idx <- which(grp == g)
      if (length(g_idx) > 1) {
        for (a in g_idx) {
          for (b in g_idx) {
            if (a != b) V[a, b] <- rho * sei[a] * sei[b]
          }
        }
      }
    }
    yb <- numeric(0)
    vb <- numeric(0)
    for (g in sort(unique(grp))) {
      g_idx <- which(grp == g)
      Vg <- V[g_idx, g_idx, drop = FALSE]
      Vg_inv <- solve(Vg)
      ones <- rep(1, length(g_idx))
      wg <- as.numeric(t(ones) %*% Vg_inv %*% ones)
      yb <- c(yb, as.numeric(t(ones) %*% Vg_inv %*% yi[g_idx]) / wg)
      vb <- c(vb, 1 / wg)
    }
    fit_uni <- rma.uni(yb, vb, method = toupper(method), control = list(threshold = 1e-10, maxiter = 1000))
    if (toupper(method) == "REML") {
      fit_mv <- rma.mv(yi, V, random = ~ 1 | grp, method = "REML")
      mv_beta <- as.numeric(fit_mv$beta[1])
      mv_se <- as.numeric(fit_mv$se[1])
      mv_tau2 <- as.numeric(fit_mv$sigma2[1])
      mv_pval <- as.numeric(fit_mv$pval[1])
    } else {
      mv_beta <- NA_real_
      mv_se <- NA_real_
      mv_tau2 <- NA_real_
      mv_pval <- NA_real_
    }
    rows[[idx]] <- data.frame(
      case_id = cid,
      r_pooled_effect = as.numeric(fit_uni$beta[1]),
      r_pooled_se = as.numeric(fit_uni$se[1]),
      r_p_value = as.numeric(fit_uni$pval[1]),
      r_q = as.numeric(fit_uni$QE),
      r_i2_percent = as.numeric(fit_uni$I2),
      r_tau2 = as.numeric(fit_uni$tau2),
      r_mv_pooled_effect = mv_beta,
      r_mv_pooled_se = mv_se,
      r_mv_tau2 = mv_tau2,
      r_mv_p_value = mv_pval,
      stringsAsFactors = FALSE
    )
  } else {
    fit <- rma.uni(yi, sei^2, method = toupper(method), control = list(threshold = 1e-10, maxiter = 1000))
    rows[[idx]] <- data.frame(
      case_id = cid,
      r_pooled_effect = as.numeric(fit$beta[1]),
      r_pooled_se = as.numeric(fit$se[1]),
      r_p_value = as.numeric(fit$pval[1]),
      r_q = as.numeric(fit$QE),
      r_i2_percent = as.numeric(fit$I2),
      r_tau2 = as.numeric(fit$tau2),
      r_mv_pooled_effect = NA_real_,
      r_mv_pooled_se = NA_real_,
      r_mv_tau2 = NA_real_,
      r_mv_p_value = NA_real_,
      stringsAsFactors = FALSE
    )
  }
}

out <- do.call(rbind, rows)
write.table(out, output_path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
cat("metafor cross-check completed for", nrow(out), "cases\n")
