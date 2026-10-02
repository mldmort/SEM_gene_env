#!/usr/bin/env Rscript
# Structural equation model of the parental age -> ASD association (Figure 2 model).
#
#   parent_age ~ rp*rv + ep*edu + pp*pgs            (genetic confounding of parental age)
#   edu        ~ re*rv + pe*pgs
#   dnm        ~ g*parent_age                       (gamma: age -> damaging DNMs)
#   asd        ~ c*parent_age + d*dnm + ra*rv + ea*edu + pa*pgs   (asd binary, probit)
#
# Covariates in EVERY equation: sex, PC1..10, five sequencing covariates (child coverage,
# mother/father coverage, mother/father CRAM mapped %) and, for SPARK WES, genotype-source
# indicators. Birth order (capped at 4; bo_* columns, linear or categorical) enters the parental-age equation (parental age at a
# child's birth = age at first birth + spacing) and the ASD equation (sibling sampling), not DNM.
# One model per parent (father_age, mother_age), fitted per cohort with WLSMV
# (asd declared ordered) and family-clustered standard errors.
#
# Controls-only check: the same model without the asd equation, fitted in
# unaffected children only, to see whether conditioning on recruitment through an
# ASD child shifts the upstream paths (gamma, rp, ep, pp, re, pe).
#
# Inputs : results/sem/input_<cohort>.tsv           (01_prepare_sem_data.py)
# Outputs: results/sem/sem_estimates.tsv            every parameter, all fits
#          results/sem/sem_fit.tsv                  n, families, fit indices per fit
#          results/sem/sem_meta.tsv                 fixed-effect IVW meta across cohorts
#
# Every model is fitted once per inherited RV definition (rv_def column in all outputs).
#
# Usage: Rscript 02_sem_lavaan.R [repo_root]

suppressMessages({
  library(lavaan)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
repo <- if (length(args) >= 1) args[1] else Sys.getenv("GXE_REPO", ".")
res_dir <- file.path(repo, "results", "sem")
cohorts <- c("spark_wgs", "ssc", "spark_wes")
parents <- c("father_age", "mother_age")
n_pcs <- 10
pcs <- paste0("PC", seq_len(n_pcs))
seq_covs <- c("coverage", "mother_coverage", "father_coverage", "mother_mapped_pct", "father_mapped_pct")
base_covs <- c("sex_male", pcs, seq_covs)
# Per cohort, genotype-source indicators (src_*, SPARK WES only) join base_covs in every
# equation; birth_order enters the parental-age and ASD equations.

model_string <- function(parent, with_asd, covs, bo) {
  cov_main <- paste(covs, collapse = " + ")
  bo_terms <- paste(bo, collapse = " + ")
  m <- c(
    sprintf("%s ~ rp*rv + ep*edu + pp*pgs + %s + %s", parent, bo_terms, cov_main),
    sprintf("edu ~ re*rv + pe*pgs + %s", cov_main),
    sprintf("dnm ~ g*%s + %s", parent, cov_main)
  )
  if (with_asd) {
    m <- c(m,
      sprintf("asd ~ c*%s + d*dnm + ra*rv + ea*edu + pa*pgs + %s + %s", parent, bo_terms, cov_main),
      # Indirect causal pathway: parental age -> DNM -> ASD
      "ind_dnm := g*d",
      # Non-causal (confounding) pathways linking parental age and ASD through shared genetics
      "conf_rv  := rp*ra",
      "conf_edu := ep*ea",
      "conf_pgs := pp*pa",
      "conf_rv_edu  := rp*re*ea + ep*re*ra",
      "conf_pgs_edu := pp*pe*ea + ep*pe*pa",
      "conf_total := rp*ra + ep*ea + pp*pa + rp*re*ea + ep*re*ra + pp*pe*ea + ep*pe*pa",
      "total_causal := c + g*d",
      "prop_mediated_dnm := g*d / (c + g*d)"
    )
  }
  paste(m, collapse = "\n")
}

fit_one <- function(dat, parent, with_asd) {
  covs <- c(base_covs, grep("^src_", names(dat), value = TRUE))
  bo <- grep("^bo_", names(dat), value = TRUE)   # bo_linear, or bo_2 / bo_3 / bo_4plus
  vars <- c("FID", "asd", parent, "dnm", "rv", "edu", "pgs", covs, bo)
  d <- dat[complete.cases(dat[, ..vars])]
  if (!with_asd) d <- d[asd == 0]
  args <- list(model = model_string(parent, with_asd, covs, bo), data = as.data.frame(d),
               estimator = if (with_asd) "WLSMV" else "MLR")
  # Cluster on family unless every family contributes one child (e.g. SSC controls:
  # one designated sibling per family), where clustering is a no-op and lavaan errors.
  if (max(table(d$FID)) > 1) args$cluster <- "FID"
  if (with_asd) args$ordered <- "asd"
  fit <- do.call(sem, args)
  list(fit = fit, n = nrow(d), families = uniqueN(d$FID),
       cases = sum(d$asd == 1), controls = sum(d$asd == 0))
}

est_rows <- list(); fit_rows <- list()
for (cohort in cohorts) {
  dat_all <- fread(file.path(res_dir, sprintf("input_%s.tsv", cohort)), na.strings = "NA")
  # One SEM per inherited rare-variant definition (input columns rv_<name>, config sem.rv_definitions).
  rv_defs <- sub("^rv_", "", grep("^rv_", names(dat_all), value = TRUE))
  for (rv_def in rv_defs) {
  dat <- copy(dat_all)
  dat[, rv := get(paste0("rv_", rv_def))]
  for (parent in parents) {
    for (sample in c("all", "controls")) {
      with_asd <- sample == "all"
      cat(sprintf("== %s | rv=%s | %s | %s\n", cohort, rv_def, parent, sample))
      r <- tryCatch(fit_one(dat, parent, with_asd), error = function(e) {
        cat("   FAILED:", conditionMessage(e), "\n"); NULL })
      if (is.null(r)) next
      pe <- as.data.table(parameterEstimates(r$fit, ci = TRUE))
      st <- as.data.table(standardizedSolution(r$fit))[, .(lhs, op, rhs, std.all = est.std)]
      pe <- merge(pe, st, by = c("lhs", "op", "rhs"), all.x = TRUE, sort = FALSE)
      # Keep every regression coefficient (covariates included) and the defined effects;
      # unlabelled covariate paths get "<lhs>~<rhs>" as their label.
      pe <- pe[op %in% c("~", ":=")]
      pe[label == "", label := paste0(lhs, "~", rhs)]
      pe[, `:=`(cohort = cohort, rv_def = rv_def, parent = parent, sample = sample)]
      est_rows[[length(est_rows) + 1]] <- pe
      fm <- fitMeasures(r$fit, c("chisq.scaled", "df", "pvalue.scaled", "cfi.scaled",
                                 "rmsea.scaled", "srmr"))
      fm <- fm[!is.na(names(fm))]
      fit_rows[[length(fit_rows) + 1]] <- data.table(
        cohort = cohort, rv_def = rv_def, parent = parent, sample = sample, n = r$n, families = r$families,
        cases = r$cases, controls = r$controls, t(fm), converged = lavInspect(r$fit, "converged"))
    }
  }
  }
}

est <- rbindlist(est_rows, fill = TRUE)
setcolorder(est, c("cohort", "rv_def", "parent", "sample", "label", "lhs", "op", "rhs"))
fwrite(est, file.path(res_dir, "sem_estimates.tsv"), sep = "\t")
fwrite(rbindlist(fit_rows, fill = TRUE), file.path(res_dir, "sem_fit.tsv"), sep = "\t")

# Fixed-effect inverse-variance meta-analysis across cohorts, per parameter.
meta <- est[is.finite(se) & se > 0, {
  w <- 1 / se^2
  b <- sum(w * est) / sum(w); s <- sqrt(1 / sum(w))
  q <- sum(w * (est - b)^2); k <- .N
  .(est = b, se = s, z = b / s, pvalue = 2 * pnorm(-abs(b / s)),
    ci.lower = b - 1.96 * s, ci.upper = b + 1.96 * s,
    k = k, cohorts = paste(cohort, collapse = ","),
    het_q = q, het_p = if (k > 1) pchisq(q, k - 1, lower.tail = FALSE) else NA_real_)
}, by = .(rv_def, parent, sample, label)]
fwrite(meta, file.path(res_dir, "sem_meta.tsv"), sep = "\t")
cat("\nWrote sem_estimates.tsv, sem_fit.tsv, sem_meta.tsv to", res_dir, "\n")
