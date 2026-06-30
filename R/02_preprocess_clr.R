#!/usr/bin/env Rscript
# =============================================================================
# 02_preprocess_clr.R
# -----------------------------------------------------------------------------
# Compositionally-correct preprocessing of MetaPhlAn species matrices.
# Implements REPORT.md sec 6b: prevalence filter -> zero handling -> CLR.
#
# WHY CLR (centered log-ratio), not TSS or rarefaction:
#   * Microbiome relative abundances are COMPOSITIONAL (sum-constrained);
#     correlations on raw proportions are spurious (Gloor 2017, doi:10/gcz9sj).
#   * Rarefaction throws away data and power -> "inadmissible" (McMurdie & Holmes
#     2014, doi:10.1371/journal.pcbi.1003531).
#   * CLR maps the simplex to real space: clr(x)_i = ln( x_i / g(x) ),
#     g(x) = geometric mean of the sample. Suitable for ENet and tree models.
#
# Inputs : data/processed/<cohort>_species.csv  (from 01_*)
# Outputs: data/processed/<cohort>_species_clr.csv
#
# Requirements: BiocManager::install("zCompositions"); install.packages("compositions")
# =============================================================================

suppressPackageStartupMessages({
  library(zCompositions)   # Bayesian-multiplicative zero replacement (cmultRepl)
})

OUT_DIR        <- file.path("data", "processed")
COHORTS        <- c("QinJ_2012", "KarlssonFH_2013")
PREVALENCE_MIN <- 0.10    # keep species present in >= 10% of samples
DETECTION_MIN  <- 1e-4    # relative-abundance threshold counting as "present"

# ---- CLR with zero replacement ----------------------------------------------
clr_transform <- function(mat) {
  # mat: samples x species, relative abundances (rows need not sum to exactly 1)
  # 1) prevalence filter
  present  <- colMeans(mat > DETECTION_MIN) >= PREVALENCE_MIN
  mat      <- mat[, present, drop = FALSE]
  message(sprintf("  prevalence filter (>=%.0f%%): %d species retained",
                  100 * PREVALENCE_MIN, ncol(mat)))

  # 2) zero replacement (CLR is undefined at zero). Bayesian-multiplicative
  #    replacement preserves the relative structure better than a flat pseudocount.
  mat_nozero <- tryCatch(
    cmultRepl(mat, method = "GBM", output = "prop", suppress.print = TRUE),
    error = function(e) {
      message("  cmultRepl failed (", conditionMessage(e),
              "); falling back to pseudocount.")
      pc <- min(mat[mat > 0]) / 2
      sweep(mat + pc, 1, rowSums(mat + pc), "/")
    }
  )

  # 3) CLR
  gm  <- exp(rowMeans(log(mat_nozero)))         # per-sample geometric mean
  clr <- log(sweep(mat_nozero, 1, gm, "/"))
  clr
}

# ---- main -------------------------------------------------------------------
for (study in COHORTS) {
  f <- file.path(OUT_DIR, paste0(study, "_species.csv"))
  if (!file.exists(f)) { message("skip (missing): ", f); next }

  mat <- as.matrix(read.csv(f, row.names = 1, check.names = FALSE))
  message("CLR for ", study, " (", nrow(mat), " samples)")
  clr <- clr_transform(mat)

  out <- file.path(OUT_DIR, paste0(study, "_species_clr.csv"))
  write.csv(clr, out)
  message("  -> wrote ", out, "\n")
}

# NOTE on harmonization for cross-cohort transfer (REPORT.md sec 6f):
#   When training on QinJ_2012 and testing on KarlssonFH_2013, align feature
#   sets to the INTERSECTION of species (Python side handles this), and consider
#   batch correction (MMUPHin / ConQuR) before pooling cohorts.
message("Next: python python/03_train_rf_shap.py")
