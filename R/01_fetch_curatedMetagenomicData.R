y#!/usr/bin/env Rscript
# =============================================================================
# 01_fetch_curatedMetagenomicData.R
# -----------------------------------------------------------------------------
# Pull harmonized shotgun-metagenomic profiles for T2DM cohorts from
# curatedMetagenomicData (Pasolli et al., Nat Methods 2017; doi:10.1038/nmeth.4468)
# and export tidy matrices for the downstream Python ML pipeline.
#
#   Discovery cohort        : QinJ_2012        (T2DM, China; Qin et al. Nature 2012)
#   External validation     : KarlssonFH_2013  (T2DM/IGT, Sweden women; Nature 2013)
#
# Outputs (data/processed/):
#   <cohort>_species.csv     samples x species (MetaPhlAn relative abundance, %)
#   <cohort>_pathways.csv    samples x pathways (HUMAnN, optional)
#   <cohort>_metadata.csv    sample metadata incl. study_condition, metformin*
#
# Requirements: R >= 4.3, Bioconductor
#   BiocManager::install(c("curatedMetagenomicData", "mia", "SummarizedExperiment"))
# =============================================================================

suppressPackageStartupMessages({
  library(curatedMetagenomicData)
  library(SummarizedExperiment)
  library(mia)            # TreeSummarizedExperiment helpers
})

# ---- config -----------------------------------------------------------------
OUT_DIR  <- file.path("data", "processed")
COHORTS  <- c("QinJ_2012", "KarlssonFH_2013")
DATATYPE <- "relative_abundance"   # MetaPhlAn species-level relative abundance
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ---- helper: fetch one cohort as a (Tree)SummarizedExperiment ---------------
fetch_cohort <- function(study, datatype = DATATYPE) {
  message("Fetching: ", study, " [", datatype, "]")
  # curatedMetagenomicData() matches by a regex on study_name + dataType.
  # dryrun = FALSE actually downloads from ExperimentHub (cached locally).
  pattern <- paste0(study, ".+", datatype)
  se_list <- curatedMetagenomicData(pattern, dryrun = FALSE, counts = FALSE)
  if (length(se_list) == 0) stop("No resource matched: ", pattern)
  se_list[[1]]
}

# ---- helper: write species matrix + metadata --------------------------------
export_cohort <- function(se, study) {
  # assay: features (species) x samples  ->  transpose to samples x features
  species <- t(as.matrix(assay(se)))
  meta    <- as.data.frame(colData(se))

  # Keep only columns we actually use downstream; tolerate missing ones.
  keep <- intersect(
    c("study_name", "subject_id", "study_condition", "disease",
      "age", "gender", "BMI", "country",
      # medication metadata is sparse but ESSENTIAL for the metformin analysis:
      "antibiotics_current_use", "treatment", "disease_subtype"),
    colnames(meta)
  )
  meta <- meta[, keep, drop = FALSE]

  # Best-effort metformin flag from whatever metadata the cohort exposes.
  # NOTE: curatedMetagenomicData metformin annotation is incomplete; verify per
  # cohort and, if absent, treat metformin status as a known unmeasured confounder
  # (see REPORT.md sec 6g, pitfall #1).
  meta$metformin <- NA
  if ("treatment" %in% colnames(meta)) {
    meta$metformin <- grepl("metformin", tolower(meta$treatment))
  }

  write.csv(species, file.path(OUT_DIR, paste0(study, "_species.csv")))
  write.csv(meta,    file.path(OUT_DIR, paste0(study, "_metadata.csv")))
  message(sprintf("  -> %s: %d samples x %d species",
                  study, nrow(species), ncol(species)))
  invisible(NULL)
}

# ---- main -------------------------------------------------------------------
for (study in COHORTS) {
  se <- tryCatch(fetch_cohort(study), error = function(e) {
    message("  !! failed for ", study, ": ", conditionMessage(e)); NULL
  })
  if (!is.null(se)) export_cohort(se, study)
}

message("\nDone. Inspect available studies with: sampleMetadata |> dplyr::distinct(study_name)")
message("Next: source('R/02_preprocess_clr.R')")
