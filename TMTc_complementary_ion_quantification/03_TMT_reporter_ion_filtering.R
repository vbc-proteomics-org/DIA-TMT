# ===========================================================================
# Script:  03_TMT_reporter_ion_filtering.R
# Purpose: Filter and log2-transform the long TMT reporter-ion protein tables
#          for the baseline and standard data sets.
# Inputs:  long_prot_baseline, long_prot_standard and
#          filter_and_log2_transform_long_prot (from script 02).
# Outputs: prot_results (named list of filtered protein data frames).
# Depends: 02_TMT_reporter_ion_import.R
# ===========================================================================

# --- Filter and log2-transform the protein data frames ---
filter_prot <- list(

  # baseline
  list(name = "long_prot_baseline_log2_quan_values", df = long_prot_baseline, quan = TRUE,  log2 = TRUE),
  list(name = "long_prot_baseline_log2", df = long_prot_baseline, quan = F,  log2 = TRUE),

  # standard
  list(name = "long_prot_standard_log2_quan_values", df = long_prot_standard, quan = TRUE,  log2 = TRUE),
  list(name = "long_prot_standard_log2", df = long_prot_standard, quan = F,  log2 = TRUE)

)

prot_results <- lapply(filter_prot, function(r) {
  filter_and_log2_transform_long_prot(
    df = r$df,
    quan_values_filter = r$quan,
    unambiguous_filter = TRUE,
    log2_transform = r$log2
  )
})

names(prot_results) <- sapply(filter_prot, `[[`, "name")

rm(filter_prot)
