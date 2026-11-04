# =============================================================================
# Script 03: Filter, transform and count PSM/protein data
# Description: Applies purity / signal-to-noise / unambiguity filters and log2
#              transformation to the imported PSM and protein data frames, then
#              counts identified and quantified rows per dataset and replicate.
# Depends on: 02_constants_functions_and_data_import.R
# =============================================================================

# --- Filter PSM data frames --------------------------------------------------

filter_psm <- list(

  list(name = "long_psm_baseline", df = long_psm_baseline, purity = 0, SN = 0, log2 = FALSE),

  list(name = "long_psm_standard_StoN_grouping_20", df = long_psm_standard_StoN_grouping_20, purity = 0, SN = 0, log2 = FALSE)
)

psm_results <- lapply(filter_psm, function(r) {
  filter_psm_long_dataframe(
    df = r$df,
    purity_filter = r$purity,
    chimerys_coeff_filter = if (!is.null(r$chimerys_coeff_filter)) r$chimerys_coeff_filter else 0,
    unambiguous_filter = TRUE,
    SN_filter = r$SN,
    log2_transform = r$log2,
    resolution_filter = if (!is.null(r$resolution_filter)) r$resolution_filter else NULL
  )
})

names(psm_results) <- sapply(filter_psm, `[[`, "name")

rm(filter_psm)

# --- Filter and log2-transform protein data frames ---------------------------

# Each entry is filtered and log2-transformed; `quan = TRUE` keeps only rows
# with quantitative values, `quan = FALSE` keeps all identified rows.
filter_prot <- list(

  # baseline
  list(name = "long_prot_baseline_log2_quan_values", df = long_prot_baseline, quan = TRUE,  log2 = TRUE),
  list(name = "long_prot_baseline_log2", df = long_prot_baseline, quan = F,  log2 = TRUE),

  # standard S/N grouping with empty / occupied
  list(name = "long_prot_standard_StoN_grouping_20_log2_quan_values", df = long_prot_standard_StoN_grouping_20, quan = TRUE,  log2 = TRUE),
  list(name = "long_prot_standard_StoN_grouping_20_log2", df = long_prot_standard_StoN_grouping_20, quan = F,  log2 = TRUE)

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

# --- Count rows per dataset and summarise over replicates --------------------

# Count rows for all protein datasets.
prot_row_counts <- lapply(prot_results, count_rows_per_dataset)
names(prot_row_counts) <- names(prot_results)

# Summarise replicates for all protein datasets.
prot_row_counts_avg <- lapply(prot_row_counts, summarise_replicates)
names(prot_row_counts_avg) <- paste0(names(prot_results), "_avg")

# Count rows for all PSM datasets.
psm_row_counts <- lapply(psm_results, count_rows_per_dataset)
names(psm_row_counts) <- names(psm_results)

# Summarise replicates for all PSM datasets.
psm_row_counts_avg <- lapply(psm_row_counts, summarise_replicates)
names(psm_row_counts_avg) <- paste0(names(psm_results), "_avg")

# --- Order the label factor by m/z and injection time ------------------------

reorder_labels <- function(df) {
  df <- df %>%
    mutate(
      mz = str_extract(label, "^m/z [0-9\\.]+"),
      time = str_extract(label, "[0-9]+ms$"),
      time_num = as.numeric(str_remove(time, "ms"))
    ) %>%
    arrange(mz, time_num) %>%
    mutate(label = factor(label, levels = unique(label))) %>%
    select(-mz, -time, -time_num)

  return(df)
}

psm_row_counts_avg <- lapply(
  psm_row_counts_avg,
  reorder_labels
)

prot_row_counts_avg <- lapply(
  prot_row_counts_avg,
  reorder_labels
)

