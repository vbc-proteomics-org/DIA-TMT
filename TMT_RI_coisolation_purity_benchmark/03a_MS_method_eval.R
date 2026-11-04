# ===========================================================================
# Script:  03a_MS_method_eval.R
# Purpose: Assemble the master prot_results and psm_results lists (one entry
#          per filtering scheme), build top20_prot and compute averaged row
#          counts.
# Inputs:  long_prot_* / long_psm_* objects and functions from 02b.
# Outputs: prot_results, psm_results, top20_prot, prot_row_counts_avg,
#          psm_row_counts_avg, reorder_labels.
# Depends: Run after 02b.
# ===========================================================================

rm(data_list)
# Apply the function to filter the PSM data frames

filter_psm <- list(

  # our dataset
  list(name = "long_psm_baseline", df = long_psm_baseline, purity = 0, SN = 0, log2 = FALSE),
  list(name = "long_psm_baseline_70perc", df = long_psm_baseline_70perc, purity = 0, SN = 0, log2 = FALSE),

  list(name = "long_psm_standard_StoN_grouping_20_nonlog", df = long_psm_standard_StoN_grouping_20, purity = 0, SN = 0, log2 = FALSE),

  # S/N titration with exemplary DDA & DIA data, always purity >= 70% and resolution >= 45k
  list(name = "long_psm_StoN_grouping_3_nonlog", df = long_psm_StoN_grouping_3, purity = 0, SN = 0, log2 = FALSE),
  list(name = "long_psm_StoN_grouping_6_nonlog", df = long_psm_StoN_grouping_6, purity = 0, SN = 0, log2 = FALSE),
  list(name = "long_psm_StoN_grouping_9_nonlog", df = long_psm_StoN_grouping_9, purity = 0, SN = 0, log2 = FALSE),
  list(name = "long_psm_StoN_grouping_12_nonlog", df = long_psm_StoN_grouping_12, purity = 0, SN = 0, log2 = FALSE),
  list(name = "long_psm_StoN_grouping_15_nonlog", df = long_psm_StoN_grouping_15, purity = 0, SN = 0, log2 = FALSE),
  list(name = "long_psm_StoN_grouping_20_nonlog", df = long_psm_StoN_grouping_20, purity = 0, SN = 0, log2 = FALSE),
  list(name = "long_psm_StoN_grouping_30_nonlog", df = long_psm_StoN_grouping_30, purity = 0, SN = 0, log2 = FALSE),
  list(name = "long_psm_StoN_grouping_40_nonlog", df = long_psm_StoN_grouping_40, purity = 0, SN = 0, log2 = FALSE),
  list(name = "long_psm_StoN_grouping_80_nonlog", df = long_psm_StoN_grouping_80, purity = 0, SN = 0, log2 = FALSE),

  # DIANN and SN data

  list(name = "long_psm_SNaut_baseline_nonlog", df = long_psm_SNaut_baseline, purity = 0, SN = 0, log2 = FALSE),
  list(name = "long_psm_DIANN_baseline_nonlog", df = long_psm_DIANN_baseline, purity = 0, SN = 0, log2 = FALSE),

  list(name = "long_psm_SNaut_StoN_grouping_20_nonlog", df = long_psm_SNaut_StoN_grouping_20, purity = 0, SN = 0, log2 = FALSE),
  list(name = "long_psm_DIANN_StoN_grouping_20_nonlog", df = long_psm_DIANN_StoN_grouping_20, purity = 0, SN = 0, log2 = FALSE)

 

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

# apply the function filter_protein_data_frames on protein_data_frames based on empty Abundance rows

# filter the datasets, log2 transform and exclude specific datasets

filter_prot <- list(

  # baseline
  list(name = "long_prot_baseline_log2_quan_values", df = long_prot_baseline, quan = TRUE,  log2 = TRUE),
  list(name = "long_prot_baseline_log2", df = long_prot_baseline, quan = F,  log2 = TRUE),

  # baseline 70% purity
  list(name = "long_prot_baseline_70perc_log2_quan_values", df = long_prot_baseline_70perc, quan = TRUE,  log2 = TRUE),
  list(name = "long_prot_baseline_70perc_log2", df = long_prot_baseline_70perc, quan = F,  log2 = TRUE),

  
  # standard S/N grouping (occupied & empty) + resolution + 70% purity
  list(name = "long_prot_standard_StoN_grouping_20_log2_quan_values", df = long_prot_standard_StoN_grouping_20, quan = TRUE,  log2 = TRUE),
  list(name = "long_prot_standard_StoN_grouping_20_log2", df = long_prot_standard_StoN_grouping_20, quan = F,  log2 = TRUE),

  # S/N filter titration for occupied and empty channels as groups + >= 70% purity + >= 45k resolution
  list(name = "long_prot_StoN_grouping_3_log2_quan_values", df = long_prot_StoN_grouping_3, quan = TRUE,  log2 = TRUE),
  list(name = "long_prot_StoN_grouping_3_log2", df = long_prot_StoN_grouping_3, quan = F,  log2 = TRUE),

  list(name = "long_prot_StoN_grouping_6_log2_quan_values", df = long_prot_StoN_grouping_6, quan = TRUE,  log2 = TRUE),
  list(name = "long_prot_StoN_grouping_6_log2", df = long_prot_StoN_grouping_6, quan = F,  log2 = TRUE),

  list(name = "long_prot_StoN_grouping_9_log2_quan_values", df = long_prot_StoN_grouping_9, quan = TRUE,  log2 = TRUE),
  list(name = "long_prot_StoN_grouping_9_log2", df = long_prot_StoN_grouping_9, quan = F,  log2 = TRUE),

  list(name = "long_prot_StoN_grouping_12_log2_quan_values", df = long_prot_StoN_grouping_12, quan = TRUE,  log2 = TRUE),
  list(name = "long_prot_StoN_grouping_12_log2", df = long_prot_StoN_grouping_12, quan = F,  log2 = TRUE),

  list(name = "long_prot_StoN_grouping_15_log2_quan_values", df = long_prot_StoN_grouping_15, quan = TRUE,  log2 = TRUE),
  list(name = "long_prot_StoN_grouping_15_log2", df = long_prot_StoN_grouping_15, quan = F,  log2 = TRUE),

  list(name = "long_prot_StoN_grouping_20_log2_quan_values", df = long_prot_StoN_grouping_20, quan = TRUE,  log2 = TRUE),
  list(name = "long_prot_StoN_grouping_20_log2", df = long_prot_StoN_grouping_20, quan = F,  log2 = TRUE),

  list(name = "long_prot_StoN_grouping_30_log2_quan_values", df = long_prot_StoN_grouping_30, quan = TRUE,  log2 = TRUE),
  list(name = "long_prot_StoN_grouping_30_log2", df = long_prot_StoN_grouping_30, quan = F,  log2 = TRUE),

  list(name = "long_prot_StoN_grouping_40_log2_quan_values", df = long_prot_StoN_grouping_40, quan = TRUE,  log2 = TRUE),
  list(name = "long_prot_StoN_grouping_40_log2", df = long_prot_StoN_grouping_40, quan = F,  log2 = TRUE),

  list(name = "long_prot_StoN_grouping_80_log2_quan_values", df = long_prot_StoN_grouping_80, quan = TRUE,  log2 = TRUE),
  list(name = "long_prot_StoN_grouping_80_log2", df = long_prot_StoN_grouping_80, quan = F,  log2 = TRUE),

  # Spectronaut baseline
  list(name = "long_prot_SNaut_baseline_log2_quan_values", df = long_prot_SNaut_baseline, quan = TRUE,  log2 = TRUE),
  list(name = "long_prot_SNaut_baseline_log2", df = long_prot_SNaut_baseline, quan = F,  log2 = TRUE),

  # DIANN baseline
  list(name = "long_prot_DIANN_baseline_log2_quan_values", df = long_prot_DIANN_baseline, quan = TRUE,  log2 = TRUE),
  list(name = "long_prot_DIANN_baseline_log2", df = long_prot_DIANN_baseline, quan = F,  log2 = TRUE),

  # Spectronaut StoN grouping
  list(name = "long_prot_SNaut_StoN_grouping_20_log2_quan_values", df = long_prot_SNaut_StoN_grouping_20, quan = TRUE,  log2 = TRUE),
  list(name = "long_prot_SNaut_StoN_grouping_20_log2", df = long_prot_SNaut_StoN_grouping_20, quan = F,  log2 = TRUE),

  # DIANN StoN grouping
  list(name = "long_prot_DIANN_StoN_grouping_20_log2_quan_values", df = long_prot_DIANN_StoN_grouping_20, quan = TRUE,  log2 = TRUE),
  list(name = "long_prot_DIANN_StoN_grouping_20_log2", df = long_prot_DIANN_StoN_grouping_20, quan = F,  log2 = TRUE)

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


# filter proteins with a threshold applied to each "label" and each maxExp column is processed individually 
# Receive one data frame per comparison
# threshold = 0.2 meaning that the top20% abundant proteins are kept for each organsims-maxExp combination

top20_prot <- filter_prot_abundance(prot_results$long_prot_standard_StoN_grouping_20_log2, 
                                     threshold = 0.2)




# ---- Get number of peptides, PSMs and Proteins for all data frames stored in the lists ----

# Count rows for all Protein datasets
prot_row_counts <- lapply(prot_results, count_rows_per_dataset)
names(prot_row_counts) <- names(prot_results)

# Summarise replicates for all Protein datasets
prot_row_counts_avg <- lapply(prot_row_counts, summarise_replicates)
names(prot_row_counts_avg) <- paste0(names(prot_results), "_avg")

# Count rows for all PSM datasets
psm_row_counts <- lapply(psm_results, count_rows_per_dataset)
names(psm_row_counts) <- names(psm_results)

# Summarise replicates for all PSM datasets
psm_row_counts_avg <- lapply(psm_row_counts, summarise_replicates)
names(psm_row_counts_avg) <- paste0(names(psm_results), "_avg")


# Reorder factor levels of the label column
reorder_labels <- function(df) {
  # Extract m/z and time parts
  df <- df %>%
    mutate(
      mz = str_extract(label, "^m/z [0-9\\.]+"),
      time = str_extract(label, "[0-9]+ms$"),
      # convert time to numeric for ordering
      time_num = as.numeric(str_remove(time, "ms"))
    ) %>%
    # Arrange by mz and time_num
    arrange(mz, time_num) %>%
    # Create a new ordered factor for label
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

