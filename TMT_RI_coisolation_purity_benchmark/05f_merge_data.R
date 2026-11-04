# ===========================================================================
# Script:  05f_merge_data.R
# Purpose: Merge protein ratio and CV data into a single long table for the
#          faceted plots.
# Inputs:  prot_results (03a), shared functions.
# Outputs: all_ratios_and_cv (feeds 05g and 05h) and
#          ratio_counts_mean_by_label.
# Depends: Run after 03a.
# ===========================================================================

# Merge the data from the different sources into a single long data frame
# This enables faceted CV% and ratio boxplots that show the influence of the filter strategies

# First, each source is filtered to the relevant columns


# columns to select
select_cols <- c("unnorm_H100Y000vsH100Y010",
                 "unnorm_H100Y000vsH100Y020",
                 "unnorm_H100Y000vsH100Y100",
                 "unnorm_H100Y100vsH100Y010",
                 "unnorm_H100Y100vsH100Y020",
                 "unnorm_H100Y020vsH100Y010",
                 "cv_H100Y100",
                 "cv_H100Y020",
                 "cv_H100Y010",
                 "label", "Organism", "dataset")

# Define the mapping between variable names and prot_results elements
data_mapping <- list(
  baseline = "long_prot_baseline_log2",
  baseline_70perc = "long_prot_baseline_log2",
  
  standard_StoN_grouping_20 = "long_prot_standard_StoN_grouping_20_log2",
  StoN_grouping_3 = "long_prot_StoN_grouping_3_log2",
  StoN_grouping_6 = "long_prot_StoN_grouping_6_log2",
  StoN_grouping_9 = "long_prot_StoN_grouping_9_log2",
  StoN_grouping_12 = "long_prot_StoN_grouping_12_log2",
  StoN_grouping_15 = "long_prot_StoN_grouping_15_log2",
  StoN_grouping_20 = "long_prot_StoN_grouping_20_log2",
  StoN_grouping_30 = "long_prot_StoN_grouping_30_log2",
  StoN_grouping_40 = "long_prot_StoN_grouping_40_log2",
  StoN_grouping_80 = "long_prot_StoN_grouping_80_log2",
  
  DIANN_baseline = "long_prot_DIANN_baseline_log2",
  SNaut_baseline = "long_prot_SNaut_baseline_log2",
  DIANN_StoN_grouping_20 = "long_prot_DIANN_StoN_grouping_20_log2",
  SNaut_StoN_grouping_20 = "long_prot_SNaut_StoN_grouping_20_log2"
  
  )


# Apply the selection to all datasets using lapply
data_list <- lapply(data_mapping, function(col_name) {
  prot_results[[col_name]] %>% select(all_of(select_cols))
})

# then these dataframes are merged and a id column is added based on the name of the dataframes
all_ratios_and_cv <- bind_rows(data_list, .id = "id")


# ---- Count the #proteins in each comparison --------------------------------

# Step 1: Count non-NA values per id, organism, and dataset
calculate_ratio_counts <- function(data) {
  # Get all unnorm_ columns
  unnorm_cols <- grep("^unnorm_", names(data), value = TRUE)

  # Count non-NA values for each unnorm_ column
  counts <- data %>%
    group_by(id, Organism, dataset, label) %>%  # Include label for step 2
    summarise(
      across(
        all_of(unnorm_cols),
        ~sum(!is.na(.x)),
        .names = "{.col}_count"
      ),
      .groups = "drop"
    )

  return(counts)
}

# Step 2: Calculate mean counts across datasets with the same label
calculate_mean_counts_by_label <- function(counts_data) {
  # Get all count columns
  count_cols <- grep("_count$", names(counts_data), value = TRUE)

  # Calculate mean across datasets for each label
  mean_counts <- counts_data %>%
    group_by(id, Organism, label) %>%
    summarise(
      across(
        all_of(count_cols),
        ~round(mean(.x, na.rm = TRUE),1),
        .names = "{.col}_mean"
      ),
      n_datasets = n(),  # Number of datasets averaged
      .groups = "drop"
    )

  return(mean_counts)
}

# Apply both steps
ratio_counts <- calculate_ratio_counts(all_ratios_and_cv)
ratio_counts_mean_by_label <- calculate_mean_counts_by_label(ratio_counts)
