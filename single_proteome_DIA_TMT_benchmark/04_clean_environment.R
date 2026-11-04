# =============================================================================
# Script 04: Clean up the environment
# Description: Removes the large per-group PSM/protein data frame lists and
#              import helpers once the downstream long_* and *_results objects
#              have been created, then triggers garbage collection.
# Depends on: 03_filter_transform_and_counts.R
# =============================================================================

# --- Define groups to remove -------------------------------------------------

# Group suffixes must match those used in the per-group variable names.

all_groups <- c(
  "baseline", 
  "standard_StoN_grouping_20"
)

# --- Remove per-group data frame lists ---------------------------------------

rm(list = c(
  paste0("psm_data_frames_",     all_groups),
  paste0("protein_data_frames_", all_groups)
))

rm(list = c(
  paste0("data_list_", all_groups)
))

gc()
