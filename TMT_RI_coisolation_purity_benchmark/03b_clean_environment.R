# ===========================================================================
# Script:  03b_clean_environment.R
# Purpose: Remove large intermediate objects from the environment to free
#          memory.
# Inputs:  Objects created in 02b and 03a.
# Outputs: Side effect only: a smaller global environment.
# Depends: Run after 03a.
# ===========================================================================

# clean up space after all relevant objects are created

# Remove all individual psm/protein data frame lists — they've been rbind-ed into long_ format
rm(list = c(
  paste0("psm_data_frames_",     all_groups),
  paste0("protein_data_frames_", all_groups)
))

# Also clean up helper objects
rm(psm_data_frames, protein_data_frames, psm_lists, protein_lists, all_lists)

# Remove the long_ flat dataframes — they're now inside psm_results / prot_results
rm(list = c(
  paste0("long_psm_", all_groups),
  paste0("long_prot_", all_groups)
))

gc()
