# ===========================================================================
# Script:  03d_concat_prot_counts.R
# Purpose: Concatenate per-dataset protein counts and build the loss
#          annotations.
# Inputs:  prot_row_counts_avg (from 03a).
# Outputs: long_prot_counts, loss_annotations, long_prot_counts_all,
#          loss_annotations_all.
# Depends: Run after 03a; required by 04 and the 05* plotting scripts.
# ===========================================================================

# Identify all "keys" (dataset bases) that have both Identified and Quantified summaries from prot_row_counts_avg
# detect pairs by stripping the trailing suffixes "_log2_avg" and "_log2_quan_values_avg".
names_avg <- names(prot_row_counts_avg)

identified_suffix <- "_log2_avg$"
quantified_suffix <- "_log2_quan_values_avg$"

identified_keys <- names_avg[str_detect(names_avg, identified_suffix)] %>%
  str_replace(identified_suffix, "")  # remove suffix -> base key

quantified_keys <- names_avg[str_detect(names_avg, quantified_suffix)] %>%
  str_replace(quantified_suffix, "")  # remove suffix -> base key

# Intersect to ensure to only handle complete pairs
common_keys <- intersect(identified_keys, quantified_keys)

# Helper to build the full list names back from a base key
name_identified <- function(key) paste0(key, "_log2_avg")
name_quantified <- function(key) paste0(key, "_log2_quan_values_avg")

# Function that reproduces per-dataset pipeline in a generic way
build_pair_outputs <- function(key) {
  # Extract the two data frames
  df_ident <- prot_row_counts_avg[[name_identified(key)]]
  df_quant <- prot_row_counts_avg[[name_quantified(key)]]

  # Apply renaming
  df_ident2 <- df_ident %>%
    rename(
      mean_noF = Num_Rows_mean,
      SD_noF   = Num_Rows_sd,
      Rep1_noF = Rep1,
      Rep2_noF = Rep2,
      Rep3_noF = Rep3
    )
  df_quant2 <- df_quant %>%
    rename(
      mean_f = Num_Rows_mean,
      SD_f   = Num_Rows_sd,
      Rep1_F = Rep1,
      Rep2_F = Rep2,
      Rep3_F = Rep3
    )

  # Join by label
  comb <- df_ident2 %>%
    left_join(df_quant2, by = "label")

  # Long-format for plotting
  long_counts <- comb %>%
    select(label, mean_noF, SD_noF, mean_f, SD_f,
           Rep1_noF, Rep2_noF, Rep3_noF, Rep1_F, Rep2_F, Rep3_F) %>%
    pivot_longer(
      cols = c(mean_noF, mean_f),
      names_to = "Type",
      values_to = "Value"
    ) %>%
    mutate(
      SD_value = case_when(
        Type == "mean_noF" ~ SD_noF,
        Type == "mean_f"   ~ SD_f,
        TRUE ~ NA_real_
      ),
      Label_Pos = Value + SD_value + 100,
      Type = recode(Type,
                    mean_noF = "Identified",
                    mean_f   = "Quantified")
    )

  # Loss annotations
  loss_ann <- long_counts %>%
    filter(Type %in% c("Identified", "Quantified")) %>%
    select(label, Type, Value,
           Rep1_noF, Rep2_noF, Rep3_noF, Rep1_F, Rep2_F, Rep3_F, SD_f) %>%
    pivot_wider(names_from = Type, values_from = Value) %>%
    mutate(
      perc_loss = round((1 - Quantified / Identified) * 100, 1),
      loss = paste0("▼ ", perc_loss, "%")
    )

  list(long_counts = long_counts, loss_annotations = loss_ann)
}

# Build outputs for all pairs and collect into named lists
pair_outputs <- map(common_keys, build_pair_outputs)
names(pair_outputs) <- common_keys

# the two lists with all the respective dfs
long_prot_counts <- map(pair_outputs, "long_counts")
loss_annotations  <- map(pair_outputs, "loss_annotations")

# concatenate all datsets of each list into one df with the identifier "id"
long_prot_counts_all <- imap_dfr(long_prot_counts, ~ mutate(.x, id = .y))
loss_annotations_all <- imap_dfr(loss_annotations,  ~ mutate(.x, id = .y))
