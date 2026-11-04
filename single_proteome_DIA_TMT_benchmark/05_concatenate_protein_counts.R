# =============================================================================
# Script 05: Concatenate identified vs. quantified protein counts
# Description: Pairs the identified and quantified per-replicate count summaries
#              for each dataset, reshapes them to long format and computes the
#              quantification loss per label, producing long_prot_counts and
#              loss_annotations used by the downstream plots.
# Depends on: 03_filter_transform_and_counts.R
# =============================================================================

# --- Detect identified / quantified pairs ------------------------------------

# Pairs are detected by stripping the trailing "_log2_avg" (identified) and
# "_log2_quan_values_avg" (quantified) suffixes to a common base key.
names_avg <- names(prot_row_counts_avg)

identified_suffix <- "_log2_avg$"
quantified_suffix <- "_log2_quan_values_avg$"

identified_keys <- names_avg[str_detect(names_avg, identified_suffix)] %>%
  str_replace(identified_suffix, "")

quantified_keys <- names_avg[str_detect(names_avg, quantified_suffix)] %>%
  str_replace(quantified_suffix, "")

# Keep only base keys that have both an identified and a quantified summary.
common_keys <- intersect(identified_keys, quantified_keys)

# --- Helpers to rebuild the full list names from a base key ------------------

name_identified <- function(key) paste0(key, "_log2_avg")
name_quantified <- function(key) paste0(key, "_log2_quan_values_avg")

# --- Build long-format counts and loss annotations per pair ------------------

build_pair_outputs <- function(key) {
  df_ident <- prot_row_counts_avg[[name_identified(key)]]
  df_quant <- prot_row_counts_avg[[name_quantified(key)]]

  # Rename to distinguish non-filtered (identified) and filtered (quantified).
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

  comb <- df_ident2 %>%
    left_join(df_quant2, by = "label")

  # Long format for plotting.
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

  # Per-label quantification loss (identified -> quantified).
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

# --- Collect outputs across all pairs ----------------------------------------

pair_outputs <- map(common_keys, build_pair_outputs)
names(pair_outputs) <- common_keys

long_prot_counts <- map(pair_outputs, "long_counts")
loss_annotations  <- map(pair_outputs, "loss_annotations")

# Concatenate all datasets of each list into one data frame keyed by "id".
long_prot_counts_all <- imap_dfr(long_prot_counts, ~ mutate(.x, id = .y))
loss_annotations_all <- imap_dfr(loss_annotations,  ~ mutate(.x, id = .y))
