# ===========================================================================
# Script:  04_protein_counts_TMT_vs_TMTc.R
# Purpose: Count quantified Human/Yeast proteins for standard TMT reporter-ion
#          (baseline/standard) and TMTc quantification, then build the
#          three-way comparison bar plot (unfiltered RI / filtered RI / TMTc).
# Inputs:  prot_results (script 03), counts_by_df_TMTc_subtr (script 01),
#          all_cond (defined below).
# Outputs: p_three_way_total and the saved plots/p_three_way_total.svg.
# Depends: 01_TMTc_import_and_noise_subtraction.R,
#          03_TMT_reporter_ion_filtering.R
# ===========================================================================

# --- Condition order used for plotting ---
# Entries must match the label_patterns / experiment-label output exactly.
all_cond <- c("DDA 50ng",
              "DDA 100ng",
              "DIA 50ng",
              "DIA 100ng",
              "HM-DIA 50ng",
              "HM-DIA 100ng",
              "DIA120 50ng",
              "DIA120 100ng")

# --- Count Human/Yeast proteins per replicate (normal TMT) ---
# One row per dataset x label x Organism x filter
filter_group_map <- c(
  baseline = "long_prot_baseline_log2_quan_values",
  standard = "long_prot_standard_log2_quan_values"
)

counts_by_df_normal <- imap_dfr(filter_group_map, function(result_key, filter_name) {
  prot_results[[result_key]] %>%
    group_by(dataset, label, Organism) %>%
    summarise(count = n(), .groups = "drop") %>%
    filter(Organism %in% c("Human", "Yeast")) %>%
    mutate(filter = filter_name)
})

# --- Count Human/Yeast proteins per replicate (TMTc) ---
# counts_by_df_TMTc_subtr was built earlier in the pipeline (script 01).
# Verify it exists and has the expected columns.
stopifnot(all(c("filter", "df_name", "human", "yeast", "experiment") %in%
                names(counts_by_df_TMTc_subtr)))

# --- Total (Human + Yeast) summary for normal TMT ---
# Sum species within each replicate first, THEN compute mean + SD.
counts_total_normal <- counts_by_df_normal %>%
  # sum Human + Yeast within each replicate
  group_by(filter, dataset, label) %>%
  summarise(total = sum(count), .groups = "drop") %>%
  # average across replicates
  group_by(filter, label) %>%
  summarise(
    count    = mean(total),
    count_sd = sd(total),
    .groups  = "drop"
  ) %>%
  rename(experiment = label) %>%
  mutate(
    experiment = factor(experiment, levels = c(
      "DDA 50ng",    "DDA 100ng",
      "DIA 50ng",    "DIA 100ng",
      "HM-DIA 50ng", "HM-DIA 100ng",
      "DIA120 50ng", "DIA120 100ng"
    )),
    filter = factor(filter, levels = c(
      "baseline", "standard"
    ))
  )

# --- Total (Human + Yeast) summary for TMTc ---
# Human + yeast already summed at replicate level -> mean + SD directly
counts_total_TMTc <- counts_by_df_TMTc_subtr %>%
  mutate(total = human + yeast) %>%
  group_by(filter, experiment) %>%
  summarise(
    count    = mean(total),
    count_sd = sd(total),
    .groups  = "drop"
  ) %>%
  mutate(
    experiment = factor(experiment, levels = c(
      "DDA 50ng",    "DDA 100ng",
      "DIA 50ng",    "DIA 100ng",
      "HM-DIA 50ng", "HM-DIA 100ng",
      "DIA120 50ng", "DIA120 100ng"
    )),
    filter = factor(filter, levels = c("F2500"))
  )

# --- Build combined summary dataframe ---
# Mean + SD of total quantified proteins per source x experiment
counts_combined_total <- bind_rows(
  # Normal TMT unfiltered
  counts_total_normal %>%
    filter(filter == "baseline") %>%
    mutate(source = "noF"),

  # Normal TMT filtered
  counts_total_normal %>%
    filter(filter == "standard") %>%
    mutate(source = "F"),

  # TMTc (F2500 filter)
  counts_total_TMTc %>%
    filter(filter == "F2500") %>%
    mutate(source = "TMTc")
) %>%
  mutate(
    source     = factor(source,     levels = c("noF", "F", "TMTc")),
    experiment = factor(experiment, levels = c(
      "DDA 50ng",    "DDA 100ng",
      "DIA 50ng",    "DIA 100ng",
      "HM-DIA 50ng", "HM-DIA 100ng",
      "DIA120 50ng", "DIA120 100ng"
    ))
  ) %>%
  select(source, experiment, count, count_sd)

# Verify: count_sd should be non-NA for all rows
counts_combined_total %>% print(n = 40)

# --- Build replicate-level dataframe for dot overlays ---
# One row per replicate file x source, total proteins
reps_noF <- counts_by_df_normal %>%
  filter(filter == "baseline") %>%
  group_by(dataset, label) %>%
  summarise(count = sum(count), .groups = "drop") %>%
  mutate(
    source     = "noF",
    experiment = factor(label, levels = levels(counts_combined_total$experiment))
  ) %>%
  select(source, experiment, dataset, count)

reps_F <- counts_by_df_normal %>%
  filter(filter == "standard") %>%
  group_by(dataset, label) %>%
  summarise(count = sum(count), .groups = "drop") %>%
  mutate(
    source     = "F",
    experiment = factor(label, levels = levels(counts_combined_total$experiment))
  ) %>%
  select(source, experiment, dataset, count)

reps_TMTc <- counts_by_df_TMTc_subtr %>%
  filter(filter == "F2500") %>%
  mutate(
    count      = human + yeast,
    source     = "TMTc",
    experiment = factor(experiment, levels = levels(counts_combined_total$experiment))
  ) %>%
  select(source, experiment, dataset = df_name, count)

# Sanity check before binding
cat("reps_noF:  class =", class(reps_noF),  "| nrow =", nrow(reps_noF),  "\n")
cat("reps_F:    class =", class(reps_F),    "| nrow =", nrow(reps_F),    "\n")
cat("reps_TMTc: class =", class(reps_TMTc), "| nrow =", nrow(reps_TMTc), "\n")

counts_combined_replicates <- bind_rows(reps_noF, reps_F, reps_TMTc) %>%
  mutate(source = factor(source, levels = c("noF", "F", "TMTc")))

# --- Plot function ---
# df already contains totals (Human + Yeast combined) with correct SD.
plot_TMT_three_way_total <- function(df,
                    df_replicates         = NULL,
                    condition_experiments = NULL,
                    y_limit               = 5500,
                    x_title               = "",
                    title_suffix          = "") {

  all_experiments <- if (!is.null(condition_experiments)) {
    condition_experiments
  } else {
    levels(df$experiment)
  }

  data_agg <- df %>%
    filter(experiment %in% all_experiments) %>%
    mutate(
      experiment = factor(experiment, levels = all_experiments),
      count_sd   = replace_na(count_sd, 0)
    )

  # Build x_key: noF > F > TMTc within each experiment group
  x_key_levels <- unlist(lapply(all_experiments, function(exp) {
    present_sources <- data_agg %>%
      filter(experiment == exp) %>%
      pull(source) %>% unique() %>% as.character()
    ordered_sources <- intersect(c("noF", "F", "TMTc"), present_sources)
    paste(exp, ordered_sources, sep = "__")
  }))

  # Insert spacer between experiment groups for visual separation
  spacer_levels <- unlist(lapply(seq_along(all_experiments), function(i) {
    exp      <- all_experiments[i]
    exp_keys <- x_key_levels[grepl(paste0("^", exp, "__"), x_key_levels)]
    if (i < length(all_experiments)) c(exp_keys, paste0("spacer_", i))
    else exp_keys
  }))

  # One tick per experiment group, placed on the center bar of the group
  center_keys <- sapply(all_experiments, function(exp) {
    keys <- x_key_levels[grepl(paste0("^", exp, "__"), x_key_levels)]
    keys[ceiling(length(keys) / 2)]
  })
  group_labels <- setNames(all_experiments, center_keys)

  # Colour scheme: RdBu palette
  source_colours <- c(
    "noF"  = "#ca0020",
    "F"    = "#f4a582",
    "TMTc" = "#92c5de"
  )

  # Add x_key and bar styling columns
  data_agg <- data_agg %>%
    mutate(
      x_key     = factor(paste(experiment, source, sep = "__"), levels = spacer_levels),
      fill_col  = source_colours[as.character(source)],
      alpha_val = case_when(
        source == "noF"  ~ 0.5,
        source == "F"    ~ 1,
        source == "TMTc" ~ 1,
        TRUE             ~ 1
      )
    )

  # Add x_key and styling to replicate dots
  if (!is.null(df_replicates)) {
    reps <- df_replicates %>%
      filter(experiment %in% all_experiments) %>%
      mutate(
        experiment = factor(experiment, levels = all_experiments),
        x_key      = factor(paste(experiment, source, sep = "__"), levels = spacer_levels),
        fill_col   = source_colours[as.character(source)],
        alpha_val  = case_when(
          source == "noF"  ~ 0.6,
          source == "F"    ~ 0.8,
          source == "TMTc" ~ 0.8,
          TRUE             ~ 0.7
        )
      )
  }

  p <- ggplot(data_agg, aes(x = x_key, y = count)) +
    # Bars — outline always black
    geom_col(
      aes(fill = source, alpha = alpha_val),
      color = "black", linewidth = 0.25, width = 0.7
    ) +
    # Error bars — always black
    geom_errorbar(aes(ymin = pmax(0, count - count_sd),
                    ymax = count + count_sd),
                  color = "black", width = 0.2, linewidth = 0.2)

  # Optional replicate dots — excluded from legend
  if (!is.null(df_replicates)) {
    p <- p +
      geom_point(
        data = reps,
        aes(x = x_key, y = count, fill = source, alpha = alpha_val),
        shape = 21, color = "black", size = 0.5,
        show.legend = FALSE,
        position = position_jitter(width = 0.08, seed = 42)
      )
  }

  p <- p +
    # Count labels above error bar cap
    geom_text(aes(y     = count + count_sd + y_limit * 0.015,
                  label = round(count, 0)),
              angle = 90, hjust = -0.1, vjust = 0.4,
              size = 1.5, color = "black") +
    # Manual fill legend — bars only, no dot shape in keys
    scale_fill_manual(
      name = "Quantification",
      values = source_colours,
      breaks = c("noF", "F", "TMTc"),
      labels = c(
        "noF"  = "TMT RI (unfiltered)",
        "F"    = "TMT RI (filtered)",
        "TMTc" = "TMTc (filtered)"
      ),
      guide = guide_legend(
        override.aes = list(
          alpha = c(0.5, 1, 1),
          color = "black"
        )
      )
    ) +
    scale_alpha_identity(guide = "none") +
    # One tick per experiment group, centered on the middle bar
    scale_x_discrete(
      breaks = center_keys,
      labels = group_labels,
      drop   = FALSE
    ) +
    scale_y_continuous(limits = c(0, y_limit), expand = c(0, 0)) +
    labs(y = "#quantified proteins", x = x_title, title = title_suffix) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title       = element_text(face = "plain", hjust = 0.5, size = 6),
      axis.text.x      = element_text(angle = 45, hjust = 1, vjust = 1, size = 6, colour = "black"),
      axis.text.y      = element_text(size = 6, colour = "black"),
      axis.title       = element_text(size = 6, colour = "black"),
      legend.title     = element_text(size = 4),
      legend.text      = element_text(size = 4),
      legend.key.size  = unit(0.2, "cm"),
      legend.position  = "right",
      panel.grid       = element_blank(),
      panel.background = element_blank(),
      axis.line        = element_line(colour = "black"),
      axis.ticks       = element_line(color = "black", linewidth = 0.25)
    )

  p
}

# --- Generate and save plot ---
p_three_way_total <- plot_TMT_three_way_total(
  df                    = counts_combined_total,
  df_replicates         = counts_combined_replicates,
  condition_experiments = all_cond,
  y_limit               = 7500,
  title_suffix          = "Quantified proteins across methods before and after filtering\n for TMT RI and TMTc quantification"
)
p_three_way_total

ggsave("plots/p_three_way_total.svg", device = "svg", p_three_way_total,
       width = 12.7, height = 8.9, units = "cm", dpi = 300)
