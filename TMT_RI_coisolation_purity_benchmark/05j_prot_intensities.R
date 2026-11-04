# ===========================================================================
# Script:  05j_prot_intensities.R
# Purpose: Summed protein reporter intensities per channel; defines
#          sum_prot_intensities_with_mix and a linear-scale plotting function.
# Inputs:  prot_results (03a), channel_names_groups (01), okabe_ito_colors
#          (02b).
# Outputs: p_prot_intensities_DDA_IT8_54ng (Fig1).
# Depends: Run after 02b and 03a.
# ===========================================================================

sum_prot_intensities_with_mix <- function(data, channel_info, backtransform = TRUE) {
  # 1. Identify valid abundance columns
  abundance_cols <- names(data)[names(data) %in% c(
    channel_info$label,
    paste0("Pure-", channel_info$label)
  )]

  if (length(abundance_cols) == 0) {
    stop("No matching abundance columns found in the dataset.")
  }

  # 2. Pivot to long format
  long_data <- data %>%
    pivot_longer(
      cols = all_of(abundance_cols),
      names_to = "channel_name",
      values_to = "raw_value"
    ) %>%
    mutate(
      channel_name = sub("^Pure-", "", channel_name),
      # Back-transform from log2 to linear if requested, otherwise use as-is
      intensity = if (backtransform) {
        ifelse(!is.na(raw_value), 2^raw_value, NA_real_)
      } else {
        raw_value
      }
    )

  # 3. Sum intensities per dataset label, replicate, and channel
  replicate_summary <- long_data %>%
    group_by(label, replicate, channel_name) %>%
    summarise(
      sum_intensity_repl = if (all(is.na(intensity))) NA_real_ else sum(intensity, na.rm = TRUE),
      .groups = "drop"
    )

  # 4. Reshape so rep1, rep2, rep3 are separate columns
  replicates_wide <- replicate_summary %>%
    select(label, channel_name, replicate, sum_intensity_repl) %>%
    pivot_wider(
      names_from = replicate,
      values_from = sum_intensity_repl,
      names_prefix = "rep"
    )

  # Dynamic replicate column names
  replicate_cols <- paste0("rep", sort(unique(replicate_summary$replicate)))
  log10_replicate_cols <- paste0("log10_", replicate_cols)

  # 5. Calculate linear and log10 summaries across replicates
  replicates_wide <- replicates_wide %>%
    rowwise() %>%
    mutate(
      # log10 of each replicate sum
      across(
        all_of(replicate_cols),
        ~ ifelse(!is.na(.x) & .x > 0, log10(.x), NA_real_),
        .names = "log10_{col}"
      ),

      # linear-space mean and SD across replicate sums
      mean_intensity = if (all(is.na(c_across(all_of(replicate_cols))))) {
        NA_real_
      } else {
        mean(c_across(all_of(replicate_cols)), na.rm = TRUE)
      },

      sd_intensity = if (sum(!is.na(c_across(all_of(replicate_cols)))) >= 2) {
        sd(c_across(all_of(replicate_cols)), na.rm = TRUE)
      } else {
        NA_real_
      },

      # log10-space mean and SD across log10 replicate sums
      log10_mean_intensity = if (all(is.na(c_across(all_of(log10_replicate_cols))))) {
        NA_real_
      } else {
        mean(c_across(all_of(log10_replicate_cols)), na.rm = TRUE)
      },

      log10_sd_intensity = if (sum(!is.na(c_across(all_of(log10_replicate_cols)))) >= 2) {
        sd(c_across(all_of(log10_replicate_cols)), na.rm = TRUE)
      } else {
        NA_real_
      }
    ) %>%
    ungroup()

  # 6. Add experiment metadata and clean up labels
  final_output <- replicates_wide %>%
    left_join(
      channel_info %>% select(label, experiment),
      by = c("channel_name" = "label")
    ) %>%
    mutate(dataset_label = label) %>%
    select(
      dataset_label,
      channel_name,
      all_of(replicate_cols),
      all_of(log10_replicate_cols),
      mean_intensity,
      sd_intensity,
      log10_mean_intensity,
      log10_sd_intensity,
      experiment
    )

  return(final_output)
}

all_prot_intensities_baseline <- sum_prot_intensities_with_mix(prot_results$long_prot_baseline_log2,
                                                               channel_info = channel_names_groups,
                                                               backtransform = TRUE)


condition_order <- c("H000Y000",  "H100Y000", "H100Y010", "H100Y020", "H100Y100")
channel_order <- c("126","127N","127C","128N","128C","129N","129C","130N","130C",
                   "131N","131C","132N","132C","133N","133C","134N","134C","135N")

# Reorder experiment based on predefined condition_order
all_prot_intensities_baseline <- all_prot_intensities_baseline %>%
  mutate(
    experiment = factor(experiment, levels = condition_order),
    channel_name = factor(channel_name, levels = channel_order)
  ) %>%
  arrange(experiment, channel_name)


plot_reporter_intensities_linear <- function(intensity_df,
                                             label_to_plot,
                                             y_limit = NULL,
                                             fill_colors = okabe_ito_colors,
                                             title = NULL) {

  # Subset dataset
  plot_data <- intensity_df %>%
    filter(dataset_label == label_to_plot)

  if (nrow(plot_data) == 0) {
    stop("No data found for dataset_label: ", label_to_plot)
  }

  # Use linear mean intensity
  plot_data <- plot_data %>%
    mutate(y_value = mean_intensity)

  # Ensure correct ordering by joining with channel_names_groups
  plot_data <- plot_data %>%
    left_join(
      channel_names_groups %>% select(label, label),
      by = c("channel_name" = "label")
    ) %>%
    mutate(
      channel_name = factor(channel_name, levels = channel_names_groups$label)
    )

  # Pivot linear replicate columns to long format
  replicate_cols <- grep("^rep[0-9]+$", names(plot_data), value = TRUE)

  plot_data_long <- plot_data %>%
    pivot_longer(
      cols = all_of(replicate_cols),
      names_to = "replicate",
      values_to = "Replicate_Value"
    )

  dodge <- position_dodge(width = 0.6)

  p <- ggplot(plot_data, aes(x = channel_name, y = y_value, fill = experiment)) +
    geom_col(
      width = 0.6,
      position = dodge,
      color = "black",
      linewidth = 0.25
    ) +
    geom_errorbar(
      aes(
        ymin = mean_intensity - sd_intensity,
        ymax = mean_intensity + sd_intensity
      ),
      position = dodge,
      width = 0.25,
      size = 0.25
    ) +
    geom_point(
      data = plot_data_long,
      aes(y = Replicate_Value, fill = experiment),
      shape = 21,
      color = "black",
      size = 0.8,
      alpha = 0.7,
      position = position_jitter(width = 0.15, height = 0)
    ) +
    geom_text(
      aes(label = scales::scientific(y_value, digits = 2)),
      position = dodge,
      vjust = 0.4,
      hjust = -1,
      angle = 90,
      size = 1.3
    ) +
    scale_fill_manual(values = fill_colors) +
    scale_y_continuous(  
      limits = y_limit,  
      breaks = seq(0, y_limit[2], by = 0.5e8),  
      expand = expansion(mult = c(0, 0))  
    ) +
    labs(
      title = if (!is.null(title)) paste0("Summed protein intensities over all replicates: ", title)  
      else "Summed protein intensities over all replicates",
      x = "Channel",
      y = "Summed protein intensities"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6, colour = "black"),
      plot.title   = element_text(lineheight = .8, face = "plain", hjust = 0.5, size = 6),
      axis.text.y = element_text(size = 6, colour = "black"),
      panel.grid = element_blank(),
      legend.position = "none",
      axis.title.y = element_text(size = 6, colour = "black"),
      axis.title.x = element_blank(),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black", size = 0.3)
    )

  return(p)
}

p_prot_intensities_DDA_IT8_54ng <- plot_reporter_intensities_linear(
  all_prot_intensities_baseline,
  label_to_plot = "DDA 0.5 m/z \nIT8, 53.75ng",
  y_limit = c(0, 3e8),
  title = "DDA IT 8 ms"
)

p_prot_intensities_DDA_IT8_54ng

ggsave(
  "plots/p_prot_intensities_DDA_IT8_54ng_linear.svg",
  p_prot_intensities_DDA_IT8_54ng,
  width = 12.7,
  height = 8.9,
  units = "cm",
  dpi = 300
)
