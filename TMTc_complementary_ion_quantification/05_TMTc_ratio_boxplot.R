# ===========================================================================
# Script:  05_TMTc_ratio_boxplot.R
# Purpose: Boxplot of background-subtracted (126-median subtracted) TMTc log2
#          ratios for the 131N and 135N channels, for a single experiment
#          label (DIA120 50ng, 2 replicates).
# Inputs:  ratio_data_TMTc_subtr (built in script 01).
# Outputs: p_boxplot_TMTc_subtr_F2500_myexp, p4_multi007.rds and the saved
#          plots/p_boxplot_TMTc_subtr_F2500_myexp.svg.
# Depends: 01_TMTc_import_and_noise_subtraction.R
# ===========================================================================

# Columns used: filter, experiment, species, ratio131N, ratio135N,
#               log2_ratio131N, log2_ratio135N

plot_boxplot_TMTc_subtr_by_experiment <- function(
    ratio_df,                    # e.g. ratio_data_TMTc_subtr
    experiment_label,               # e.g. "DDA 50ng" - flexible, single condition
    filter_label = "F2500",         # which filter level to use
    channel_cols = c("log2_ratio131N", "log2_ratio135N"),
    channel_labels = c("131N", "135N"),
    ratio_human = NULL,             # optional named vector, names = channel_cols
    ratio_yeast = NULL,             # optional named vector, names = channel_cols
    y_limits = NULL
) {

  # --- Filter to one filter-level + one experiment ---
  df <- ratio_df %>%
    filter(filter == filter_label, experiment == experiment_label)

  if (nrow(df) == 0) {
    stop("No data found for filter = '", filter_label,
         "' and experiment = '", experiment_label, "'. Check spelling/levels.")
  }

  # --- Reshape wide -> long over the two channel columns ---
  df <- df %>%
    select(species, replicate, all_of(channel_cols)) %>%
    pivot_longer(
      cols = all_of(channel_cols),
      names_to = "comparison",
      values_to = "value"
    ) %>%
    filter(is.finite(value))

  df$comparison <- factor(df$comparison, levels = channel_cols)
  x_labels <- channel_labels
  df$x_num <- as.numeric(df$comparison)

  # --- Organism directly from existing species column ---
  df <- df %>%
    rename(Organism = species) %>%
    filter(Organism %in% c("Human", "Yeast"))

  # --- Medians per channel & organism ---
  medians <- df %>%
    group_by(comparison, x_num, Organism) %>%
    summarise(median_val = median(value, na.rm = TRUE), .groups = "drop")

  color_values <- c("Human" = "#56B4E9", "Yeast" = "#E69F00")
  dodge_width <- 0.7

  # --- Y-axis range ---
  if (is.null(y_limits)) {
    y_min <- min(df$value, unlist(ratio_human), unlist(ratio_yeast), na.rm = TRUE)
    y_max <- max(df$value, unlist(ratio_human), unlist(ratio_yeast), na.rm = TRUE)
  } else {
    y_min <- y_limits[1]
    y_max <- y_limits[2]
  }
  annotation_y <- y_max * 0.95

  # --- Reference-line segments (optional) ---
  ref_lines <- NULL
  if (!is.null(ratio_human) || !is.null(ratio_yeast)) {
    ref_lines <- tibble(
      x_num = seq_along(channel_cols),
      ratio_human = if (!is.null(ratio_human)) ratio_human[channel_cols] else NA_real_,
      ratio_yeast = if (!is.null(ratio_yeast)) ratio_yeast[channel_cols] else NA_real_
    ) %>%
      mutate(xmin = x_num - 0.45, xmax = x_num + 0.45)
  }

  # --- Plot ---
  p <- ggplot(df, aes(x = x_num, y = value, fill = Organism,
                    group = interaction(x_num, Organism)))

  if (!is.null(ref_lines)) {
    p <- p +
      geom_segment(data = ref_lines,
                   aes(x = xmin, xend = xmax, y = ratio_human, yend = ratio_human),
                   inherit.aes = FALSE, linetype = "dashed", color = "#2E6A87", linewidth = 0.5) +
      geom_segment(data = ref_lines,
                   aes(x = xmin, xend = xmax, y = ratio_yeast, yend = ratio_yeast),
                   inherit.aes = FALSE, linetype = "dashed", color = "#B87F00", linewidth = 0.5)
  }

  p <- p +
    geom_boxplot(width = 0.6,
                 position = position_dodge2(width = dodge_width, padding = 0.1),
                 outlier.shape = NA,
                 box.linewidth = 0.25,
                 whisker.linewidth = 0.25,
                 median.linewidth = 0.25) +
    geom_text(
      data = medians,
      aes(x = x_num, y = annotation_y, label = round(median_val, 2), group = Organism, fill = NULL),
      position = position_dodge(width = dodge_width),
      vjust = 0, size = 1.7, color = "black"
    ) +
    scale_fill_manual(values = color_values) +
    scale_x_continuous(breaks = seq_along(channel_cols), labels = x_labels,
                    expand = expansion(add = 0.6)) +
    scale_y_continuous(limits = c(y_min, y_max), expand = c(0, 0), breaks = seq(floor(y_min), ceiling(y_max), by = 1)) +
    theme_minimal(base_size = 10) +
    labs(
      y = expression(Log[2]~"[Ratio]"),
      fill = "Organism",
      title = paste0(experiment_label, " | ", filter_label, " (TMTc, subtr.)")
    ) +
    theme(
      plot.title = element_text(size = 6, face = "bold", hjust = 0.5, colour = "black"),
      axis.title.x = element_blank(),
      axis.line = element_line(color = "black"),
      axis.title.y = element_text(size = 6, face = "bold", colour = "black"),
      axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, size = 6, colour = "black"),
      axis.text.y = element_text(size = 6, colour = "black"),
      legend.position = "inside",
      legend.position.inside = c(0.9, 0.8),
      legend.text = element_text(size = 5.5),
      legend.key.size = unit(0.5, "cm"),
      legend.title = element_blank(),
      panel.grid = element_blank(),
      axis.ticks = element_line(color = "black", size = 0.25)
    )

  return(p)
}


# --- Build the figure-panel boxplot (DIA120 50ng, F2500) ---
p_boxplot_TMTc_subtr_F2500_myexp <- plot_boxplot_TMTc_subtr_by_experiment(
  ratio_df         = ratio_data_TMTc_subtr,
  experiment_label = "DIA120 50ng",      # experiment with 2 replicates used in the figure
  filter_label     = "F2500",
  channel_cols     = c("log2_ratio131N", "log2_ratio135N"),
  channel_labels   = c("131N", "135N"),
  ratio_human      = c(log2_ratio131N = 0, log2_ratio135N = 0),
  ratio_yeast      = c(log2_ratio131N = log2(5), log2_ratio135N = log2(5)),
  y_limits         = c(-2, 6)
)

p_boxplot_TMTc_subtr_F2500_myexp

saveRDS(p_boxplot_TMTc_subtr_F2500_myexp, "p4_multi007.rds")

ggsave("plots/p_boxplot_TMTc_subtr_F2500_myexp.svg", p_boxplot_TMTc_subtr_F2500_myexp,
       width = 12.7, height = 8.9, units = "cm", dpi = 300
)
