# =============================================================================
# Script 09: Ratio boxplots and CV violin plots
# Description: Defines plot_boxplot_vertical_prot() and plot_violin_cv_prot() for
#              per-method log2 fold-change and %CV distributions, and builds the
#              H100/H010 ratio boxplot and H100 CV violin figures for the paper.
# Depends on: 03_filter_transform_and_counts.R, 06_condition_subsets_and_colors.R
# =============================================================================

# --- Fold-change boxplot function --------------------------------------------

plot_boxplot_vertical_prot <- function(
    long_df,
    comparison_col,
    labels_to_plot = NULL,
    ratio_human = 0,
    ratio_yeast = 2.32,
    y_limits = NULL,
    title = NULL,
    plot_yeast = TRUE  # Controls both the reference line and the plotted data.
) {

  # --- Filter data by labels ---
  df <- long_df
  if (!is.null(labels_to_plot)) {
    df <- df %>% filter(label %in% labels_to_plot)
    df$label <- factor(df$label, levels = labels_to_plot)
  }

  # --- Determine organism ---
  df <- df %>%
    mutate(Organism = case_when(
      str_detect(Description, "OS=Homo sapiens") ~ "Human",
      str_detect(Description, "OS=Saccharomyces") ~ "Yeast",
      TRUE ~ NA_character_
    )) %>%
    filter(!is.na(Organism))

  # --- Drop yeast entirely if not requested ---
  if (!plot_yeast) {
    df <- df %>% filter(Organism != "Yeast")
  }

  # --- Compute medians per label & organism ---
  medians <- df %>%
    group_by(label, Organism) %>%
    summarise(median_val = median(!!sym(comparison_col), na.rm = TRUE), .groups = "drop")

  # --- Define colors ---
  color_values <- c("Human" = "#56B4E9", "Yeast" = "#E69F00")
  dodge_width <- 0.7

  # --- Determine y-axis range ---
  if (is.null(y_limits)) {
    # Only include ratio_yeast in the range when yeast is actually plotted.
    y_vals <- c(df[[comparison_col]], ratio_human)
    if (plot_yeast) y_vals <- c(y_vals, ratio_yeast)

    y_min <- min(y_vals, na.rm = TRUE)
    y_max <- max(y_vals, na.rm = TRUE)
  } else {
    y_min <- y_limits[1]
    y_max <- y_limits[2]
  }

  annotation_y <- y_max * 0.95

  # --- Plot ---
  p <- ggplot(df, aes(x = label, y = !!sym(comparison_col), fill = Organism)) +
    geom_hline(yintercept = ratio_human, linetype = "dashed", color = "#2E6A87", size = 0.5)

  # Only add the yeast reference line when requested.
  if (plot_yeast) {
    p <- p + geom_hline(yintercept = ratio_yeast, linetype = "dashed", color = "#B87F00", size = 0.5)
  }

  p <- p +
    geom_boxplot(width = 0.6,
                 position = position_dodge(width = dodge_width),
                 outlier.shape = NA,
                 box.linewidth = 0.25,
                 whisker.linewidth = 0.25,
                 median.linewidth = 0.25) +
    geom_text(
      data = medians,
      aes(
        x = label,
        y = annotation_y,
        label = round(median_val, 2),
        group = Organism
      ),
      position = position_dodge(width = dodge_width),
      vjust = 0,
      size = 1.7,
      color = "black"
    ) +
    scale_fill_manual(values = color_values) +
    scale_y_continuous(limits = c(y_min, y_max), expand = c(0, 0)) +
    theme_minimal(base_size = 10) +
    labs(
      title = if (!is.null(title)) paste0("Log2 FC across methods at 54 ng input: ", title)
      else "Log2 FC across methods at 54 ng input",
      y = expression(Log[2]~"[FC]"),
      fill = "Organism"
    ) +
    theme(
      plot.title   = element_text(lineheight = .8, face = "plain", hjust = 0.5, size = 6),
      axis.title.x = element_blank(),
      axis.line = element_line(color = "black"),
      axis.title.y = element_text(size = 6, face = "bold", colour = "black"),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6, colour = "black"),
      axis.text.y = element_text(size = 6, colour = "black"),
      legend.position = "inside",
      legend.position.inside = c(0.9,0.8),
      legend.text = element_text(size = 5.5),
      legend.key.size = unit(0.5, "cm"),
      legend.title = element_blank(),
      panel.grid = element_blank(),
      axis.ticks = element_line(color = "black", size = 0.25)
    )

  return(p)
}

# --- CV violin plot function -------------------------------------------------

plot_violin_cv_prot <- function(
    long_df,
    cv_col,
    labels_to_plot = NULL,
    y_limits = NULL,
    title = NULL,
    plot_yeast = TRUE  # Controls whether yeast data is plotted.
) {

  # --- Filter data ---
  df <- long_df
  if (!is.null(labels_to_plot)) {
    df <- df %>% filter(label %in% labels_to_plot)
    df$label <- factor(df$label, levels = labels_to_plot)
  }

  # --- Determine organism ---
  df <- df %>%
    mutate(Organism = case_when(
      str_detect(Description, "OS=Homo sapiens") ~ "Human",
      str_detect(Description, "OS=Saccharomyces") ~ "Yeast",
      TRUE ~ NA_character_
    )) %>%
    filter(!is.na(Organism))

  # --- Drop yeast entirely if not requested ---
  if (!plot_yeast) {
    df <- df %>% filter(Organism != "Yeast")
  }

  # --- Compute medians per label & organism ---
  medians <- df %>%
    group_by(label, Organism) %>%
    summarise(median_val = median(!!sym(cv_col), na.rm = TRUE), .groups = "drop")

  # --- Determine y-axis range ---
  if (is.null(y_limits)) {
    y_min <- min(df[[cv_col]], na.rm = TRUE)
    y_max <- max(df[[cv_col]], na.rm = TRUE)
  } else {
    y_min <- y_limits[1]
    y_max <- y_limits[2]
  }

  annotation_y <- y_max * 0.96

  color_values <- c("Human" = "#56B4E9", "Yeast" = "#E69F00")
  dodge_width <- 0.8

  # --- Plot ---
  p <- ggplot(df, aes(x = label, y = !!sym(cv_col), fill = Organism)) +
    geom_violin(
      scale = "width",
      trim = TRUE,
      width = 0.7,
      position = position_dodge(width = dodge_width),
      alpha = 0.6,
      color = "white"
    ) +
    geom_boxplot(
      aes(fill = Organism),
      width = 0.12,
      box.linewidth = 0.25,
      whisker.linewidth = 0.25,
      median.linewidth = 0.25,
      position = position_dodge(width = dodge_width),
      outlier.shape = NA,
      color = "black",
      alpha = 0.8
    ) +
    geom_text(
      data = medians,
      aes(
        x = label,
        y = annotation_y,
        label = sprintf("%.1f", median_val),
        group = Organism
      ),
      position = position_dodge(width = dodge_width),
      vjust = 0,
      size = 1.7,
      color = "black"
    ) +
    scale_fill_manual(values = color_values) +
    scale_y_continuous(limits = c(y_min, y_max), expand = c(0, 0)) +
    theme_minimal(base_size = 10) +
    labs(
      title = if (!is.null(title)) paste0("%CV across methods at 54 ng input: ", title)
      else "%CV across methods at 54 ng input",
      y = "CV [%]",
      fill = "Organism"
    ) +
    theme(
      plot.title   = element_text(lineheight = .8, face = "plain", hjust = 0.5, size = 6),
      axis.title.x = element_blank(),
      axis.line = element_line(color = "black"),
      axis.title.y = element_text(size = 6, colour = "black"),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6, colour = "black"),
      axis.text.y = element_text(size = 6, colour = "black"),
      legend.position = "inside",
      legend.text = element_text(size = 5.5),
      legend.position.inside = c(0.8,0.8),
      legend.title = element_blank(),
      panel.grid = element_blank(),
      axis.ticks = element_line(color = "black", size = 0.3)
    )

  return(p)
}

# --- Paper figure: H100/H010 ratio boxplot -----------------------------------

p_prot_ratio_boxplots_H100H010_54ng_StoN_grouping <- plot_boxplot_vertical_prot(long_df = prot_results$long_prot_standard_StoN_grouping_20_log2,
                    comparison_col = "unnorm_H100vsH010",
                    labels_to_plot = methods_54ng,
                    ratio_human = 3.32,
                    y_limits = c(-1,5),
                    plot_yeast = F,
                    title = "H100/H010")
p_prot_ratio_boxplots_H100H010_54ng_StoN_grouping

p_prot_ratio_boxplots_H100H010_54ng_StoN_grouping <- p_prot_ratio_boxplots_H100H010_54ng_StoN_grouping+
  scale_x_discrete(labels = c("DDA", "DIA", "DIA HCD 29%", "DIA no loop control")) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 6)
  )
p_prot_ratio_boxplots_H100H010_54ng_StoN_grouping

ggsave("plots/p_prot_ratio_boxplots_H100H010_54ng_StoN_grouping.svg", p_prot_ratio_boxplots_H100H010_54ng_StoN_grouping,
       width = 12.7, height = 8.9, units = "cm", dpi = 300, scale = 1)

# --- Paper figure: H100 CV violin --------------------------------------------

p_violin_methods_cv_H100_54ng_StoN_grouping <- plot_violin_cv_prot(long_df = prot_results$long_prot_standard_StoN_grouping_20_log2,
                    cv_col = "cv_H100",
                    labels_to_plot = methods_54ng,
                    y_limits = c(0,100),
                    plot_yeast = F,
                    title = "H100")
p_violin_methods_cv_H100_54ng_StoN_grouping

p_violin_methods_cv_H100_54ng_StoN_grouping <- p_violin_methods_cv_H100_54ng_StoN_grouping +
  scale_x_discrete(labels = c("DDA", "DIA", "DIA HCD 29%", "DIA no loop control")) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 6)
  )

p_violin_methods_cv_H100_54ng_StoN_grouping

ggsave("plots/p_violin_methods_cv_H100_54ng_StoN_grouping.svg", p_violin_methods_cv_H100_54ng_StoN_grouping,
       width = 12.7, height = 8.9, units = "cm", dpi = 300, scale = 1)
