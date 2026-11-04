# ===========================================================================
# Script:  05i_boxplots_method_wise_all_ratios.R
# Purpose: Method-wise boxplots of all ratios; defines
#          plot_boxplot_vertical_prot_by_method.
# Inputs:  prot_results (03a), shared functions.
# Outputs: p_boxplot_isol_w_IT_53ng_all_ratios_DDA (re-saved via saveRDS in 07)
#          and ..._DIA.
# Depends: Run after 03a.
# ===========================================================================

plot_boxplot_vertical_prot_by_method <- function(
    long_df,
    method_label,
    comparison_cols,
    ratio_human = NULL,
    ratio_yeast = NULL,
    comparison_labels = NULL,
    y_limits = NULL
) {

# ---- Filter to a single method ---------------------------------------------
  df <- long_df %>% filter(label == method_label)

# ---- Reshape wide -> long over the ratio columns ---------------------------
  df <- df %>%
    pivot_longer(
      cols = all_of(comparison_cols),
      names_to = "comparison",
      values_to = "value"
    )

  df$comparison <- factor(df$comparison, levels = comparison_cols)
  x_labels <- if (!is.null(comparison_labels)) comparison_labels else comparison_cols

# ---- Numeric x position (needed to combine with geom_segment) --------------
  df$x_num <- as.numeric(df$comparison)

# ---- Determine organism ----------------------------------------------------
  df <- df %>%
    mutate(Organism = case_when(
      str_detect(Description, "OS=Homo sapiens") ~ "Human",
      str_detect(Description, "OS=Saccharomyces") ~ "Yeast",
      TRUE ~ NA_character_
    )) %>%
    filter(!is.na(Organism))

# ---- Compute medians per comparison & organism -----------------------------
  medians <- df %>%
    group_by(comparison, x_num, Organism) %>%
    summarise(median_val = median(value, na.rm = TRUE), .groups = "drop")

# ---- Colors / dodge --------------------------------------------------------
  color_values <- c("Human" = "#56B4E9", "Yeast" = "#E69F00")
  dodge_width <- 0.7

# ---- Y-axis range ----------------------------------------------------------
  if (is.null(y_limits)) {
    y_min <- min(df$value, unlist(ratio_human), unlist(ratio_yeast), na.rm = TRUE)
    y_max <- max(df$value, unlist(ratio_human), unlist(ratio_yeast), na.rm = TRUE)
  } else {
    y_min <- y_limits[1]
    y_max <- y_limits[2]
  }
  annotation_y <- y_max * 0.95

# ---- Build reference-line segments (one per comparison category) -----------
  ref_lines <- NULL
  if (!is.null(ratio_human) || !is.null(ratio_yeast)) {
    ref_lines <- tibble(
      x_num = seq_along(comparison_cols),
      ratio_human = if (!is.null(ratio_human)) ratio_human[comparison_cols] else NA_real_,
      ratio_yeast = if (!is.null(ratio_yeast)) ratio_yeast[comparison_cols] else NA_real_
    ) %>%
      mutate(xmin = x_num - 0.45, xmax = x_num + 0.45)
  }

# ---- Plot ------------------------------------------------------------------
  p <- ggplot(df, aes(x = x_num, y = value, fill = Organism, group = interaction(x_num, Organism)))

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
    scale_x_continuous(breaks = seq_along(comparison_cols), labels = x_labels,
                       expand = expansion(add = 0.6)) +
    scale_y_continuous(limits = c(y_min, y_max), expand = c(0, 0), breaks = seq(floor(y_min), ceiling(y_max), by = 1)) +
    theme_minimal(base_size = 10) +
    labs(
      y = expression(Log[2]~"[FC]"),
      fill = "Organism",
      title = method_label
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


# ---- Plot ------------------------------------------------------------------

p_boxplot_isol_w_IT_53ng_all_ratios_DIA <- plot_boxplot_vertical_prot_by_method(
  long_df = prot_results$long_prot_standard_StoN_grouping_20_log2,
  method_label = "DIA 1.0 m/z \nIT8, 53.75ng",
  comparison_cols = c("unnorm_H100Y100vsH100Y010", "unnorm_H100Y100vsH100Y020", "unnorm_H100Y020vsH100Y010"),
  comparison_labels = c("H100Y100/H100Y010", "H100Y100/H100Y020", "H100Y020/H100Y010"),
  ratio_human = c(unnorm_H100Y100vsH100Y010 = 0, unnorm_H100Y100vsH100Y020 = 0, unnorm_H100Y020vsH100Y010 = 0),
  ratio_yeast = c(unnorm_H100Y100vsH100Y010 = 3.32, unnorm_H100Y100vsH100Y020 = 2.32, unnorm_H100Y020vsH100Y010 = 1),
  y_limits = c(-2, 6)
)
p_boxplot_isol_w_IT_53ng_all_ratios_DIA <- p_boxplot_isol_w_IT_53ng_all_ratios_DIA + labs(title = "DIA 1.0 m/z, IT8, 53.75ng")
p_boxplot_isol_w_IT_53ng_all_ratios_DIA

ggsave("plots/p_boxplot_isol_w_IT_53ng_all_ratios_DIA.svg", p_boxplot_isol_w_IT_53ng_all_ratios_DIA, 
       width = 12.7, height = 8.9, units = "cm", dpi = 300#, scale = 1 
)

saveRDS(p_boxplot_isol_w_IT_53ng_all_ratios_DIA, "p2_multi005_DIA.rds")  


p_boxplot_isol_w_IT_53ng_all_ratios_DDA <- plot_boxplot_vertical_prot_by_method(
  long_df = prot_results$long_prot_standard_StoN_grouping_20_log2,
  method_label = "DDA 0.5 m/z \nIT8, 53.75ng",
  comparison_cols = c("unnorm_H100Y100vsH100Y010", "unnorm_H100Y100vsH100Y020", "unnorm_H100Y020vsH100Y010"),
  comparison_labels = c("H100Y100/H100Y010", "H100Y100/H100Y020", "H100Y020/H100Y010"),
  ratio_human = c(unnorm_H100Y100vsH100Y010 = 0, unnorm_H100Y100vsH100Y020 = 0, unnorm_H100Y020vsH100Y010 = 0),
  ratio_yeast = c(unnorm_H100Y100vsH100Y010 = 3.32, unnorm_H100Y100vsH100Y020 = 2.32, unnorm_H100Y020vsH100Y010 = 1),
  y_limits = c(-2, 6)
)
p_boxplot_isol_w_IT_53ng_all_ratios_DDA <- p_boxplot_isol_w_IT_53ng_all_ratios_DDA + labs(title = "DDA 0.5 m/z, IT8, 53.75ng")
p_boxplot_isol_w_IT_53ng_all_ratios_DDA

ggsave("plots/p_boxplot_isol_w_IT_53ng_all_ratios_DDA.svg", p_boxplot_isol_w_IT_53ng_all_ratios_DDA, 
       width = 12.7, height = 8.9, units = "cm", dpi = 300#, scale = 1 
)

saveRDS(p_boxplot_isol_w_IT_53ng_all_ratios_DDA, "p1_multi005_DDA.rds")  
