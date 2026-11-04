# ===========================================================================
# Script:  05e_prot_ratio_overview.R
# Purpose: Protein ratio boxplots and CV violin plots per condition.
# Inputs:  prot_results, top20_prot (03a), condition subsets (04).
# Outputs: p_boxplot_prot_isol_w_IT_53ng_* and p_violin_prot_isol_w_IT_53ng_*
#          (Fig3).
# Depends: Run after 03a and 04.
# ===========================================================================

plot_boxplot_vertical_prot <- function(
    long_df,
    comparison_col,
    labels_to_plot = NULL,
    ratio_human = 0,
    ratio_yeast = 2.32,
    y_limits = NULL   
) {

# ---- Filter data -----------------------------------------------------------
  df <- long_df
  if (!is.null(labels_to_plot)) {
    df <- df %>% filter(label %in% labels_to_plot)
    df$label <- factor(df$label, levels = labels_to_plot)
  }

# ---- Determine organism ----------------------------------------------------
  df <- df %>%
    mutate(Organism = case_when(
      str_detect(Description, "OS=Homo sapiens") ~ "Human",
      str_detect(Description, "OS=Saccharomyces") ~ "Yeast",
      TRUE ~ NA_character_
    )) %>%
    filter(!is.na(Organism))

# ---- Compute medians per label & organism ----------------------------------
  medians <- df %>%
    group_by(label, Organism) %>%
    summarise(median_val = median(!!sym(comparison_col), na.rm = TRUE), .groups = "drop")

# ---- Define colors ---------------------------------------------------------
  color_values <- c("Human" = "#56B4E9", "Yeast" = "#E69F00")
  dodge_width <- 0.7

# ---- Determine y-axis range ------------------------------------------------
  if (is.null(y_limits)) {
    y_min <- min(df[[comparison_col]], ratio_human, ratio_yeast, na.rm = TRUE)
    y_max <- max(df[[comparison_col]], ratio_human, ratio_yeast, na.rm = TRUE)
  } else {
    y_min <- y_limits[1]
    y_max <- y_limits[2]
  }

# ---- Place annotations at a consistent height (e.g., 95% of y_max) ---------
  annotation_y <- y_max * 0.95

# ---- Plot ------------------------------------------------------------------
  p <- ggplot(df, aes(x = label, y = !!sym(comparison_col), fill = Organism)) +
    geom_hline(yintercept = ratio_human, linetype = "dashed", color = "#2E6A87", size = 0.5) +
    geom_hline(yintercept = ratio_yeast, linetype = "dashed", color = "#B87F00", size = 0.5) +
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
        y = annotation_y,           # fixed relative position
        label = round(median_val, 1),
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
      y = expression(Log[2]~"[FC]"),
      fill = "Organism"
    ) +
    theme(
      plot.title = element_text(size = 6, face = "bold", hjust = 0.5, colour = "black"),
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

plot_violin_cv_prot <- function(
    long_df,
    cv_col,
    labels_to_plot = NULL,
    y_limits = NULL    
) {

# ---- Filter data -----------------------------------------------------------
  df <- long_df
  if (!is.null(labels_to_plot)) {
    df <- df %>% filter(label %in% labels_to_plot)
    df$label <- factor(df$label, levels = labels_to_plot)
  }

# ---- Determine organism ----------------------------------------------------
  df <- df %>%
    mutate(Organism = case_when(
      str_detect(Description, "OS=Homo sapiens") ~ "Human",
      str_detect(Description, "OS=Saccharomyces") ~ "Yeast",
      TRUE ~ NA_character_
    )) %>%
    filter(!is.na(Organism))

# ---- Compute medians per label & organism ----------------------------------
  medians <- df %>%
    group_by(label, Organism) %>%
    summarise(median_val = median(!!sym(cv_col), na.rm = TRUE), .groups = "drop")

# ---- Determine y-axis range (consistent across plots if provided) ----------
  if (is.null(y_limits)) {
    y_min <- min(df[[cv_col]], na.rm = TRUE)
    y_max <- max(df[[cv_col]], na.rm = TRUE)
  } else {
    y_min <- y_limits[1]
    y_max <- y_limits[2]
  }

# ---- Fixed annotation height (e.g., 95% of y_max) --------------------------
  annotation_y <- y_max * 0.96

# ---- Define colors and dodge width -----------------------------------------
  color_values <- c("Human" = "#56B4E9", "Yeast" = "#E69F00")
  dodge_width <- 0.8

# ---- Plot ------------------------------------------------------------------
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
# ---- Add median text (consistent height) -----------------------------------
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
      y = "CV [%]",
      fill = "Organism"
    ) +
    theme(
      plot.title = element_text(size = 6, face = "bold", hjust = 0.5, colour = "black"),
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


# ---- Figure drafts paper ---------------------------------------------------

# boxplots
# "normal" boxplot (H100Y100 vs H100Y020)
p_boxplot_prot_isol_w_IT_53ng_H100Y100vsH100Y020 <- plot_boxplot_vertical_prot(long_df = prot_results$long_prot_standard_StoN_grouping_20_log2,
                                                                               comparison_col = "unnorm_H100Y100vsH100Y020",
                                                                               labels_to_plot = isol_w_IT_53ng,
                                                                               ratio_human = 0,
                                                                               ratio_yeast = 2.32,
                                                                               y_limits = c(-1,5))
p_boxplot_prot_isol_w_IT_53ng_H100Y100vsH100Y020

p_boxplot_prot_isol_w_IT_53ng_H100Y100vsH100Y020 <- p_boxplot_prot_isol_w_IT_53ng_H100Y100vsH100Y020+
  scale_x_discrete(labels = c("IT 8ms", "IT 20ms", "IT 8ms", "IT 8ms", "IT 20ms", "IT 40ms", "IT 8ms")) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 6)#,axis.title.x = element_text(size = 6)
  )
p_boxplot_prot_isol_w_IT_53ng_H100Y100vsH100Y020

ggsave("plots/p_boxplot_prot_isol_w_IT_53ng_H100Y100vsH100Y020.svg", p_boxplot_prot_isol_w_IT_53ng_H100Y100vsH100Y020, 
       width = 12.7, height = 8.9, units = "cm", dpi = 300, scale = 1)

# "normal" boxplot (H100Y020 vs H100Y010)

p_boxplot_prot_isol_w_IT_53ng_H100Y020vsH100Y010 <- plot_boxplot_vertical_prot(long_df = prot_results$long_prot_standard_StoN_grouping_20_log2,
                                                                               comparison_col = "unnorm_H100Y020vsH100Y010",
                                                                               labels_to_plot = isol_w_IT_53ng,
                                                                               ratio_human = 0,
                                                                               ratio_yeast = 1,
                                                                               y_limits = c(-1,2))
p_boxplot_prot_isol_w_IT_53ng_H100Y020vsH100Y010

p_boxplot_prot_isol_w_IT_53ng_H100Y020vsH100Y010 <- p_boxplot_prot_isol_w_IT_53ng_H100Y020vsH100Y010+
  scale_x_discrete(labels = c("IT 8ms", "IT 20ms", "IT 8ms", "IT 8ms", "IT 20ms", "IT 40ms", "IT 8ms")) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 6)#,axis.title.x = element_text(size = 6)
  )
p_boxplot_prot_isol_w_IT_53ng_H100Y020vsH100Y010

ggsave("plots/p_boxplot_prot_isol_w_IT_53ng_H100Y020vsH100Y010.svg", p_boxplot_prot_isol_w_IT_53ng_H100Y020vsH100Y010, 
       width = 12.7, height = 8.9, units = "cm", dpi = 300, scale = 1)


# same boxplots but with only the 20% highest abundant proteins (H100Y100 vs H100Y020)
p_boxplot_prot_isol_w_IT_53ng_H100Y100vsH100Y020_top20 <- plot_boxplot_vertical_prot(long_df = top20_prot$unnorm_meanExp_H100Y100vsH100Y020,
                                                                               comparison_col = "unnorm_H100Y100vsH100Y020",
                                                                               labels_to_plot = isol_w_IT_53ng,
                                                                               ratio_human = 0,
                                                                               ratio_yeast = 2.32,
                                                                               y_limits = c(-1,5))
p_boxplot_prot_isol_w_IT_53ng_H100Y100vsH100Y020_top20

p_boxplot_prot_isol_w_IT_53ng_H100Y100vsH100Y020_top20 <- p_boxplot_prot_isol_w_IT_53ng_H100Y100vsH100Y020_top20+
  scale_x_discrete(labels =  c("IT 8ms", "IT 20ms", "IT 8ms", "IT 8ms", "IT 20ms", "IT 40ms", "IT 8ms")) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 6),#axis.title.x = element_text(size = 6)
  )
p_boxplot_prot_isol_w_IT_53ng_H100Y100vsH100Y020_top20

ggsave("plots/p_boxplot_prot_isol_w_IT_53ng_H100Y100vsH100Y020_top20.svg", p_boxplot_prot_isol_w_IT_53ng_H100Y100vsH100Y020_top20, 
       width = 12.7, height = 8.9, units = "cm", dpi = 300, scale = 1)


# same boxplots but with only the 20% highest abundant proteins (H100Y020 vs H100Y010)
p_boxplot_prot_isol_w_IT_53ng_H100Y020vsH100Y010_top20 <- plot_boxplot_vertical_prot(long_df = top20_prot$unnorm_meanExp_H100Y020vsH100Y010,
                                                                                     comparison_col = "unnorm_H100Y020vsH100Y010",
                                                                                     labels_to_plot = isol_w_IT_53ng,
                                                                                     ratio_human = 0,
                                                                                     ratio_yeast = 1,
                                                                                     y_limits = c(-1,2))
p_boxplot_prot_isol_w_IT_53ng_H100Y020vsH100Y010_top20

p_boxplot_prot_isol_w_IT_53ng_H100Y020vsH100Y010_top20 <- p_boxplot_prot_isol_w_IT_53ng_H100Y020vsH100Y010_top20+
  scale_x_discrete(labels =  c("IT 8ms", "IT 20ms", "IT 8ms", "IT 8ms", "IT 20ms", "IT 40ms", "IT 8ms")) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 6),#axis.title.x = element_text(size = 6)
  )
p_boxplot_prot_isol_w_IT_53ng_H100Y020vsH100Y010_top20

ggsave("plots/p_boxplot_prot_isol_w_IT_53ng_H100Y020vsH100Y010_top20.svg", p_boxplot_prot_isol_w_IT_53ng_H100Y020vsH100Y010_top20, 
       width = 12.7, height = 8.9, units = "cm", dpi = 300, scale = 1)


# violin plots
# cv_H100Y100 
p_violin_prot_isol_w_IT_53ng_H100Y100 <- plot_violin_cv_prot(long_df = prot_results$long_prot_standard_StoN_grouping_20_log2,
                                                     cv_col = "cv_H100Y100",
                                                     labels_to_plot = isol_w_IT_53ng,
                                                     y_limits = c(0,120))
p_violin_prot_isol_w_IT_53ng_H100Y100

p_violin_prot_isol_w_IT_53ng_H100Y100 <- p_violin_prot_isol_w_IT_53ng_H100Y100 +
  scale_x_discrete(labels = c("IT 8ms", "IT 20ms", "IT 8ms", "IT 8ms", "IT 20ms", "IT 40ms", "IT 8ms")) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 6) #,axis.title.x = element_text(size = 6)
  )

p_violin_prot_isol_w_IT_53ng_H100Y100

ggsave("plots/p_violin_prot_isol_w_IT_53ng_H100Y100.svg", p_violin_prot_isol_w_IT_53ng_H100Y100, 
       width = 12.7, height = 8.9, units = "cm", dpi = 300, scale = 1)


p_violin_prot_isol_w_IT_53ng_H100Y020 <- plot_violin_cv_prot(long_df = prot_results$long_prot_standard_StoN_grouping_20_log2,
                                                             cv_col = "cv_H100Y020",
                                                             labels_to_plot = isol_w_IT_53ng,
                                                             y_limits = c(0,120))
p_violin_prot_isol_w_IT_53ng_H100Y020

p_violin_prot_isol_w_IT_53ng_H100Y020 <- p_violin_prot_isol_w_IT_53ng_H100Y020 +
  scale_x_discrete(labels = c("IT 8ms", "IT 20ms", "IT 8ms", "IT 8ms", "IT 20ms", "IT 40ms", "IT 8ms")) +
  labs(x = "IT [ms]") +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 6)#,axis.title.x = element_text(size = 6)
  )

p_violin_prot_isol_w_IT_53ng_H100Y020

ggsave("plots/p_violin_prot_isol_w_IT_53ng_H100Y020.svg", p_violin_prot_isol_w_IT_53ng_H100Y020, 
       width = 12.7, height = 8.9, units = "cm", dpi = 300, scale = 1)

p_violin_prot_isol_w_IT_53ng_H100Y010 <- plot_violin_cv_prot(long_df = prot_results$long_prot_standard_StoN_grouping_20_log2,
                                                             cv_col = "cv_H100Y010",
                                                             labels_to_plot = isol_w_IT_53ng,
                                                             y_limits = c(0,120))
p_violin_prot_isol_w_IT_53ng_H100Y010

p_violin_prot_isol_w_IT_53ng_H100Y010 <- p_violin_prot_isol_w_IT_53ng_H100Y010 +
  scale_x_discrete(labels = c("IT 8ms", "IT 20ms", "IT 8ms", "IT 8ms", "IT 20ms", "IT 40ms", "IT 8ms")) +
  labs(x = "IT [ms]") +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 6)#,axis.title.x = element_text(size = 6)
  )

p_violin_prot_isol_w_IT_53ng_H100Y010

ggsave("plots/p_violin_prot_isol_w_IT_53ng_H100Y010.svg", p_violin_prot_isol_w_IT_53ng_H100Y010, 
       width = 12.7, height = 8.9, units = "cm", dpi = 300, scale = 1)
