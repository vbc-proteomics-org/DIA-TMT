# ===========================================================================
# Script:  05h_violinplots_cv_faceted.R
# Purpose: Faceted CV violin plots across the S/N titration and the software
#          comparison.
# Inputs:  all_ratios_and_cv (05f); ratio_StoN_titration /
#          ratio_software_comparison vectors (05g); condition subsets (04).
# Outputs: p_StoN_titration_cv_H100Y020 (sup_fig2),
#          p_software_comparison_cv_H100Y020 /
#          p_software_comparison_cv_H100Y010 (sup_fig4).
# Depends: Run after 05f and 05g.
# ===========================================================================

plot_violin_cv_faceted <- function(
    data,
    cv_col,
    ids = NULL,
    labels = NULL, 
    y_limit = NULL,
    facet_var = "id",
    facet_labels = NULL,
    nrcols = NULL
) {

# ---- Filter and prep -------------------------------------------------------
  plot_data <- data
  if (!is.null(ids)) plot_data <- plot_data %>% filter(id %in% ids)
  if (!is.null(labels)) plot_data <- plot_data %>% filter(label %in% labels)
  plot_data <- plot_data %>% filter(!is.na(.data[[cv_col]]))

# ---- Force facet order according to 'ids' argument -------------------------
  if (!is.null(ids)) {
    plot_data[[facet_var]] <- factor(plot_data[[facet_var]], levels = ids)
  } else {
    plot_data[[facet_var]] <- factor(plot_data[[facet_var]])
  }

# ---- Determine organism ----------------------------------------------------
  df <- plot_data %>%
    mutate(
      Organism = case_when(
        str_detect(Organism, "Human") ~ "Human",
        str_detect(Organism, "Yeast") ~ "Yeast",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(Organism))

# ---- Compute medians per label & organism & id -----------------------------
  medians <- df %>%
    group_by(id, label, Organism) %>%     
    summarise(median_val = median(.data[[cv_col]], na.rm = TRUE), .groups = "drop")

# ---- Determine y-axis range (consistent across plots if provided) ----------
  if (is.null(y_limit)) {
    y_min <- min(df[[cv_col]], na.rm = TRUE)
    y_max <- max(df[[cv_col]], na.rm = TRUE)
  } else {
    y_min <- y_limit[1]
    y_max <- y_limit[2]
  }

# ---- Fixed annotation height (e.g., 95% of y_max) --------------------------
  annotation_y <- y_max * 0.96

# ---- Define colors and dodge width -----------------------------------------
  color_values <- c("Human" = "#56B4E9", "Yeast" = "#E69F00")
  dodge_width <- 0.8

# ---- Choose facet labeler --------------------------------------------------
  if (is.null(facet_labels)) {
    labeller_used <- label_value
  } else {
    labeller_used <- labeller(!!facet_var := facet_labels)
  }

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
        label = round(median_val, 1),
        group = Organism,
        facet_var = id          
      ),
      position = position_dodge(width = dodge_width),
      angle = 45,
      vjust = 0,                   
      size = 1.7,
      color = "black",
      inherit.aes = FALSE 
    ) +
    scale_fill_manual(values = color_values) +
    scale_y_continuous(limits = c(y_min, y_max), expand = c(0, 0)) +
    theme_minimal(base_size = 10) +
    labs(
      y = "CV [%]",
      fill = "Organism"
    ) +
    facet_wrap(vars(.data[[facet_var]]),
               scales = "free_x",
               ncol = nrcols,
               labeller = labeller_used) +
    theme(
      strip.text = element_text(size = 6),
      plot.title = element_text(size = 6, face = "bold", hjust = 0.5, colour = "black"),
      axis.title.x = element_blank(),
      axis.line = element_line(color = "black"),
      axis.title.y = element_text(size = 6, colour = "black"),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6, colour = "black"),
      axis.text.y = element_text(size = 6, colour = "black"),
      legend.position = "inside",
      legend.text = element_text(size = 5.5),
      legend.key.size = unit(0.5, "cm"),
      legend.position.inside = c(0.8,0.8),
      legend.title = element_blank(),
      panel.grid = element_blank(),
      axis.ticks = element_line(color = "black", size = 0.3)
    )

  return(p)
}


# ---- Plot the S/N titration ------------------------------------------------



p_StoN_titration_cv_H100Y020 <- plot_violin_cv_faceted(
  data = all_ratios_and_cv,
  ids = ratio_StoN_titration,            
  labels = test_conditions,
  cv_col = "cv_H100Y020",
  y_limit = c(0,110),
  facet_var = "id",
  facet_labels = ratio_StoN_titration_labels,
  nrcols = 11
)

p_StoN_titration_cv_H100Y020

p_StoN_titration_cv_H100Y020 <- p_StoN_titration_cv_H100Y020 +
  scale_x_discrete(labels = rep(c("", ""), times = 7 )) +
  theme(
    strip.text = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 6)
  )

p_StoN_titration_cv_H100Y020


ggsave("plots/p_StoN_titration_cv_H100Y020.svg", device = "svg", p_StoN_titration_cv_H100Y020, 
       width = 18.3, height = 8.9, units = "cm", dpi = 300, scale = 1)


# ---- Violin plots for comparison to SN & DIANN -----------------------------


p_software_comparison_cv_H100Y020 <- plot_violin_cv_faceted(
  data = all_ratios_and_cv,
  ids = ratio_software_comparison,            
  labels = test_condition_DIA,
  cv_col = "cv_H100Y020",
  y_limit = c(0,110),
  facet_var = "id",
  facet_labels = ratio_software_comparison_labels,
  nrcols = 3
)

p_software_comparison_cv_H100Y020

p_software_comparison_cv_H100Y020 <- p_software_comparison_cv_H100Y020 +
  scale_x_discrete(labels = rep(c("DDA", "DIA"), times = 3 )) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 6)
  )

p_software_comparison_cv_H100Y020


ggsave("plots/p_software_comparison_cv_H100Y020_new.svg", device = "svg", p_software_comparison_cv_H100Y020, 
       width = 18.3, height = 8.9, units = "cm", dpi = 300, scale = 1)

p_software_comparison_cv_H100Y010 <- plot_violin_cv_faceted(
  data = all_ratios_and_cv,
  ids = ratio_software_comparison,            
  labels = test_condition_DIA,
  cv_col = "cv_H100Y010",
  y_limit = c(0,110),
  facet_var = "id",
  facet_labels = ratio_software_comparison_labels,
  nrcols = 3
)

p_software_comparison_cv_H100Y010

p_software_comparison_cv_H100Y010 <- p_software_comparison_cv_H100Y010 +
  scale_x_discrete(labels = rep(c("DDA", "DIA"), times = 3 )) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 6)
  )

p_software_comparison_cv_H100Y010


ggsave("plots/p_software_comparison_cv_H100Y010_new.svg", device = "svg", p_software_comparison_cv_H100Y010, 
       width = 18.3, height = 8.9, units = "cm", dpi = 300, scale = 1)