# ===========================================================================
# Script:  05g_boxplots_faceted.R
# Purpose: Faceted ratio boxplots across the S/N titration and the software
#          comparison.
# Inputs:  all_ratios_and_cv (05f), condition subsets (04).
# Outputs: ratio_StoN_titration / ratio_software_comparison selector vectors
#          and p_StoN_titration_boxplot_H100Y020vsH100Y010 (sup_fig2),
#          p_software_comparison_boxplot_H100Y020vsH100Y010 (sup_fig4).
# Depends: Run after 05f. NOTE: also defines ratio_StoN_titration(_labels) and
#          ratio_software_comparison(_labels), which 05h reuses.
# ===========================================================================

plot_boxplot_faceted <- function(
    data,
    ids = NULL,
    labels = NULL, 
    comparison_col,
    y_limit = NULL,
    facet_var = "id",
    facet_labels = NULL,
    nrcols = NULL,
    ratio_human = 0,
    ratio_yeast = 2.32) {

# ---- Filter and prep -------------------------------------------------------
  plot_data <- data
  if (!is.null(ids)) plot_data <- plot_data %>% filter(id %in% ids)
  if (!is.null(labels)) plot_data <- plot_data %>% filter(label %in% labels)
  plot_data <- plot_data %>% filter(!is.na(.data[[comparison_col]]))

# ---- Force facet order according to 'ids' argument -------------------------
  if (!is.null(ids)) {
    plot_data[[facet_var]] <- factor(plot_data[[facet_var]], levels = ids)
  } else {
    plot_data[[facet_var]] <- factor(plot_data[[facet_var]])
  }

# ---- Determine organism ----
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
    summarise(median_val = median(.data[[comparison_col]], na.rm = TRUE), .groups = "drop")


# ---- Define colors ---------------------------------------------------------
  color_values <- c("Human" = "#56B4E9", "Yeast" = "#E69F00")
  dodge_width <- 0.7

# ---- Determine y-axis range ----
  if (is.null(y_limit)) {
    y_min <- min(df[[comparison_col]], ratio_human, ratio_yeast, na.rm = TRUE)
    y_max <- max(df[[comparison_col]], ratio_human, ratio_yeast, na.rm = TRUE)
  } else {
    y_min <- y_limit[1]
    y_max <- y_limit[2]
  }

  annotation_y <- y_max * 0.95

# ---- Choose facet labeler --------------------------------------------------
  if (is.null(facet_labels)) {
    labeller_used <- label_value
  } else {
    labeller_used <- labeller(!!facet_var := facet_labels)
  }

# ---- Plot ------------------------------------------------------------------
  p <- ggplot(df, aes(x = label, y = .data[[comparison_col]], fill = Organism)) +
    geom_hline(yintercept = ratio_human, linetype = "dashed", color = "#2E6A87", size = 0.5) +
    geom_hline(yintercept = ratio_yeast, linetype = "dashed", color = "#B87F00", size = 0.5) +
    geom_boxplot(width = 0.5, 
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
      y = expression(Log[2]~"[FC]"),
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
      axis.title.y = element_text(size = 6, face = "bold", colour = "black"),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6, colour = "black"),
      axis.text.y = element_text(size = 6, colour = "black"),
      legend.position = "inside", 
      legend.position.inside = c(0.1,0.1),
      legend.text = element_text(size = 5.5),
      legend.key.size = unit(0.5, "cm"),
      legend.title = element_blank(),
      panel.grid = element_blank(),
      axis.ticks = element_line(color = "black", size = 0.25)
    )

  return(p)
}

# Define the data to plot together (using the 'id' column from all_ratios_and_cv)
# These conditions are also reused for the faceted violin CV plots (05h)


# purity titration
ratio_StoN_titration <- c(
  "baseline",
  "baseline_70perc",
  "StoN_grouping_3",
  "StoN_grouping_6",
  "StoN_grouping_9",
  "StoN_grouping_12",
  "StoN_grouping_15",
  "StoN_grouping_20",
  "StoN_grouping_30",
  "StoN_grouping_40",
  "StoN_grouping_80"
)

ratio_StoN_titration_labels <- c(
  "baseline" = "No Filter",
  "baseline_70perc" = "70% purity",
  "StoN_grouping_3" = "S/N=3",
  "StoN_grouping_6" = "S/N=6",
  "StoN_grouping_9" = "S/N=9",
  "StoN_grouping_12" = "S/N=12",
  "StoN_grouping_15" = "S/N=15",
  "StoN_grouping_20" = "S/N=20",
  "StoN_grouping_30" = "S/N=30",
  "StoN_grouping_40" = "S/N=40",
  "StoN_grouping_80" = "S/N=80"
)

# comparison SN & DIANN

ratio_software_comparison <- c(
  "StoN_grouping_20",
  "DIANN_StoN_grouping_20",
  "SNaut_StoN_grouping_20"
)

ratio_software_comparison_labels <- c(
  "StoN_grouping_20" = "PD",
  "DIANN_StoN_grouping_20" = "DIA-NN",
  "SNaut_StoN_grouping_20" = "SN"
)


# ---- Plot the S/N titration in 11 facets -----------------------------------

# H100Y020 vs. H100Y010
p_StoN_titration_boxplot_H100Y020vsH100Y010 <- plot_boxplot_faceted(
  data = all_ratios_and_cv,
  ids = ratio_StoN_titration,            
  labels = test_conditions,
  comparison_col = "unnorm_H100Y020vsH100Y010",
  y_limit = c(-1,1.5),
  facet_var = "id",
  facet_labels = ratio_StoN_titration_labels,
  nrcols = 11,
  ratio_human = 0,
  ratio_yeast = 1
)

p_StoN_titration_boxplot_H100Y020vsH100Y010

p_StoN_titration_boxplot_H100Y020vsH100Y010 <- p_StoN_titration_boxplot_H100Y020vsH100Y010 +
  scale_x_discrete(labels = rep(c("DDA", "DIA"), times = 7 )) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 6)
  )

p_StoN_titration_boxplot_H100Y020vsH100Y010

ggsave("plots/p_StoN_titration_boxplot_H100Y020vsH100Y010.svg", device = "svg", 
       p_StoN_titration_boxplot_H100Y020vsH100Y010, 
       width = 18.3, height = 8.9, units = "cm", dpi = 300, scale = 1)

# ---- Plot the comparison to SN & DIANN in 3 facets -------------------------


# H100Y020 vs. H100Y010
p_software_comparison_boxplot_H100Y020vsH100Y010 <- plot_boxplot_faceted(
  data = all_ratios_and_cv,
  ids = ratio_software_comparison,            
  labels = test_condition_DIA,
  comparison_col = "unnorm_H100Y020vsH100Y010",
  y_limit = c(-1,2),
  facet_var = "id",
  facet_labels = ratio_software_comparison_labels,
  nrcols = 3,
  ratio_human = 0,
  ratio_yeast = 1
)

p_software_comparison_boxplot_H100Y020vsH100Y010


p_software_comparison_boxplot_H100Y020vsH100Y010 <- p_software_comparison_boxplot_H100Y020vsH100Y010 +
  scale_x_discrete(labels = c("PD", "DIA-NN", "Spectronaut")) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 6)
  )

p_software_comparison_boxplot_H100Y020vsH100Y010


ggsave("plots/p_software_comparison_boxplot_H100Y020vsH100Y010_new.svg", device = "svg", 
       p_software_comparison_boxplot_H100Y020vsH100Y010, 
       width = 18.3, height = 8.9, units = "cm", dpi = 300, scale = 1)