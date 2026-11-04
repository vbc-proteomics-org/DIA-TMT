# ===========================================================================
# Script:  06_plots_psms_coisolation.R
# Purpose: Co-isolation purity distribution ridge plots annotated with the
#          percentage of precursors above 0.7.
# Inputs:  psm_results$long_psm_baseline (03a), label_colors (04), condition
#          subsets (04).
# Outputs: p_pur_dist_input_tit, p_pur_dist_isol_w_IT_53ng (Fig2).
# Depends: Run after 04. BUG: desired_order is used before it is defined (see
#          inline BUG flag).
# ===========================================================================

# ---- Plot the distribution of purity values for the datasets ---------------

data_coisolation <- psm_results$long_psm_baseline

# force the order of datasets
desired_order <- c(
  "DIA 0.5 m/z \nIT8, 3ng",
  "DIA 0.5 m/z \nIT8, 53.75ng",
  "DIA 0.5 m/z \nIT8, 100ng",
  
  "DIA 1.0 m/z \nIT8, 3ng",
  "DIA 1.0 m/z \nIT8, 53.75ng",
  "DIA 1.0 m/z \nIT8, 100ng",
  
  "DIA 1.0 m/z \nIT20, 3ng",
  "DIA 1.0 m/z \nIT20, 53.75ng",
  "DIA 1.0 m/z \nIT20, 100ng",
  
  "DIA 1.0 m/z \nIT40, 3ng",
  "DIA 1.0 m/z \nIT40, 5ng",
  "DIA 1.0 m/z \nIT40, 10ng",
  "DIA 1.0 m/z \nIT40, 20ng",
  "DIA 1.0 m/z \nIT40, 25ng",
  "DIA 1.0 m/z \nIT40, 30ng",
  "DIA 1.0 m/z \nIT40, 40ng",
  "DIA 1.0 m/z \nIT40, 53.75ng",
  "DIA 1.0 m/z \nIT40, 100ng",
  
  "DIA 2.0 m/z \nIT8, 3ng",
  "DIA 2.0 m/z \nIT8, 53.75ng",
  "DIA 2.0 m/z \nIT8, 100ng",
  
  "DDA 0.5 m/z \nIT8, 3ng",
  "DDA 0.5 m/z \nIT8, 53.75ng",
  "DDA 0.5 m/z \nIT8, 100ng",
  
  "DDA 0.5 m/z \nIT20, 3ng",
  "DDA 0.5 m/z \nIT20, 53.75ng",
  "DDA 0.5 m/z \nIT20, 100ng"
)



data_coisolation <- data_coisolation %>%
  mutate(label = factor(label, levels = desired_order))



# ---- Purity distribution with the added % of precursors above 0.7 ----------

# compute % above 0.7 for each label
perc_over0c7 <- data_coisolation %>%
  group_by(label, replicate) %>%
  summarise(
    perc_above_0.7_rep = mean(`Co-Isolation Purity` >= 0.7, na.rm = TRUE) * 100,
    total_psm_rep = n(),
    psm_over70_rep = sum(`Co-Isolation Purity` >= 0.7, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(label) %>%
  summarise(
    perc_above_0.7 = mean(perc_above_0.7_rep, na.rm = TRUE),
    total_psm = mean(total_psm_rep, na.rm = TRUE),  # or sum(), depending on what you want
    psm_over70 = mean(psm_over70_rep, na.rm = TRUE)
  )


# Plotting function for the purity distribution

plot_coisolation_with_percent <- function(data_coisolation, 
                                          perc_over0c7,
                                          selected_labels = NULL, 
                                          label_colors) {
  # Subset data if requested
  plot_data <- if (!is.null(selected_labels)) {
    data_coisolation %>%
      dplyr::filter(label %in% selected_labels)
  } else {
    data_coisolation
  }

  if (nrow(plot_data) == 0) {
    stop("No data found for the selected labels.")
  }

  # Subset perc_over0c7 to match selected labels
  perc_data <- if (!is.null(selected_labels)) {
    perc_over0c7 %>%
      dplyr::filter(label %in% selected_labels)
  } else {
    perc_over0c7
  }

# ---- Force consistent label order ------------------------------------------
  if (!is.null(selected_labels)) {
    plot_data <- plot_data %>%
      dplyr::mutate(label = factor(label, levels = selected_labels))

    perc_data <- perc_data %>%
      dplyr::mutate(label = factor(label, levels = selected_labels))
  } else {
    # fallback: preserve the order from the data
    plot_data <- plot_data %>%
      dplyr::mutate(label = factor(label, levels = unique(label)))

    perc_data <- perc_data %>%
      dplyr::mutate(label = factor(label, levels = levels(plot_data$label)))
  }

# ---- Build the plot --------------------------------------------------------
  p <- ggplot(plot_data, aes(x = `Co-Isolation Purity`, y = label, fill = label)) + 
    annotate(
      "rect",
      xmin = 0.7, xmax = 1.2,
      ymin = -Inf, ymax = Inf,
      fill = "#0073c2ff", alpha = 0.2
    ) +
    ggridges::geom_density_ridges(rel_min_height = 0.009,
                                  scale = 0.99,
                                  linewidth = 0.15
                                  ) +
    scale_fill_manual(values = label_colors) +
    labs(
      x = "Co-Isolation Purity",
      y = "Input [ng]"
    ) +
    scale_x_continuous(breaks = seq(0, 1, by = 0.1)) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 6, face = "bold", colour = "black"),
      axis.line = element_line(color = "black"),
      axis.text.x = element_text(size = 6, colour = "black"),
      axis.text.y = element_text(size = 6, colour = "black"),
      axis.title.x = element_text(size = 6, colour = "black"),
      axis.title.y = element_text(size = 6, colour = "black"),
      legend.position = "none",
      plot.margin = margin(5.5, 5.5, 5.5, 5.5),
      panel.grid = element_blank(),
      axis.ticks = element_line(color = "black", size = 0.3)
    ) +
    geom_text(
      data = perc_data,
      aes(
        x = 0.25,  # text inside the plot
        y = label,
        label = sprintf("total precursors: %.0f",
                        total_psm)
      ),
      hjust = 0,
      vjust = -1,
      size = 1.5,
      color = "black"
    ) +
    geom_text(
      data = perc_data,
      aes(
        x = 0.72,  # text inside the plot
        y = label,
        label = sprintf("%.0f (%.1f%%)",
                        psm_over70, perc_above_0.7)
      ),
      hjust = 0,
      vjust = -1,
      size = 1.5,
      color = "black"
    ) +
    geom_vline(xintercept = 0.7,
               linetype = "dotted") +
    coord_cartesian(xlim = c(0, 1), clip = "off")  # fully inside 


  return(p)
}


# ---- Figures paper draft ---------------------------------------------------

p_pur_dist_input_tit <- plot_coisolation_with_percent(data_coisolation,
                                                           perc_over0c7, 
                                                           selected_labels = input_titration,
                                                           label_colors = label_colors)

p_pur_dist_input_tit

p_pur_dist_input_tit <- p_pur_dist_input_tit +
  scale_y_discrete(labels = c("3", "5", "10", "20", "25", "30","40", "54", "100" )) 
p_pur_dist_input_tit

p_pur_dist_input_tit <- p_pur_dist_input_tit +
  scale_fill_manual(values = c(
    "DIA 1.0 m/z \nIT40, 3ng" = "#FFF5EB",
    "DIA 1.0 m/z \nIT40, 5ng" = "#FEE6CE",
    "DIA 1.0 m/z \nIT40, 10ng" = "#FDD0A2",
    "DIA 1.0 m/z \nIT40, 20ng" = "#FDAE6B",
    "DIA 1.0 m/z \nIT40, 25ng" = "#FD8D3C",
    "DIA 1.0 m/z \nIT40, 30ng" = "#F16913",
    "DIA 1.0 m/z \nIT40, 40ng" = "#D94801",
    "DIA 1.0 m/z \nIT40, 53.75ng" = "#A63603",
    "DIA 1.0 m/z \nIT40, 100ng" = "#7F2704"
  )) +
  theme(legend.position = "none")

p_pur_dist_input_tit

ggsave("plots/p_pur_dist_input_titration.svg", p_pur_dist_input_tit, 
       width = 12.7, height = 8.9, units = "cm", dpi = 300#, scale = 1
       )

p_pur_dist_isol_w_IT_53ng <- plot_coisolation_with_percent(data_coisolation,
                                                           perc_over0c7, 
                                                           selected_labels = isol_w_IT_53ng,
                                                           label_colors = label_colors)

p_pur_dist_isol_w_IT_53ng

p_pur_dist_isol_w_IT_53ng <- p_pur_dist_isol_w_IT_53ng +
  scale_y_discrete(labels = c("IT 8ms", "IT 20ms", "IT 8ms", "IT 8ms", "IT 20ms", "IT 40ms", "IT 8ms")) +
  theme(axis.title.y = element_blank())
p_pur_dist_isol_w_IT_53ng

p_pur_dist_isol_w_IT_53ng <- p_pur_dist_isol_w_IT_53ng +
  scale_fill_manual(values = c(
    "DDA 0.5 m/z \nIT8, 53.75ng" = "#5B859EFF",
    "DDA 0.5 m/z \nIT20, 53.75ng" = "#1E395FFF",
    "DIA 0.5 m/z \nIT8, 53.75ng" = "#1E5A46FF",
    "DIA 1.0 m/z \nIT8, 53.75ng" = "#D48F90FF",
    "DIA 1.0 m/z \nIT20, 53.75ng" = "#DF8D71FF",
    "DIA 1.0 m/z \nIT40, 53.75ng" = "#AB84A5FF",
    "DIA 2.0 m/z \nIT8, 53.75ng" = "#D8B847FF"
  )) +
  theme(legend.position = "none")

p_pur_dist_isol_w_IT_53ng

ggsave("plots/p_pur_dist_isol_w_IT_53ng.svg", p_pur_dist_isol_w_IT_53ng, 
       width = 12.7, height = 8.9, units = "cm", dpi = 300, scale = 1)