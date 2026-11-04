# ===========================================================================
# Script:  05c_plots_prot_faceted.R
# Purpose: Faceted protein-count plots across the S/N titration series.
# Inputs:  loss_annotations_all (03d), condition subsets (04).
# Outputs: p_prot_StoN_titration (sup_fig2).
# Depends: Run after 03d and 04.
# ===========================================================================

plot_protcount_faceted <- function(data,
                                   ids = NULL,
                                   labels = NULL,
                                   replicate_cols = c("Rep1_F", "Rep2_F", "Rep3_F"),
                                   y_limit = NULL,
                                   facet_var = "id",
                                   facet_labels = NULL,
                                   nrcols = NULL) {   

# ---- Filter and prep -------------------------------------------------------
  plot_data <- data
  if (!is.null(ids)) plot_data <- plot_data %>% filter(id %in% ids)
  if (!is.null(labels)) plot_data <- plot_data %>% filter(label %in% labels)
  plot_data <- plot_data %>% filter(!is.na(Quantified))

# ---- Force facet order according to 'ids' argument -------------------------
  if (!is.null(ids)) {
    plot_data[[facet_var]] <- factor(plot_data[[facet_var]], levels = ids)
  } else {
    plot_data[[facet_var]] <- factor(plot_data[[facet_var]])
  }

# ---- Handle replicates -----------------------------------------------------
  replicate_cols_present <- intersect(replicate_cols, names(plot_data))
  plot_data_long <- plot_data %>%
    pivot_longer(cols = all_of(replicate_cols_present),
                 names_to = "Replicate",
                 values_to = "Replicate_Value") %>%
    filter(!is.na(Replicate_Value))

# ---- Colors ----------------------------------------------------------------
  all_labels <- unique(plot_data$label)
  fallback_colors <- RColorBrewer::brewer.pal(12, "Set3")
  label_colors <- rep(fallback_colors, length.out = length(all_labels))
  names(label_colors) <- all_labels

# ---- Y limit ---------------------------------------------------------------
  if (is.null(y_limit)) y_limit <- max(plot_data$Quantified, na.rm = TRUE) * 1.2

# ---- Facet layout ----------------------------------------------------------
  n_ids <- length(unique(plot_data[[facet_var]]))
  ncol_auto <- ifelse(n_ids <= 3, n_ids, 3)

# ---- Labeller setup --------------------------------------------------------
  labeller_used <- if (!is.null(facet_labels)) {
    if (is.function(facet_labels)) {
      facet_labels
    } else {
      as_labeller(facet_labels)
    }
  } else {
    "label_value"
  }

# ---- Plot ------------------------------------------------------------------
  p <- ggplot(plot_data, aes(x = label, y = Quantified, fill = label)) +
    geom_bar(stat = "identity", color = "black", width = 0.5, linewidth = 0.2) +
    geom_errorbar(aes(ymin = pmax(0, Quantified - SD_f), ymax = Quantified + SD_f),
                  width = 0.2, color = "black", size = 0.2) +
    geom_point(
      data = plot_data_long,
      aes(y = Replicate_Value, fill = label),
      shape = 21, color = "black", size = 0.7, alpha = 0.7,
      position = position_jitter(width = 0.1),
      show.legend = FALSE
    ) +
    geom_text(aes(y = Quantified + (y_limit * 0.02), label = round(Quantified, 0)),
              size = 1.7, angle = 90, vjust = 0.5, hjust = -0.2) +
    scale_fill_manual(values = label_colors#, guide = "none"
                      ) +
    scale_y_continuous(limits = c(0, y_limit), expand = c(0,0)) +
    facet_wrap(vars(.data[[facet_var]]),
               scales = "free_x",
               ncol = nrcols,
               labeller = labeller_used) +
    labs(x = NULL, y = "#quantified Proteins") +
    theme_minimal(base_size = 10) +
    theme(
      strip.text = element_text(size = 5#, margin = margin(0.4, 0.4, 0.4, 0.4, "pt")
                                ),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6, color = "black"),
      axis.text.y = element_text(size = 6, color = "black"),
      axis.title.y = element_text(size = 6, color = "black"),
      panel.grid = element_blank(),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(size = 0.3, color = "black")
    )

  return(p)
}

# Define the data to plot together and the facet headings for the resulting plot
# S/N titration
IDs_StoN_titration <- c(
  "long_prot_baseline",
  "long_prot_baseline_70perc",
  "long_prot_StoN_grouping_3",
  "long_prot_StoN_grouping_6",
  "long_prot_StoN_grouping_9",
  "long_prot_StoN_grouping_12",
  "long_prot_StoN_grouping_15",
  "long_prot_StoN_grouping_20",
  "long_prot_StoN_grouping_30",
  "long_prot_StoN_grouping_40",
  "long_prot_StoN_grouping_80"
)

IDs_StoN_titration_labels <- c(
  "long_prot_baseline" = "No Filter",
  "long_prot_baseline_70perc" = "70% purity",
  "long_prot_StoN_grouping_3" = "S/N=3",
  "long_prot_StoN_grouping_6" = "S/N=6",
  "long_prot_StoN_grouping_9" = "S/N=9",
  "long_prot_StoN_grouping_12" = "S/N=12",
  "long_prot_StoN_grouping_15" = "S/N=15",
  "long_prot_StoN_grouping_20" = "S/N=20",
  "long_prot_StoN_grouping_30" = "S/N=30",
  "long_prot_StoN_grouping_40" = "S/N=40",
  "long_prot_StoN_grouping_80" = "S/N=80"
)

# ---- Plot the S/N titration in 11 facets -----------------------------------

p_prot_StoN_titration <- plot_protcount_faceted(
  data = loss_annotations_all,
  ids = IDs_StoN_titration,            
  labels = test_conditions,
  replicate_cols = c("Rep1_F", "Rep2_F", "Rep3_F"),
  y_limit = 12000,
  facet_labels = IDs_StoN_titration_labels,
  nrcols = 11
)

p_prot_StoN_titration

p_prot_StoN_titration <- p_prot_StoN_titration +
  scale_x_discrete(labels = rep("", 22 )) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 6),
    axis.title.x = element_text(size = 6)
  )

p_prot_StoN_titration

p_prot_StoN_titration <- p_prot_StoN_titration +
  scale_fill_manual(
    values = c(
      "DDA 0.5 m/z \nIT8, 53.75ng" = "#5B859EFF",
      "DIA 1.0 m/z \nIT8, 53.75ng" = "#D48F90FF"
    ),
    labels = c(
      "DDA 0.5 m/z \nIT8, 53.75ng" = "DDA 0.5 m/z IT 8 ms",
      "DIA 1.0 m/z \nIT8, 53.75ng" = "DIA 1.0 m/z IT 8 ms"
    )
  ) +
  theme(
    strip.text = element_blank(),
    legend.position = "inside",
    legend.text = element_text(size = 6),
    legend.position.inside = c(0.8, 0.8),
    legend.title = element_blank(),
    legend.key.size = unit(5, "pt"),
    legend.key.spacing.y = unit(1, "pt")
  )
p_prot_StoN_titration

ggsave("plots/p_prot_StoN_titration.svg", device = "svg", p_prot_StoN_titration, 
       width = 8, height = 6, units = "cm", dpi = 300, scale = 1)
