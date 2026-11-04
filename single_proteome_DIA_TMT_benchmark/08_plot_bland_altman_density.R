# =============================================================================
# Script 08: Bland-Altman plots with density panels
# Description: Defines plot_blant_altmann_prot(), which draws a Bland-Altman
#              scatter of log2 fold-change vs. log2 abundance with a side density
#              panel per organism, and builds the per-method BA figures. The two
#              "_StoN_grouping" figures feed the final arranged figure.
# Depends on: 03_filter_transform_and_counts.R
# =============================================================================

# --- Bland-Altman + density plot function ------------------------------------

plot_blant_altmann_prot <- function(
    long_df,
    comparison_col,
    dataset_filter = NULL,
    ratio_hs = NULL,
    ratio_sc = NULL,
    show_yeast = TRUE,
    y_limits = c(-4, 5),
    x_limits = c(5, 20),
    title = NULL,
    hline_hs = 0,               # Fold-change reference line(s) for Homo sapiens.
    hline_sc = 2.32,            # Fold-change reference line(s) for Saccharomyces.
    hline_color_hs = "#2E6A87",
    hline_color_sc = "#B87F00"
) {

  # ---- Filter dataset if requested ----
  df <- long_df
  if (!is.null(dataset_filter)) {
    df <- df %>% filter(dataset == dataset_filter)
  }

  # ---- Determine organism from description ----
  df <- df %>%
    mutate(Organism = case_when(
      str_detect(`Organism`, "Human") ~ "Homo sapiens",
      str_detect(`Organism`, "Yeast") ~ "Saccharomyces",
      TRUE ~ "Other"
    ))

  # ---- Optionally remove yeast ----
  if (!show_yeast) {
    df <- df %>% filter(Organism != "Saccharomyces")
  }

  # ---- Define columns to plot ----
  ratio_col <- comparison_col
  meanexp_col <- if (paste0("meanExp_", comparison_col) %in% names(df)) {
    paste0("meanExp_", comparison_col)
  } else if (paste0("unnorm_meanExp_", gsub("^unnorm_", "", comparison_col)) %in% names(df)) {
    paste0("unnorm_meanExp_", gsub("^unnorm_", "", comparison_col))
  } else {
    stop("No matching meanExp column found for ", comparison_col)
  }

  # ---- Keep only rows with valid numeric data ----
  df <- df %>%
    filter(!is.na(.data[[ratio_col]]), !is.na(.data[[meanexp_col]]))

  # ---- Count number of actually plotted proteins per organism ----
  counts <- df %>%
    group_by(Organism) %>%
    summarise(n = n(), .groups = "drop")

  n_hs <- counts$n[counts$Organism == "Homo sapiens"]
  if (length(n_hs) == 0) n_hs <- 0
  n_sc <- counts$n[counts$Organism == "Saccharomyces"]
  if (length(n_sc) == 0) n_sc <- 0

  # ---- Define color scale ----
  color_values <- c("Homo sapiens" = "#56B4E9")
  if (show_yeast) color_values["Saccharomyces"] <- "#E69F00"

  # ---- Define dynamic legend labels ----
  labels_values <- c(
    "Homo sapiens" = paste0("Human: ", n_hs,
                    if (!is.null(ratio_hs)) paste0(" (", ratio_hs, ")") else "")
  )
  if (show_yeast) {
    labels_values["Saccharomyces"] <- paste0("Yeast: ", n_sc,
                    if (!is.null(ratio_sc)) paste0(" (", ratio_sc, ")") else "")
  }

  # ---- Base scatter plot ----
  p_scatter <- ggplot(df, aes_string(x = meanexp_col, y = ratio_col, color = "Organism")) +
    geom_point(alpha = 0.4, size = 0.05) +
    scale_color_manual(values = color_values, labels = labels_values, name = NULL) +
    labs(title = if (!is.null(title)) paste0("BA + density plot: ", title)
         else "BA + density plot",
      y = expression(Log[2]~"[FC]"),
      x = expression(Log[2]~"[abundance]")
    ) +
    scale_y_continuous(limits = y_limits, breaks = seq(y_limits[1], y_limits[2], by = 1)) +
    scale_x_continuous(limits = x_limits, breaks = seq(x_limits[1], x_limits[2], by = 2)) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title   = element_text(lineheight = .8, face = "plain", hjust = 0.5, size = 6),
      axis.line = element_line(color = "black"),
      axis.text.x = element_text(size = 6, colour = "black"),
      axis.text.y = element_text(size = 6, colour = "black"),
      axis.title.x = element_text(size = 6, face = "bold", colour = "black"),
      axis.title.y = element_text(size = 6, face = "bold", colour = "black"),
      legend.position = "inside",
      legend.justification = c(0,1),
      legend.direction = "vertical",
      legend.text = element_text(size = 6, colour = "black"),
      panel.grid = element_blank(),
      axis.ticks = element_line(color = "black", size = 0.25)
    ) +
    guides(color = guide_legend(override.aes = list(size = 1.7)))

  # ---- Add horizontal lines for Homo sapiens if provided ----
  if (!is.null(hline_hs) && length(hline_hs) > 0) {
    for (v in as.numeric(hline_hs)) {
      p_scatter <- p_scatter + geom_hline(
        yintercept = v,
        linetype = "dashed",
        color = hline_color_hs,
        size = 0.5
      )
    }
  }

  # ---- Add horizontal lines for Saccharomyces if provided and shown ----
  if (show_yeast && !is.null(hline_sc) && length(hline_sc) > 0) {
    for (v in as.numeric(hline_sc)) {
      p_scatter <- p_scatter + geom_hline(
        yintercept = v,
        linetype = "dashed",
        color = hline_color_sc,
        size = 0.5
      )
    }
  }

  # ---- Calculate mean & SD per organism ----
  sd_summary <- df %>%
    group_by(Organism) %>%
    summarise(
      mean = mean(.data[[ratio_col]], na.rm = TRUE),
      sd = sd(.data[[ratio_col]], na.rm = TRUE),
      median = median(.data[[ratio_col]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      FC_linear = 2^mean,
      FC_log2   = mean
    )

  # ---- Density panel (right) ----
  p_density <- ggplot(df, aes_string(x = ratio_col, fill = "Organism", color = "Organism")) +
    geom_density(alpha = 0.3, fill = NA, linewidth = 0.25) +
    geom_vline(data = sd_summary, aes(xintercept = mean, color = Organism),
               linetype = "dotted", linewidth = 0.3) +
    geom_text(
      data = sd_summary,
      aes(
        x = mean + 0.3,
        y = 1.25,
        label = paste0("FC = ", round(FC_linear, 2), "\nSD = ", round(sd, 2))
      ),
      colour = "black",
      size = 1.7,
      hjust = 0.5,
      vjust = 0,
      angle = 0,
      inherit.aes = FALSE
    ) +
    labs(
      y = "Density") +
    scale_fill_manual(values = color_values) +
    scale_color_manual(values = color_values) +
    coord_flip() +
    scale_x_continuous(
      limits = y_limits,
      breaks = seq(y_limits[1], y_limits[2], by = 1)  # Match BA y-axis breaks.
    ) +
    scale_y_continuous(expand = c(0,0.02), limits = c(0,4), breaks = seq(0,4, 1)) +
    theme_minimal(base_size = 10) +
    theme(
      legend.position = "none",
      axis.title = element_blank(),
      axis.ticks = element_line(color = "black", size = 0.25),
      panel.grid = element_blank(),
      axis.text.x = element_text(size = 6, colour = "black"),
      axis.text.y = element_blank(),
      axis.title.x = element_text(size = 6, colour = "black"),
      axis.line = element_line(color = "black")
    )

  # ---- Combine scatter + density ----
  p_final <- p_scatter + p_density + plot_layout(widths = c(4, 1))

  # ---- Return both plot and SD table ----
  return(list(
    plot = p_final,
    sd_summary = sd_summary
  ))
}

# --- Figures for the paper (feed script 11) ----------------------------------

p_prot_ratio_DDA_IT20_54ng_rep1_H100vsH010_StoN_grouping <- plot_blant_altmann_prot(prot_results$long_prot_standard_StoN_grouping_20_log2,
                    comparison_col = "unnorm_H100vsH010",
                    dataset_filter = "20260214_Neo6A1_TH136_TT_THIDmulti008A_DDA_mz0c5_ij20_53c75ng_1_Proteins_purity_tmt_quant",
                    ratio_hs = "100/10",
                    ratio_sc = "0/10",
                    show_yeast = F,
                    x_limits = c(4,20),
                    y_limits = c(0,6),
                    hline_hs = 3.32,
                    hline_sc = NULL,
                    title = "DDA")

p_prot_ratio_DDA_IT20_54ng_rep1_H100vsH010_StoN_grouping


ggsave("plots/p_prot_ratio_DDA_IT20_54ng_rep1_H100vsH010_StoN_grouping.svg", p_prot_ratio_DDA_IT20_54ng_rep1_H100vsH010_StoN_grouping$plot,
       width = 12.7, height = 8.9, units = "cm", dpi = 300
)

p_prot_ratio_DIA_IT8_HCD29_54ng_rep1_H100vsH010_StoN_grouping <- plot_blant_altmann_prot(prot_results$long_prot_standard_StoN_grouping_20_log2,
                    comparison_col = "unnorm_H100vsH010",
                    dataset_filter = "20260214_Neo6A1_TH136_TT_THIDmulti008A_DIA_mz1_ij8_HCD29_53c75ng_1_Proteins_purity_tmt_quant",
                    ratio_hs = "100/10",
                    ratio_sc = "0/10",
                    show_yeast = F,
                    x_limits = c(4,20),
                    y_limits = c(0,6),
                    hline_hs = 3.32,
                    hline_sc = NULL,
                    title = "DIA HCD 29%")

p_prot_ratio_DIA_IT8_HCD29_54ng_rep1_H100vsH010_StoN_grouping


ggsave("plots/p_prot_ratio_DIA_IT8_HCD29_54ng_rep1_H100vsH010_StoN_grouping.svg", p_prot_ratio_DIA_IT8_HCD29_54ng_rep1_H100vsH010_StoN_grouping$plot,
       width = 12.7, height = 8.9, units = "cm", dpi = 300
)
