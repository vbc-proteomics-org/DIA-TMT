# ===========================================================================
# Script:  05d_plots_prot_BA_plots_with_density.R
# Purpose: Bland-Altman (ratio vs mean) protein plots with marginal density;
#          defines plot_blant_altmann_prot and batch-saves the BA plots.
# Inputs:  prot_results (03a), condition subsets (04).
# Outputs: p_prot_ratio_DDA_IT8_rep1_H100Y000vsH100Y010_wo_yeast and the
#          p_prot_BA_* objects (Fig1 / Fig3 / sup_fig3), plus PNGs written to
#          disk.
# Depends: Run after 03a and 04.
# ===========================================================================

plot_blant_altmann_prot <- function(
    long_df,
    comparison_col,
    dataset_filter = NULL,
    ratio_hs = NULL,
    ratio_sc = NULL,
    show_yeast = TRUE,
    y_limits = c(-4, 5),
    x_limits = c(5, 20),
    density_limits = c(0, 2),   # <-- new: controls density (y) axis range in density panel
    density_breaks = 1,         # <-- new: step size for breaks
    hline_hs = 0,
    hline_sc = NULL,
    hline_color_hs = "#2E6A87",
    hline_color_sc = "#B87F00",
    title_suffix = NULL
) {

# ---- Filter dataset if requested -------------------------------------------
  df <- long_df
  if (!is.null(dataset_filter)) {
    df <- df %>% filter(dataset == dataset_filter)
  }

# ---- Determine organism from description -----------------------------------
  df <- df %>%
    mutate(Organism = case_when(
      str_detect(`Organism`, "Human") ~ "Homo sapiens",
      str_detect(`Organism`, "Yeast") ~ "Saccharomyces",
      TRUE ~ "Other"
    ))

# ---- Optionally remove yeast -----------------------------------------------
  if (!show_yeast) {
    df <- df %>% filter(Organism != "Saccharomyces")
  }

# ---- Define columns to plot ------------------------------------------------
  ratio_col <- comparison_col
  meanexp_col <- if (paste0("meanExp_", comparison_col) %in% names(df)) {
    paste0("meanExp_", comparison_col)
  } else if (paste0("unnorm_meanExp_", gsub("^unnorm_", "", comparison_col)) %in% names(df)) {
    paste0("unnorm_meanExp_", gsub("^unnorm_", "", comparison_col))
  } else {
    stop("No matching meanExp column found for ", comparison_col)
  }

# ---- Keep only rows with valid numeric data --------------------------------
  df <- df %>%
    filter(!is.na(.data[[ratio_col]]), !is.na(.data[[meanexp_col]]))

# ---- Count number of actually plotted proteins per organism ----------------
  counts <- df %>%
    group_by(Organism) %>%
    summarise(n = n(), .groups = "drop")

  n_hs <- counts$n[counts$Organism == "Homo sapiens"]
  if (length(n_hs) == 0) n_hs <- 0
  n_sc <- counts$n[counts$Organism == "Saccharomyces"]
  if (length(n_sc) == 0) n_sc <- 0

# ---- Define color scale ----------------------------------------------------
  color_values <- c("Homo sapiens" = "#56B4E9")
  if (show_yeast) color_values["Saccharomyces"] <- "#E69F00"

# ---- Define dynamic legend labels ------------------------------------------
  labels_values <- c(
    "Homo sapiens" = paste0("Human: ", n_hs,
                            if (!is.null(ratio_hs)) paste0(" (", ratio_hs, ")") else "")
  )
  if (show_yeast) {
    labels_values["Saccharomyces"] <- paste0("Yeast: ", n_sc,
                                             if (!is.null(ratio_sc)) paste0(" (", ratio_sc, ")") else "")
  }

# ---- Base scatter plot -----------------------------------------------------
  p_scatter <- ggplot(df, aes_string(x = meanexp_col, y = ratio_col, color = "Organism")) +
    geom_point(alpha = 0.4, size = 0.05) +
    scale_color_manual(values = color_values, labels = labels_values, name = NULL) +
    labs(
      title = title_suffix,
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

# ---- Add horizontal lines for Homo sapiens if provided ---------------------
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

# ---- Add horizontal lines for Saccharomyces if provided and shown ----------
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

# ---- Calculate mean & SD per organism --------------------------------------
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

# ---- Density panel (right) -------------------------------------------------
  p_density <- ggplot(df, aes_string(x = ratio_col, fill = "Organism", color = "Organism")) +  
    geom_density(alpha = 0.3, fill = NA, linewidth = 0.25) +  
    geom_vline(data = sd_summary, aes(xintercept = mean, color = Organism),  
               linetype = "dotted", linewidth = 0.3) +  
    geom_text(  
      data = sd_summary,  
      aes(  
        x = mean + 0.3,  
        y = 1.25,                      # place near baseline of density  
        label = paste0("FC = ", round(FC_linear, 2), "\nSD = ", round(sd, 2))  
      ),  
      colour = "black",  
      size = 1.3,  
      hjust = 0.5,  
      vjust = 0,  
      angle = 0,  
      inherit.aes = F  
    ) +
    labs(
      y = "Density") +
    scale_fill_manual(values = color_values) +
    scale_color_manual(values = color_values) +
    coord_flip() +
    scale_x_continuous(
      limits = y_limits,
      breaks = seq(y_limits[1], y_limits[2], by = 1)  # same breaks as BA y-axis
    ) +
    scale_y_continuous(
      expand = c(0, 0.02),
      limits = density_limits,
      breaks = seq(density_limits[1], density_limits[2], by = density_breaks)
    ) +
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

# ---- Combine scatter + density ---------------------------------------------
  p_final <- p_scatter + p_density + plot_layout(widths = c(4, 1))

# ---- Return both plot and SD table -----------------------------------------
  return(list(
    plot = p_final,
    sd_summary = sd_summary
  ))
}


p_prot_ratio_DDA_IT8_rep1_H100Y000vsH100Y010_wo_yeast <- plot_blant_altmann_prot(prot_results$long_prot_baseline_log2, 
                                                                                 comparison_col = "unnorm_H100Y000vsH100Y010", 
                                                                                 dataset_filter = "20250920_Neo6A1_TH073_TT_THIDmulti005_DDA_mz0c5_ij8_53c75ng_1_Proteins_purity_tmt_quant",
                                                                                 ratio_hs = "100/100",
                                                                                 ratio_sc = "0/10",
                                                                                 show_yeast = F,
                                                                                 x_limits = c(4,20),
                                                                                 y_limits = c(-3,3),
                                                                                 density_limits = c(0, 2.1),
                                                                                 hline_hs = 0, 
                                                                                 hline_sc = NULL,
                                                                                 title_suffix = "DDA 0.5 m/z, IT 8 ms, 54 ng")
p_prot_ratio_DDA_IT8_rep1_H100Y000vsH100Y010_wo_yeast

ggsave("plots/p_prot_ratio_DDA_IT8_rep1_H100Y000vsH100Y010_wo_yeast.svg", p_prot_ratio_DDA_IT8_rep1_H100Y000vsH100Y010_wo_yeast$plot, 
       width = 12.7, height = 8.9, units = "cm", dpi = 300#, scale = 1 
)

# ---- Plot all possible combinations ----------------------------------------

# ---- Parameters ------------------------------------------------------------

comparisons <- c(
  "unnorm_H100Y000vsH100Y010",
  "unnorm_H100Y000vsH100Y020",
  "unnorm_H100Y000vsH100Y100",
  "unnorm_H100Y100vsH100Y010",
  "unnorm_H100Y100vsH100Y020",
  "unnorm_H100Y020vsH100Y010"
)

yeast_hlines <- c(
  "unnorm_H100Y100vsH100Y010" = 3.32,
  "unnorm_H100Y100vsH100Y020" = 2.32,
  "unnorm_H100Y020vsH100Y010" = 1
)

ratio_hs_list <- c(
  "unnorm_H100Y000vsH100Y010" = "100/100",
  "unnorm_H100Y000vsH100Y020" = "100/100",
  "unnorm_H100Y000vsH100Y100" = "100/100",
  "unnorm_H100Y100vsH100Y010" = "100/100",
  "unnorm_H100Y100vsH100Y020" = "100/100",
  "unnorm_H100Y020vsH100Y010" = "100/100"
)

ratio_sc_list <- c(
  "unnorm_H100Y100vsH100Y010" = "100/10",
  "unnorm_H100Y100vsH100Y020" = "100/20",
  "unnorm_H100Y020vsH100Y010" = "20/10"
)

show_yeast_list <- c(
  "unnorm_H100Y000vsH100Y010" = FALSE,
  "unnorm_H100Y000vsH100Y020" = FALSE,
  "unnorm_H100Y000vsH100Y100" = FALSE,
  "unnorm_H100Y100vsH100Y010" = TRUE,
  "unnorm_H100Y100vsH100Y020" = TRUE,
  "unnorm_H100Y020vsH100Y010" = TRUE
)

# ---- Save function (generate + save + discard each plot immediately) -------

save_all_blant_altmann_plots <- function(
    prot_results,
    comparisons,
    yeast_hlines,
    base_dir        = file.path(getwd(), "plots/prot_ratios"),
    ratio_hs_list   = NULL,
    ratio_sc_list   = NULL,
    show_yeast_list = NULL,
    result_filter   = NULL
) {

  if (!dir.exists(base_dir)) dir.create(base_dir, recursive = TRUE)

  if (!is.null(result_filter)) {
    prot_results <- prot_results[intersect(result_filter, names(prot_results))]
    if (length(prot_results) == 0) stop("No matching result names found in prot_results.")
  }

  total   <- length(comparisons) * sum(sapply(prot_results, function(df) length(unique(df$dataset))))
  counter <- 0L

  for (cmp in comparisons) {
    comparison_dir <- file.path(base_dir, cmp)
    if (!dir.exists(comparison_dir)) dir.create(comparison_dir, recursive = TRUE)

    for (result_name in names(prot_results)) {
      df_full      <- prot_results[[result_name]]
      raw_datasets <- unique(df_full$dataset)

      for (raw_ds in raw_datasets) {

        counter <- counter + 1L

        dataset_label <- sub(".*multi005_(.*?)_Proteins.*", "\\1", raw_ds)

        yline      <- if (cmp %in% names(yeast_hlines))    yeast_hlines[[cmp]]    else NA
        ratio_hs   <- if (!is.null(ratio_hs_list)   && cmp %in% names(ratio_hs_list))   ratio_hs_list[[cmp]]   else NULL
        ratio_sc   <- if (!is.null(ratio_sc_list)   && cmp %in% names(ratio_sc_list))   ratio_sc_list[[cmp]]   else NULL
        show_yeast <- if (!is.null(show_yeast_list) && cmp %in% names(show_yeast_list)) show_yeast_list[[cmp]] else TRUE

        # Horizontal lines are passed into the function, not added to the patchwork afterwards
        res <- plot_blant_altmann_prot(
          long_df        = df_full,
          comparison_col = cmp,
          dataset_filter = raw_ds,
          ratio_hs       = ratio_hs,
          ratio_sc       = ratio_sc,
          show_yeast     = show_yeast,
          hline_hs       = 0,
          hline_sc       = if (!is.na(yline)) yline else NULL
        )

        p <- res$plot +
          ggtitle(paste0(cmp, " — ", result_name, " — ", dataset_label)) +
          theme(plot.title = element_text(hjust = 1, size = 6, colour = "black"))

        clean_label <- stringr::str_replace_all(dataset_label, "[^A-Za-z0-9_-]", "_")
        file_path   <- file.path(comparison_dir, paste0(result_name, "_", clean_label, ".png"))

        ggsave(
          filename = file_path,
          plot     = p,
          device   = ragg::agg_png,
          width    = 12.7,
          height   = 8.9,
          units    = "cm",
          dpi      = 300,
          scale    = 1
        )

        message(sprintf("[%d/%d] Saved: %s", counter, total, basename(file_path)))

        rm(p, res)
        if (counter %% 10L == 0L) gc()
      }
    }
  }

  message("Done — ", total, " plots saved to: ", base_dir)
}

# ---- Run -------------------------------------------------------------------

save_all_blant_altmann_plots(
  prot_results    = prot_results,
  comparisons     = comparisons,
  yeast_hlines    = yeast_hlines,
  ratio_hs_list   = ratio_hs_list,
  ratio_sc_list   = ratio_sc_list,
  show_yeast_list = show_yeast_list,
  result_filter   = c("long_prot_baseline_log2", "long_prot_standard_StoN_grouping_20_log2", 
                      "long_prot_StoN_grouping_40_log2", "long_prot_StoN_grouping_80_log2",
        # NOTE: 'long_prot_S_filter_all_2000_log2' is not present in prot_results;
        #       result_filter intersects with existing names, so it is silently ignored.
                      "long_prot_S_filter_all_2000_log2")
)

# ---- Selected conditions for the draft figures -----------------------------

# DDA H100Y100 vs. H100Y020
p_prot_BA_100vs20_DDA_54ng_IT8_rep1_standard <- plot_blant_altmann_prot(prot_results$long_prot_standard_StoN_grouping_20_log2,
                                                               comparison_col = "unnorm_H100Y100vsH100Y020", 
                                                               dataset_filter = "20250920_Neo6A1_TH073_TT_THIDmulti005_DDA_mz0c5_ij8_53c75ng_1_Proteins_purity_tmt_quant",
                                                               ratio_hs = "100/100",
                                                               ratio_sc = "100/20",
                                                               show_yeast = TRUE,
                                                               x_limits = c(4,20),
                                                               y_limits = c(-2,5),
                                                               hline_hs = 0, 
                                                               hline_sc = 2.32)

p_prot_BA_100vs20_DDA_54ng_IT8_rep1_standard

ggsave("plots/p_prot_BA_100vs20_DDA_54ng_IT8_rep1_standard.svg", p_prot_BA_100vs20_DDA_54ng_IT8_rep1_standard$plot, 
       width = 12.7, height = 8.9, units = "cm", dpi = 300#, scale = 1
)


# DDA H100Y100 vs. H100Y010
p_prot_BA_100vs10_DDA_54ng_IT8_rep1_standard <- plot_blant_altmann_prot(prot_results$long_prot_standard_StoN_grouping_20_log2,
                                                                        comparison_col = "unnorm_H100Y100vsH100Y010", 
                                                                        dataset_filter = "20250920_Neo6A1_TH073_TT_THIDmulti005_DDA_mz0c5_ij8_53c75ng_1_Proteins_purity_tmt_quant",
                                                                        ratio_hs = "100/100",
                                                                        ratio_sc = "100/10",
                                                                        show_yeast = TRUE,
                                                                        x_limits = c(4,20),
                                                                        y_limits = c(-2,5),
                                                                        hline_hs = 0, 
                                                                        hline_sc = 3.32)

p_prot_BA_100vs10_DDA_54ng_IT8_rep1_standard

ggsave("plots/p_prot_BA_100vs10_DDA_54ng_IT8_rep1_standard.svg", p_prot_BA_100vs10_DDA_54ng_IT8_rep1_standard$plot, 
       width = 12.7, height = 8.9, units = "cm", dpi = 300#, scale = 1
)

# DDA H100Y020 vs. H100Y010
p_prot_BA_20vs10_DDA_54ng_IT8_rep1_standard <- plot_blant_altmann_prot(prot_results$long_prot_standard_StoN_grouping_20_log2,
                                                                        comparison_col = "unnorm_H100Y020vsH100Y010", 
                                                                        dataset_filter = "20250920_Neo6A1_TH073_TT_THIDmulti005_DDA_mz0c5_ij8_53c75ng_1_Proteins_purity_tmt_quant",
                                                                        ratio_hs = "100/100",
                                                                        ratio_sc = "20/10",
                                                                        show_yeast = TRUE,
                                                                        x_limits = c(8,20),
                                                                        y_limits = c(-2,3),
                                                                       density_limits = c(0, 2.2),
                                                                        hline_hs = 0, 
                                                                        hline_sc = 1)

p_prot_BA_20vs10_DDA_54ng_IT8_rep1_standard

ggsave("plots/p_prot_BA_20vs10_DDA_54ng_IT8_rep1_standard.svg", p_prot_BA_20vs10_DDA_54ng_IT8_rep1_standard$plot, 
       width = 12.7, height = 8.9, units = "cm", dpi = 300#, scale = 1
)


# DIA H100Y100 vs. H100Y020
p_prot_BA_100vs20_DIA_mz1_54ng_IT8_rep1_standard <- plot_blant_altmann_prot(prot_results$long_prot_standard_StoN_grouping_20_log2,
                                                                   comparison_col = "unnorm_H100Y100vsH100Y020", 
                                                                   dataset_filter = "20250920_Neo6A1_TH073_TT_THIDmulti005_DIA_mz1_ij8_53c75ng_1_Proteins_purity_tmt_quant",
                                                                   ratio_hs = "100/100",
                                                                   ratio_sc = "100/20",
                                                                   show_yeast = T,
                                                                   x_limits = c(4,20),
                                                                   y_limits = c(-2,5),
                                                                   hline_hs = 0, 
                                                                   hline_sc = 2.32)


p_prot_BA_100vs20_DIA_mz1_54ng_IT8_rep1_standard

ggsave("plots/p_prot_BA_100vs20_DIA_mz1_54ng_IT8_rep1_standard.svg", p_prot_BA_100vs20_DIA_mz1_54ng_IT8_rep1_standard$plot, 
       width = 12.7, height = 8.9, units = "cm", dpi = 300#, scale = 1
)

# DIA H100Y100 vs. H100Y010
p_prot_BA_100vs10_DIA_mz1_54ng_IT8_rep1_standard <- plot_blant_altmann_prot(prot_results$long_prot_standard_StoN_grouping_20_log2,
                                                                            comparison_col = "unnorm_H100Y100vsH100Y010", 
                                                                            dataset_filter = "20250920_Neo6A1_TH073_TT_THIDmulti005_DIA_mz1_ij8_53c75ng_1_Proteins_purity_tmt_quant",
                                                                            ratio_hs = "100/100",
                                                                            ratio_sc = "100/10",
                                                                            show_yeast = T,
                                                                            x_limits = c(4,20),
                                                                            y_limits = c(-2,5),
                                                                            hline_hs = 0, 
                                                                            hline_sc = 3.32)


p_prot_BA_100vs10_DIA_mz1_54ng_IT8_rep1_standard

ggsave("plots/p_prot_BA_100vs10_DIA_mz1_54ng_IT8_rep1_standard.svg", p_prot_BA_100vs10_DIA_mz1_54ng_IT8_rep1_standard$plot, 
       width = 12.7, height = 8.9, units = "cm", dpi = 300#, scale = 1
)

# DIA H100Y020 vs. H100Y010
p_prot_BA_20vs10_DIA_mz1_54ng_IT8_rep1_standard <- plot_blant_altmann_prot(prot_results$long_prot_standard_StoN_grouping_20_log2,
                                                                            comparison_col = "unnorm_H100Y020vsH100Y010", 
                                                                            dataset_filter = "20250920_Neo6A1_TH073_TT_THIDmulti005_DIA_mz1_ij8_53c75ng_1_Proteins_purity_tmt_quant",
                                                                            ratio_hs = "100/100",
                                                                            ratio_sc = "20/10",
                                                                            show_yeast = T,
                                                                            x_limits = c(8,20),
                                                                            y_limits = c(-2,3),
                                                                           density_limits = c(0, 2.2),
                                                                            hline_hs = 0, 
                                                                            hline_sc = 1)


p_prot_BA_20vs10_DIA_mz1_54ng_IT8_rep1_standard

ggsave("plots/p_prot_BA_20vs10_DIA_mz1_54ng_IT8_rep1_standard.svg", p_prot_BA_20vs10_DIA_mz1_54ng_IT8_rep1_standard$plot, 
       width = 12.7, height = 8.9, units = "cm", dpi = 300#, scale = 1
)


