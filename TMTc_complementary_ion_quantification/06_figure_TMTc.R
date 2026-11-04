# ===========================================================================
# Script:  06_figure_TMTc.R
# Purpose: Assemble the final TMTc figure from the three panels built earlier:
#          (A) summed protein intensities, (B) three-way protein counts,
#          (C) log2 ratio boxplots.
# Inputs:  plot_list_species_sums_subtr, plots_ratio_by_filter_subtr (script 01),
#          p_three_way_total (script 04).
# Outputs: figure_TMTc and the saved plots/figure_TMTc.svg.
# Depends: 01_TMTc_import_and_noise_subtraction.R,
#          04_protein_counts_TMT_vs_TMTc.R
# ===========================================================================

# --- Panel A: summed protein intensities (DIA120, 50ng) ---
p_int_subtr_F2500_DIA120_50ng <- plot_list_species_sums_subtr$p_int_subtr_F2500_DIA120_50ng +
  scale_y_continuous(limits = c(0,8E6), expand = expansion(mult = c(0, 0))) +
  labs(title = "Summed protein intensities over both replicates: DIA120, 50ng") +
  theme(plot.title = element_text(size = 4),
        axis.text.x      = element_text(size = 4, colour = "black"),
        axis.text.y      = element_text(size = 4, colour = "black")
        )

p_int_subtr_F2500_DIA120_50ng

# --- Panel C: log2 fold-change boxplots across methods (F2500) ---
prots_ratio_F2500 <- plots_ratio_by_filter_subtr$F2500 +
  labs(title = "Log2 FC across methods") +
  theme(plot.title = element_text(size = 4),
        axis.text.x      = element_text(size = 4, colour = "black"),
        axis.text.y      = element_text(size = 4, colour = "black"),
        legend.position = "none")

prots_ratio_F2500

# --- Panel B: three-way protein counts ---
p_three_way_total <- p_three_way_total +
  theme(plot.title = element_text(size = 4),
        axis.text.x      = element_text(size = 4, colour = "black"),
        axis.text.y      = element_text(size = 4, colour = "black"),
        legend.position = "inside")

p_three_way_total

# --- Assemble final figure (patchwork) ---
figure5_TMTc <-
  (p_int_subtr_F2500_DIA120_50ng |
     p_three_way_total) /
  prots_ratio_F2500 +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 10, face = "bold"),
        plot.margin = margin(1, 1, 1, 1, "pt")   # top, right, bottom, left
        )

figure5_TMTc

ggsave("plots/figure5_TMTc.svg", figure5_TMTc,
       width = 17, height = 13, units = "cm", dpi = 300
)
