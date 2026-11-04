# =============================================================================
# Script 11: Arranged multi-panel figure (final output)
# Description: Assembles the individual protein-level plots into the final
#              six-panel paper figure and saves it to SVG. This is the final
#              "customer" of the pipeline.
# Depends on: 07_plot_protein_overlay.R, 08_plot_bland_altman_density.R,
#             09_plot_ratio_boxplots_and_cv.R, 10_plot_protein_intensities.R
# =============================================================================

# --- Assemble the six-panel figure -------------------------------------------

fig4_multi008A_StoN_grouping <- plot_grid(p_prot_intensities_DDA_IT20_54ng,
                    p_overlay_54ng_StoN_grouping,
                    p_prot_ratio_boxplots_H100H010_54ng_StoN_grouping,
                    p_violin_methods_cv_H100_54ng_StoN_grouping,
                    p_prot_ratio_DDA_IT20_54ng_rep1_H100vsH010_StoN_grouping$plot,
                    p_prot_ratio_DIA_IT8_HCD29_54ng_rep1_H100vsH010_StoN_grouping$plot,
                    nrow = 3,
                    ncol = 2,
                    align = "vh",
                    axis = "lb",
                    labels = "AUTO",
                    label_size = 10,
                    vjust = 1.2)

fig4_multi008A_StoN_grouping

ggsave("plots/fig4_multi008A_StoN_grouping.svg", fig4_multi008A_StoN_grouping,
       width = 18.3, height = 16, units = "cm", dpi = 300, scale = 1
)
