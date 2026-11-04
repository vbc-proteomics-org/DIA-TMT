# ===========================================================================
# Script:  07_arranged_plots.R
# Purpose: Final figure assembly: arrange the upstream plot objects into Fig1-3
#          and Supplemental Figures 2-4 via cowplot::plot_grid and write the
#          SVGs. Terminal deliverable of the pipeline.
# Inputs:  All p_* plot objects produced by the 05* and 06 scripts.
# Outputs: plots/Fig1.svg ... plots/sup_fig4_correct_FASTA.svg and
#          p1_multi005_DDA.rds.
# Depends: Run last, after every upstream script.
# ===========================================================================

# Figure 1
# plots: A: exp_design_plot, B: method_design, 
# C: p_prot_intensities_DDA_IT8_54ng, D: p_prot_ratio_DDA_IT8_rep1_H100Y000vsH100Y010_wo_yeast$plot 

fig1 <- plot_grid(exp_design_plot,
                  method_design,
                  p_prot_intensities_DDA_IT8_54ng,
                  p_prot_ratio_DDA_IT8_rep1_H100Y000vsH100Y010_wo_yeast$plot,
                  nrow = 2,
                  ncol = 2,
                  rel_heights = c(1,1),
                  rel_widths = c(1,1),
                  align = "vh",
                  axis = "lb",
                  labels = "AUTO",
                  label_size = 8,
                  vjust = 1.2)

fig1

ggsave("plots/Fig1.svg", fig1,
       width = 17, height = 11, units = "cm", dpi = 300, scale = 1
       )


# Figure 2 
# plots: A: p_quan_prot_input_titration, B: p_pur_dist_input_tit, C: p_quan_prot_isol_w_IT_53ng_noF, 
# D: p_pur_dist_isol_w_IT_53ng, E: p_quan_prot_isol_w_IT_53ng_F, F: p_SN_isol_w_IT_53ng


fig2 <- plot_grid(p_overlay_input,
                  p_overlay_isol_w_54_ng,
                  p_pur_dist_input_tit,
                  p_pur_dist_isol_w_IT_53ng,
                  nrow = 2,
                  ncol = 2,
                  rel_heights = c(1,1),
                  rel_widths = c(1,1),
                  align = "vh",
                  labels = "AUTO",
                  label_size = 8,
                  axis = "lb",
                  vjust = 1.2)

fig2

ggsave("plots/Fig2.svg", fig2,
       width = 18.3, height = 14, units = "cm", dpi = 300, scale = 1)

# Figure 3
# plots: A: , B: p_prot_BA_100vs20_DIA_mz1_54ng_IT8_rep1, 
# C: p_boxplot_prot_isol_w_IT_53ng_H100Y020vsH100Y010 , D: p_boxplot_prot_isol_w_IT_53ng_H100Y020vsH100Y010_top20
# E: p_violin_prot_isol_w_IT_53ng_H100Y020, F: p_violin_prot_isol_w_IT_53ng_H100Y010 
fig3 <- plot_grid(p_prot_BA_20vs10_DDA_54ng_IT8_rep1_standard$plot,
                  p_prot_BA_20vs10_DIA_mz1_54ng_IT8_rep1_standard$plot,
                  p_boxplot_prot_isol_w_IT_53ng_H100Y020vsH100Y010,
                  p_boxplot_prot_isol_w_IT_53ng_H100Y020vsH100Y010_top20,
                  p_violin_prot_isol_w_IT_53ng_H100Y020,
                  p_violin_prot_isol_w_IT_53ng_H100Y010,
                  nrow = 3,
                  ncol = 2,
                  rel_heights = c(1,1,1),
                  rel_widths = c(1,1,1),
                  align = "vh",
                  labels = "AUTO",
                  label_size = 10,
                  vjust = 1,
                  axis = "lb")

fig3
ggsave("plots/Fig3.svg", fig3,
       width = 18.3, height = 17, units = "cm", dpi = 300, scale = 1)

# Supplemental Figure 2
# plots: A: p_prot_StoN_titration, B: p_StoN_titration_cv_H100Y020, 
# C: p_StoN_titration_boxplot_H100Y100vsH100Y020

sup_fig2 <- plot_grid(p_prot_StoN_titration,
                  p_StoN_titration_cv_H100Y020,
                  p_StoN_titration_boxplot_H100Y020vsH100Y010,
                  nrow = 3,
                  ncol = 1,
                  rel_heights = c(1,1,1),
                  rel_widths = c(1,1,1),
                  align = "vh",
                  labels = "AUTO",
                  label_size = 10,
                  vjust = 1,
                  axis = "lb")

sup_fig2
ggsave("plots/sup_fig2.svg", sup_fig2,
       width = 18.3, height = 17, units = "cm", dpi = 300, scale = 1)


# Supplemental Figure 3
# plots: A: p_prot_StoN_titration, B: p_StoN_titration_cv_H100Y020, 
# C: p_StoN_titration_boxplot_H100Y100vsH100Y020

sup_fig3 <- plot_grid(p_prot_BA_100vs10_DDA_54ng_IT8_rep1_standard$plot,
                      p_prot_BA_100vs10_DIA_mz1_54ng_IT8_rep1_standard$plot,
                      p_prot_BA_100vs20_DDA_54ng_IT8_rep1_standard$plot,
                      p_prot_BA_100vs20_DIA_mz1_54ng_IT8_rep1_standard$plot,
                      nrow = 2,
                      ncol = 2,
                      rel_heights = c(1,1,1),
                      rel_widths = c(1,1,1),
                      align = "vh",
                      labels = "AUTO",
                      label_size = 10,
                      vjust = 1,
                      axis = "lb")

sup_fig3
ggsave("plots/sup_fig3.svg", sup_fig3,
       width = 18.3, height = 14, units = "cm", dpi = 300, scale = 1)

# Supplemental Figure 4
# plots: A: p_prot_StoN_titration, B: p_StoN_titration_cv_H100Y020, 
# C: p_StoN_titration_boxplot_H100Y100vsH100Y020

sup_fig4 <- plot_grid(p_overlay_software_comp,
                      p_software_comparison_boxplot_H100Y020vsH100Y010,
                      p_software_comparison_cv_H100Y020,
                      p_software_comparison_cv_H100Y010,
                      nrow = 2,
                      ncol = 2,
                      rel_heights = c(1,1,1),
                      rel_widths = c(1,1,1),
                      align = "vh",
                      labels = "AUTO",
                      label_size = 10,
                      vjust = 1,
                      axis = "lb")

sup_fig4
ggsave("plots/sup_fig4.svg", sup_fig4,
       width = 16, height = 12, units = "cm", dpi = 300, scale = 1)