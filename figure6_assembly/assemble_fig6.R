# ── assemble_fig6.R ────
# Figure assembly for the manuscript (Fig. 6 and Supplementary Fig. 1).
#
# SCOPE
#   This script belongs in a SEPARATE figure-assembly project, NOT in the
#   multi009 analysis project. It combines panels produced by four independent
#   R projects:
#       multi005  (TMT RI, DDA + DIA)
#       multi007  (TMTc)
#       multi009  (DIA LFQ)
#
# PREREQUISITES
#   Each analysis project must be run FIRST. Every project publishes its panel
#   as an .rds file (e.g. multi009 writes p3_multi009.rds via 02_DIA_boxplot_ratios.R).
#   This script only reads those pre-built panels and arranges them.
#
#   Required packages (ggplot2 and cowplot) must be loaded in the session.
#
# SUGGESTED DIRECTORY LAYOUT (relative paths below assume this)
#   figure_assembly/
#     assemble_fig6.R        
#   ../Fig6/p1_multi005_DDA.rds
#   ../Fig6/p2_multi005_DIA.rds
#   ../Fig6/p3_multi009.rds
#   ../Fig6/p4_multi007.rds

# ────

# Selected boxplots: multi009 LFQ, multi005 TMT RI, multi008A TMT RI, multi007 TMTc
p1 <- readRDS("../Fig6/p1_multi005_DDA.rds")
p2 <- readRDS("../Fig6/p2_multi005_DIA.rds")
p3 <- readRDS("../Fig6/p3_multi009.rds")
p4 <- readRDS("../Fig6/p4_multi007.rds")

# Harmonise theme differences at the assembly stage
common_theme <- theme(
  plot.title = element_text(size = 6, face = "bold", hjust = 0.5, colour = "black"),
  axis.title.x = element_blank(),
  axis.line = element_line(color = "black"),
  axis.title.y = element_text(size = 6, face = "bold", colour = "black"),
  axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, size = 6, colour = "black"),
  axis.text.y = element_text(size = 6, colour = "black"),
  legend.position = "inside",
  legend.position.inside = c(0.9, 0.8),
  legend.text = element_text(size = 6),
)

p1 <- p1 + common_theme
p2 <- p2 + common_theme
p3 <- p3 + common_theme
p4 <- p4 + common_theme

fig6 <- plot_grid(
  p1, p2, p3, p4,
  nrow = 2,
  ncol = 2,
  align = "vh",
  labels = "AUTO",
  label_size = 10,
  vjust = 1,
  axis = "lb"
)

fig6

ggsave("plots/Fig6.svg", fig6,
       width = 18.3, height = 14, units = "cm", dpi = 300, scale = 1)



# ── Supplementary Fig. 1: experimental design for multi008A & multi007 ────

p5 <- readRDS("../multi008A/exp_design_multi008A.rds")
p6 <- readRDS("../multi007/exp_design_multi007.rds")

common_theme2 <- theme(
  plot.title = element_text(size = 6, face = "bold", hjust = 0.5, colour = "black"),
  axis.title.x = element_blank(),
  axis.line = element_line(color = "black"),
  axis.title.y = element_text(size = 6, face = "bold", colour = "black"),
  axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6, colour = "black"),
  axis.text.y = element_text(size = 6, colour = "black"),
  legend.position = "inside",
  legend.position.inside = c(0.9, 0.8),
  legend.text = element_text(size = 6),
)


p5 <- p5 + common_theme2
p6 <- p6 + common_theme2

sup_fig1 <- plot_grid(
  p5, p6,
  ncol = 2,
  align = "vh",
  labels = "AUTO",
  label_size = 10,
  vjust = 1,
  axis = "lb"
)

sup_fig1

ggsave("plots/sup_fig1.svg", sup_fig1,
       width = 18.3, height = 8, units = "cm", dpi = 300, scale = 1)
