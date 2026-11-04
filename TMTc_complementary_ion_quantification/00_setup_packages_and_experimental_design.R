# ===========================================================================
# Script:  00_setup_packages_and_experimental_design.R
# Purpose: Load all required R packages, define the plot output directory and
#          the TMT channel metadata, and build the experimental-design
#          schematic plot used in the final figure.
# Inputs:  None (packages loaded via library(); design defined inline).
# Outputs: plots_dir, channel_names_groups, exp_design_plot and the saved
#          files exp_design_multi007.rds, plots/exp_design_TMTc.svg/.png.
# Depends: First script in the pipeline; run before all other scripts.
# ===========================================================================

library(writexl)
library(cowplot)
library(patchwork)
library(gridExtra)
library(tidyverse)
library(svglite)
library(Hmisc)
library(RColorBrewer)
library(zoo)
library(reshape2)
library(viridis)
library(ggpattern)
library(colorspace)
library(grid)
library(ggpubr)
library(data.table)
library(ggridges)
library(ggdist)
library(DescTools)
library(matrixStats)
library(ggtext)
library(paletteer)
library(readxl)

# --- Plot output directory ---
plots_dir <- "plots"
dir.create(plots_dir, showWarnings = FALSE)

# --- Define TMT channels and group assignments ---
mix_levels <- c("Label only", "Mix1", "Mix2")

channel_names_groups <- tibble(
  label = c("126", "131N", "135N"),
  mix = c("Label only", "Mix1", "Mix2"),
  experiment = c("H000Y000", "H100Y005", "H100Y025")
) %>%
  mutate(
    mix = factor(mix, levels = mix_levels),
    total = if_else(mix == "Label only", 0, 1),
    yeast_ratio = case_when(
      mix == "Mix1" ~ 0.05,
      mix == "Mix2" ~ 0.25,
      TRUE ~ 0.0
    ),
    hela_ratio = case_when(
      mix %in% c("Mix1", "Mix2") ~ 1.0,
      TRUE ~ 0.0
    ),
    yeast = yeast_ratio * total,
    hela  = hela_ratio  * total,
    label = factor(label, levels = unique(label))
  )

# --- Prepare long data; suppress Label only bars ---
exp_design <- channel_names_groups %>%
  select(label, mix, hela, yeast) %>%
  pivot_longer(cols = c(hela, yeast),
               names_to = "protein", values_to = "abundance") %>%
  mutate(
    protein   = recode(protein, hela = "HeLa", yeast = "Yeast"),
    abundance = if_else(as.character(mix) == "Label only", NA_real_, abundance)
  )

# --- X marks for Label only channels ---
x_marks <- channel_names_groups %>%
  filter(mix == "Label only") %>%
  transmute(
    label     = factor(label, levels = levels(channel_names_groups$label)),
    abundance = 0.5
  )

# --- Group labels above the plot (non-Label-only channels only) ---
group_labels <- channel_names_groups %>%
  mutate(x = seq_along(label)) %>%
  filter(mix != "Label only") %>%
  group_by(mix, experiment) %>%
  summarise(
    x           = mean(x),
    group_label = paste0(as.character(first(mix)), "\n", first(experiment)),
    .groups     = "drop"
  )

# --- Final plot ---
exp_design_plot <- ggplot(exp_design, aes(x = label, y = abundance, fill = protein)) +
  geom_col(
    position = position_dodge(width = 0.3),
    width    = 0.3,
    na.rm    = TRUE
  ) +
  geom_point(
    data        = x_marks,
    aes(x = label, y = abundance),
    inherit.aes = FALSE,
    shape       = 4,
    size        = 3,
    stroke      = 0.8,
    color       = "grey50"
  ) +
  geom_text(
    data        = group_labels,
    aes(x = x, y = 1.05, label = group_label),
    inherit.aes = FALSE,
    size        = 1.8,
    color       = "black",
    lineheight  = 0.9
  ) +
  scale_x_discrete(drop = FALSE) +
  scale_fill_manual(values = c("HeLa" = "#56B4E9", "Yeast" = "#E69F00")) +
  scale_y_continuous(
    expand = c(0, 0),
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    oob    = scales::oob_keep
  ) +
  labs(y = "Ratio", fill = "") +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x      = element_text(angle = 90, hjust = 1, vjust = 0.3, size = 8, colour = "black"),
    axis.text.y      = element_text(size = 8, colour = "black"),
    axis.title.y     = element_text(size = 8, colour = "black"),
    axis.title.x     = element_blank(),
    legend.position  = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line        = element_line(colour = "black"),
    axis.ticks       = element_line(color = "black", linewidth = 0.25),
    plot.margin      = margin(25, 10, 10, 10)   # extra top margin for group labels
  ) +
  geom_hline(yintercept = 0, color = "black") +
  coord_cartesian(clip = "off")

exp_design_plot

# --- Save design plot and RDS ---
ggsave(file.path(plots_dir, "exp_design_TMTc.svg"),
       plot = exp_design_plot,
       device = "svg",
       width = 8.9,
       height = 7,
       unit = "cm")

saveRDS(exp_design_plot, "exp_design_multi007.rds")

ggsave(file.path(plots_dir, "exp_design_TMTc.png"),
       plot = exp_design_plot,
       device = "png",
       width = 8.9,
       height = 7,
       dpi = 600,
       unit = "cm")