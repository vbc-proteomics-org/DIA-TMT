# =============================================================================
# Script 01: Experimental design and package loading
# Description: Loads all required packages, defines the TMT channel / group
#              layout (channel_names_groups) used throughout the pipeline, and
#              builds the experimental-design overview plot.
# Depends on: none (run first)
# =============================================================================

# --- Packages ----------------------------------------------------------------
library(writexl)
library(cowplot)
library(patchwork)
library(gridExtra)
library(tidyverse)
library(ggplot2)
library(svglite)
library(Hmisc)
library(RColorBrewer)
library(purrr)
library(zoo)
library(gridExtra)
library(reshape2)
library(viridis)
library(ggpattern)
library(colorspace)
library(grid)
library(ggpubr)
library(readr)
library(data.table)
library(ggridges)
library(scales)
library(ggdist)
library(DescTools)
library(matrixStats)
library(ggtext)
library(ggrepel)

# --- Output directory --------------------------------------------------------

# Subdirectory where plots are written.
plots_dir <- "plots"

# --- TMT channels and group assignments --------------------------------------

mix_levels <- c("Label only", "H1", "H2", "H3")

channel_names_groups <- tibble(
  label = c("126", "127N", "127C", "128N", "128C", "129N", 
            "129C", "130N", "130C", "131N", "131C", "132N",
            "132C", "133N", "133C", "134N", "134C", "135N"),
  mix = c("Label only", "H1", "H1", "H1", "H1", "Label only",
          "H2", "H2", "H2", "H2", "Label only", "Label only",
          "H3", "H3", "H3", "H3", "Label only", "Label only"),
  experiment = c("H000", "H010", "H010", "H010", "H010", "H000",
                 "H050", "H050", "H050", "H050", "H000", "H000",
                 "H100", "H100", "H100", "H100", "H000", "H000")
) %>%
  mutate(
    mix = factor(mix, levels = mix_levels),
    total = if_else(mix == "Label only", 0, 1),
    hela_ratio = case_when(
      mix == "H1" ~ 0.1,
      mix == "H2" ~ 0.5,
      mix == "H3" ~ 1.0,
      TRUE ~ 0.0
    ),
    hela = hela_ratio * total,
    label = factor(label, levels = unique(label))  # Preserve channel order for plotting.
  )

# --- Prepare data for plotting -----------------------------------------------

exp_design <- channel_names_groups %>%
  select(label, mix, experiment, hela) %>%
  pivot_longer(
    cols = c(hela),
    names_to = "protein",
    values_to = "abundance"
  ) %>%
  mutate(
    protein = recode(protein, hela = "HeLa"),
    abundance = if_else(as.character(mix) == "Label only", NA_real_, abundance),
    label = factor(label, levels = levels(channel_names_groups$label))
  )

# X marks for label-only channels.
x_marks <- channel_names_groups %>%
  filter(mix == "Label only") %>%
  transmute(
    label = factor(label, levels = levels(channel_names_groups$label)),
    abundance = 0.5
  )

# Group labels drawn above the plot.
group_labels <- channel_names_groups %>%
  mutate(x = seq_along(label)) %>%
  filter(mix != "Label only") %>%
  group_by(mix, experiment) %>%
  summarise(
    x = mean(x),
    group_label = paste0(as.character(first(mix)), "\n", first(experiment)),
    .groups = "drop"
  )

# --- Experimental-design plot ------------------------------------------------

exp_design_plot <- ggplot(exp_design, aes(x = label, y = abundance, fill = protein)) +
  geom_col(
    position = position_dodge(width = 0.3),
    width = 0.3,
    na.rm = TRUE
  ) +
  geom_point(
    data = x_marks,
    aes(x = label, y = abundance),
    inherit.aes = FALSE,
    shape = 4,
    size = 3,
    stroke = 0.8,
    color = "grey50"
  ) +
  geom_text(
    data = group_labels,
    aes(x = x, y = 1.05, label = group_label),
    inherit.aes = FALSE,
    size = 1.3,
    color = "black",
    lineheight = 0.9
  ) +
  scale_x_discrete(drop = FALSE) +
  scale_fill_manual(values = c("HeLa" = "#56B4E9")) +
  scale_y_continuous(
    expand = c(0, 0),
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    oob = scales::oob_keep
  ) +
  labs(x = "TMT Channel", y = "Ratio", fill = "") +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x  = element_text(angle = 90, hjust = 1, vjust = 0.3, size = 6, colour = "black"),
    axis.text.y  = element_text(size = 6, colour = "black"),
    axis.title.y = element_text(size = 6, colour = "black"),
    axis.title.x = element_blank(),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line  = element_line(colour = "black"),
    axis.ticks = element_line(color = "black", linewidth = 0.25),
    plot.margin = margin(25, 10, 10, 10)
  ) +
  geom_hline(yintercept = 0, color = "black") +
  coord_cartesian(clip = "off")

exp_design_plot

ggsave(file.path(plots_dir, "exp_design_multi008A.svg"), 
       plot = exp_design_plot, 
       device = "svg",
       width = 8.9, 
       height = 7,
       unit = "cm")