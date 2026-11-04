# ===========================================================================
# Script:  01_experimental_design_and_packages.R
# Purpose: Load all required R packages and build the experimental- and
#          method-design schematic plots plus the TMT channel metadata used
#          throughout the pipeline.
# Inputs:  None (packages loaded via library(); design defined inline).
# Outputs: exp_design_plot, method_design, channel_names_groups and related
#          design objects.
# Depends: First script in the pipeline; run before all other scripts.
# ===========================================================================

# Load packages

library(tidyverse)
library(patchwork)
library(svglite)
library(RColorBrewer)
library(viridis)
library(grid)
library(data.table)
library(scales)
library(matrixStats)
library(cowplot)
library(ragg)
library(arrow)


# Define the subdirectory where plots should be saved "plots"
plots_dir <- "plots"

# ---- Define TMT channels and group assignments -----------------------------
mix_levels <- c("Label only", "HeLa control", "Mix1", "Mix2", "Mix3")

channel_names_groups <- tibble(
  label = c("126", "127N", "127C", "128N", "128C", "129N", 
            "129C", "130N", "130C", "131N", "131C", "132N",
            "132C", "133N", "133C", "134N", "134C", "135N"),
  mix = c("Label only", "Label only", "Label only", "HeLa control", "HeLa control", "HeLa control",
          "Label only", "Mix1", "Mix1", "Mix1", "Label only", "Mix2",
          "Mix2", "Mix2", "Label only", "Mix3", "Mix3", "Mix3"),
  experiment = c("H000Y000", "H000Y000", "H000Y000", "H100Y000", "H100Y000", "H100Y000",
                 "H000Y000", "H100Y010", "H100Y010", "H100Y010", "H000Y000", "H100Y020",
                 "H100Y020", "H100Y020", "H000Y000", "H100Y100", "H100Y100", "H100Y100")
) %>%
  mutate(
    mix = factor(mix, levels = mix_levels),
    total = if_else(mix == "Label only", 0, 1),
    yeast_ratio = case_when(
      mix == "Mix1" ~ 0.1,
      mix == "Mix2" ~ 0.2,
      mix == "Mix3" ~ 1.0,
      TRUE ~ 0.0
    ),
    hela_ratio = case_when(
      mix %in% c("Mix1", "Mix2", "Mix3") ~ 1.0,
      mix == "HeLa control" ~ 1.0,
      TRUE ~ 0.0
    ),
    yeast = yeast_ratio * total,
    hela = hela_ratio * total,
    label = factor(label, levels = unique(label))  # Maintain order for plotting
  )

# ---- Add labels for facet strip: "Type\nRatio: x:y" ------------------------
mix_labels <- c(
  "Label only"     = "Label only\nH000Y000",
  "HeLa control"  = "HeLa control\nH100Y000",
  "Mix1"          = "Mix1\nH100Y010",
  "Mix2"          = "Mix2\nH100Y020",
  "Mix3"          = "Mix3\nH100Y100"
)

# ---- Prepare data for plotting ---------------------------------------------
exp_design <- channel_names_groups %>%  
  select(label, mix, experiment, hela, yeast) %>%  
  pivot_longer(  
    cols = c(hela, yeast),  
    names_to = "protein",  
    values_to = "abundance"  
  ) %>%  
  mutate(  
    protein = recode(protein, hela = "HeLa", yeast = "Yeast"),  
    abundance = if_else(as.character(mix) == "Label only", NA_real_, abundance),  
    label = factor(label, levels = levels(channel_names_groups$label))  
  )  

# ---- X marks for label-only channels ---------------------------------------
x_marks <- channel_names_groups %>%  
  filter(mix == "Label only") %>%  
  transmute(  
    label = factor(label, levels = levels(channel_names_groups$label)),  
    abundance = 0.5  
  )  

# ---- Group labels above the plot (excluding Label only) --------------------
group_labels <- channel_names_groups %>%  
  mutate(x = seq_along(label)) %>%  
  filter(mix != "Label only") %>%  
  group_by(mix, experiment) %>%  
  summarise(  
    x = mean(x),  
    group_label = paste0(as.character(first(mix)), "\n", first(experiment)),  
    .groups = "drop"  
  )  

# ---- Final plot ------------------------------------------------------------
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
  scale_fill_manual(values = c("HeLa" = "#56B4E9", "Yeast" = "#E69F00")) +  
  scale_y_continuous(  
    expand = c(0, 0),  
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    oob = scales::oob_keep        # keeps out-of-bounds values instead of dropping them
  ) +  
  labs(x = "TMT Channel", y = "Ratio", fill = ""#, title ="Experimental design"
       ) +  
  theme_minimal(base_size = 10) +  
  theme(
    plot.title   = element_text(lineheight = .8, face = "plain", hjust = 0.5, size = 6),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3, size = 6, colour = "black"),  
    axis.text.y = element_text(size = 6, colour = "black"),  
    axis.title.y = element_text(size = 6, colour = "black"),  
    axis.title.x = element_blank(),  
    legend.position = "none",  
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    panel.background = element_blank(),  
    axis.line = element_line(colour = "black"),  
    axis.ticks = element_line(color = "black", linewidth = 0.25),  
    plot.margin = margin(25, 10, 10, 10)  # increased top margin for labels
  ) +  
  geom_hline(yintercept = 0, color = "black") +  
  coord_cartesian(clip = "off")  

exp_design_plot

ggsave(file.path(plots_dir, "exp_design.svg"), 
       plot = exp_design_plot, 
       device = "svg",
       width = 8.9, 
       height = 7,
       unit = "cm")


# Define all possible inputs
inputs <- c(3, 5, 10, 20, 25, 30, 40, 54, 100)

# Define IT per method
IT <- tibble(
  Method = c("DDA 0.5 m/z IT8ms", 
             "DDA 0.5 m/z IT20ms",
             "DIA 0.5 m/z IT8ms", 
             "DIA 1.0 m/z IT8ms", 
             "DIA 1.0 m/z IT20ms", 
             "DIA 1.0 m/z IT40ms",  
             "DIA 2.0 m/z IT8ms"),
  IT = c(8,20,8,8,20,40,8)
)

# Define tested inputs for each method
tested_inputs <- tribble(
  ~Method, ~Input,
  "DDA 0.5 m/z IT8ms", 3,
  "DDA 0.5 m/z IT8ms", 54,
  "DDA 0.5 m/z IT8ms", 100,
  "DDA 0.5 m/z IT20ms", 3,
  "DDA 0.5 m/z IT20ms", 54,
  "DDA 0.5 m/z IT20ms", 100,
  "DIA 0.5 m/z IT8ms", 3,
  "DIA 0.5 m/z IT8ms", 54,
  "DIA 0.5 m/z IT8ms", 100,
  "DIA 1.0 m/z IT8ms", 3,
  "DIA 1.0 m/z IT8ms", 54,
  "DIA 1.0 m/z IT8ms", 100,
  "DIA 1.0 m/z IT20ms", 3,
  "DIA 1.0 m/z IT20ms", 54,
  "DIA 1.0 m/z IT20ms", 100,
  "DIA 1.0 m/z IT40ms", 3,
  "DIA 1.0 m/z IT40ms", 5,
  "DIA 1.0 m/z IT40ms", 10,
  "DIA 1.0 m/z IT40ms", 20,
  "DIA 1.0 m/z IT40ms", 25,
  "DIA 1.0 m/z IT40ms", 30,
  "DIA 1.0 m/z IT40ms", 40,
  "DIA 1.0 m/z IT40ms", 54,
  "DIA 1.0 m/z IT40ms", 100,
  "DIA 2.0 m/z IT8ms", 3,
  "DIA 2.0 m/z IT8ms", 54,
  "DIA 2.0 m/z IT8ms", 100
)

full <- expand_grid(Method = IT$Method, Input = inputs) %>%
  # Join and mark if tested
  left_join(tested_inputs %>% mutate(Tested = TRUE), by = c("Method", "Input")) %>%
  mutate(Tested = if_else(is.na(Tested), FALSE, TRUE)) %>%
  left_join(IT, by = "Method") 

# Assign IT groups for shapes
full <- full %>%
  mutate(IT_group = case_when(
    IT <= 8 ~ "8",
    IT <= 20 ~ "20",
    TRUE ~ "40"
  ))

full <- full %>%
  mutate(Tested_label = if_else(Tested, "Tested", "Not Tested"))


full$Method <- factor(full$Method, levels = c("DDA 0.5 m/z IT8ms", 
                                              "DDA 0.5 m/z IT20ms",
                                              "DIA 0.5 m/z IT8ms", 
                                              "DIA 1.0 m/z IT8ms", 
                                              "DIA 1.0 m/z IT20ms", 
                                              "DIA 1.0 m/z IT40ms",  
                                              "DIA 2.0 m/z IT8ms"))


# method_design plot
method_design <- ggplot(full, aes(x = factor(Input), y = Method, 
                 shape = IT_group, color = Tested_label)) +
  geom_point(size = 1.7) +
  scale_color_manual(values = c("Not Tested" = "white", "Tested" = "#0073C2FF")) +
  scale_shape_manual(values = c("8" = 15, "20" = 16, "40" = 17)) +  # triangle, circle, square
  theme_minimal(base_size = 10) +
  labs(
    title = "Method design",
    x = "Input amount [ng]",
    y = "Isolation window [m/z]",
    color = "Tested",
    shape = "IT [ms]"
  ) +
  theme(
    plot.title   = element_text(lineheight = .8, face = "plain", hjust = 0.5, size = 6),
    axis.text.x = element_text(size = 6, colour = "black"),
    axis.text.y = element_text(size = 6, colour = "black"),
    axis.title.y = element_text(size = 6, colour = "black"),
    axis.title.x = element_text(size = 6), colour = "black",
    legend.title = element_text(size = 5, face = "bold", colour = "black"),
    legend.text = element_text(size = 5, colour = "black"),
    legend.key.size = unit(0.3, "cm"),    
    legend.position = "inside",
    panel.grid = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.ticks.y = element_line(color = "black", size = 0.25),
    axis.ticks.x = element_line(color = "black", size = 0.25)
  ) +
  scale_y_discrete(labels = c("0.5", "0.5","0.5", "1.0", "1.0", "1.0", "2.0"))


method_design


ggsave(filename = file.path(plots_dir, "method_design.svg"), 
       plot = method_design,
       device = "svg",
       width = 8.9, 
       height = 7,
       units = "cm",
       )
