# ===========================================================================
# Script:  05a_plots_prot_overlayed.R
# Purpose: Overlayed protein-count plots (input titration and isolation-window
#          comparisons).
# Inputs:  long_prot_counts (03d), condition subsets (04).
# Outputs: p_overlay_input, p_overlay_isol_w_54_ng (Fig2).
# Depends: Run after 03d and 04.
# ===========================================================================

color_vals <- c(
  "DIA 1.0 m/z \nIT40, 3ng"      = "#FFF5EB",
  "DIA 1.0 m/z \nIT40, 5ng"      = "#FEE6CE",
  "DIA 1.0 m/z \nIT40, 10ng"     = "#FDD0A2",
  "DIA 1.0 m/z \nIT40, 20ng"     = "#FDAE6B",
  "DIA 1.0 m/z \nIT40, 25ng"     = "#FD8D3C",
  "DIA 1.0 m/z \nIT40, 30ng"     = "#F16913",
  "DIA 1.0 m/z \nIT40, 40ng"     = "#D94801",
  "DIA 1.0 m/z \nIT40, 53.75ng"  = "#A63603",
  "DIA 1.0 m/z \nIT40, 100ng"    = "#7F2704"
)

colors_method_comparison <- c(
  "DDA 0.5 m/z \nIT8, 53.75ng" = "#5B859EFF",
  "DDA 0.5 m/z \nIT20, 53.75ng" = "#1E395FFF",
  "DIA 0.5 m/z \nIT8, 53.75ng" = "#1E5A46FF",
  "DIA 1.0 m/z \nIT8, 53.75ng" = "#D48F90FF",
  "DIA 1.0 m/z \nIT20, 53.75ng" = "#DF8D71FF",
  "DIA 1.0 m/z \nIT40, 53.75ng" = "#AB84A5FF",
  "DIA 2.0 m/z \nIT8, 53.75ng" = "#D8B847FF"
)


plot_quan_protein_overlay <- function(long_prot_noF,
                                      long_prot_F,
                                      condition_proteins = NULL,
                                      replicate_cols = c("Rep1_F", "Rep2_F", "Rep3_F"),
                                      y_limit = 12000,
                                      color_values,
                                      x_labels = NULL,      # NEW: custom x-axis labels
                                      x_title = "Input [ng]", # NEW: custom x-axis title
                                      title_suffix = "") {

# ---- Prepare both datasets -------------------------------------------------
  prep_data <- function(df, source_label) {
    d <- df
    if (!is.null(condition_proteins)) d <- d %>% filter(label %in% condition_proteins)
    d <- d %>%
      filter(Type == "Quantified") %>%
      mutate(
        source = source_label,
        SD_value = replace_na(SD_value, 0),
        Label_Pos = Value + SD_value + 100
      )
    if (!is.null(condition_proteins)) d <- d %>% mutate(label = factor(label, levels = condition_proteins))
    d
  }

  data_noF <- prep_data(long_prot_noF, "noF")
  data_F   <- prep_data(long_prot_F,   "F")

# ---- Combined for arrows + delta -------------------------------------------
  combined <- data_F %>%
    select(label, Value_F = Value, SD_F = SD_value) %>%
    left_join(data_noF %>% select(label, Value_noF = Value, SD_noF = SD_value), by = "label") %>%
    mutate(
      label_num     = as.numeric(label),
      arrow_x       = label_num + 0.3,
      Delta         = round(Value_noF - Value_F, 0),
      Label_Pos_F   = Value_F   + SD_F   + 100,
      Label_Pos_noF = Value_noF + SD_noF + 100,
      Arrow_mid_y   = (Value_F + SD_F + Value_noF - SD_noF) / 2
    )

# ---- Replicate long format -------------------------------------------------
  make_long <- function(df) {
    cols_present <- intersect(replicate_cols, names(df))
    if (length(cols_present) > 0) {
      df %>%
        pivot_longer(cols = all_of(cols_present), names_to = "Replicate", values_to = "Replicate_Value") %>%
        filter(!is.na(Replicate_Value))
    } else {
      df %>% mutate(Replicate_Value = NA)
    }
  }

  long_noF <- make_long(data_noF)
  long_F   <- make_long(data_F)

# ---- X-axis labels: use custom if provided, else use condition_proteins levels ----
  x_scale <- if (!is.null(x_labels)) {
    scale_x_discrete(labels = x_labels)
  } else {
    scale_x_discrete()
  }

  ggplot() +
    # noF bars — faded, drawn first (behind)
    geom_bar(data = data_noF,
             aes(x = label, y = Value, fill = label),
             stat = "identity", color = "black", width = 0.5, linewidth = 0.2, alpha = 0.5) +
    geom_errorbar(data = data_noF,
                  aes(x = label, ymin = pmax(0, Value - SD_value), ymax = Value + SD_value),
                  width = 0.2, color = "black", size = 0.2, alpha = 0.3) +
    # F bars — full color, on top
    geom_bar(data = data_F,
             aes(x = label, y = Value, fill = label),
             stat = "identity", color = "black", width = 0.45, linewidth = 0.2, alpha = 1) +
    geom_errorbar(data = data_F,
                  aes(x = label, ymin = pmax(0, Value - SD_value), ymax = Value + SD_value),
                  width = 0.15, color = "black", size = 0.2) +
    # Replicate points
    geom_point(data = long_noF,
               aes(x = label, y = Replicate_Value, fill = label),
               shape = 21, color = "black", size = 0.3, alpha = 0.3,
               position = position_jitter(width = 0.1)) +
    geom_point(data = long_F,
               aes(x = label, y = Replicate_Value, fill = label),
               shape = 21, color = "black", size = 0.3, alpha = 0.7,
               position = position_jitter(width = 0.1)) +
    # Label on noF bar
    geom_text(data = combined,
              aes(x = as.numeric(label), y = Label_Pos_noF, label = round(Value_noF, 0)),
              size = 1.5, vjust = 0.4, angle = 90, hjust = -0.2, color = "grey50") +
    # Label on F bar
    geom_text(data = combined,
              aes(x = as.numeric(label), y = Value_F - SD_F - 400, label = round(Value_F, 0)),
              size = 1.5, vjust = 0.4, angle = 90, hjust = -0.2, color = "black") +
    scale_fill_manual(values = color_values, guide = "none") +
    x_scale +
    scale_y_continuous(limits = c(0, y_limit), expand = c(0, 0)) +
    labs(y = expression("#quantified proteins"), x = x_title, title = title_suffix) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title   = element_text(lineheight = .8, face = "plain", hjust = 0.5, size = 6),
      axis.text.x  = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 6, colour = "black"),
      axis.text.y  = element_text(size = 6, colour = "black"),
      axis.title   = element_text(size = 6, colour = "black"),
      panel.grid   = element_blank(),
      panel.background = element_blank(),
      axis.line    = element_line(colour = "black"),
      axis.ticks   = element_line(color = "black", size = 0.25)
    )
}


# same plots but with the S/N filter for empty and occupied as two groups

p_overlay_input <- plot_quan_protein_overlay(  
  long_prot_noF      = long_prot_counts$long_prot_baseline,  
  long_prot_F        = long_prot_counts$long_prot_standard_StoN_grouping_20,  
  condition_proteins = input_titration,  
  y_limit            = 12000,  
  color_values       = color_vals,  
  x_labels           = c("3", "5", "10", "20", "25", "30", "40", "54", "100"),  
  x_title            = "Input [ng]",  
  title_suffix       = "Input titration: DIA 1 m/z IT 40 ms"  
)
p_overlay_input

ggsave("plots/p_overlay_input.svg", device = "svg", p_overlay_input, 
       width = 12.7, height = 8.9, units = "cm", dpi = 300, scale = 1)

p_overlay_isol_w_54_ng <- plot_quan_protein_overlay(  
  long_prot_noF      = long_prot_counts$long_prot_baseline,  
  long_prot_F        = long_prot_counts$long_prot_standard_StoN_grouping_20,  
  condition_proteins = isol_w_IT_53ng,  
  y_limit            = 12000,  
  color_values       = colors_method_comparison,  
  x_labels           = c("IT 8ms", "IT 20ms", "IT 8ms", "IT 8ms", "IT 20ms", "IT 40ms", "IT 8ms"),  
  x_title            = "",  
  title_suffix       = "Method comparison at 54 ng input"  
)

p_overlay_isol_w_54_ng

ggsave("plots/p_overlay_isol_w_54_ng.svg", device = "svg", p_overlay_isol_w_54_ng, 
       width = 12.7, height = 8.9, units = "cm", dpi = 300, scale = 1)