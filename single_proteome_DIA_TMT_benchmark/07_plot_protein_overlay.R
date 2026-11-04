# =============================================================================
# Script 07: Overlay plot of quantified proteins
# Description: Defines plot_quan_protein_overlay(), which overlays identified
#              ("noF") and quantified ("F") protein counts per method with error
#              bars, replicate points and delta labels, and builds the 54 ng
#              method-comparison overlay figure.
# Depends on: 05_concatenate_protein_counts.R, 06_condition_subsets_and_colors.R
# =============================================================================

# --- Colour values per method ------------------------------------------------

color_vals <- c(
  "DDA 0.5 m/z \nIT20, 53.75ng"     = "#E78AC3",
  "DIA 1.0 m/z \nIT8, 53.75ng"      = "#66C2A5",
  "DIA 1.0 m/z \nHCD29, 53.75ng"      = "#FC8D62",
  "DIA 1.0 m/z \nALL, 53.75ng"    = "#8DA0CB"
  )

# --- Overlay plot function ---------------------------------------------------

plot_quan_protein_overlay <- function(long_prot_noF,
                    long_prot_F,
                    condition_proteins = NULL,
                    replicate_cols = c("Rep1_F", "Rep2_F", "Rep3_F"),
                    y_limit = 12000,
                    color_values,
                    x_labels = NULL,
                    x_title = "Input [ng]",
                    title_suffix = "") {

  # --- Prepare both datasets ---
  prep_data <- function(df, source_label) {
    d <- df
    if (!is.null(condition_proteins)) d <- d %>% filter(label %in% condition_proteins)
    d <- d %>%
      filter(Type == "Quantified") %>%
      mutate(
        source = source_label,
        SD_value = replace_na(SD_value, 0),
        Label_Pos = Value + SD_value + 100,
        label = factor(label, levels = names(color_values))  # Force order from color_values.
      )
    if (!is.null(condition_proteins)) d <- d %>% mutate(label = factor(label, levels = condition_proteins))
    d
  }

  data_noF <- prep_data(long_prot_noF, "noF")
  data_F   <- prep_data(long_prot_F,   "F")

  # --- Combined for arrows + delta ---
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

  # --- Replicate long format ---
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

  # --- x-axis labels: use custom if provided, else default discrete scale ---
  x_scale <- if (!is.null(x_labels)) {
    scale_x_discrete(labels = x_labels)
  } else {
    scale_x_discrete()
  }

  ggplot() +
    # noF bars - faded, drawn first (behind).
    geom_bar(data = data_noF,
             aes(x = label, y = Value, fill = label),
             stat = "identity", color = "black", width = 0.5, linewidth = 0.2, alpha = 0.5) +
    geom_errorbar(data = data_noF,
                  aes(x = label, ymin = pmax(0, Value - SD_value), ymax = Value + SD_value),
                  width = 0.2, color = "black", size = 0.2, alpha = 0.3) +
    # F bars - full colour, on top.
    geom_bar(data = data_F,
             aes(x = label, y = Value, fill = label),
             stat = "identity", color = "black", width = 0.45, linewidth = 0.2, alpha = 1) +
    geom_errorbar(data = data_F,
                  aes(x = label, ymin = pmax(0, Value - SD_value), ymax = Value + SD_value),
                  width = 0.15, color = "black", size = 0.2) +
    # Replicate points.
    geom_point(data = long_noF,
               aes(x = label, y = Replicate_Value, fill = label),
               shape = 21, color = "black", size = 0.3, alpha = 0.3,
               position = position_jitter(width = 0.1)) +
    geom_point(data = long_F,
               aes(x = label, y = Replicate_Value, fill = label),
               shape = 21, color = "black", size = 0.3, alpha = 0.7,
               position = position_jitter(width = 0.1)) +
    # Count label on noF bar.
    geom_text(data = combined,
              aes(x = as.numeric(label), y = Label_Pos_noF, label = round(Value_noF, 0)),
              size = 1.5, vjust = 0.4, angle = 90, hjust = -0.2, color = "grey50") +
    # Count label on F bar.
    geom_text(data = combined,
              aes(x = as.numeric(label), y = Value_F - SD_F - 700, label = round(Value_F, 0)),
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

# --- Build the 54 ng method-comparison overlay -------------------------------

p_overlay_54ng_StoN_grouping <- plot_quan_protein_overlay(
  long_prot_noF      = long_prot_counts$long_prot_baseline,
  long_prot_F        = long_prot_counts$long_prot_standard_StoN_grouping_20,
  condition_proteins = methods_54ng,
  y_limit            = 12000,
  color_values       = color_vals,
  x_labels           = c("DDA", "DIA", "DIA HCD 29%", "DIA no loop control"),
  x_title            = "Input [ng]",
  title_suffix       = "Method comparison at 53 ng input"
)
p_overlay_54ng_StoN_grouping

ggsave("plots/p_overlay_54ng_StoN_grouping.svg", device = "svg", p_overlay_54ng_StoN_grouping,
       width = 12.7, height = 8.9, units = "cm", dpi = 300, scale = 1)
