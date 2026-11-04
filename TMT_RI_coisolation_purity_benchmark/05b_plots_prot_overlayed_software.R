# ===========================================================================
# Script:  05b_plots_prot_overlayed_software.R
# Purpose: Overlayed protein counts comparing search software (PD / DIA-NN /
#          Spectronaut).
# Inputs:  long_prot_counts (03d), condition subsets (04).
# Outputs: p_overlay_software_comp (sup_fig4).
# Depends: Run after 03d and 04.
# ===========================================================================

plot_quan_protein_overlay_multi <- function(software_pairs,
                                            condition_proteins = NULL,
                                            replicate_cols = c("Rep1_F", "Rep2_F", "Rep3_F"),
                                            y_limit = 12000,
                                            color_values,        # named by software name
                                            x_labels = NULL,     # custom labels for software x-axis
                                            x_title = "",
                                            title_suffix = "",
                                            bar_width = 0.25,
                                            dodge_width = 0.8) {

# ---- Prepare one noF + F pair per software ---------------------------------
  prep_data <- function(df, source_label, sw_name) {
    d <- df
    if (!is.null(condition_proteins)) d <- d %>% filter(label %in% condition_proteins)
    d <- d %>%
      filter(Type == "Quantified") %>%
      mutate(
        source   = source_label,
        software = sw_name,
        SD_value = replace_na(SD_value, 0)
      )
    if (!is.null(condition_proteins))
      d <- d %>% mutate(label = factor(label, levels = condition_proteins))
    d
  }

  all_noF <- purrr::imap_dfr(software_pairs, ~prep_data(.x$noF, "noF", .y))
  all_F   <- purrr::imap_dfr(software_pairs, ~prep_data(.x$F,   "F",   .y))

# ---- Software as x-axis ----
  sw_levels <- names(software_pairs)
  all_noF$software <- factor(all_noF$software, levels = sw_levels)
  all_F$software   <- factor(all_F$software,   levels = sw_levels)

# ---- Replicate long format -------------------------------------------------
  make_long <- function(df) {
    cols_present <- intersect(replicate_cols, names(df))
    if (length(cols_present) > 0) {
      df %>%
        pivot_longer(cols = all_of(cols_present),
                     names_to = "Replicate", values_to = "Replicate_Value") %>%
        filter(!is.na(Replicate_Value))
    } else {
      df %>% mutate(Replicate_Value = NA_real_)
    }
  }

  long_noF <- make_long(all_noF)
  long_F   <- make_long(all_F)

# ---- Value labels ----------------------------------------------------------
  combined <- all_F %>%
    select(software, label, Value_F = Value, SD_F = SD_value) %>%
    left_join(
      all_noF %>% select(software, label, Value_noF = Value, SD_noF = SD_value),
      by = c("software", "label")
    ) %>%
    mutate(
      Label_Pos_noF = Value_noF + SD_noF + 100,
      Label_Pos_F   = pmax(0, Value_F - SD_F - 700)
    )

# ---- X-axis labels: custom if provided, else software names ----------------
  x_scale <- if (!is.null(x_labels)) {
    scale_x_discrete(labels = setNames(x_labels, sw_levels))
  } else {
    scale_x_discrete()
  }

  dodge <- position_dodge(width = dodge_width)

  ggplot() +
    # noF bars — wide, faded (behind)
    geom_bar(data = all_noF,
             aes(x = software, y = Value, fill = software, group = label),
             stat = "identity", color = "black",
             width = bar_width + 0.05, linewidth = 0.2, alpha = 0.4,
             position = dodge) +
    geom_errorbar(data = all_noF,
                  aes(x = software, ymin = pmax(0, Value - SD_value),
                      ymax = Value + SD_value, group = label),
                  width = 0.08, color = "black", size = 0.2, alpha = 0.3,
                  position = dodge) +
    # F bars — solid, narrower (on top)
    geom_bar(data = all_F,
             aes(x = software, y = Value, fill = software, group = label),
             stat = "identity", color = "black",
             width = bar_width, linewidth = 0.2, alpha = 1,
             position = dodge) +
    geom_errorbar(data = all_F,
                  aes(x = software, ymin = pmax(0, Value - SD_value),
                      ymax = Value + SD_value, group = label),
                  width = 0.06, color = "black", size = 0.2,
                  position = dodge) +
    # Replicate dots — noF
    geom_point(data = long_noF,
               aes(x = software, y = Replicate_Value, fill = software, group = label),
               shape = 21, color = "black", size = 0.3, alpha = 0.3,
               position = position_jitterdodge(jitter.width = 0.05, dodge.width = dodge_width)) +
    # Replicate dots — F
    geom_point(data = long_F,
               aes(x = software, y = Replicate_Value, fill = software, group = label),
               shape = 21, color = "black", size = 0.3, alpha = 0.7,
               position = position_jitterdodge(jitter.width = 0.05, dodge.width = dodge_width)) +
    # Value labels — noF (grey, above bar)
    geom_text(data = combined,
              aes(x = software, y = Label_Pos_noF, label = round(Value_noF, 0), group = label),
              size = 1.5, angle = 90, vjust = 0.4, hjust = -0.2, color = "grey50",
              position = dodge) +
    # Value labels — F (black, inside bar)
    geom_text(data = combined,
              aes(x = software, y = Label_Pos_F, label = round(Value_F, 0), group = label),
              size = 1.5, angle = 90, vjust = 0.4, hjust = -0.2, color = "black",
              position = dodge) +
    scale_fill_manual(values = color_values, guide = "none") +
    x_scale +
    scale_y_continuous(limits = c(0, y_limit), expand = c(0, 0)) +
    labs(y = expression("#quantified proteins"), x = x_title, title = title_suffix) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title       = element_text(lineheight = .8, face = "plain", hjust = 0.5, size = 6),
      axis.text.x      = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 6, colour = "black"),
      axis.text.y      = element_text(size = 6, colour = "black"),
      axis.title       = element_text(size = 6, colour = "black"),
      panel.grid       = element_blank(),
      panel.background = element_blank(),
      axis.line        = element_line(colour = "black"),
      axis.ticks       = element_line(color = "black", size = 0.25)
    )
}

colors_software_comp <- c(
  "PD"    = "#fee08b",
  "DIANN" = "#e6f598",
  "SN"    = "#99d594"
)

p_overlay_software_comp <- plot_quan_protein_overlay_multi(
  software_pairs = list(
    "PD"    = list(noF = long_prot_counts$long_prot_baseline,      F = long_prot_counts$long_prot_standard_StoN_grouping_20),
    "DIANN" = list(noF = long_prot_counts$long_prot_DIANN_baseline, F = long_prot_counts$long_prot_DIANN_StoN_grouping_20),
    "SN"    = list(noF = long_prot_counts$long_prot_SNaut_baseline, F = long_prot_counts$long_prot_SNaut_StoN_grouping_20)
  ),
  condition_proteins = test_condition_DIA,
  y_limit      = 12000,
  color_values = colors_software_comp,
  x_labels     = c("PD", "DIANN", "SN"),
  x_title      = "",
  title_suffix = "Software comparison: DIA 1.0 m/z IT 8 ms, 53.75 ng"
)

p_overlay_software_comp

ggsave("plots/p_overlay_software_comp_new.svg", device = "svg", p_overlay_software_comp,
       width = 12.7, height = 8.9, units = "cm", dpi = 300, scale = 1)
