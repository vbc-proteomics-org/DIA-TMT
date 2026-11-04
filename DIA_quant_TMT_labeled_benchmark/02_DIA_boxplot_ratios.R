# ── 02_DIA_boxplot_ratios.R ────
# Combined boxplot of 3 LFQ log2 ratios for the selected dataset (by tag),
# in the style of the multi005 plot_boxplot_vertical_prot_by_method() template.
#
# This script is self-contained and can be run independently: it sources
# 00_shared_functions.R and loads/processes the dataset itself. It does NOT
# depend on 01_DIA_analysis.R.
#
# USAGE
#   1. Set `selected_tag` below to the dataset tag to plot.
#   2. source("02_DIA_boxplot_ratios.R")
# ────

source("00_shared_functions.R")

# ════
# USER CONFIGURATION
# ════

selected_tag <- "MS2_Norm"

# Dataset definition for the selected tag.
datasets <- list(
  list(
    path  = "data/20261008_083751_TT_multi009_single_channels_DIA_LFQ_ME_grouping_correct_FASTA_correct_parsing_Report_non_pivot.tsv",
    label = "MS2 + cross-run normalisation",
    tag   = "MS2_Norm"
  )
)

plots_dir <- "plots"
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

# Comparisons to combine into one plot
comparison_cols <- c(
  "H100Y100vsH100Y010",
  "H100Y100vsH100Y020",
  "H100Y020vsH100Y010"
)
comparison_labels <- c(
  "H100Y100/H100Y010",
  "H100Y100/H100Y020",
  "H100Y020/H100Y010"
)

# ════
# LOAD + PROCESS THE SELECTED DATASET
# ════

ds <- datasets[[which(sapply(datasets, `[[`, "tag") == selected_tag)]]
if (length(ds) == 0) stop("selected_tag not found in datasets: ", selected_tag)

cat("Selected dataset -> tag:", ds$tag, "| label:", ds$label, "\n")

res <- process_sn_file(ds$path, ds$label)
df  <- res$df   # log2 scale, Human + Yeast, with log2ratio_* columns

# ════
# PLOT FUNCTION (adapted from multi005 plot_boxplot_vertical_prot_by_method)
# ════

plot_boxplot_ratios_selected <- function(df,
                    comparison_cols,
                    comparison_labels = NULL,
                    title = NULL,
                    y_limits = NULL) {

  ratio_cols <- paste0("log2ratio_", comparison_cols)
  missing_cols <- ratio_cols[!ratio_cols %in% names(df)]
  if (length(missing_cols) > 0)
    stop("Missing columns in df: ", paste(missing_cols, collapse = ", "))

  # Wide -> long over the requested ratio columns
  long_df <- df %>%
    select(Accession, Species, all_of(ratio_cols)) %>%
    pivot_longer(
      cols      = all_of(ratio_cols),
      names_to  = "comparison",
      values_to = "value"
    ) %>%
    mutate(comparison = str_remove(comparison, "^log2ratio_")) %>%
    filter(!is.na(value), Species %in% c("Human", "Yeast")) %>%
    rename(Organism = Species) %>%
    mutate(Organism = recode(Organism, Human = "Human", Yeast = "Yeast"))

  long_df$comparison <- factor(long_df$comparison, levels = comparison_cols)
  x_labels <- if (!is.null(comparison_labels)) comparison_labels else comparison_cols
  long_df$x_num <- as.numeric(long_df$comparison)

  # Expected ratios (Human always 0; Yeast via expected_log2())
  ratio_human <- setNames(rep(0, length(comparison_cols)), comparison_cols)
  ratio_yeast <- setNames(
    sapply(comparison_cols, function(cn) {
      gA <- str_extract(cn, "^.+(?=vs)")
      gB <- str_extract(cn, "(?<=vs).+")
      expected_log2(gA, gB, "Yeast")
    }),
    comparison_cols
  )

  # Medians per comparison & organism
  medians <- long_df %>%
    group_by(comparison, x_num, Organism) %>%
    summarise(median_val = median(value, na.rm = TRUE), .groups = "drop")

  color_values <- c("Human" = "#56B4E9", "Yeast" = "#E69F00")
  dodge_width  <- 0.7

  # Y-axis range
  if (is.null(y_limits)) {
    y_min <- min(long_df$value, unlist(ratio_human), unlist(ratio_yeast), na.rm = TRUE) - 0.3
    y_max <- max(long_df$value, unlist(ratio_human), unlist(ratio_yeast), na.rm = TRUE) + 0.5
  } else {
    y_min <- y_limits[1]
    y_max <- y_limits[2]
  }
  annotation_y <- y_max * 0.95

  # Reference-line segments (one per comparison)
  ref_lines <- tibble(
    x_num       = seq_along(comparison_cols),
    ratio_human = ratio_human[comparison_cols],
    ratio_yeast = ratio_yeast[comparison_cols]
  ) %>%
    mutate(xmin = x_num - 0.45, xmax = x_num + 0.45)

  # Plot
  p <- ggplot(long_df, aes(x = x_num, y = value, fill = Organism,
                    group = interaction(x_num, Organism))) +
    geom_segment(data = ref_lines,
                 aes(x = xmin, xend = xmax, y = ratio_human, yend = ratio_human),
                 inherit.aes = FALSE, linetype = "dashed", color = "#2E6A87", linewidth = 0.5) +
    geom_segment(data = ref_lines %>% filter(!is.na(ratio_yeast)),
                 aes(x = xmin, xend = xmax, y = ratio_yeast, yend = ratio_yeast),
                 inherit.aes = FALSE, linetype = "dashed", color = "#B87F00", linewidth = 0.5) +
    geom_boxplot(width = 0.6,
                 position = position_dodge2(width = dodge_width, padding = 0.1),
                 outlier.shape = NA,
                 box.linewidth = 0.25,
                 whisker.linewidth = 0.25,
                 median.linewidth = 0.25) +
    geom_text(
      data = medians,
      aes(x = x_num, y = annotation_y, label = round(median_val, 2), group = Organism, fill = NULL),
      position = position_dodge(width = dodge_width),
      vjust = 0, size = 1.7, color = "black"
    ) +
    scale_fill_manual(values = color_values) +
    scale_x_continuous(breaks = seq_along(comparison_cols), labels = x_labels,
                    expand = expansion(add = 0.6)) +
    scale_y_continuous(limits = c(y_min, y_max), expand = c(0, 0), breaks = seq(floor(y_min), ceiling(y_max), by = 1)) +
    theme_minimal(base_size = 10) +
    labs(
      y     = expression(Log[2]~"[FC]"),
      fill  = "Organism",
      title = title
    ) +
    theme(
      plot.title              = element_text(size = 6, face = "bold", hjust = 0.5, colour = "black"),
      axis.title.x            = element_blank(),
      axis.line               = element_line(color = "black"),
      axis.title.y            = element_text(size = 6, face = "bold", colour = "black"),
      axis.text.x             = element_text(angle = 0, hjust = 0.5, vjust = 1, size = 6, colour = "black"),
      axis.text.y             = element_text(size = 6, colour = "black"),
      legend.position         = "inside",
      legend.position.inside  = c(0.9, 0.85),
      legend.text             = element_text(size = 5.5),
      legend.key.size         = unit(0.5, "cm"),
      legend.title            = element_blank(),
      panel.grid              = element_blank(),
      axis.ticks              = element_line(color = "black", linewidth = 0.25)
    )

  return(p)
}

# ════
# BUILD + SAVE
# ════

p_boxplot_selected <- plot_boxplot_ratios_selected(
  df                 = df,
  comparison_cols    = comparison_cols,
  comparison_labels  = comparison_labels,
  title              = ds$label,
  y_limits = c(-2,6)
)

p_boxplot_selected

ggsave("plots/p_boxplot_DIA_LFQ_multi009.svg", p_boxplot_selected,
       width = 12.7, height = 8.9, units = "cm", dpi = 300)

# Publish this panel for the separate figure-assembly project (assemble_fig6.R).
saveRDS(p_boxplot_selected, "p3_multi009.rds")
