# ── 00_shared_functions.R ────
# Shared constants, helpers, and plotting functions for the DIA/SN pipeline.
# Source this file once; it has no side-effects (no I/O, no plotting).
#
#   source("00_shared_functions.R")
# ────

suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
  library(cowplot)
  library(patchwork)
})

# ════
# 1.  CONSTANTS
# ════

channel_map <- c(
  "128C" = "H100Y000", "128N" = "H100Y000", "129N" = "H100Y000",
  "130C" = "H100Y010", "130N" = "H100Y010", "131N" = "H100Y010",
  "132C" = "H100Y020", "132N" = "H100Y020", "133N" = "H100Y020",
  "134C" = "H100Y100", "134N" = "H100Y100", "135N" = "H100Y100"
)

channels         <- names(channel_map)
channel_order    <- c("128C","128N","129N","130C","130N","131N",
                    "132C","132N","133N","134C","134N","135N")
groups           <- c("H100Y000","H100Y010","H100Y020","H100Y100")
experiment_order <- groups

group_info <- tibble(channel = names(channel_map), experiment = unname(channel_map))

# Yeast spike-in amounts (relative, arbitrary units)
yeast_amount <- c(H100Y000 = 0, H100Y010 = 1, H100Y020 = 2, H100Y100 = 10)
human_label  <- c(H100Y000 = "1", H100Y010 = "1", H100Y020 = "1", H100Y100 = "1")
yeast_label  <- c(H100Y000 = "0", H100Y010 = "1", H100Y020 = "2", H100Y100 = "10")

group_channels_map <- list(
  H100Y000 = c("128C","128N","129N"),
  H100Y010 = c("130C","130N","131N"),
  H100Y020 = c("132C","132N","133N"),
  H100Y100 = c("134C","134N","135N")
)

comparison_pairs <- list(
  c("H100Y000","H100Y010"), c("H100Y000","H100Y020"), c("H100Y000","H100Y100"),
  c("H100Y100","H100Y010"), c("H100Y100","H100Y020"), c("H100Y020","H100Y010")
)

okabe_ito_colors <- c(
  "H100Y000" = "#56B4E9", "H100Y010" = "#009E73",
  "H100Y020" = "#E69F00", "H100Y100" = "#D55E00"
)
species_colors <- c("Human" = "#56B4E9", "Yeast" = "#E69F00")

# ════
# 2.  EXPECTED LOG2 RATIO
# ════

#' Expected log2 ratio for a given group pair and organism.
#' Human is always 0 (equal HeLa input volume across all channels).
#' Yeast returns NA when either group is H100Y000 (no yeast in ctrl).
expected_log2 <- function(groupA, groupB, organism) {
  if (organism == "Human") return(0)
  yA <- yeast_amount[[groupA]]; yB <- yeast_amount[[groupB]]
  if (yA == 0 || yB == 0) return(NA_real_)
  log2(yA / yB)
}

# ════
# 3.  DATA PROCESSING
# ════

#' Read, pivot, classify, and compute group statistics for one SN non-pivot TSV.
#'
#' @param path   Full path to the TSV.
#' @param label  Display label for console messages.
#' @return List: $df_raw (linear-scale wide, all species),
#'               $df     (log2-scale wide, Human+Yeast, with group stats).
process_sn_file <- function(path, label = basename(path)) {

  cat("\n", strrep("─", 60), "\n", sep = "")
  cat("Processing:", label, "\n")

  df_long <- read_tsv(path, show_col_types = FALSE) %>%
    mutate(channel     = str_extract(R.FileName, "[0-9]+[NC]$"),
           PG.Quantity = as.numeric(PG.Quantity))

  ch_found <- sort(unique(df_long$channel))
  cat("  Rows:", nrow(df_long), "| Channels:", paste(ch_found, collapse = ", "), "\n")

  # ── Pivot ────
  df_wide <- df_long %>%
    select(Accession = PG.ProteinAccessions,
           Description = PG.ProteinDescriptions,
           Organisms = PG.Organisms,
           channel, PG.Quantity) %>%
    group_by(Accession, Description, Organisms, channel) %>%
    summarise(PG.Quantity = mean(PG.Quantity, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = channel, values_from = PG.Quantity)

  # ── Species classification ────
  df_wide <- df_wide %>%
    mutate(Species = case_when(
      str_detect(Organisms, regex("Saccharomyces cerevisiae", ignore_case = TRUE)) &
        !str_detect(Organisms, regex("Homo sapiens",           ignore_case = TRUE)) ~ "Yeast",
      str_detect(Organisms, regex("Homo sapiens",             ignore_case = TRUE)) &
        !str_detect(Organisms, regex("Saccharomyces cerevisiae", ignore_case = TRUE)) ~ "Human",
      TRUE ~ "Other"
    ))

  df_raw <- df_wide                    # linear, all species
  df     <- df_wide %>% filter(Species %in% c("Human","Yeast"))

  cat("  Human:", sum(df$Species == "Human"),
      "| Yeast:", sum(df$Species == "Yeast"),
      "| Other (removed):", sum(df_wide$Species == "Other"), "\n")

  # ── Only work on channels actually present in this file ────
  ch_present <- intersect(channels, names(df))

  # ── Log2 transform ────
  df <- df %>%
    mutate(across(all_of(ch_present), as.numeric)) %>%
    mutate(across(all_of(ch_present), ~ log2(if_else(. > 0, ., NA_real_))))

  cat("  Valid quants per channel:\n  ",
      paste(sprintf("%s:%d", ch_present,
                    sapply(ch_present, \(c) sum(!is.na(df[[c]])))),
            collapse = "  "), "\n")

  # ── Group means, CVs, pairwise log2 ratios ────
  gc <- group_info %>%
    filter(channel %in% ch_present) %>%
    split(.$experiment) %>%
    lapply(\(x) x$channel)

  safe_sub  <- function(a, b) ifelse(!is.na(a) & !is.na(b), a - b,       NA_real_)
  safe_mean <- function(a, b) ifelse(!is.na(a) & !is.na(b), (a + b) / 2, NA_real_)

  for (grp in names(gc)) {
    cols <- gc[[grp]]
    mat  <- as.matrix(df[, cols, drop = FALSE])
    ok   <- !apply(mat, 1, anyNA)
    means <- cvs <- rep(NA_real_, nrow(mat))
    if (any(ok)) {
      m          <- mat[ok, , drop = FALSE]
      means[ok]  <- rowMeans(m)
      if (ncol(m) > 1) cvs[ok] <- 100 * apply(m, 1, sd) / abs(rowMeans(m))
    }
    df[[paste0("mean_", grp)]] <- means
    df[[paste0("cv_",   grp)]] <- cvs
  }

  for (p in comparison_pairs) {
    pn <- paste0(p[1], "vs", p[2])
    ma <- paste0("mean_", p[1]); mb <- paste0("mean_", p[2])
    if (ma %in% names(df) && mb %in% names(df)) {
      df[[paste0("log2ratio_", pn)]] <- safe_sub( df[[ma]], df[[mb]])
      df[[paste0("meanExpr_",  pn)]] <- safe_mean(df[[ma]], df[[mb]])
    }
  }

  list(df_raw = df_raw, df = df, ch_present = ch_present)
}

# ════
# 4.  SHARED THEME
# ════

theme_dia <- function(base_size = 10) {
  theme_minimal(base_size = base_size) %+replace% theme(
    axis.line  = element_line(color = "black"),
    axis.ticks = element_line(color = "black", linewidth = 0.3),
    panel.grid = element_blank(),
    axis.text  = element_text(size = 6, colour = "black"),
    axis.title = element_text(size = 6, colour = "black"),
    plot.title = element_text(face = "plain", hjust = 0.5, size = 7)
  )
}

# ════
# 5.  PLOT FUNCTIONS
# ════

# ── 5a. Summed channel intensities (total) ────
make_channel_intensities <- function(df_raw, group_info, ch_present,
                    channel_order, experiment_order) {
  ch_use <- intersect(ch_present, channel_order)
  df_raw %>%
    filter(Species %in% c("Human","Yeast")) %>%
    mutate(across(all_of(ch_use), as.numeric)) %>%
    summarise(across(all_of(ch_use), ~ sum(., na.rm = TRUE))) %>%
    pivot_longer(everything(), names_to = "channel", values_to = "sum_intensity") %>%
    left_join(group_info, by = "channel") %>%
    mutate(channel    = factor(channel,    levels = channel_order),
           experiment = factor(experiment, levels = experiment_order)) %>%
    arrange(channel)
}

plot_channel_intensities <- function(intensity_df, title = NULL,
                    fill_colors = okabe_ito_colors) {
  # Compute y_max explicitly so geom_text labels never push the axis out of range.
  y_max <- max(intensity_df$sum_intensity, na.rm = TRUE) * 1.22

  ggplot(intensity_df, aes(x = channel, y = sum_intensity, fill = experiment)) +
    geom_col(width = 0.7, color = "black", linewidth = 0.25) +
    geom_text(aes(label = scales::scientific(sum_intensity, digits = 2)),
              vjust = -0.3, size = 1.5) +
    scale_fill_manual(values = fill_colors, name = "Group") +
    scale_y_continuous(limits = c(0, y_max), expand = c(0, 0)) +
    labs(title = title, x = NULL, y = "Summed protein intensity (linear)") +
    theme_dia() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          legend.text  = element_text(size = 6),
          legend.title = element_text(size = 6))
}

# ── 5b. Summed channel intensities by species ────
make_channel_intensities_by_species <- function(df_raw, group_info, ch_present,
                    channel_order, experiment_order) {
  ch_use <- intersect(ch_present, channel_order)
  df_raw %>%
    filter(Species %in% c("Human","Yeast")) %>%
    mutate(across(all_of(ch_use), as.numeric)) %>%
    group_by(Species) %>%
    summarise(across(all_of(ch_use), ~ sum(., na.rm = TRUE)), .groups = "drop") %>%
    pivot_longer(-Species, names_to = "channel", values_to = "sum_intensity") %>%
    left_join(group_info, by = "channel") %>%
    mutate(channel    = factor(channel,    levels = channel_order),
           experiment = factor(experiment, levels = experiment_order),
           Species    = factor(Species,    levels = c("Human","Yeast"))) %>%
    arrange(channel, Species)
}

plot_channel_intensities_by_species <- function(intensity_df, title = NULL,
                    fill_colors = species_colors) {
  # Hard y_max so rotated text labels (hjust = -0.1, angle = 90) never push
  # the axis beyond the data range.
  y_max <- max(intensity_df$sum_intensity, na.rm = TRUE) * 1.30

  ggplot(intensity_df, aes(x = channel, y = sum_intensity, fill = Species)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.65,
             color = "black", linewidth = 0.25) +
    geom_text(aes(label = scales::scientific(sum_intensity, digits = 2)),
              position = position_dodge(width = 0.7),
              vjust = 0.4, hjust = -0.1, angle = 90, size = 1.3) +
    scale_fill_manual(values = fill_colors, name = NULL) +
    scale_y_continuous(limits = c(0, y_max), expand = c(0, 0)) +
    labs(title = title, x = NULL, y = "Summed protein intensity (linear)") +
    theme_dia() +
    theme(axis.text.x          = element_text(angle = 90, hjust = 1, vjust = 0.5),
          legend.text          = element_text(size = 6),
          legend.position      = "inside",
          legend.justification = c(1, 1))
}

# ── 5c. Group-level BA plot (scatter + density) ────
#
# y_limits and x_limits are both NULL by default → computed from this dataset's data.
# Pass explicit values only to force a specific range.
plot_ba_dia <- function(df, comparison_col,
                    ratio_hs = NULL, ratio_sc = NULL,
                    show_yeast     = TRUE,
                    y_limits       = NULL,   # NULL = auto from data
                    x_limits       = NULL,   # NULL = auto from data
                    hline_hs       = 0,
                    hline_sc       = NULL,
                    hline_color_hs = "#2E6A87",
                    hline_color_sc = "#B87F00",
                    title_suffix   = NULL) {

  df <- df %>% mutate(Organism = case_when(
    Species == "Human" ~ "Homo sapiens",
    Species == "Yeast" ~ "Saccharomyces",
    TRUE               ~ "Other"))
  if (!show_yeast) df <- df %>% filter(Organism != "Saccharomyces")

  ratio_col   <- paste0("log2ratio_", comparison_col)
  meanexp_col <- paste0("meanExpr_",  comparison_col)
  if (!ratio_col   %in% names(df)) { message("Skipping (column missing): ", ratio_col);   return(invisible(NULL)) }
  if (!meanexp_col %in% names(df)) { message("Skipping (column missing): ", meanexp_col); return(invisible(NULL)) }

  df <- df %>% filter(!is.na(.data[[ratio_col]]), !is.na(.data[[meanexp_col]]))
  if (nrow(df) == 0) { message("Skipping (no data): ", comparison_col); return(invisible(NULL)) }

  # ── Auto-compute limits from this dataset ────
  if (is.null(x_limits)) {
    x_lo     <- quantile(df[[meanexp_col]], 0.005, na.rm = TRUE) - 0.5
    x_hi     <- quantile(df[[meanexp_col]], 0.995, na.rm = TRUE) + 0.5
    x_limits <- c(x_lo, x_hi)
  }

  if (is.null(y_limits)) {
    # Use all species being plotted (human always, yeast if show_yeast)
    y_vals   <- df[[ratio_col]]
    y_lo     <- floor(  quantile(y_vals, 0.005, na.rm = TRUE) - 0.3)
    y_hi     <- ceiling(quantile(y_vals, 0.995, na.rm = TRUE) + 0.3)
    # Always include 0 and any expected lines so they are never cut off
    ref_vals <- c(0, hline_hs)
    if (!is.null(hline_sc) && !is.na(hline_sc)) ref_vals <- c(ref_vals, hline_sc)
    y_lo     <- min(y_lo, floor(  min(ref_vals) - 0.3))
    y_hi     <- max(y_hi, ceiling(max(ref_vals) + 0.3))
    y_limits <- c(y_lo, y_hi)
  }

  counts <- df %>% group_by(Organism) %>% summarise(n = n(), .groups = "drop")
  n_hs   <- { v <- counts$n[counts$Organism == "Homo sapiens"];  if (!length(v)) 0L else v }
  n_sc   <- { v <- counts$n[counts$Organism == "Saccharomyces"]; if (!length(v)) 0L else v }

  color_values  <- c("Homo sapiens" = "#56B4E9")
  if (show_yeast) color_values["Saccharomyces"] <- "#E69F00"

  labels_values <- c("Homo sapiens" = paste0("Human: ", n_hs,
                    if (!is.null(ratio_hs)) paste0(" (", ratio_hs, ")") else ""))
  if (show_yeast)
    labels_values["Saccharomyces"] <- paste0("Yeast: ", n_sc,
                    if (!is.null(ratio_sc)) paste0(" (", ratio_sc, ")") else "")

  p_scatter <- ggplot(df, aes(.data[[meanexp_col]], .data[[ratio_col]], color = Organism)) +
    geom_point(alpha = 0.4, size = 0.05) +
    scale_color_manual(values = color_values, labels = labels_values, name = NULL) +
    labs(title = title_suffix,
         y = expression(Log[2]~"[FC]"), x = expression(Log[2]~"[abundance]")) +
    scale_y_continuous(limits = y_limits, breaks = seq(y_limits[1], y_limits[2], by = 1)) +
    scale_x_continuous(limits = x_limits, breaks = scales::pretty_breaks(n = 5)) +
    theme_dia() +
    theme(plot.title           = element_text(lineheight = .8, size = 6),
          axis.title           = element_text(size = 6, face = "bold"),
          legend.position      = "inside",
          legend.justification = c(0, 1),
          legend.direction     = "vertical",
          legend.text          = element_text(size = 6)) +
    guides(color = guide_legend(override.aes = list(size = 1.7)))

  if (!is.null(hline_hs))
    for (v in as.numeric(hline_hs))
      p_scatter <- p_scatter +
        geom_hline(yintercept = v, linetype = "dashed", color = hline_color_hs, linewidth = 0.5)

  if (show_yeast && !is.null(hline_sc) && !is.na(hline_sc))
    for (v in as.numeric(hline_sc))
      p_scatter <- p_scatter +
        geom_hline(yintercept = v, linetype = "dashed", color = hline_color_sc, linewidth = 0.5)

  sd_summary <- df %>%
    group_by(Organism) %>%
    summarise(mean   = mean(  .data[[ratio_col]], na.rm = TRUE),
              sd     = sd(    .data[[ratio_col]], na.rm = TRUE),
              median = median(.data[[ratio_col]], na.rm = TRUE),
              .groups = "drop")

  p_density <- ggplot(df, aes(.data[[ratio_col]], fill = Organism, color = Organism)) +
    geom_density(alpha = 0.3, fill = NA, linewidth = 0.25) +
    geom_vline(data = sd_summary, aes(xintercept = mean, color = Organism),
               linetype = "dotted", linewidth = 0.3) +
    geom_text(data = sd_summary,
              aes(x = mean + 0.3, y = 1.25, label = paste0("SD=", round(sd, 2))),
              colour = "black", size = 1.3, hjust = 0.5, vjust = 0, inherit.aes = FALSE) +
    labs(y = "Density") +
    scale_fill_manual(values = color_values) +
    scale_color_manual(values = color_values) +
    coord_flip() +
    scale_x_continuous(limits = y_limits, breaks = seq(y_limits[1], y_limits[2], by = 1)) +
    scale_y_continuous(expand = c(0, 0.02), limits = c(0, 2), breaks = 0:2) +
    theme_dia() +
    theme(legend.position = "none", axis.title = element_blank(), axis.text.y = element_blank())

  list(plot = p_scatter + p_density + plot_layout(widths = c(4, 1)),
       sd_summary = sd_summary,
       y_limits   = y_limits,
       x_limits   = x_limits)
}

# ── 5d. Single channel-vs-channel BA panel ────
# Axis limits are always computed from the data of this specific channel pair.
draw_ba_panel <- function(df, cA, cB, h_offset, y_offset,
                    ctrl_involved, filter_yeast = FALSE) {

  if (!cA %in% names(df) || !cB %in% names(df)) {
    return(ggplot() + labs(title = paste0(cA, " vs ", cB, " — missing")) +
             theme_void() + theme(plot.title = element_text(size = 6)))
  }

  plot_df <- df %>%
    select(Accession, Species, all_of(c(cA, cB))) %>%
    filter(!is.na(.data[[cA]]), !is.na(.data[[cB]]),
           Species %in% c("Human","Yeast")) %>%
    mutate(mean_ab = (.data[[cA]] + .data[[cB]]) / 2,
           diff_ab =  .data[[cA]] - .data[[cB]])

  if (filter_yeast) plot_df <- plot_df %>% filter(Species == "Human")
  if (nrow(plot_df) == 0) return(ggplot() + labs(title = paste0(cA," vs ",cB," — no data")) + theme_void())

  stats <- plot_df %>%
    group_by(Species) %>%
    summarise(n = n(), median = median(diff_ab, na.rm = TRUE),
              sd = sd(diff_ab, na.rm = TRUE), .groups = "drop")

  stat_label <- stats %>%
    mutate(lbl = paste0(Species, " med=", sprintf("%+.2f", median),
                    " SD=", sprintf("%.2f", sd), " n=", n)) %>%
    pull(lbl) %>% paste(collapse = "\n")

  # Limits from this pair's data, always including reference lines
  x_lo <- quantile(plot_df$mean_ab, 0.005, na.rm = TRUE) - 0.3
  x_hi <- quantile(plot_df$mean_ab, 0.995, na.rm = TRUE) + 0.3
  y_lo <- quantile(plot_df$diff_ab, 0.005, na.rm = TRUE) - 0.3
  y_hi <- quantile(plot_df$diff_ab, 0.995, na.rm = TRUE) + 0.3
  refs <- c(0, h_offset)
  if (!is.na(y_offset) && !ctrl_involved) refs <- c(refs, y_offset)
  y_lo <- min(y_lo, min(refs) - 0.2)
  y_hi <- max(y_hi, max(refs) + 0.2)

  p <- ggplot(plot_df, aes(x = mean_ab, y = diff_ab, color = Species)) +
    geom_point(data = ~ filter(., Species == "Human"), size = 0.15, alpha = 0.25) +
    geom_point(data = ~ filter(., Species == "Yeast"), size = 0.4,  alpha = 0.55) +
    geom_hline(yintercept = 0,        linetype = "dashed", color = "black",
               linewidth = 0.4, alpha = 0.6) +
    geom_hline(yintercept = h_offset, color = "#0072B2", linewidth = 0.2, alpha = 0.85) +
    scale_color_manual(values = species_colors) +
    coord_cartesian(xlim = c(x_lo, x_hi), ylim = c(y_lo, y_hi)) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 4)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
    annotate("label", x = x_lo + 0.02*(x_hi-x_lo), y = y_hi - 0.02*(y_hi-y_lo),
             label = stat_label, hjust = 0, vjust = 1,
             size = 1.6, label.size = 0.15, fill = "white", alpha = 0.8) +
    labs(title = paste0(cA, " vs ", cB),
         x = paste0("mean log\u2082(", cA, ", ", cB, ")"),
         y = paste0("log\u2082(", cA, "/", cB, ")")) +
    theme_minimal(base_size = 7) +
    theme(plot.title      = element_text(size = 6, face = "bold", hjust = 0.5),
          axis.title      = element_text(size = 5.5, face = "bold", colour = "black"),
          axis.text       = element_text(size = 5, colour = "black"),
          axis.line       = element_line(color = "black", linewidth = 0.3),
          axis.ticks      = element_line(color = "black", linewidth = 0.2),
          panel.grid      = element_blank(),
          legend.position = "none")

  if (!ctrl_involved && !is.na(y_offset))
    p <- p + geom_hline(yintercept = y_offset, color = "#D55E00",
                    linewidth = 0.2, alpha = 0.85)
  else if (ctrl_involved && !filter_yeast)
    p <- p + annotate("label",
                    x = x_hi - 0.02*(x_hi-x_lo), y = y_lo + 0.05*(y_hi-y_lo),
                    label = "Yeast \u2248 0\nin ctrl",
                    hjust = 1, vjust = 0, size = 1.6, label.size = 0.15,
                    color = "#E69F00", fill = "white", alpha = 0.8)
  p
}

# ── 5e. Ratio overview boxplot ────
# y_limits = NULL → computed from this dataset's data.
plot_ratio_overview <- function(df, comparisons = NULL,
                    filter_yeast_ctrl = TRUE,
                    y_limits = NULL, title = NULL) {

  ratio_cols <- names(df)[startsWith(names(df), "log2ratio_")]
  if (length(ratio_cols) == 0) { message("No log2ratio columns found."); return(invisible(NULL)) }

  long_df <- df %>%
    select(Accession, Species, all_of(ratio_cols)) %>%
    pivot_longer(all_of(ratio_cols), names_to = "comparison", values_to = "log2ratio") %>%
    mutate(comparison = str_remove(comparison, "^log2ratio_")) %>%
    filter(!is.na(log2ratio), Species %in% c("Human","Yeast"))

  if (!is.null(comparisons)) {
    long_df <- long_df %>% filter(comparison %in% comparisons)
    long_df$comparison <- factor(long_df$comparison, levels = comparisons)
  } else {
    long_df$comparison <- factor(long_df$comparison)
  }

  if (filter_yeast_ctrl)
    long_df <- long_df %>%
      filter(!(Species == "Yeast" & str_detect(as.character(comparison), "H100Y000")))

  comp_levels <- levels(long_df$comparison)

  expected_long <- tibble(comparison = comp_levels) %>%
    mutate(gA = str_extract(comparison, "^.+(?=vs)"),
           gB = str_extract(comparison, "(?<=vs).+"),
           Human = map2_dbl(gA, gB, \(a,b) expected_log2(a, b, "Human")),
           Yeast = map2_dbl(gA, gB, \(a,b) expected_log2(a, b, "Yeast"))) %>%
    select(comparison, Human, Yeast) %>%
    pivot_longer(c(Human, Yeast), names_to = "Species", values_to = "expected") %>%
    filter(!is.na(expected),
           !(filter_yeast_ctrl & Species == "Yeast" & str_detect(comparison, "H100Y000"))) %>%
    mutate(comparison = factor(comparison, levels = comp_levels))

  medians <- long_df %>%
    group_by(comparison, Species) %>%
    summarise(median_val = median(log2ratio, na.rm = TRUE), .groups = "drop")

  # Auto y-limits from this dataset, always covering 0 and expected values
  if (is.null(y_limits)) {
    ref_vals <- c(0, na.omit(expected_long$expected))
    y_min <- min(quantile(long_df$log2ratio, 0.005, na.rm = TRUE) - 0.3, min(ref_vals) - 0.3)
    y_max <- max(quantile(long_df$log2ratio, 0.995, na.rm = TRUE) + 0.5, max(ref_vals) + 0.5)
  } else { y_min <- y_limits[1]; y_max <- y_limits[2] }

  annotation_y <- y_max - 0.03*(y_max - y_min)
  dodge_width  <- 0.7

  ggplot(long_df, aes(x = comparison, y = log2ratio, fill = Species)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.3) +
    geom_boxplot(width = 0.6, position = position_dodge(width = dodge_width),
                 outlier.shape = NA, linewidth = 0.25) +
    geom_errorbar(data = expected_long,
                  aes(x = comparison, ymin = expected, ymax = expected, color = Species),
                  position = position_dodge(width = dodge_width),
                  width = 0.52, linewidth = 0.6, inherit.aes = FALSE) +
    geom_text(data = expected_long,
              aes(x = comparison, y = expected, label = sprintf("%.2f", expected), group = Species),
              position = position_dodge(width = dodge_width),
              hjust = 2, vjust = 0.5, size = 1.6, color = "black", inherit.aes = FALSE) +
    geom_text(data = medians,
              aes(x = comparison, y = annotation_y,
                  label = sprintf("%.2f", median_val), group = Species),
              position = position_dodge(width = dodge_width),
              vjust = 0, size = 1.7, color = "black") +
    scale_fill_manual(values  = c(Human = "#56B4E9", Yeast = "#E69F00")) +
    scale_color_manual(values = c(Human = "#0072B2", Yeast = "#D55E00")) +
    scale_y_continuous(limits = c(y_min, y_max), expand = c(0, 0)) +
    labs(title = title, y = expression(log[2]~"ratio"), fill = NULL, color = NULL) +
    theme_dia(base_size = 7) +
    theme(axis.title.x           = element_blank(),
          axis.text.x            = element_text(angle = 45, hjust = 1, size = 6),
          legend.position        = "inside",
          legend.position.inside = c(0.92, 0.85),
          legend.text            = element_text(size = 5.5),
          legend.key.size        = unit(0.4, "cm"))
}

cat("00_shared_functions.R loaded.\n")
