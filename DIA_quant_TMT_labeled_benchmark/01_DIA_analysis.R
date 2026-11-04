# ── 01_DIA_analysis.R ────
# Master orchestration script for the DIA / Spectronaut pipeline (multi009).
#
# USAGE
#   1. Open R / RStudio with working directory = this folder.
#   2. Edit the "USER CONFIGURATION" block below (path, label, tag).
#   3. source("01_DIA_analysis.R")
#
# The single dataset in `datasets` is processed independently:
#   • y-axis limits are auto-computed from its own data percentiles
#     (no shared axes, no hard-coded ranges).
#   • Specific limits may be forced via optional per-dataset fields.
#
# The experimental design figure is produced ONCE (it is the same for all
# datasets because the sample composition is always the same).
# ────

source("00_shared_functions.R")

# ════
# USER CONFIGURATION  ← edit this section only
# ════
#
# The dataset entry must have at minimum:
#   path  – full or relative path to the Spectronaut non-pivot TSV
#   label – display name shown in plot titles and console output
#   tag   – short identifier used in output file names (no spaces)
#
# Optional per-dataset overrides:
#   ba_y_limits    = c(lo, hi)  – force y-range on group-level BA plots
#   ratio_y_limits = c(lo, hi)  – force y-range on ratio overview boxplot
# Leave out (or set to NULL) to let each plot auto-scale from its own data.
# ────

datasets <- list(

  list(
    path  = "data/20261008_083751_TT_multi009_single_channels_DIA_LFQ_ME_grouping_correct_FASTA_correct_parsing_Report_non_pivot.tsv",
    label = "MS2 + cross-run normalisation",
    tag   = "MS2_Norm"
  )
)

# Output directory (relative to working directory)
plots_dir <- "plots"
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

# ════
# EXPERIMENTAL DESIGN FIGURE  (produced once — identical across datasets)
# ════

cat("\n", strrep("═", 60), "\n", sep = "")
cat("Generating experimental design figure (once)...\n")

# Relative protein amounts per group
# HeLa base = 10 au (equal in every channel); Yeast spike = 0, 1, 2, 10 au
amounts_df <- tibble(
  Group   = factor(c("H100Y000","H100Y010","H100Y020","H100Y100"),
                   levels = c("H100Y000","H100Y010","H100Y020","H100Y100")),
  Human   = c(10, 10, 10, 10),
  Yeast   = c( 0,  1,  2, 10)
) %>%
  pivot_longer(c(Human, Yeast), names_to = "Species", values_to = "Amount") %>%
  mutate(Species = factor(Species, levels = c("Yeast","Human")))

p_comp <- ggplot(amounts_df, aes(x = Group, y = Amount, fill = Species)) +
  geom_col(position = position_dodge(width = 0.65), width = 0.6,
           color = "black", linewidth = 0.2) +
  geom_text(data = filter(amounts_df, Amount > 0),
            aes(label = Amount),
            position = position_dodge(width = 0.65),
            vjust = -0.4, size = 2.0, color = "black") +
  scale_fill_manual(values = c(Human = "#56B4E9", Yeast = "#E69F00")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12)), limits = c(0, NA)) +
  scale_x_discrete(labels = c(
    H100Y000 = "H100Y000\n(no yeast)",
    H100Y010 = "H100Y010\n(yeast \u00d71)",
    H100Y020 = "H100Y020\n(yeast \u00d72)",
    H100Y100 = "H100Y100\n(yeast \u00d710)")) +
  labs(y = "Relative protein amount (a.u.)", fill = NULL, title = "A") +
  theme_minimal(base_size = 7) +
  theme(plot.title      = element_text(size = 7, face = "bold", colour = "black"),
        axis.title.x    = element_blank(),
        axis.title.y    = element_text(size = 6, face = "bold", colour = "black"),
        axis.text       = element_text(size = 6, colour = "black"),
        axis.line       = element_line(color = "black", linewidth = 0.3),
        axis.ticks      = element_line(color = "black", linewidth = 0.25),
        panel.grid      = element_blank(),
        legend.position = "top",
        legend.text     = element_text(size = 5.5),
        legend.key.size = unit(0.35, "cm"),
        legend.margin   = margin(0, 0, 2, 0))

comp_order_exp <- c(
  "H100Y000vsH100Y010", "H100Y000vsH100Y020", "H100Y000vsH100Y100",
  "H100Y100vsH100Y010", "H100Y100vsH100Y020", "H100Y020vsH100Y010"
)

exp_df <- tibble(comparison = comp_order_exp) %>%
  mutate(gA    = str_extract(comparison, "^.+(?=vs)"),
         gB    = str_extract(comparison, "(?<=vs).+"),
         Human = map2_dbl(gA, gB, \(a, b) expected_log2(a, b, "Human")),
         Yeast = map2_dbl(gA, gB, \(a, b) expected_log2(a, b, "Yeast"))) %>%
  select(comparison, Human, Yeast) %>%
  pivot_longer(c(Human, Yeast), names_to = "Species", values_to = "expected") %>%
  mutate(comparison = factor(comparison, levels = rev(comp_order_exp)),
         Species    = factor(Species, levels = c("Human","Yeast")))

# ════
# PER-DATASET FUNCTION
# ════

run_dataset <- function(ds, plots_dir) {

  tag   <- ds$tag
  label <- ds$label

  # Optional y-limit overrides; NULL = auto-computed from this dataset
  ba_y_limits    <- ds[["ba_y_limits"]]
  ratio_y_limits <- ds[["ratio_y_limits"]]

  cat("\n", strrep("═", 60), "\n", sep = "")
  cat("Dataset : ", label, "\n", sep = "")
  cat("Tag     : ", tag,   "\n", sep = "")
  cat("Path    : ", ds$path, "\n", sep = "")
  cat(strrep("─", 60), "\n")

  # ── Per-dataset output subfolder ────
  # All figures specific to this dataset go here.
  # Cross-dataset comparison figures stay in plots_dir directly.
  ds_dir <- file.path(plots_dir, tag)
  dir.create(ds_dir, showWarnings = FALSE, recursive = TRUE)
  cat("  Output folder: ", ds_dir, "\n", sep = "")

  # ── Load & process ────
  res        <- process_sn_file(ds$path, label)
  df_raw     <- res$df_raw   # linear scale, all species (for intensity figs)
  df         <- res$df       # log2 scale, Human + Yeast (for ratio/BA figs)
  ch_present <- res$ch_present

  # ────
  # FIG 2  Channel intensities – total (summed across all channels)
  # ────
  cat("\n  [Fig 2] Channel intensities (total)...\n")

  int_df <- make_channel_intensities(df_raw, group_info, ch_present,
                    channel_order, experiment_order)
  p2 <- plot_channel_intensities(int_df, title = label)
  print(p2)
  out2 <- file.path(ds_dir, paste0("channel_intensities_", tag, ".svg"))
  ggsave(out2, p2, width = 12.7, height = 8.9, units = "cm", dpi = 300)
  cat("    Saved:", file.path(tag, basename(out2)), "\n")

  # ────
  # FIG 3  Channel intensities – by species
  # ────
  cat("\n  [Fig 3] Channel intensities by species...\n")

  int_sp_df <- make_channel_intensities_by_species(df_raw, group_info, ch_present,
                    channel_order, experiment_order)
  p3 <- plot_channel_intensities_by_species(int_sp_df, title = label)
  print(p3)
  out3 <- file.path(ds_dir, paste0("channel_intensities_by_species_", tag, ".svg"))
  ggsave(out3, p3, width = 12.7, height = 8.9, units = "cm", dpi = 300)
  cat("    Saved:", file.path(tag, basename(out3)), "\n")

  # ────
  # FIG 4  Group-level BA plots – one panel per comparison pair
  # ────
  cat("\n  [Fig 4] Group-level BA plots...\n")
  cat("  ── SD summaries ────\n")

  for (pair in comparison_pairs) {

    gA <- pair[1]; gB <- pair[2]
    pn <- paste0(gA, "vs", gB)

    r_hs <- 0                    # Human always 0
    r_sc <- expected_log2(gA, gB, "Yeast")             # NA if ctrl involved

    res4 <- plot_ba_dia(
      df             = df,
      comparison_col = pn,
      ratio_hs       = sprintf("%.3f", r_hs),
      ratio_sc       = if (!is.na(r_sc)) sprintf("%.3f", r_sc) else NA_character_,
      show_yeast     = TRUE,
      y_limits       = ba_y_limits,
      hline_hs       = r_hs,
      hline_sc       = if (!is.na(r_sc)) r_sc else NULL,
      title_suffix   = paste0(gA, " / ", gB, "\n(", label, ")")
    )

    if (is.null(res4)) {
      cat("    SKIP:", pn, "(no data)\n")
      next
    }

    out4 <- file.path(ds_dir, paste0("BA_", pn, "_", tag, ".svg"))
    ggsave(out4, res4$plot, width = 120, height = 80, units = "mm", dpi = 300)
    cat("    Saved:", file.path(tag, basename(out4)), "\n")

    # SD summary
    tbl <- res4$sd_summary %>%
      mutate(comparison = pn, .before = 1) %>%
      mutate(across(where(is.numeric), \(x) round(x, 3)))
    cat("    ", pn, ":\n", sep = "")
    print(as.data.frame(tbl))
    cat("\n")
  }

  # ────
  # FIG 5  Channel-by-channel BA 3×3 grids
  # ────
  cat("\n  [Fig 5] Channel-by-channel BA grids (", nrow(df), " proteins)...\n",
      sep = "")

  for (pair in comparison_pairs) {

    gA <- pair[1]; gB <- pair[2]
    pn <- paste0(gA, "vs", gB)

    # Channels for each group that are actually present in this file
    chA <- intersect(group_channels_map[[gA]], ch_present)
    chB <- intersect(group_channels_map[[gB]], ch_present)

    if (length(chA) == 0 || length(chB) == 0) {
      cat("    SKIP:", pn, "(missing channel group in this dataset)\n")
      next
    }

    ctrl_involved <- (gA == "H100Y000" || gB == "H100Y000")

    h_offset <- 0                    # Human log2 offset = 0
    yA <- yeast_amount[[gA]]; yB <- yeast_amount[[gB]]
    y_offset <- if (yA == 0 || yB == 0) NA_real_ else log2(yA / yB)

    # Build panels row = chA, col = chB
    panels <- vector("list", length(chA) * length(chB))
    idx <- 1L
    for (i in seq_along(chA)) {
      for (j in seq_along(chB)) {
        panels[[idx]] <- draw_ba_panel(
          df            = df,
          cA            = chA[i],
          cB            = chB[j],
          h_offset      = h_offset,
          y_offset      = y_offset,
          ctrl_involved = ctrl_involved
        )
        idx <- idx + 1L
      }
    }

    # Replace any NULL slots with blank panels
    panels <- lapply(panels, function(p) if (is.null(p)) ggplot() + theme_void() else p)

    subtitle <- if (!ctrl_involved && !is.na(y_offset))
      paste0("Expected: Human = 0.000, Yeast = ", sprintf("%.3f", y_offset))
    else
      "Expected: Human = 0.000  (Yeast not comparable \u2014 ctrl has 0 yeast)"

    grid_plot <- wrap_plots(panels, nrow = length(chA), ncol = length(chB)) +
      plot_annotation(
        title    = paste0(gA, " vs ", gB, "  \u2014  ", label),
        subtitle = subtitle,
        theme    = theme(
          plot.title    = element_text(size = 7,  face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 6,  hjust = 0.5, colour = "grey30"),
          plot.margin   = margin(4, 4, 4, 4)
        )
      )

    out5_svg <- file.path(ds_dir, paste0("BA_channels_", pn, "_", tag, ".svg"))
    out5_png <- file.path(ds_dir, paste0("BA_channels_", pn, "_", tag, ".png"))
    ggsave(out5_svg, grid_plot, width = 180, height = 180, units = "mm", dpi = 300)
    ggsave(out5_png, grid_plot, width = 180, height = 180, units = "mm", dpi = 300)
    cat("    Saved:", file.path(tag, basename(out5_svg)),
        "/", basename(out5_png), "\n")
  }

  # ────
  # FIG 6  Ratio overview boxplot
  # ────
  cat("\n  [Fig 6] Ratio overview...\n")

  comp_order_full <- c(
    "H100Y000vsH100Y010", "H100Y000vsH100Y020", "H100Y000vsH100Y100",
    "H100Y100vsH100Y010", "H100Y100vsH100Y020", "H100Y020vsH100Y010"
  )

  p6 <- plot_ratio_overview(
    df                = df,
    comparisons       = comp_order_full,
    filter_yeast_ctrl = TRUE,
    y_limits          = ratio_y_limits,
    title             = label
  )
  print(p6)
  out6 <- file.path(ds_dir, paste0("ratio_overview_", tag, ".svg"))
  ggsave(out6, p6, width = 12, height = 7, units = "cm", dpi = 300)
  cat("    Saved:", file.path(tag, basename(out6)), "\n")

  cat("\n  ", strrep("\u2500", 55), "\n", sep = "")
  cat("  Dataset '", tag, "' complete.\n", sep = "")
  cat("  All figures in: plots/", tag, "/\n", sep = "")
  invisible(NULL)
}

# ════
# MAIN LOOP  ← no edits needed beyond the USER CONFIGURATION block above
# ════

cat("\n", strrep("═", 60), "\n", sep = "")
cat("Running ", length(datasets), " dataset(s)...\n", sep = "")

for (ds in datasets) run_dataset(ds, plots_dir)

cat("\n", strrep("═", 60), "\n", sep = "")
cat("All ", length(datasets), " dataset(s) processed.\n", sep = "")
cat("Plots saved to: ", normalizePath(plots_dir, mustWork = FALSE), "\n", sep = "")
cat(strrep("═", 60), "\n")
