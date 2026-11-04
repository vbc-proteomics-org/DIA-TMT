# ===========================================================================
# Script:  01_TMTc_import_and_noise_subtraction.R
# Purpose: Import the TMTc (complementary-ion) protein tables for the F2500
#          filter, subtract the 126 background per species, and build the
#          summed-intensity plots, protein counts and log2-ratio boxplots.
# Inputs:  data/prot_TMTc/F2500/*.csv ; plots_dir (from script 00).
# Outputs: data_list_prot_TMTc_subtr, plot_list_species_sums_subtr,
#          counts_by_df_TMTc_subtr, ratio_data_TMTc_subtr,
#          plots_ratio_by_filter_subtr and the saved SVG plots.
# Depends: 00_setup_packages_and_experimental_design.R
# ===========================================================================

# --- Data directory ---
data_dir_prot_TMTc_subtr <- file.path("data/prot_TMTc")

# --- Import helpers ---
clean_name <- function(path) tools::file_path_sans_ext(basename(path))

read_one_file <- function(file) {
  ext <- tolower(tools::file_ext(file))

  tryCatch(
    {
      if (ext == "csv") {
        readr::read_csv(file, show_col_types = FALSE)
      } else if (ext %in% c("xlsx", "xls")) {
        readxl::read_excel(file)
      } else {
        warning("Skipping unsupported file type: ", file)
        return(NULL)
      }
    },
    error = function(e) {
      warning("Could not read file: ", file, "\nReason: ", e$message)
      NULL
    }
  )
}

import_data_from_dir <- function(dir_path) {
  files <- list.files(
    path = dir_path,
    pattern = "\\.(csv|xlsx|xls)$",
    full.names = TRUE,
    ignore.case = TRUE
  )

  out <- list()

  for (file in files) {
    nm <- clean_name(file)
    message("Reading and storing: ", basename(dir_path), " / ", nm)

    dat <- read_one_file(file)
    if (!is.null(dat)) out[[nm]] <- dat
  }

  out
}

# --- Experiment-label helpers (dataset name -> readable label) ---
experiment_map <- c(
  "DDA_50ng"    = "DDA 50ng",
  "DDA_100ng"   = "DDA 100ng",
  "DIA_paper_50ng"    = "DIA 50ng",
  "DIA_paper_100ng"   = "DIA 100ng",
  "DIA_high_mz_50ng" = "HM-DIA 50ng",
  "DIA_high_mz_100ng"= "HM-DIA 100ng",
  "DIA_long_grad_50ng"  = "DIA120 50ng",
  "DIA_long_grad_100ng" = "DIA120 100ng"
)

get_experiment_label <- function(dataset_name) {
  matched <- names(experiment_map)[sapply(names(experiment_map), function(pat) grepl(pat, dataset_name))]
  if (length(matched) == 0) return(NA_character_)
  experiment_map[[matched[1]]]
}

# --- Plot function for the protein intensities ---
plot_species_averaged <- function(df_avg, df_raw = NULL, title = "", y_limit = NULL, log10_scale = FALSE) {
  df_long <- df_avg %>%
    pivot_longer(
      cols = c(mean_135N, mean_131N, mean_126),
      names_to = "Channel", values_to = "Intensity"
    ) %>%
    mutate(
      SD = case_when(
        Channel == "mean_135N" ~ df_avg$sd_135N[match(species, df_avg$species)],
        Channel == "mean_131N" ~ df_avg$sd_131N[match(species, df_avg$species)],
        Channel == "mean_126"  ~ df_avg$sd_126[match(species, df_avg$species)]
      ),
      Channel = recode(Channel,
                    mean_135N = "135N",
                    mean_131N = "131N",
                    mean_126  = "126")
    )

  # transform if log10
  if (log10_scale) {
    df_long <- df_long %>%
      mutate(
        SD        = SD / (Intensity * log(10)),  # delta method for SD on log10 scale
        Intensity = log10(Intensity)
      )
    y_label <- "Log10 summed protein intensities"
  } else {
    y_label <- "Summed protein intensities"
  }

  # prepare raw replicate data if provided
  if (!is.null(df_raw)) {
    df_raw_long <- df_raw %>%
      pivot_longer(
        cols = c(sum_135N, sum_131N, sum_126),
        names_to = "Channel", values_to = "Replicate_Value"
      ) %>%
      mutate(Channel = recode(Channel,
                    sum_135N = "135N",
                    sum_131N = "131N",
                    sum_126  = "126"))

    if (log10_scale) {
      df_raw_long <- df_raw_long %>%
        mutate(Replicate_Value = log10(Replicate_Value))
    }
  }

  ggplot(df_long, aes(x = Channel, y = Intensity, fill = species)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7,
             color = "black", linewidth = 0.25) +
    geom_errorbar(aes(ymin = Intensity - SD, ymax = Intensity + SD),
                  position = position_dodge(width = 0.8), width = 0.25,
                  linewidth = 0.3) +
    # individual replicate points
    { if (!is.null(df_raw))
      geom_point(
        data = df_raw_long,
        aes(x = Channel, y = Replicate_Value, fill = species),
        shape = 21,
        color = "black",
        size = 0.8,
        alpha = 0.7,
        position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8)
      )
    } +
    geom_text(
      aes(y     = Intensity + SD,
          label = if (log10_scale) round(Intensity, 2) else scales::scientific(Intensity, digits = 3)),
      position = position_dodge(width = 0.8),
      vjust = -0.5, size = 1.5
    ) +
    labs(title = title, x = NULL, y = y_label, fill = "Species") +
    scale_y_continuous(limits = y_limit, expand = expansion(mult = c(0, 0))) +
    scale_fill_manual(values = c("Human" = "#56B4E9", "Yeast" = "#E69F00")) +
    theme_minimal() +
    theme(
      plot.title   = element_text(hjust = 0.5, size = 6),
      axis.text    = element_text(size = 6, colour = "black"),
      axis.title   = element_text(size = 6, colour = "black"),
      panel.grid   = element_blank(),
      axis.line    = element_line(colour = "black"),
      axis.ticks   = element_line(color = "black", size = 0.25),
      legend.position = "none"
    ) +
    expand_limits(y = max(df_long$Intensity, na.rm = TRUE) * 1.15)
}

# --- Per-protein processing: subtract 126 background per species ---
process_one_df_subtr <- function(df, use_sd = FALSE, min_intensity = 1) {

  # 1. Rename and identify species
  df_work <- df %>%
    dplyr::filter(!grepl("contam", Protein, fixed = TRUE)) %>%
    dplyr::rename(
      `135N_raw` = comp_ion1_intensity,
      `131N_raw` = comp_ion2_intensity,
      `126_raw`  = comp_ion3_intensity
    ) %>%
    dplyr::mutate(is_yeast = grepl("YEAST", Protein, ignore.case = TRUE))

  # 2. Calculate subtrahends per species (median, optionally + SD)
  calc_sub <- function(df_sub) {
    m <- median(df_sub$`126_raw`, na.rm = TRUE)
    s <- if (use_sd) sd(df_sub$`126_raw`, na.rm = TRUE) else 0
    val <- m + s
    if (is.na(val)) 0 else val
  }

  sub_y <- calc_sub(df_work %>% dplyr::filter( is_yeast))
  sub_h <- calc_sub(df_work %>% dplyr::filter(!is_yeast))

  # 3. Apply corrections and calculate ratios
  df_work %>%
    dplyr::mutate(
      current_sub = dplyr::if_else(is_yeast, sub_y, sub_h),

      `135N_corr` = pmax(`135N_raw` - current_sub, 0),
      `131N_corr` = pmax(`131N_raw` - current_sub, 0),
      `126_corr`  = pmax(`126_raw`  - current_sub, 0),

      is_quantifiable = (`135N_corr` > min_intensity) & (`131N_corr` > min_intensity),

      weighted_avg = dplyr::if_else(
        is_yeast,
        (`135N_corr` + `131N_corr`) / (5 + 1),
        (`135N_corr` + `131N_corr`) / (1 + 1)
      ),

      ratio135N = dplyr::if_else(is_quantifiable, `135N_corr` / weighted_avg, NA_real_),
      ratio131N = dplyr::if_else(is_quantifiable, `131N_corr` / weighted_avg, NA_real_)
    ) %>%
    dplyr::select(-current_sub)
}

# --- Import and process (F2500; median + SD subtraction) ---
data_list_prot_TMTc_subtr <- list(
  F2500 = map(
    import_data_from_dir(file.path(data_dir_prot_TMTc_subtr, "F2500")),
    ~ process_one_df_subtr(.x, use_sd = T)
  )
)

# --- Species sums ---

species_sums_TMTc_subtr <- lapply(names(data_list_prot_TMTc_subtr), function(name) {
  df_list <- data_list_prot_TMTc_subtr[[name]]
  result <- lapply(names(df_list), function(df_name) {
    df <- df_list[[df_name]]
    df %>%
      mutate(species = case_when(
        grepl("HUMAN", Protein, ignore.case = TRUE) ~ "Human",
        grepl("YEAST", Protein, ignore.case = TRUE) ~ "Yeast",
        TRUE ~ "Other"
      )) %>%
      filter(species != "Other") %>%
      group_by(species) %>%
      summarise(
        sum_135N = sum(`135N_corr`, na.rm = TRUE),
        sum_131N = sum(`131N_corr`, na.rm = TRUE),
        sum_126  = sum(`126_corr`,  na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(
        log10_135N = log10(sum_135N),
        log10_131N = log10(sum_131N),
        log10_126  = log10(sum_126)
      )
  })
  names(result) <- names(df_list)
  return(result)
})
names(species_sums_TMTc_subtr) <- names(data_list_prot_TMTc_subtr)

# flatten and average across replicates
species_sums_flat_subtr <- imap_dfr(species_sums_TMTc_subtr, function(filter_list, filter_name) {
  imap_dfr(filter_list, function(df, dataset_name) {
    df %>% mutate(
      filter     = filter_name,
      dataset    = dataset_name,
      experiment = get_experiment_label(dataset_name)
    )
  })
})

species_sums_avg_subtr <- species_sums_flat_subtr %>%
  group_by(filter, experiment, species) %>%
  summarise(
    mean_135N = mean(sum_135N, na.rm = TRUE),
    sd_135N   = sd(sum_135N,   na.rm = TRUE),
    mean_131N = mean(sum_131N, na.rm = TRUE),
    sd_131N   = sd(sum_131N,   na.rm = TRUE),
    mean_126  = mean(sum_126,  na.rm = TRUE),
    sd_126    = sd(sum_126,    na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    experiment = factor(experiment, levels = c(
      "DDA 50ng", "DDA 100ng",
      "DIA 50ng", "DIA 100ng",
      "HM-DIA 50ng", "HM-DIA 100ng",
      "DIA120 50ng", "DIA120 100ng"
    )),
    filter = factor(filter, levels = c("F2500"))
  )

# global y max excluding HM-DIA 100ng
global_y_max_subtr <- species_sums_avg_subtr %>%
  filter(species != "Other", experiment != "HM-DIA 100ng") %>%
  summarise(max_val = max(mean_135N + sd_135N,
                    mean_131N + sd_131N,
                    mean_126  + sd_126,
                    na.rm = TRUE)) %>%
  pull(max_val) * 1.05

hm_dia_y_max_subtr <- species_sums_avg_subtr %>%
  filter(species != "Other", experiment == "HM-DIA 100ng") %>%
  summarise(max_val = max(mean_135N + sd_135N,
                    mean_131N + sd_131N,
                    mean_126  + sd_126,
                    na.rm = TRUE)) %>%
  pull(max_val) * 1.05

dir.create(file.path(plots_dir, "species_sums_subtr"), showWarnings = FALSE)

plot_list_species_sums_subtr <- list()

species_sums_avg_subtr %>%
  group_by(filter, experiment) %>%
  group_walk(function(df_slice, keys) {

    y_lim <- if (keys$experiment == "HM-DIA 100ng") {
      c(0, hm_dia_y_max_subtr)
    } else {
      c(0, global_y_max_subtr)
    }

    df_raw_slice <- species_sums_flat_subtr %>%
      filter(filter == keys$filter, experiment == keys$experiment)

    plot_title <- paste0(keys$experiment, " | ", keys$filter)
    file_name  <- paste0("p_int_subtr_", keys$filter, "_", gsub("[ -]", "_", keys$experiment), ".svg")
    list_key   <- paste0("p_int_subtr_", keys$filter, "_", gsub("[ -]", "_", keys$experiment))

    p <- plot_species_averaged(
      df_avg  = df_slice,
      df_raw  = df_raw_slice,
      title   = plot_title,
      y_limit = y_lim,
      log10_scale = FALSE
    )

    plot_list_species_sums_subtr[[list_key]] <<- p

    ggsave(
      filename = file.path(plots_dir, "species_sums_subtr", file_name),
      plot     = p,
      device   = "svg",
      width    = 12.7, height = 10, units = "cm"
    )
    message("Saved: ", file_name)
  })


# --- Protein counts ---

counts_by_df_TMTc_subtr <- imap_dfr(data_list_prot_TMTc_subtr, function(filter_group, filter_name) {
  imap_dfr(filter_group, function(df, file_name) {
    df_quantified <- df %>%
      filter(
        (!is.na(`135N_corr`) & `135N_corr` > 0) |
          (!is.na(`131N_corr`) & `131N_corr` > 0)
      )
    tibble(
      filter   = filter_name,
      df_name  = file_name,
      human    = sum(grepl("HUMAN", df_quantified$Protein, ignore.case = TRUE)),
      yeast    = sum(grepl("YEAST", df_quantified$Protein, ignore.case = TRUE))
    )
  })
}) %>%
  mutate(
    clean_name = str_remove(df_name, "_prot_level_hyperscore$"),
    replicate  = as.integer(str_extract(clean_name, "\\d+$")),
    group      = str_remove(clean_name, "_\\d+$"),
    experiment = sapply(group, get_experiment_label)
  )

counts_by_group_TMTc_subtr <- counts_by_df_TMTc_subtr %>%
  group_by(filter, experiment) %>%
  summarise(
    human_mean   = mean(human),
    human_sd     = sd(human),
    yeast_mean   = mean(yeast),
    yeast_sd     = sd(yeast),
    n_replicates = n(),
    .groups = "drop"
  )

counts_long_TMTc_subtr <- counts_by_group_TMTc_subtr %>%
  select(filter, experiment, human_mean, yeast_mean) %>%
  pivot_longer(
    cols      = c(human_mean, yeast_mean),
    names_to  = "species",
    values_to = "count"
  ) %>%
  mutate(
    species    = recode(species, human_mean = "Human", yeast_mean = "Yeast"),
    species    = factor(species, levels = c("Yeast", "Human")),
    label      = as.character(round(count)),
    experiment = factor(experiment, levels = c(
      "DDA 50ng", "DDA 100ng",
      "DIA 50ng", "DIA 100ng",
      "HM-DIA 50ng", "HM-DIA 100ng",
      "DIA120 50ng", "DIA120 100ng"
    )),
    filter = factor(filter, levels = c("F2500"))
  )

dir.create(file.path(plots_dir, "prot_counts_subtr"), showWarnings = FALSE)

plots_by_filter_subtr <- map(levels(counts_long_TMTc_subtr$filter), function(filt) {
  plot_data <- counts_long_TMTc_subtr %>% filter(filter == filt)

  ggplot(plot_data, aes(x = experiment, y = count, fill = species)) +
    geom_col(position = "stack", color = "black", linewidth = 0.25) +
    geom_text(
      aes(label = label),
      position = position_stack(vjust = 0.5),
      color = "black", size = 1.5
    ) +
    scale_y_continuous(limits = c(0, 4000), expand = c(0, 0)) +
    scale_fill_manual(values = c("Human" = "#56B4E9", "Yeast" = "#E69F00")) +
    labs(title = filt, y = "#quantified proteins", fill = NULL) +
    theme_minimal() +
    theme(
      axis.text.x     = element_text(angle = 45, hjust = 1, vjust = 1, size = 6, color = "black"),
      axis.text.y     = element_text(size = 6, color = "black"),
      axis.title.y    = element_text(size = 6, color = "black"),
      axis.title.x    = element_blank(),
      panel.grid      = element_blank(),
      axis.line       = element_line(color = "black"),
      axis.ticks      = element_line(size = 0.3, color = "black"),
      legend.text     = element_text(size = 5),
      legend.key.size = unit(0.3, "cm"),
      plot.title      = element_text(size = 6, hjust = 0.5)
    )
})

names(plots_by_filter_subtr) <- levels(counts_long_TMTc_subtr$filter)

iwalk(plots_by_filter_subtr, function(p, filter_name) {
  ggsave(
    filename = file.path(plots_dir, "prot_counts_subtr", paste0("prot_counts_TMTc_subtr_", filter_name, ".svg")),
    plot     = p,
    device   = "svg",
    width    = 12.7, height = 10, units = "cm"
  )
  message("Saved: prot_counts_TMTc_subtr_", filter_name)
})

# --- Ratios ---

ratio_data_TMTc_subtr <- imap_dfr(data_list_prot_TMTc_subtr, function(filter_group, filter_name) {
  imap_dfr(filter_group, function(df, file_name) {
    df %>%
      filter(!is.na(ratio135N) & !is.na(ratio131N)) %>%
      select(Protein, ratio135N, ratio131N) %>%
      mutate(
        species = case_when(
          grepl("HUMAN", Protein, ignore.case = TRUE) ~ "Human",
          grepl("YEAST", Protein, ignore.case = TRUE) ~ "Yeast",
          TRUE ~ "Other"
        )
      ) %>%
      filter(species %in% c("Human", "Yeast")) %>%
      select(species, ratio135N, ratio131N) %>%
      mutate(filter = filter_name, df_name = file_name)
  })
}) %>%
  mutate(
    clean_name     = str_remove(df_name, "_prot_level_hyperscore$"),
    replicate      = as.integer(str_extract(clean_name, "_(\\d+)$", group = 1)),
    group          = str_remove(clean_name, "_\\d+$"),
    experiment     = sapply(group, get_experiment_label),
    log2_ratio135N = log2(ratio135N),
    log2_ratio131N = log2(ratio131N)
  ) %>%
  filter(is.finite(log2_ratio135N) & is.finite(log2_ratio131N)) %>%
  mutate(
    experiment = factor(experiment, levels = c(
      "DDA 50ng", "DDA 100ng",
      "DIA 50ng", "DIA 100ng",
      "HM-DIA 50ng", "HM-DIA 100ng",
      "DIA120 50ng", "DIA120 100ng"
    )),
    filter = factor(filter, levels = c("F2500"))
  )

median_df_TMTc_subtr <- ratio_data_TMTc_subtr %>%
  group_by(filter, experiment, species) %>%
  summarise(
    median_log2_ratio135N = median(log2_ratio135N, na.rm = TRUE),
    median_log2_ratio131N = median(log2_ratio131N, na.rm = TRUE),
    .groups = "drop"
  )

ratio_data_TMTc_long_subtr <- ratio_data_TMTc_subtr %>%
  pivot_longer(
    cols      = c(log2_ratio135N, log2_ratio131N, ratio135N, ratio131N),
    names_to  = "metric",
    values_to = "value"
  ) %>%
  mutate(
    channel = ifelse(grepl("135N", metric), "135N", "131N"),
    type    = ifelse(grepl("log2", metric), "log2_ratio", "normal_ratio")
  )

median_df_TMTc_long_subtr <- median_df_TMTc_subtr %>%
  pivot_longer(
    cols      = c(median_log2_ratio135N, median_log2_ratio131N),
    names_to  = "channel",
    values_to = "median_value"
  ) %>%
  mutate(channel = recode(channel,
                    median_log2_ratio135N = "135N",
                    median_log2_ratio131N = "131N"))

dir.create(file.path(plots_dir, "ratios_subtr"), showWarnings = FALSE)

filters_to_plot_subtr <- levels(ratio_data_TMTc_long_subtr$filter)
if (is.null(filters_to_plot_subtr)) filters_to_plot_subtr <- unique(ratio_data_TMTc_long_subtr$filter)

plots_ratio_by_filter_subtr <- map(filters_to_plot_subtr, function(filt) {

  df_plot <- ratio_data_TMTc_long_subtr %>%
    filter(filter == filt, type == "log2_ratio")

  df_meds <- median_df_TMTc_long_subtr %>%
    filter(filter == filt) %>%
    mutate(label = sprintf("%.2f", median_value))

  ggplot(df_plot, aes(x = experiment, y = value, fill = species)) +
    geom_boxplot(
      position      = position_dodge(width = 0.85),
      outlier.shape = NA,
      linewidth     = 0.2,
      median.linewidth = 0.3
    ) +
    geom_text(
      data = df_meds,
      aes(x = experiment, y = median_value, label = label, group = species),
      position    = position_dodge(width = 0.85),
      vjust       = -0.6, size = 1.3,
      inherit.aes = FALSE
    ) +
    facet_wrap(~ channel, nrow = 1) +
    geom_hline(yintercept = log2(5), linetype = "dashed", color = "#E69F00", linewidth = 0.4) +
    geom_hline(yintercept = 0,       linetype = "dashed", color = "#56B4E9", linewidth = 0.4) +
    scale_fill_manual(values = c("Human" = "#56B4E9", "Yeast" = "#E69F00")) +
    scale_y_continuous(limits = c(-2, 3)) +
    labs(title = filt, x = NULL, y = "log2(Ratio)", fill = NULL) +
    theme_minimal(base_size = 6) +
    theme(
      plot.title  = element_text(size = 6, hjust = 0.5, color = "black"),
      axis.text.x = element_text(size = 6, angle = 45, hjust = 1, color = "black"),
      panel.grid  = element_blank(),
      axis.line   = element_line(color = "black"),
      axis.ticks  = element_line(linewidth = 0.3, color = "black")
    )
})

names(plots_ratio_by_filter_subtr) <- filters_to_plot_subtr

iwalk(plots_ratio_by_filter_subtr, function(p, filt) {
  ggsave(
    filename = file.path(plots_dir, "ratios_subtr", paste0("ratios_TMTc_subtr_", filt, ".svg")),
    plot     = p,
    device   = "svg",
    width    = 12.7, height = 10, units = "cm"
  )
  message("Saved: ratios_TMTc_subtr_", filt)
})

plots_ratio_by_filter_subtr$F2500
