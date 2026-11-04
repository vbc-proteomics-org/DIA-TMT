# =============================================================================
# Script 02: Constants, functions and data import
# Description: Imports the PSM and protein result tables (baseline and
#              standard_StoN_grouping_20), cleans and renames columns, removes
#              contaminants, annotates organism, computes per-group means / CVs /
#              ratios, and builds the long-format tables consumed downstream.
#              Also defines the filtering, log2-transform and counting functions
#              and constant parameters used by later scripts.
# Depends on: 01_experimental_design_and_packages.R (channel_names_groups)
# =============================================================================

# --- Data and plot directories -----------------------------------------------
# NOTE: Paths are relative to the working directory. Update to match your local
# environment.

data_dir_baseline <- file.path("data", "baseline")
data_dir_standard_StoN_grouping_20 <- file.path("data", "standard_StoN_grouping_20")

plots_dir <- "plots"

# List all .txt result files in each data subdirectory.
files_baseline <- list.files(path = data_dir_baseline, pattern = "\\.txt$", full.names = TRUE)
files_standard_StoN_grouping_20 <- list.files(path = data_dir_standard_StoN_grouping_20, pattern = "\\.txt$", full.names = TRUE)

# --- Data import -------------------------------------------------------------

# Import each file into a named list; names are the file base names.
import_data <- function(files, existing_data = list()) {
  existing_names <- names(existing_data)
  
  # Keep the full base name (no extension, no shortening).
  clean_name <- function(file) {
    tools::file_path_sans_ext(basename(file))
  }
  
  file_names <- sapply(files, clean_name, USE.NAMES = FALSE)
  new_files <- files[!file_names %in% existing_names]
  new_names <- file_names[!file_names %in% existing_names]
  
  for (i in seq_along(new_files)) {
    file <- new_files[i]
    full_name <- new_names[i]
    
    message("Reading and storing: ", full_name)
    
    data <- tryCatch(
      read_tsv(file, show_col_types = FALSE),
      error = function(e) {
        warning("Could not read file: ", file, "\nReason: ", e$message)
        return(NULL)
      }
    )
    
    if (!is.null(data)) {
      existing_data[[full_name]] <- data
    }
  }
  
  return(existing_data)
}

# First load: all data frames stored in a list per condition group.
data_list_baseline <- import_data(files_baseline)
data_list_standard_StoN_grouping_20 <- import_data(files_standard_StoN_grouping_20)

# Split each list into PSM and protein data frames.
psm_data_frames_baseline <- data_list_baseline[names(data_list_baseline)[grepl("PSM", names(data_list_baseline))]]
protein_data_frames_baseline <- data_list_baseline[names(data_list_baseline)[grepl("Proteins", names(data_list_baseline))]]

psm_data_frames_standard_StoN_grouping_20 <- data_list_standard_StoN_grouping_20[names(data_list_standard_StoN_grouping_20)[grepl("PSM", names(data_list_standard_StoN_grouping_20))]]
protein_data_frames_standard_StoN_grouping_20 <- data_list_standard_StoN_grouping_20[names(data_list_standard_StoN_grouping_20)[grepl("Proteins", names(data_list_standard_StoN_grouping_20))]]

psm_names <- c(
  "psm_data_frames_baseline", "psm_data_frames_standard_StoN_grouping_20"
)

protein_names <- c(
  "protein_data_frames_baseline", "protein_data_frames_standard_StoN_grouping_20"
)

all_names <- c(psm_names, protein_names)

# Fetch all list objects from the environment.
all_lists <- mget(all_names, envir = .GlobalEnv)

# --- Column renaming ---------------------------------------------------------

# Rename abundance, purity, resolution, ratio, injection-time and S/N columns to
# short, consistent names.
rename_abundance_columns <- function(df) {
  # Define TMT channels 
  tmt_channels <- c("126", "127N", "127C", "128N", "128C", "129N", 
                    "129C", "130N", "130C", "131N", "131C", "132N",
                    "132C", "133N", "133C", "134N", "134C", "135N")
  
  # 1. Rename "Abundance: " columns (PSMs) 
  abundance_cols_grouped <- paste0("Abundance: ", tmt_channels)
  matched_grouped <- abundance_cols_grouped[abundance_cols_grouped %in% names(df)]
  new_names_grouped <- tmt_channels[abundance_cols_grouped %in% names(df)]
  names(df)[match(matched_grouped, names(df))] <- new_names_grouped
  
  # 2. Rename "Abundance: ..." columns (Proteins, possibly with F2 and extra text)
  abundance_cols_psm <- grep("^Abundance:.*\\d+[NC]?", names(df), value = TRUE)
  for (col in abundance_cols_psm) {
    match <- stringr::str_match(col, ".*?(\\d{3}[NC]?)[, ]")[, 2]
    if (!is.na(match) && match %in% tmt_channels) {
      names(df)[names(df) == col] <- match
    }
  }
  
  # 3. Rename "Annotated TMTpro-<channel>" -> "Pure-<channel>" (purity PSM frames).
  annotated_psm_cols <- grep("^Annotated TMTpro-\\d{3}[NC]?$", names(df), value = TRUE)
  for (col in annotated_psm_cols) {
    match <- stringr::str_match(col, "Annotated TMTpro-(\\d{3}[NC]?)")[, 2]
    if (!is.na(match) && match %in% tmt_channels) {
      names(df)[names(df) == col] <- paste0("Pure-", match)
    }
  }
  
  # 4. Rename "Annotated protein-level TMTpro-<channel>" -> "Pure-<channel>" (purity protein frames).
  annotated_protein_cols <- grep("^Annotated protein-level TMTpro-\\d{3}[NC]?$", names(df), value = TRUE)
  for (col in annotated_protein_cols) {
    match <- stringr::str_match(col, "Annotated protein-level TMTpro-(\\d{3}[NC]?)")[, 2]
    if (!is.na(match) && match %in% tmt_channels) {
      names(df)[names(df) == col] <- paste0("Pure-", match)
    }
  }
  
  # 5. Rename "Abundance Ratio: (1, mix2) / (1, mix1)" -> "mix2/mix1".
  # Only relevant when using the ratios pre-calculated by Proteome Discoverer.
  ratio_cols <- grep("^Abundance Ratio: \\(\\d+, mix\\d+\\) / \\(\\d+, mix\\d+\\)$", names(df), value = TRUE)
  for (col in ratio_cols) {
    matches <- stringr::str_match(col, "^Abundance Ratio: \\(\\d+, (mix\\d+)\\) / \\(\\d+, (mix\\d+)\\)$")
    if (!is.na(matches[1, 2]) && !is.na(matches[1, 3])) {
      new_name <- paste0(matches[1, 2], "/", matches[1, 3])
      names(df)[names(df) == col] <- new_name
    }
  }
  
  # 6. Rename "RESGUI_<channel> Resolution" -> "res_<channel>".
  res_cols <- grep("^RESGUI_\\d{3}[NC]? Resolution$", names(df), value = TRUE)
  for (col in res_cols) {
    matches <- stringr::str_match(col, "^RESGUI_(\\d{3}[NC]?) Resolution$")
    if (!is.na(matches[1, 2])) {
      new_name <- paste0("res_", matches[1, 2])
      names(df)[names(df) == col] <- new_name
    }
  }
  
  # 7. Rename "Ion Inject Time [ms]" -> "IT[ms]".
  if ("Ion Inject Time [ms]" %in% names(df)) {
    names(df)[names(df) == "Ion Inject Time [ms]"] <- "IT[ms]"
  }
  
  # 8. Rename "Condition_SN_<suffix>" -> "SN_<suffix>".
  StoN_cols <- grep("^Condition_SN_", names(df), value = TRUE)  
  
  for (col in StoN_cols) {
    matches <- str_match(col, "^Condition_SN_(.*)$")
    suffix <- matches[1, 2] 
    
    if (!is.na(suffix)) {  
      new_name <- paste0("SN_", suffix)  
      names(df)[names(df) == col] <- new_name  
    }
  }
  
  return(df)
}

all_lists <- lapply(all_lists, function(lst) {
  lapply(lst, rename_abundance_columns)
})

list2env(all_lists, envir = .GlobalEnv)

# --- Column definitions and subsetting ---------------------------------------

# Abundance columns.
abundance_cols <- c("126", "127N", "127C", "128N", "128C", "129N", 
                    "129C", "130N", "130C", "131N", "131C", "132N",
                    "132C", "133N", "133C", "134N", "134C", "135N")

pure_abund_cols <- c("Pure-126", "Pure-127N", "Pure-127C", "Pure-128N", "Pure-128C", "Pure-129N", 
                    "Pure-129C", "Pure-130N", "Pure-130C", "Pure-131N", "Pure-131C", "Pure-132N",
                    "Pure-132C", "Pure-133N", "Pure-133C", "Pure-134N", "Pure-134C", "Pure-135N")

# Resolution columns (detected from an example data frame).
resolution_cols <- grep("^res_", names(psm_data_frames_baseline[[1]]), value = TRUE)

StoN_cols <- c("SN_H000","SN_H100", "SN_H050", "SN_H010")

psm_purity_columns_to_keep <- c("Annotated Sequence", "First Scan", "Master Protein Accessions", "Ion Inject Time ms", 
                    "PSM Ambiguity", "Average Reporter S/N", "Normalized CHIMERYS Coefficient", "Co-Isolation Purity",
                    "Master Protein Descriptions", 
                    resolution_cols, pure_abund_cols, StoN_cols)

protein_purity_columns_to_keep <- c("Accession", "Description",
                    pure_abund_cols)

# Keep only the requested columns that exist in df.
keep_columns <- function(df, columns_to_keep) {
  df[, intersect(names(df), columns_to_keep), drop = FALSE]
}

# Apply column subsetting to PSM lists.
psm_lists <- mget(psm_names, envir = .GlobalEnv)
psm_lists <- lapply(psm_lists, function(lst) {
  lapply(lst, keep_columns, columns_to_keep = psm_purity_columns_to_keep)
})
list2env(psm_lists, envir = .GlobalEnv)

# Apply column subsetting to protein lists.
protein_lists <- mget(protein_names, envir = .GlobalEnv)
protein_lists <- lapply(protein_lists, function(lst) {
  lapply(lst, keep_columns, columns_to_keep = protein_purity_columns_to_keep)
})
list2env(protein_lists, envir = .GlobalEnv)

# --- Remove contaminants -----------------------------------------------------

# Drop rows whose accession column flags a contaminant.
remove_cont <- function(df) {
  if ("Master Protein Accessions" %in% names(df)) {
    df <- df[!grepl("cont", df$`Master Protein Accessions`, ignore.case = TRUE), ]
  } else if ("Accession" %in% names(df)) {
    df <- df[!grepl("cont", df$Accession, ignore.case = TRUE), ]
  }
  
  return(df)
}

all_lists <- mget(all_names, envir = .GlobalEnv)

all_lists <- lapply(all_lists, function(lst) {
  lapply(lst, remove_cont)
})

list2env(all_lists, envir = .GlobalEnv)

# --- Annotate organism -------------------------------------------------------

# Add an Organism column derived from the description column.
add_organism <- function(df) {
  desc_col <- NULL
  
  if ("Description" %in% names(df)) {
    desc_col <- "Description"
  } else if ("Master Protein Descriptions" %in% names(df)) {
    desc_col <- "Master Protein Descriptions"
  }
  
  if (!is.null(desc_col)) {
    df <- df %>%
      mutate(Organism = case_when(
        str_detect(.data[[desc_col]], "OS=Homo sapiens") ~ "Human",
        str_detect(.data[[desc_col]], "OS=Saccharomyces") ~ "Yeast",
        TRUE ~ "Other"
      ))
  } else {
    df$Organism <- NA
  }
  
  return(df)
}

all_lists <- mget(all_names, envir = .GlobalEnv)

all_lists <- lapply(all_lists, function(lst) {
  lapply(lst, add_organism)
})

list2env(all_lists, envir = .GlobalEnv)

# --- Replace zeros with NA ---------------------------------------------------
# replace 0 with NA in the abundance and resolution columns.

replace_zeros_with_na <- function(df, target_cols) {
  cols <- intersect(names(df), target_cols)  # only process existing columns
  df[cols] <- lapply(df[cols], function(x) {
    x[x == 0] <- NA
    x
  })
  df
}

all_lists <- mget(all_names, envir = .GlobalEnv)

all_lists <- lapply(all_lists, function(lst) {
  lapply(lst, replace_zeros_with_na, target_cols = c(pure_abund_cols, resolution_cols))
})

list2env(all_lists, envir = .GlobalEnv)

# --- Means, CVs and ratios ---------------------------------------------------
# Means / CVs / ratios are computed per row only when all replicate values in a
# group are present (no NA).

calculate_means_ratios_wide <- function(df, channel_info, exclude_channels = NULL, prefix = NULL) {
  
  # Convert selected columns to a numeric matrix.
  get_numeric_matrix <- function(cols) {
    mat <- df[, cols, drop = FALSE]                  # keep as data.frame
    mat <- as.data.frame(lapply(mat, as.numeric))    # convert to numeric
    as.matrix(mat)
  }
  
  # Apply channel prefix (e.g. "Pure-") if requested.
  if (!is.null(prefix)) {
    channel_info <- channel_info %>%
      mutate(label_prefixed = paste0(prefix, label))
    label_col <- "label_prefixed"
  } else {
    label_col <- "label"
  }
  
  # Channel names per experiment group.
  group_channels <- channel_info %>%
    filter(experiment %in% c("H010", "H050", "H100")) %>%
    split(.$experiment) %>%
    lapply(function(x) x[[label_col]])
  
  # Optionally exclude specific channels before computing stats.
  if (!is.null(exclude_channels)) {
    group_channels <- lapply(group_channels, function(cols) setdiff(cols, exclude_channels))
  }
  
  # Numeric matrix per experiment.
  mat_list <- lapply(group_channels, function(ch) get_numeric_matrix(intersect(ch, names(df))))
  
  # Row means and CV, only for rows with all values present.
  compute_valid_stats <- function(mat) {
    all_present <- !apply(mat, 1, function(r) any(is.na(r)))
    means <- cv <- rep(NA_real_, nrow(mat))
    if (any(all_present)) {
      means[all_present] <- rowMeans(mat[all_present, , drop = FALSE])
      cv_percent <- function(m) 100 * (matrixStats::rowSds(m) / matrixStats::rowMeans2(m))
      cv[all_present] <- cv_percent(mat[all_present, , drop = FALSE])
    }
    list(means = means, cv = cv)
  }
  
  stats <- lapply(mat_list, compute_valid_stats)
  
  for (grp in names(stats)) {
    df[[paste0("unnorm_", grp)]] <- stats[[grp]]$means
    df[[paste0("cv_", grp)]] <- stats[[grp]]$cv
  }
  
  # Ratios and mean expressions, only where both inputs are non-NA.
  safe_div <- function(a, b) ifelse(!is.na(a) & !is.na(b), a / b, NA_real_)
  safe_mean <- function(a, b) ifelse(!is.na(a) & !is.na(b), (a + b) / 2, NA_real_)
  
  df$unnorm_H100vsH010 <- safe_div(df$unnorm_H100, df$unnorm_H010)
  df$unnorm_H100vsH050 <- safe_div(df$unnorm_H100, df$unnorm_H050)
  df$unnorm_H050vsH010 <- safe_div(df$unnorm_H050, df$unnorm_H010)
 
  df$unnorm_meanExp_H100vsH010 <- safe_mean(df$unnorm_H100, df$unnorm_H010)
  df$unnorm_meanExp_H100vsH050 <- safe_mean(df$unnorm_H100, df$unnorm_H050)
  df$unnorm_meanExp_H050vsH010 <- safe_mean(df$unnorm_H050, df$unnorm_H010)
  
  
  if (!is.null(exclude_channels)) {
    message("Excluded channels from mean and CV calculations: ", paste(exclude_channels, collapse = ", "))
  }
  
  return(df)
}

all_lists <- mget(all_names, envir = .GlobalEnv)

all_lists <- lapply(all_lists, function(lst) {
  lapply(lst, calculate_means_ratios_wide,
         channel_info = channel_names_groups,
         prefix = "Pure-",
         exclude_channels = NULL)
})

list2env(all_lists, envir = .GlobalEnv)

# --- Label and replicate assignment ------------------------------------------

# Map raw dataset file names to human-readable method labels.
label_patterns <- c(
  
  "DIA_mz1_ij8_53c75ng" = "DIA 1.0 m/z \nIT8, 53.75ng",
  
  "DIA_mz1_ij8_HCD29_53c75ng" = "DIA 1.0 m/z \nHCD29, 53.75ng",
 
  "DIA_mz1_ij8_noLoop_53c75ng" = "DIA 1.0 m/z \nALL, 53.75ng",
  
  "DDA_mz0c5_ij20_53c75ng" = "DDA 0.5 m/z \nIT20, 53.75ng"
 
)

# Add dataset, label and replicate columns to each data frame.
add_labels_and_replicates <- function(df, dataset_name, label_patterns) {
  dt <- data.table::as.data.table(df)
  
  dt[, dataset := dataset_name]
  
  # Strip known suffixes (longest first) to obtain the clean dataset name.
  dataset_clean <- dataset_name
  dataset_clean <- gsub("_Proteins_purity_tmt_quant$",            "", dataset_clean)
  dataset_clean <- gsub("_PSMs_purity_tmt_quant_conditions$",     "", dataset_clean)
  dataset_clean <- gsub("_purity_tmt_quant$",                    "", dataset_clean)
  dataset_clean <- gsub("_Proteins$",                    "", dataset_clean)
  dataset_clean <- gsub("_PSMs$",                    "", dataset_clean)
  
  # Replicate extraction: first try the PSM-specific "_X_PSMs" pattern.
  replicate <- stringr::str_extract(dataset_name, "_\\d+_PSMs")
  replicate <- stringr::str_extract(replicate, "\\d+")
  replicate <- as.integer(replicate)
  
  # Fallback for proteins and DIA runs: last "_X" of the cleaned name.
  if (is.na(replicate)) {
    replicate <- stringr::str_extract(dataset_clean, "_(\\d+)$")
    replicate <- as.integer(gsub("_", "", replicate))
  }
  
  # Label matching: strip the replicate suffix, then match patterns.
  dataset_base <- sub("_\\d+$", "", dataset_clean)
  
  label <- NA_character_
  for (pattern in names(label_patterns)) {
    if (grepl(pattern, dataset_base, fixed = TRUE)) {
      label <- label_patterns[[pattern]]
      break
    }
  }
  
  dt[, label := label]
  dt[, replicate := replicate]
  
  return(dt)
}

rm(all_lists, psm_lists, protein_lists)

# Assign labels and replicate numbers to every data frame.
protein_data_frames_baseline <- Map(add_labels_and_replicates, protein_data_frames_baseline,
                    names(protein_data_frames_baseline), MoreArgs = list(label_patterns = label_patterns))
psm_data_frames_baseline <- Map(add_labels_and_replicates, psm_data_frames_baseline,
                    names(psm_data_frames_baseline), MoreArgs = list(label_patterns = label_patterns))

protein_data_frames_standard_StoN_grouping_20 <- Map(add_labels_and_replicates, protein_data_frames_standard_StoN_grouping_20, 
                    names(protein_data_frames_standard_StoN_grouping_20), MoreArgs = list(label_patterns = label_patterns))
psm_data_frames_standard_StoN_grouping_20 <- Map(add_labels_and_replicates, psm_data_frames_standard_StoN_grouping_20, 
                    names(psm_data_frames_standard_StoN_grouping_20), MoreArgs = list(label_patterns = label_patterns))

# --- Long-format tables ------------------------------------------------------

long_psm_baseline <- rbindlist(psm_data_frames_baseline, 
                    use.names = TRUE, fill = TRUE)
long_prot_baseline <- rbindlist(protein_data_frames_baseline,
                    use.names = TRUE, fill = TRUE)

long_psm_standard_StoN_grouping_20 <- rbindlist(psm_data_frames_standard_StoN_grouping_20 , 
                    use.names = TRUE, fill = TRUE)
long_prot_standard_StoN_grouping_20  <- rbindlist(protein_data_frames_standard_StoN_grouping_20 , 
                    use.names = TRUE, fill = TRUE)

# --- Filtering functions -----------------------------------------------------

# Filter a long PSM table.
# purity_filter          numeric purity threshold; applied to "Co-Isolation Purity" (also drops NA).
# chimerys_coeff_filter  numeric threshold for the normalized CHIMERYS coefficient.
# unambiguous_filter     TRUE drops rows containing both species in the descriptions.
# SN_filter              numeric S/N threshold; applied to "Average Reporter S/N".
# resolution_filter      numeric resolution threshold applied to the res_ columns.
filter_psm_long_dataframe <- function(df,
                    purity_filter = NULL,
                    chimerys_coeff_filter = NULL,
                    unambiguous_filter = TRUE,
                    SN_filter = NULL,
                    log2_transform = FALSE,
                    resolution_filter = NULL) {
  df <- copy(df)
  
  # Conditional numeric filter helper.
  apply_filter <- function(df, col, min_value = NULL) {
    if (col %in% names(df)) {
      df <- df[!is.na(df[[col]])]
      if (!is.null(min_value)) {
        df <- df[df[[col]] >= min_value]
      }
    }
    df
  }
  
  df <- apply_filter(df, "Co-Isolation Purity", purity_filter)
  df <- apply_filter(df, "Average Reporter S/N", SN_filter)
  df <- apply_filter(df, "Normalized CHIMERYS Coefficient", chimerys_coeff_filter)
  
  # Remove mixed-species PSMs.
  if (unambiguous_filter && "Master Protein Descriptions" %in% names(df)) {
    both_species <- grepl("OS=Homo sapiens", df$`Master Protein Descriptions`, ignore.case = TRUE) &
      grepl("OS=Saccharomyces cerevisiae", df$`Master Protein Descriptions`, ignore.case = TRUE)
    df <- df[!both_species]
  }
  
  # Resolution filter (excludes label-only and interference-prone channels).
  if (!is.null(resolution_filter)) {
    resolution_cols <- grep("^res_", names(df), value = TRUE)
    exclude <- c("res_126", "res_127N", "res_127C", "res_129C", "res_131C", "res_133C")
    resolution_cols <- setdiff(resolution_cols, exclude)
    
    if (length(resolution_cols) > 0) {
      df <- df[apply(df[, ..resolution_cols], 1, function(row) {
        all(!is.na(row) & row >= resolution_filter)
      })]
    }
  }
  
  # Optional log2 transform of abundance and ratio columns.
  if (log2_transform) {
    base_abundance_columns <- c("126", "127N", "127C", "128N", "128C", "129N", 
                    "129C", "130N", "130C", "131N", "131C", "132N",
                    "132C", "133N", "133C", "134N", "134C", "135N")
    prefixed_abundance_columns <- paste0("Pure-", base_abundance_columns)
    
    abundance_cols_present <- intersect(c(base_abundance_columns, prefixed_abundance_columns),
                    names(df))
    
    ratio_cols_present <- grep("^(unnorm_|H100Y|meanExp|unnorm_meanExp)", names(df), value = TRUE)
    
    all_cols <- unique(c(abundance_cols_present, ratio_cols_present))
    
    if (length(all_cols) > 0) {
      df[, (all_cols) := lapply(.SD, function(x) {
        x <- suppressWarnings(as.numeric(x))
        x[x <= 0 | is.infinite(x)] <- NA  # avoid log2(0) or Inf
        log2(x)
      }), .SDcols = all_cols]
    } else {
      warning("No quantitative columns found for log2 transformation.")
    }
  }
  
  df
}

# Filter a long protein table on empty abundance rows; optionally drop
# mixed-species proteins and log2-transform abundance / ratio columns.
filter_and_log2_transform_long_prot <- function(df,
                    datasets_to_keep = NULL,
                    datasets_to_exclude = NULL,
                    quan_values_filter = FALSE,
                    unambiguous_filter = TRUE,
                    log2_transform = TRUE) {
  df <- copy(df)
  
  # 1. Optional dataset filtering.
  if (!is.null(datasets_to_keep)) {
    df <- df[df$dataset %in% datasets_to_keep, ]
  } else if (!is.null(datasets_to_exclude)) {
    df <- df[!df$dataset %in% datasets_to_exclude, ]
  }
  
  # 2. Abundance columns (base and "Pure-" prefixed).
  base_abundance_columns <- c("126", "127N", "127C", "128N", "128C", "129N", 
                    "129C", "130N", "130C", "131N", "131C", "132N",
                    "132C", "133N", "133C", "134N", "134C", "135N")
  prefixed_abundance_columns <- paste0("Pure-", base_abundance_columns)
  
  # Channels that actually carry proteins (used for the quant filter).
  used_abund_cols <- c("127N", "127C", "128N", 
                    "128C", "129C", "130N", 
                    "130C","131N", "132C", 
                    "133N", "133C", "134N")
  used_prefixed_abund_cols <- paste0("Pure-", used_abund_cols)
  
  used_abund_cols_present <- intersect(c(used_abund_cols, used_prefixed_abund_cols),
                    colnames(df))
  
  abundance_cols_present <- intersect(c(base_abundance_columns, prefixed_abundance_columns),
                    colnames(df))

  
  # 3. Optional unambiguous (mixed-species) filter.
  if (unambiguous_filter && "Description" %in% names(df)) {
    has_both_species <- grepl("OS=Homo sapiens", df$Description, ignore.case = TRUE) &
      grepl("OS=Saccharomyces cerevisiae", df$Description, ignore.case = TRUE)
    df <- df[!has_both_species, ]
  }
  
  # 4. Optional quantitative abundance filter.
  if (quan_values_filter) {
    if (length(used_abund_cols_present) == 0) {
      warning("No abundance columns found for quantitative filtering.")
    } else {
      abundance_data <- df[, used_abund_cols_present, with = FALSE]
      
      is_empty_row <- apply(abundance_data, 1, function(row) {
        all(is.na(row) | row == 0 | row == "")
      })
      
      df <- df[!is_empty_row, ]
    }
  }
  
  # 5. Optional log2 transform.
  if (log2_transform) {
    ratio_cols_present <- grep("^(unnorm_|H100Y|meanExp|unnorm_meanExp)", names(df), value = TRUE)
    
    all_cols_to_transform <- unique(c(abundance_cols_present, ratio_cols_present))
    
    if (length(all_cols_to_transform) > 0) {
      for (col in all_cols_to_transform) {
        df[[col]] <- suppressWarnings(as.numeric(df[[col]]))
        df[[col]][df[[col]] <= 0 | is.infinite(df[[col]])] <- NA
        df[[col]] <- log2(df[[col]])
      }
    } else {
      warning("No quantitative columns found for log2 transformation.")
    }
  }
  
  return(df)
}

# --- Constant parameters -----------------------------------------------------

okabe_ito_colors <- palette.colors(palette = "Okabe-Ito")
okabe_ito_colors <- okabe_ito_colors[-(1:3)]  # drop the first three colors

# --- Row-count helper functions ----------------------------------------------

# Count rows per dataset (and label), deriving label from dataset if missing.
count_rows_per_dataset <- function(long_df) {
  if (!"label" %in% colnames(long_df)) {
    long_df <- long_df %>%
      mutate(
        label = dataset %>%
          gsub("_[0-9]+_", "_", .) %>%     # replace replicate number in the middle
          gsub("_[0-9]+$", "", .)          # remove trailing replicate number
      )
  }
  
  long_df %>%
    group_by(dataset, label) %>%
    summarise(Num_Rows = n(), .groups = "drop") %>%
    rename(Dataset = dataset)
}

# Average row counts across replicates and expand replicate values into columns.
summarise_replicates <- function(df, count_columns = c("Num_Rows")) {
  if (!"replicate" %in% names(df)) {
    df <- df %>%
      mutate(
        replicate = str_extract(Dataset, "_[0-9]+$") %>%
          str_remove_all("_")
      )
  }
  
  df <- df %>% mutate(Num_Rows = as.numeric(Num_Rows))
  
  df_summary <- df %>%
    group_by(label) %>%
    summarise(
      Num_Rows_mean = mean(Num_Rows, na.rm = TRUE),
      Num_Rows_sd   = sd(Num_Rows, na.rm = TRUE),
      Num_Replicates = n_distinct(replicate),
      replicate_values = list(Num_Rows),
      .groups = "drop"
    )
  
  # Expand replicate_values into Rep1..RepN columns, padding with NA.
  max_reps <- max(sapply(df_summary$replicate_values, length))
  
  replicate_cols <- t(sapply(df_summary$replicate_values, function(x) {
    length(x) <- max_reps
    x
  }))
  
  replicate_cols <- as_tibble(replicate_cols)
  colnames(replicate_cols) <- paste0("Rep", seq_len(ncol(replicate_cols)))
  
  df_summary <- bind_cols(df_summary %>% select(-replicate_values), replicate_cols)
  
  return(df_summary)
}