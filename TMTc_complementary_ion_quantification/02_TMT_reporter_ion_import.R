# ===========================================================================
# Script:  02_TMT_reporter_ion_import.R
# Purpose: Import and clean the standard TMT reporter-ion protein tables: 
#          rename channels, keep relevant columns, remove
#          contaminants, annotate organism, and compute means/ratios.
# Inputs:  data/baseline/*.txt, data/standard/*.txt ;
#          channel_names_groups (from script 00).
# Outputs: protein_data_frames_baseline/_standard, long_prot_baseline,
#          long_prot_standard and the function
#          filter_and_log2_transform_long_prot (used by script 03).
# Depends: 00_setup_packages_and_experimental_design.R
# ===========================================================================

# --- Import function ---
import_data <- function(files, existing_data = list()) {
  clean_name <- function(file) tools::file_path_sans_ext(basename(file))

  file_names <- sapply(files, clean_name, USE.NAMES = FALSE)
  new_files  <- files[!file_names %in% names(existing_data)]
  new_names  <- file_names[!file_names %in% names(existing_data)]

  for (i in seq_along(new_files)) {
    message("Reading and storing: ", new_names[i])
    data <- tryCatch(
      read_tsv(new_files[i], show_col_types = FALSE),
      error = function(e) {
        warning("Could not read file: ", new_files[i], "\nReason: ", e$message)
        NULL
      }
    )
    if (!is.null(data)) existing_data[[new_names[i]]] <- data
  }
  return(existing_data)
}

# --- Define data directories (baseline and standard only) ---
data_dirs <- c(
  baseline = file.path("data", "baseline"),
  standard = file.path("data", "standard")
)

# --- Import all .txt data (per folder) ---
data_list <- imap(data_dirs, function(path, name) {
  files <- list.files(path = path, pattern = "\\.txt$", full.names = TRUE)
  if (length(files) > 0) {
    message("--- Importing Group: ", name, " ---")
    import_data(files)
  } else { list() }
})

# --- Split into protein lists and push to global environment ---
protein_data_frames <- map(data_list, ~ .x[grepl("Proteins", names(.x))])

list2env(setNames(protein_data_frames, paste0("protein_data_frames_", names(protein_data_frames))), envir = .GlobalEnv)

# --- Name vectors ---
txt_groups    <- names(data_dirs)
protein_names <- paste0("protein_data_frames_", txt_groups)
all_names     <- protein_names

# --- Fetch all objects ---
all_lists <- mget(all_names, envir = .GlobalEnv)

# --- Rename the relevant abundance / resolution / S-N columns ---
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

  # 3. Rename "Annotated TMTpro-<channel>" columns -> "Pure-<channel>" (purity PSM data frames)
  annotated_psm_cols <- grep("^Annotated TMTpro-\\d{3}[NC]?$", names(df), value = TRUE)
  for (col in annotated_psm_cols) {
    match <- stringr::str_match(col, "Annotated TMTpro-(\\d{3}[NC]?)")[, 2]
    if (!is.na(match) && match %in% tmt_channels) {
      names(df)[names(df) == col] <- paste0("Pure-", match)
    }
  }

  # 4. Rename "Annotated protein-level TMTpro-<channel>" -> "Pure-<channel>" (purity protein data frames)
  annotated_protein_cols <- grep("^Annotated protein-level TMTpro-\\d{3}[NC]?$", names(df), value = TRUE)
  for (col in annotated_protein_cols) {
    match <- stringr::str_match(col, "Annotated protein-level TMTpro-(\\d{3}[NC]?)")[, 2]
    if (!is.na(match) && match %in% tmt_channels) {
      names(df)[names(df) == col] <- paste0("Pure-", match)
    }
  }

  # 5. Rename ratio columns: "Abundance Ratio: (1, mix2) / (1, mix1)" -> "mix2/mix1"
  ratio_cols <- grep("^Abundance Ratio: \\(\\d+, mix\\d+\\) / \\(\\d+, mix\\d+\\)$", names(df), value = TRUE)
  for (col in ratio_cols) {
    matches <- stringr::str_match(col, "^Abundance Ratio: \\(\\d+, (mix\\d+)\\) / \\(\\d+, (mix\\d+)\\)$")
    if (!is.na(matches[1, 2]) && !is.na(matches[1, 3])) {
      new_name <- paste0(matches[1, 2], "/", matches[1, 3])
      names(df)[names(df) == col] <- new_name
    }
  }

  # 6. Rename resolution columns: "RESGUI_126 Resolution" -> "res_126"
  res_cols <- grep("^RESGUI_\\d{3}[NC]? Resolution$", names(df), value = TRUE)
  for (col in res_cols) {
    matches <- stringr::str_match(col, "^RESGUI_(\\d{3}[NC]?) Resolution$")
    if (!is.na(matches[1, 2])) {
      new_name <- paste0("res_", matches[1, 2])
      names(df)[names(df) == col] <- new_name
    }
  }

  # 7. Rename "Ion Inject Time [ms]" to "IT[ms]" if present
  if ("Ion Inject Time [ms]" %in% names(df)) {
    names(df)[names(df) == "Ion Inject Time [ms]"] <- "IT[ms]"
  }

  # 8. Rename S/N columns: "Condition_SN_H100Y000" -> "SN_H100Y000"
  StoN_cols <- grep("^Condition_SN_", names(df), value = TRUE)

  for (col in StoN_cols) {
    matches <- str_match(col, "^Condition_SN_(.*)$")
    suffix <- matches[1, 2]  # the captured group, e.g. "H100Y000"

    if (!is.na(suffix)) {
      new_name <- paste0("SN_", suffix)  # e.g. "SN_H100Y000"
      names(df)[names(df) == col] <- new_name
    }
  }

  return(df)
}

all_lists <- lapply(all_lists, function(lst) {
  lapply(lst, rename_abundance_columns)
})

list2env(all_lists, envir = .GlobalEnv)

# --- Column specifications ---
# label abundance columns
pure_abund_cols <- c("Pure-126", "Pure-131N", "Pure-135N")

# resolution columns: defined from a protein data frame
resolution_cols <- grep("^res_", names(protein_data_frames_baseline[[1]]), value = TRUE)

# columns to keep for cleaner data frames
protein_purity_columns_to_keep <- c("Accession", "Description",
                    pure_abund_cols, "Annotated mean purity", "Annotated median purity")

keep_columns <- function(df, columns_to_keep) {
  df[, intersect(names(df), columns_to_keep), drop = FALSE]
}

# Apply to all protein dataframes
protein_lists <- mget(protein_names, envir = .GlobalEnv)
protein_lists <- lapply(protein_lists, function(lst) {
  lapply(lst, keep_columns, columns_to_keep = protein_purity_columns_to_keep)
})
list2env(protein_lists, envir = .GlobalEnv)

# --- Remove contaminants from "Accession" / "Master Protein Accessions" ---
remove_cont <- function(df) {
  # NOTE: "cont" is matched anywhere in the accession string (case-insensitive),
  # which could in theory also match non-contaminant entries containing "cont".
  if ("Master Protein Accessions" %in% names(df)) {
    df <- df[!grepl("cont", df$`Master Protein Accessions`, ignore.case = TRUE), ]
  } else if ("Accession" %in% names(df)) {
    df <- df[!grepl("cont", df$Accession, ignore.case = TRUE), ]
  }

  return(df)
}

# update all_lists
all_lists <- mget(all_names, envir = .GlobalEnv)

# Apply remove_cont() to every dataframe inside each list
all_lists <- lapply(all_lists, function(lst) {
  lapply(lst, remove_cont)
})

list2env(all_lists, envir = .GlobalEnv)

# --- Add organism column based on the description ---
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

# --- Replace zeros with NA (0 would otherwise be imputed) ---
replace_zeros_with_na <- function(df, target_cols) {
  cols <- intersect(names(df), target_cols)  # only process columns that exist in df
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

# --- Median normalization + ratio calculation ---
# Means/CV and ratios are only computed for rows where no value is NA, i.e.
# only when all replicates in one group per row are valid.

calculate_means_ratios_wide <- function(df, channel_info, exclude_channels = NULL, prefix = NULL) {

  # Helper: convert selected columns to numeric matrix safely
  get_numeric_matrix <- function(cols) {
    mat <- df[, cols, drop = FALSE]                  # keep as data.frame even if one column
    mat <- as.data.frame(lapply(mat, as.numeric))    # convert to numeric
    as.matrix(mat)                    # ensure matrix
  }

  # Add prefix if specified (such as "Pure-")
  if (!is.null(prefix)) {
    channel_info <- channel_info %>%
      mutate(label_prefixed = paste0(prefix, label))
    label_col <- "label_prefixed"
  } else {
    label_col <- "label"
  }

  # Get channel names per experiment group
  group_channels <- channel_info %>%
    filter(experiment %in% c("H000Y000", "H100Y005", "H100Y025")) %>%
    split(.$experiment) %>%
    lapply(function(x) x[[label_col]])

  # Exclude channels if requested
  if (!is.null(exclude_channels)) {
    group_channels <- lapply(group_channels, function(cols) setdiff(cols, exclude_channels))
  }

  # Compute numeric matrices for each experiment
  mat_list <- lapply(group_channels, function(ch) get_numeric_matrix(intersect(ch, names(df))))

  # Function to compute rowMeans and CV only if all values present
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

  # Compute stats for each group
  stats <- lapply(mat_list, compute_valid_stats)

  # Add to df
  for (grp in names(stats)) {
    df[[paste0("unnorm_", grp)]] <- stats[[grp]]$means
    df[[paste0("cv_", grp)]] <- stats[[grp]]$cv
  }

  # Compute ratios and mean expressions only where all relevant means are non-NA
  safe_div <- function(a, b) ifelse(!is.na(a) & !is.na(b), a / b, NA_real_)
  safe_mean <- function(a, b) ifelse(!is.na(a) & !is.na(b), (a + b) / 2, NA_real_)

  df$unnorm_H100Y025vsH100Y005 <- safe_div(df$unnorm_H100Y025, df$unnorm_H100Y005)

  df$unnorm_meanExp_H100Y025vsH100Y005 <- safe_mean(df$unnorm_H100Y025, df$unnorm_H100Y005)

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

# --- Pattern-based label mapping for source identification ---

label_patterns <- c(
  "DIA_mz1_50ng" = "DIA 50ng",
  "DIA_mz1_100ng" = "DIA 100ng",
  "DIA_mz1_120min_50ng" = "DIA120 50ng",
  "DIA_mz1_120min_100ng" = "DIA120 100ng",

  "DDA_mz0c5_ij8_50ng" = "DDA 50ng",
  "DDA_mz0c5_ij8_100ng" = "DDA 100ng"

)

# add a label column to each long-format dataframe based on label_patterns
add_labels_and_replicates <- function(df, dataset_name, label_patterns) {
  dt <- data.table::as.data.table(df)

  # Add dataset column
  dt[, dataset := dataset_name]

  # Clean dataset string (remove trailing _purity_tmt_quant if present)
  dataset_clean <- gsub("_purity_tmt_quant$", "", dataset_name)

  # Assign label: look for substring match
  label <- NA_character_
  patterns_ordered <- names(label_patterns)[order(nchar(names(label_patterns)), decreasing = TRUE)]
  for (pattern in patterns_ordered) {
    if (stringr::str_detect(dataset_clean, stringr::fixed(pattern))) {
      label <- label_patterns[[pattern]]
      break
    }
  }

  # Extract last number as replicate
  replicate <- stringr::str_extract(dataset_clean, "(\\d+)(?!.*\\d)")
  replicate <- as.integer(replicate)

  # Add columns
  dt[, label := label]
  dt[, replicate := replicate]

  return(dt)
}

rm(all_lists, protein_lists)

apply_labels_to_list <- function(obj_name) {
  lst <- get(obj_name, envir = .GlobalEnv)
  lst <- Map(
    add_labels_and_replicates,
    lst,
    names(lst),
    MoreArgs = list(label_patterns = label_patterns)
  )
  assign(obj_name, lst, envir = .GlobalEnv)
}

all_names_with_special <- c(
  all_names
)

invisible(lapply(all_names_with_special, apply_labels_to_list))

# --- Convert each protein data frame list into a long-format data frame ---
all_groups <- c(names(data_dirs))

for (grp in all_groups) {
  assign(paste0("long_prot_", grp),
         rbindlist(get(paste0("protein_data_frames_", grp)),
                   use.names = TRUE, fill = TRUE),
         envir = .GlobalEnv)
}

# --- Protein-level filter + log2 transform function (used in script 03) ---
# Filters long protein data frames on empty abundance rows, optionally removes
# rows with both species in "Description", and log2-transforms abundances/ratios.
filter_and_log2_transform_long_prot <- function(df,
                    datasets_to_keep = NULL,
                    datasets_to_exclude = NULL,
                    quan_values_filter = FALSE,
                    unambiguous_filter = TRUE,
                    log2_transform = TRUE) {
  df <- copy(df)

  # --- 1. Optional dataset filtering ---
  if (!is.null(datasets_to_keep)) {
    df <- df[df$dataset %in% datasets_to_keep, ]
  } else if (!is.null(datasets_to_exclude)) {
    df <- df[!df$dataset %in% datasets_to_exclude, ]
  }

  # --- 2. Define abundance columns (with and without "Pure-" prefix) ---
  base_abundance_columns <- c("126", "127N", "127C", "128N", "128C", "129N",
                    "129C", "130N", "130C", "131N", "131C", "132N",
                    "132C", "133N", "133C", "134N", "134C", "135N")
  prefixed_abundance_columns <- paste0("Pure-", base_abundance_columns)

  # define columns that actually contain peptides/proteins for quant filtering
  used_abund_cols <- c("128N", "128C", "129N",
                    "130N", "130C", "131N",
                    "132N","132C", "133N",
                    "134N", "134C", "135N")
  used_prefixed_abund_cols <- paste0("Pure-", used_abund_cols)

  used_abund_cols_present <- intersect(c(used_abund_cols, used_prefixed_abund_cols),
                    colnames(df))

  # detect which version exists in df
  abundance_cols_present <- intersect(c(base_abundance_columns, prefixed_abundance_columns),
                    colnames(df))

  # --- 3. Optional unambiguous filter ---
  if (unambiguous_filter && "Description" %in% names(df)) {
    has_both_species <- grepl("OS=Homo sapiens", df$Description, ignore.case = TRUE) &
      grepl("OS=Saccharomyces cerevisiae", df$Description, ignore.case = TRUE)
    df <- df[!has_both_species, ]
  }

  # --- 4. Optional quantitative abundance filter ---
  if (quan_values_filter) {
    if (length(used_abund_cols_present) == 0) {
      warning("No abundance columns found for quantitative filtering.")
    } else {
      abundance_data <- df[, used_abund_cols_present, with = FALSE]

      # Drop rows where all abundance columns are empty, 0, or NA
      is_empty_row <- apply(abundance_data, 1, function(row) {
        all(is.na(row) | row == 0 | row == "")
      })

      df <- df[!is_empty_row, ]
    }
  }

  # --- 5. Optional log2 transform ---
  if (log2_transform) {
    # Ratio and meanExp columns (both normalized and unnormalized)
    ratio_cols_present <- grep("^(unnorm_|H100Y|meanExp|unnorm_meanExp)", names(df), value = TRUE)

    # Combine abundance + ratio columns
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

