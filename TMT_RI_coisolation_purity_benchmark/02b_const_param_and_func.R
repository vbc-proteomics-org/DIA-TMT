# ===========================================================================
# Script:  02b_const_param_and_func.R
# Purpose: Import all result tables and define the constants and functions used
#          to clean, reshape, label and filter the PSM and protein data.
# Inputs:  Raw result files under data/; helper functions from 02a.
# Outputs: long_psm_* / long_prot_* data frames, all_groups, okabe_ito_colors
#          and the data-processing and filtering functions used from 03a
#          onward.
# Depends: Run after 02a.
# ===========================================================================

# ---- Read result files from the data subdirectory and clean the data ----

# ---- Import function -------------------------------------------------------
import_data <- function(files, existing_data = list()) {  
  clean_name <- function(file) tools::file_path_sans_ext(basename(file))  

  file_names <- sapply(files, clean_name, USE.NAMES = FALSE)  
  new_files  <- files[!file_names %in% names(existing_data)]  
  new_names  <- file_names[!file_names %in% names(existing_data)]  

  total_new <- length(new_files)
  message("Total files to read: ", total_new, " (", length(files) - total_new, " already loaded)")

  for (i in seq_along(new_files)) {  
    message("Reading and storing (", i, "/", total_new, "): ", new_names[i],
            "  [", total_new - i, " remaining]")  
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

# ---- Define all data directories -------------------------------------------
data_dirs <- c(
  baseline        = file.path("data", "baseline"),
  baseline_70perc = file.path("data", "baseline_70perc"),
  standard_StoN_grouping_20 = file.path("data", "standard_StoN_grouping_20"),
  StoN_grouping_3 = file.path("data", "StoN_grouping_3"),
  StoN_grouping_6 = file.path("data", "StoN_grouping_6"),
  StoN_grouping_9 = file.path("data", "StoN_grouping_9"),
  StoN_grouping_12 = file.path("data", "StoN_grouping_12"),
  StoN_grouping_15 = file.path("data", "StoN_grouping_15"),
  StoN_grouping_20 = file.path("data", "StoN_grouping_20"),
  StoN_grouping_30 = file.path("data", "StoN_grouping_30"),
  StoN_grouping_40 = file.path("data", "StoN_grouping_40"),
  StoN_grouping_80 = file.path("data", "StoN_grouping_80")

)


# ---- Import all standard .txt data (per folder) ----------------------------
data_list <- imap(data_dirs, function(path, name) {
  files <- list.files(path = path, pattern = "\\.txt$", full.names = TRUE)
  if (length(files) > 0) {
    message("--- Importing Group: ", name, " ---")
    import_data(files)
  } else { list() }
})

# ---- Import SN and DIANN with dedicated functions --------------------------

data_list[["SNaut_baseline"]]    <- import_data_SN_csv(list.files(file.path("data", "SNaut_baseline"), 
                                                                  pattern = "\\.csv$", full.names = TRUE))
data_list[["DIANN_baseline"]] <- import_data_DIANN_parquet(list.files(file.path("data", "DIANN_baseline"), 
                                                                      pattern = "\\.parquet$", full.names = TRUE))
data_list[["SNaut_StoN_grouping_20"]]    <- import_data_SN_csv(list.files(file.path("data", "SNaut_StoN_grouping_20"), 
                                                                          pattern = "\\.csv$", full.names = TRUE))
data_list[["DIANN_StoN_grouping_20"]] <- import_data_DIANN_parquet(list.files(file.path("data", "DIANN_StoN_grouping_20"), 
                                                                              pattern = "\\.parquet$", full.names = TRUE))

# ---- Helper to pick naming patterns ----------------------------------------
psm_pat     <- function(n) if (n %in% c("SNaut_baseline", "DIANN_baseline",  
                                        "SNaut_StoN_grouping_20", "DIANN_StoN_grouping_20")) "precursor" else "PSM"
protein_pat <- function(n) if (n %in% c("SNaut_baseline", "DIANN_baseline",
                                        "DIANN_StoN_grouping_20", "SNaut_StoN_grouping_20")) "PG"        else "Proteins"

# ---- Split into PSM and protein lists --------------------------------------
psm_data_frames     <- imap(data_list, ~ .x[grepl(psm_pat(.y),     names(.x))])
protein_data_frames <- imap(data_list, ~ .x[grepl(protein_pat(.y), names(.x))])

# ---- Push to global environment --------------------------------------------
list2env(setNames(psm_data_frames,     paste0("psm_data_frames_",     names(psm_data_frames))),     envir = .GlobalEnv)
list2env(setNames(protein_data_frames, paste0("protein_data_frames_", names(protein_data_frames))), envir = .GlobalEnv)

# ---- Name vectors (excluding SN and DIANN) ---------------------------------
txt_groups    <- names(data_dirs)  # only the .txt-based groups
psm_names     <- paste0("psm_data_frames_",     txt_groups)
protein_names <- paste0("protein_data_frames_", txt_groups)
all_names     <- c(psm_names, protein_names)

# ---- Fetch all objects (excludes SN and DIANN) -----------------------------
all_lists <- mget(all_names, envir = .GlobalEnv)


# clean up the data by renaming the relevant columns and removing irrelevant columns
# rename the abundance columns to have easier names 
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

  # 3. Rename "Annotated TMTpro-<channel>" columns -> "Pure-<channel>" in purity psm data frames!
  annotated_psm_cols <- grep("^Annotated TMTpro-\\d{3}[NC]?$", names(df), value = TRUE)
  for (col in annotated_psm_cols) {
    match <- stringr::str_match(col, "Annotated TMTpro-(\\d{3}[NC]?)")[, 2]
    if (!is.na(match) && match %in% tmt_channels) {
      names(df)[names(df) == col] <- paste0("Pure-", match)
    }
  }

  # 4. Rename "Annotated protein-level TMTpro-<channel>" -> "Pure-<channel>" purity Protein data frames!
  annotated_protein_cols <- grep("^Annotated protein-level TMTpro-\\d{3}[NC]?$", names(df), value = TRUE)
  for (col in annotated_protein_cols) {
    match <- stringr::str_match(col, "Annotated protein-level TMTpro-(\\d{3}[NC]?)")[, 2]
    if (!is.na(match) && match %in% tmt_channels) {
      names(df)[names(df) == col] <- paste0("Pure-", match)
    }
  }

  # 5. Rename ratio columns in format: "Abundance Ratio: (1, mix2) / (1, mix1)" -> "mix2/mix1"
  ratio_cols <- grep("^Abundance Ratio: \\(\\d+, mix\\d+\\) / \\(\\d+, mix\\d+\\)$", names(df), value = TRUE)
  for (col in ratio_cols) {
    matches <- stringr::str_match(col, "^Abundance Ratio: \\(\\d+, (mix\\d+)\\) / \\(\\d+, (mix\\d+)\\)$")
    if (!is.na(matches[1, 2]) && !is.na(matches[1, 3])) {
      new_name <- paste0(matches[1, 2], "/", matches[1, 3])
      names(df)[names(df) == col] <- new_name
    }
  }

  # 6. Rename resolution columns: "RESGUI_126 Resolution" to "resolution_126"
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

  # 8. Rename S/N columns: "Condition_SN_H100Y000" to "SN_H100Y000"
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


# specifiy the columns that contain the label abundances
abundance_cols <- c("126", "127N", "127C", "128N", "128C", "129N", 
                       "129C", "130N", "130C", "131N", "131C", "132N",
                       "132C", "133N", "133C", "134N", "134C", "135N")

pure_abund_cols <- c("Pure-126", "Pure-127N", "Pure-127C", "Pure-128N", "Pure-128C", "Pure-129N", 
                    "Pure-129C", "Pure-130N", "Pure-130C", "Pure-131N", "Pure-131C", "Pure-132N",
                    "Pure-132C", "Pure-133N", "Pure-133C", "Pure-134N", "Pure-134C", "Pure-135N")

# specifiy the columns that contain the label_resolution
resolution_cols <- grep("^res_", names(psm_data_frames_baseline[[1]]), value = TRUE)

StoN_cols <- c("SN_H000Y000","SN_H100Y000", "SN_H100Y010", "SN_H100Y020", "SN_H100Y100")


# filter for columns that are needed to make the data.frames cleaner
protein_columns_to_keep <- c("Accession", "Description",
                             abundance_cols#, ratio_cols
                             )

peptide_columns_to_keep <- c("Master Protein Accessions", "Modifications", 
                             abundance_cols#, ratio_cols
                             )

psm_columns_to_keep <- c("Annotated Sequence", "First Scan", "Master Protein Accessions", "Ion Inject Time ms", 
                         "PSM Ambiguity", "Average Reporter S/N", "Normalized CHIMERYS Coefficient",
                         "Master Protein Descriptions", abundance_cols, "m/z [Da]", "IT[ms]")

psm_purity_columns_to_keep <- c("Annotated Sequence", "First Scan", "Master Protein Accessions", "Ion Inject Time ms", 
                         "PSM Ambiguity", "Average Reporter S/N", "Normalized CHIMERYS Coefficient", "Co-Isolation Purity",
                         "Master Protein Descriptions", 
                         resolution_cols, pure_abund_cols, StoN_cols, "Annotated mean purity", "Annotated median purity")

protein_purity_columns_to_keep <- c("Accession", "Description",
                             pure_abund_cols, "Annotated mean purity", "Annotated median purity")

psm_purity_columns_to_keep_SN <- c("Co-Isolation Purity","PEP.StrippedSequence","PG.FastaFiles", "PG.ProteinDescriptions", "PG.Organisms", #abundance_cols, 
                                   resolution_cols, pure_abund_cols, StoN_cols, "Annotated mean purity", "Annotated median purity")

protein_purity_columns_to_keep_SN <- c("PG.FastaFiles", "PG.GroupLabel",
                                    pure_abund_cols, "Annotated mean purity", "Annotated median purity")

psm_purity_columns_to_keep_DIANN <- c("Protein.Group", "Stripped.Sequence", "Protein.Names", "Protein.Ids", "Co-Isolation Purity",
                                      "Master Protein Descriptions", 
                                      resolution_cols, pure_abund_cols, StoN_cols, "Annotated mean purity", "Annotated median purity")

protein_purity_columns_to_keep_DIANN <- c("Protein.Group", "Protein.Names",
                                    pure_abund_cols, "Annotated mean purity", "Annotated median purity")

keep_columns <- function(df, columns_to_keep) {
  df[, intersect(names(df), columns_to_keep), drop = FALSE]
}

# apply to all psm dataframes
psm_lists <- mget(psm_names, envir = .GlobalEnv)
psm_lists <- lapply(psm_lists, function(lst) {
  lapply(lst, keep_columns, columns_to_keep = psm_purity_columns_to_keep)
})
list2env(psm_lists, envir = .GlobalEnv)

# Apply to all Protein dataframes
protein_lists <- mget(protein_names, envir = .GlobalEnv)
protein_lists <- lapply(protein_lists, function(lst) {
  lapply(lst, keep_columns, columns_to_keep = protein_purity_columns_to_keep)
})
list2env(protein_lists, envir = .GlobalEnv)

# remove contaminants from "Accession" and "Master Protein Accessions"
remove_cont <- function(df) {  
  if ("Master Protein Accessions" %in% names(df)) {  
    df <- df[!grepl("cont", df$`Master Protein Accessions`, ignore.case = TRUE), ]  
  } else if ("Accession" %in% names(df)) {  
    df <- df[!grepl("cont", df$Accession, ignore.case = TRUE), ]  
  } else if ("PG.FastaFiles" %in% names(df)) {  
    df <- df[!grepl("TT_crap_", df$PG.FastaFiles, ignore.case = TRUE), ]  
  } else if ("Protein.Group" %in% names(df)) {  
    df <- df[!grepl("cont", df$Protein.Group, ignore.case = TRUE), ]  
  }  

  return(df)  
}


# update all_lists 
all_lists <- mget(all_names, envir = .GlobalEnv)

# Apply remove_cont() to every dataframe inside each list
all_lists <- lapply(all_lists, function(lst) {
  lapply(lst, remove_cont)
})

# Save results back into the original variable names
list2env(all_lists, envir = .GlobalEnv)

# function to add the organism in a new column

add_organism <- function(df) {
  # Determine which column to use
  desc_col <- NULL

  if ("Description" %in% names(df)) {
    desc_col <- "Description"
  } else if ("Master Protein Descriptions" %in% names(df)) {
    desc_col <- "Master Protein Descriptions"
  }

  # If a relevant column exists, assign organism based on it
  if (!is.null(desc_col)) {
    df <- df %>%
      mutate(Organism = case_when(
        str_detect(.data[[desc_col]], "OS=Homo sapiens") ~ "Human",
        str_detect(.data[[desc_col]], "OS=Saccharomyces") ~ "Yeast",
        TRUE ~ "Other"
      ))
  } else {
    # If neither column exists, just add NA column
    df$Organism <- NA
  }

  return(df)
}


all_lists <- mget(all_names, envir = .GlobalEnv)

# Apply add_organism to every dataframe inside each list
all_lists <- lapply(all_lists, function(lst) {
  lapply(lst, add_organism)
})

# Save results back into the original variable names
list2env(all_lists, envir = .GlobalEnv)

# Replace 0 with NA in all data frames 
# Function: replace zeros with NAs in selected columns
replace_zeros_with_na <- function(df, target_cols) {
  cols <- intersect(names(df), target_cols)  # only process columns that exist in df
  df[cols] <- lapply(df[cols], function(x) {
    x[x == 0] <- NA
    x
  })
  df
}

# update all_lists 
all_lists <- mget(all_names, envir = .GlobalEnv)

# Apply replace_zeros_with_na to each dataframe in each list
all_lists <- lapply(all_lists, function(lst) {
  lapply(lst, replace_zeros_with_na, target_cols = c(pure_abund_cols, resolution_cols))
})

# Assign back to environment
list2env(all_lists, envir = .GlobalEnv)

# ---- Calculations on the data: optional median normalization + ratio calculation ----

# cv and ratio calculation only for rows where no value is NA
# so calculation only when all replicates in one group per row are valid 

calculate_means_ratios_wide <- function(df, channel_info, exclude_channels = NULL, prefix = NULL) {

  get_numeric_matrix <- function(cols) {
    mat <- df[, cols, drop = FALSE]
    mat <- as.data.frame(lapply(mat, as.numeric))
    as.matrix(mat)
  }

  if (!is.null(prefix)) {
    channel_info <- channel_info %>%
      mutate(label_prefixed = paste0(prefix, label))
    label_col <- "label_prefixed"
  } else {
    label_col <- "label"
  }

  group_channels <- channel_info %>%
    filter(experiment %in% c("H100Y000", "H100Y010", "H100Y020", "H100Y100")) %>%
    split(.$experiment) %>%
    lapply(function(x) x[[label_col]])

  if (!is.null(exclude_channels)) {
    group_channels <- lapply(group_channels, function(cols) setdiff(cols, exclude_channels))
  }

  mat_list <- lapply(group_channels, function(ch) get_numeric_matrix(intersect(ch, names(df))))

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

  safe_div <- function(a, b) ifelse(!is.na(a) & !is.na(b), a / b, NA_real_)
  safe_mean <- function(a, b) ifelse(!is.na(a) & !is.na(b), (a + b) / 2, NA_real_)

  df$unnorm_H100Y000vsH100Y010 <- safe_div(df$unnorm_H100Y000, df$unnorm_H100Y010)
  df$unnorm_H100Y000vsH100Y020 <- safe_div(df$unnorm_H100Y000, df$unnorm_H100Y020)
  df$unnorm_H100Y000vsH100Y100 <- safe_div(df$unnorm_H100Y000, df$unnorm_H100Y100)
  df$unnorm_H100Y100vsH100Y010 <- safe_div(df$unnorm_H100Y100, df$unnorm_H100Y010)
  df$unnorm_H100Y100vsH100Y020 <- safe_div(df$unnorm_H100Y100, df$unnorm_H100Y020)
  df$unnorm_H100Y020vsH100Y010 <- safe_div(df$unnorm_H100Y020, df$unnorm_H100Y010)

  df$unnorm_meanExp_H100Y000vsH100Y010 <- safe_mean(df$unnorm_H100Y000, df$unnorm_H100Y010)
  df$unnorm_meanExp_H100Y000vsH100Y020 <- safe_mean(df$unnorm_H100Y000, df$unnorm_H100Y020)
  df$unnorm_meanExp_H100Y000vsH100Y100 <- safe_mean(df$unnorm_H100Y000, df$unnorm_H100Y100)
  df$unnorm_meanExp_H100Y100vsH100Y010 <- safe_mean(df$unnorm_H100Y100, df$unnorm_H100Y010)
  df$unnorm_meanExp_H100Y100vsH100Y020 <- safe_mean(df$unnorm_H100Y100, df$unnorm_H100Y020)
  df$unnorm_meanExp_H100Y020vsH100Y010 <- safe_mean(df$unnorm_H100Y020, df$unnorm_H100Y010)

  return(df)
}


# update all_lists 
all_lists <- mget(all_names, envir = .GlobalEnv)

total_outer <- length(all_lists)
outer_names <- names(all_lists)

all_lists <- lapply(seq_along(all_lists), function(i) {

  outer_name <- outer_names[i]
  inner_list <- all_lists[[i]]
  inner_names <- names(inner_list)
  total_inner <- length(inner_list)

  message("=== Outer group (", i, "/", total_outer, "): ", outer_name,
          "  [", total_outer - i, " groups remaining] ===")

  result <- lapply(seq_along(inner_list), function(j) {
    message("  -> Processing (", j, "/", total_inner, "): ", inner_names[j],
            "  [", total_inner - j, " remaining in this group]")

    calculate_means_ratios_wide(
      df = inner_list[[j]],
      channel_info = channel_names_groups,
      prefix = "Pure-",
      exclude_channels = "Pure-135N"
    )
  })
  names(result) <- inner_names
  result
})
names(all_lists) <- outer_names

# Assign modified dataframes back to the environment
list2env(all_lists, envir = .GlobalEnv)

# ---- Function for the later source identification of each datapoint. -------

# Pattern-based label mapping for the input datasets

label_patterns <- c(
  "DIA_mz0c5_ij8_3ng" = "DIA 0.5 m/z \nIT8, 3ng",
  "DIA_mz0c5_ij8_53c75ng" = "DIA 0.5 m/z \nIT8, 53.75ng",
  "DIA_mz0c5_ij8_100ng" = "DIA 0.5 m/z \nIT8, 100ng",

  "DIA_mz0c5_ij20_3ng" = "DIA 0.5 m/z \nIT20, 3ng",
  "DIA_mz0c5_ij20_53c75ng" = "DIA 0.5 m/z \nIT20, 53.75ng",
  "DIA_mz0c5_ij20_100ng" = "DIA 0.5 m/z \nIT20, 100ng",

  "DIA_mz1_ij8_3ng" = "DIA 1.0 m/z \nIT8, 3ng",
  "DIA_mz1_ij8_53c75ng" = "DIA 1.0 m/z \nIT8, 53.75ng",
  "DIA_mz1_ij8_100ng" = "DIA 1.0 m/z \nIT8, 100ng",

  "DIA_mz1_ij20_3ng" = "DIA 1.0 m/z \nIT20, 3ng",
  "DIA_mz1_ij20_53c75ng" = "DIA 1.0 m/z \nIT20, 53.75ng",
  "DIA_mz1_ij20_100ng" = "DIA 1.0 m/z \nIT20, 100ng",

  "DIA_mz1_ij40_3ng_" = "DIA 1.0 m/z \nIT40, 3ng",
  "DIA_mz1_ij40_53c75ng_" = "DIA 1.0 m/z \nIT40, 53.75ng",
  "DIA_mz1_ij40_100ng_" = "DIA 1.0 m/z \nIT40, 100ng",

  "DIA_mz1_ij40_5ng_" = "DIA 1.0 m/z \nIT40, 5ng",
  "DIA_mz1_ij40_10ng_" = "DIA 1.0 m/z \nIT40, 10ng",
  "DIA_mz1_ij40_20ng_" = "DIA 1.0 m/z \nIT40, 20ng",
  "DIA_mz1_ij40_25ng_" = "DIA 1.0 m/z \nIT40, 25ng",
  "DIA_mz1_ij40_30ng_" = "DIA 1.0 m/z \nIT40, 30ng",
  "DIA_mz1_ij40_40ng_" = "DIA 1.0 m/z \nIT40, 40ng",

  "DIA_mz2_ij8_3ng" = "DIA 2.0 m/z \nIT8, 3ng",
  "DIA_mz2_ij8_53c75ng" = "DIA 2.0 m/z \nIT8, 53.75ng",
  "DIA_mz2_ij8_100ng" = "DIA 2.0 m/z \nIT8, 100ng",

  "DDA_mz0c5_ij8_3ng" = "DDA 0.5 m/z \nIT8, 3ng",
  "DDA_mz0c5_ij8_53c75ng" = "DDA 0.5 m/z \nIT8, 53.75ng",
  "DDA_mz0c5_ij8_100ng" = "DDA 0.5 m/z \nIT8, 100ng",

  "DDA_mz0c5_ij20_3ng" = "DDA 0.5 m/z \nIT20, 3ng",
  "DDA_mz0c5_ij20_53c75ng" = "DDA 0.5 m/z \nIT20, 53.75ng",
  "DDA_mz0c5_ij20_100ng" = "DDA 0.5 m/z \nIT20, 100ng"
)

# function to add a label column to each long formated dataframe based on the above defined label_pattern

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

# ---- Registry of special (non-txt) groups ----------------------------------
special_groups <- list(  
  SNaut_baseline = list(psm_cols = psm_purity_columns_to_keep_SN,    prot_cols = protein_purity_columns_to_keep_SN,    organism_fn = add_organism_SN),  
  DIANN_baseline = list(psm_cols = psm_purity_columns_to_keep_DIANN, prot_cols = protein_purity_columns_to_keep_DIANN, organism_fn = add_organism_DIANN),
  SNaut_StoN_grouping_20 = list(psm_cols = psm_purity_columns_to_keep_SN,    prot_cols = protein_purity_columns_to_keep_SN,    organism_fn = add_organism_SN),  
  DIANN_StoN_grouping_20 = list(psm_cols = psm_purity_columns_to_keep_DIANN, prot_cols = protein_purity_columns_to_keep_DIANN, organism_fn = add_organism_DIANN)
)  

# ---- Apply all processing steps to special groups via loop -----------------
for (grp in names(special_groups)) {  
  cfg      <- special_groups[[grp]]  
  psm_obj  <- paste0("psm_data_frames_",     grp)  
  prot_obj <- paste0("protein_data_frames_", grp)  

  # rename  
  assign(psm_obj,  lapply(get(psm_obj),  rename_abundance_columns), envir = .GlobalEnv)  
  assign(prot_obj, lapply(get(prot_obj), rename_abundance_columns), envir = .GlobalEnv)  

  # keep columns  
  assign(psm_obj,  lapply(get(psm_obj),  keep_columns, columns_to_keep = cfg$psm_cols),  envir = .GlobalEnv)  
  assign(prot_obj, lapply(get(prot_obj), keep_columns, columns_to_keep = cfg$prot_cols), envir = .GlobalEnv)  

  # remove contaminants  
  assign(psm_obj,  lapply(get(psm_obj),  remove_cont), envir = .GlobalEnv)  
  assign(prot_obj, lapply(get(prot_obj), remove_cont), envir = .GlobalEnv)  

  # add organism  
  assign(psm_obj,  lapply(get(psm_obj),  cfg$organism_fn), envir = .GlobalEnv)  
  assign(prot_obj, lapply(get(prot_obj), cfg$organism_fn), envir = .GlobalEnv)  

  # replace zeros with NA  
  assign(psm_obj,  lapply(get(psm_obj),  replace_zeros_with_na, target_cols = c(pure_abund_cols, resolution_cols)), envir = .GlobalEnv)  
  assign(prot_obj, lapply(get(prot_obj), replace_zeros_with_na, target_cols = c(pure_abund_cols, resolution_cols)), envir = .GlobalEnv)  

  # calculate means/ratios  
  assign(psm_obj,  lapply(get(psm_obj),  calculate_means_ratios_wide, channel_info = channel_names_groups, prefix = "Pure-", exclude_channels = "Pure-135N"), envir = .GlobalEnv)  
  assign(prot_obj, lapply(get(prot_obj), calculate_means_ratios_wide, channel_info = channel_names_groups, prefix = "Pure-", exclude_channels = "Pure-135N"), envir = .GlobalEnv)  
}

rm(all_lists, psm_lists, protein_lists)

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
  all_names,  
  paste0("psm_data_frames_",     names(special_groups)),  
  paste0("protein_data_frames_", names(special_groups))  
)

invisible(lapply(all_names_with_special, apply_labels_to_list))

# Convert each data frame in the list to a long-format data frame


all_groups <- c(names(data_dirs), names(special_groups))

for (grp in all_groups) {
  assign(paste0("long_psm_",  grp),
         rbindlist(get(paste0("psm_data_frames_", grp)), #idcol = "dataset",
                   use.names = TRUE, fill = TRUE),
         envir = .GlobalEnv)
  assign(paste0("long_prot_", grp),
         rbindlist(get(paste0("protein_data_frames_", grp)), #idcol = "dataset",
                   use.names = TRUE, fill = TRUE),
         envir = .GlobalEnv)
}

# ---- Functions to filter datasets ------------------------------------------

# the following section contains functions to filter on psm, peptide and protein level. 
# also to log2 transform the data + filter specific data frames


#Filter PSM data frames based on "Master Protein Descriptions", "First Scan" and a numeric SN_filter based on "Average Reporter S/N"
#purity_filter -> # numeric threshold for purity
#chimerys_coeff_filter -> numeric filter for co-isolation
#unambiguous_filter = TRUE -> filter out all rows, that contain both species in "Master Protein Descriptions"
#SN_filter -> # numeric threshold for S/N; NULL = no filtering

filter_psm_long_dataframe <- function(df,
                                      purity_filter = NULL,
                                      chimerys_coeff_filter = NULL,
                                      unambiguous_filter = TRUE,
                                      SN_filter = NULL,
                                      log2_transform = FALSE,
                                      resolution_filter = NULL) {
  df <- copy(df)

# ---- Helper for conditional filtering --------------------------------------
  apply_filter <- function(df, col, min_value = NULL) {
    if (col %in% names(df)) {
      df <- df[!is.na(df[[col]])]
      if (!is.null(min_value)) {
        df <- df[df[[col]] >= min_value]
      }
    }
    df
  }

# ---- Apply numeric filters -------------------------------------------------
  df <- apply_filter(df, "Co-Isolation Purity", purity_filter)
  df <- apply_filter(df, "Average Reporter S/N", SN_filter)
  df <- apply_filter(df, "Normalized CHIMERYS Coefficient", chimerys_coeff_filter)

# ---- Remove mixed-species PSMs ---------------------------------------------
  if (unambiguous_filter && "Master Protein Descriptions" %in% names(df)) {
    both_species <- grepl("OS=Homo sapiens", df$`Master Protein Descriptions`, ignore.case = TRUE) &
      grepl("OS=Saccharomyces cerevisiae", df$`Master Protein Descriptions`, ignore.case = TRUE)
    df <- df[!both_species]
  }

# ---- Resolution filter -----------------------------------------------------
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

# ---- Log2 transform option -------------------------------------------------
  if (log2_transform) {
    # Base abundance columns and prefixed variants
    base_abundance_columns <- c("126", "127N", "127C", "128N", "128C", "129N", 
                                "129C", "130N", "130C", "131N", "131C", "132N",
                                "132C", "133N", "133C", "134N", "134C", "135N")
    prefixed_abundance_columns <- paste0("Pure-", base_abundance_columns)

    # Detect which are present
    abundance_cols_present <- intersect(c(base_abundance_columns, prefixed_abundance_columns),
                                        names(df))

    # Ratio columns (normalized and unnormalized)
    ratio_cols_present <- grep("^(unnorm_|H100Y|meanExp|unnorm_meanExp)", names(df), value = TRUE)

    # Combine all columns to log2-transform
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


#filter protein_data_frames based on empty abundance rows -> if TRUE = remove rows with all NA/0/"" in abundance cols
#also applicable on peptide_data_frames
#specific datasubset can be generated with a vector argument in datasets_to_keep
#unambiguous_filter = TRUE -> filter out all rows, that contain both species in "Protein.Descriptions
#also transforms the abundances & ratios into log2

filter_and_log2_transform_long_prot <- function(df,
                                                datasets_to_keep = NULL,
                                                datasets_to_exclude = NULL,
                                                quan_values_filter = FALSE,
                                                unambiguous_filter = TRUE,
                                                log2_transform = TRUE) {
  df <- copy(df)

# ---- 1. Optional dataset filtering -----------------------------------------
  if (!is.null(datasets_to_keep)) {
    df <- df[df$dataset %in% datasets_to_keep, ]
  } else if (!is.null(datasets_to_exclude)) {
    df <- df[!df$dataset %in% datasets_to_exclude, ]
  }

# ---- 2. Define abundance columns (with and without "Pure-" prefix) ---------
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


# ---- 3. Optional unambiguous filter ----------------------------------------
  if (unambiguous_filter && "Description" %in% names(df)) {
    has_both_species <- grepl("OS=Homo sapiens", df$Description, ignore.case = TRUE) &
      grepl("OS=Saccharomyces cerevisiae", df$Description, ignore.case = TRUE)
    df <- df[!has_both_species, ]
  }


# ---- 4. Optional quantitative abundance filter -----------------------------
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

# ---- 5. Optional log2 transform --------------------------------------------
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


filter_prot_abundance <- function(df,
                                  comparison_column = NULL,
                                  threshold = 0.2) {

# ---- Add Organism column based on Description ------------------------------
  df <- df %>%
    mutate(Organism = case_when(
      str_detect(`Description`, "OS=Homo sapiens") ~ "Homo sapiens",
      str_detect(`Description`, "OS=Saccharomyces") ~ "Saccharomyces",
      TRUE ~ "Other"
    ))

# ---- Identify unnorm_meanExp columns ---------------------------------------
  if (!is.null(comparison_column)) {
    cols_to_filter <- comparison_column
  } else {
    cols_to_filter <- names(df)[str_detect(names(df), "^unnorm_meanExp_")]
  }

# ---- Check for columns and valid threshold ---------------------------------
  if (length(cols_to_filter) == 0) {
    stop("No columns matching 'unnorm_meanExp_' found in the dataframe.")
  }

  if (threshold <= 0 | threshold > 1) {
    stop("Threshold must be between 0 and 1 (e.g., 0.2 for top 20%).")
  }

# ---- Filter per column, per label, and per organism ------------------------
  filtered_list <- lapply(cols_to_filter, function(colname) {
    message("Filtering by: ", colname)

    df %>%
      group_by(label, Organism) %>%
      group_modify(~ {
        cutoff <- quantile(.x[[colname]], probs = 1 - threshold, na.rm = TRUE)
        .x %>%
          filter(.data[[colname]] >= cutoff)
      }) %>%
      ungroup() %>%
      mutate(filtered_by = colname)
  })

# ---- Name list elements after columns --------------------------------------
  names(filtered_list) <- cols_to_filter

  return(filtered_list)
}


# ---- Establish constant parameters -----------------------------------------

okabe_ito_colors <- palette.colors(palette = "Okabe-Ito")
okabe_ito_colors <- okabe_ito_colors[-(1:3)]  # remove first 3 colors


# ---- Function to count rows for peptides, proteins and PSMs in the different conditions ----

# count all rows in the specified dataframes

count_rows_per_dataset <- function(long_df) {
  # If 'label' column doesn't exist, create it from 'dataset'
  if (!"label" %in% colnames(long_df)) {
    long_df <- long_df %>%
      mutate(
        label = dataset %>%
          gsub("_[0-9]+_", "_", .) %>%     # Replace replicate number in middle
          gsub("_[0-9]+$", "", .)          # Remove replicate number at end
      )
  }

  long_df %>%
    group_by(dataset, label) %>%
    summarise(Num_Rows = n(), .groups = "drop") %>%
    rename(Dataset = dataset)
}


# function to take the average of row counts based on identifiers e.g. mz2_2 (when there are n > 1, replicates of a condition)

summarise_replicates <- function(df, count_columns = c("Num_Rows")) {
  # Extract replicate from Dataset column if not already present
  if (!"replicate" %in% names(df)) {
    df <- df %>%
      mutate(
        replicate = str_extract(Dataset, "_[0-9]+$") %>%
          str_remove_all("_")
      )
  }

  # Make sure Num_Rows is numeric
  df <- df %>% mutate(Num_Rows = as.numeric(Num_Rows))

  # Group by label
  df_summary <- df %>%
    group_by(label) %>%
    summarise(
      # Summary statistics
      Num_Rows_mean = mean(Num_Rows, na.rm = TRUE),
      Num_Rows_sd   = sd(Num_Rows, na.rm = TRUE),
      Num_Replicates = n_distinct(replicate),
      # Replicate columns
      replicate_values = list(Num_Rows),
      .groups = "drop"
    )

  # Expand replicate_values into separate columns
  max_reps <- max(sapply(df_summary$replicate_values, length))

  replicate_cols <- t(sapply(df_summary$replicate_values, function(x) {
    # pad with NA if fewer replicates than max
    length(x) <- max_reps
    x
  }))

  replicate_cols <- as_tibble(replicate_cols)
  colnames(replicate_cols) <- paste0("Rep", seq_len(ncol(replicate_cols)))

  # Combine summary + replicate columns
  df_summary <- bind_cols(df_summary %>% select(-replicate_values), replicate_cols)

  return(df_summary)
}