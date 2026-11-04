# ===========================================================================
# Script:  02a_functions_for_SN_DIANN.R
# Purpose: Import/parsing helper functions for Spectronaut (SN) and DIA-NN
#          result tables.
# Inputs:  None (function definitions only).
# Outputs: import_data_DIANN_parquet, import_data_SN_csv, add_organism_DIANN,
#          add_organism_SN.
# Depends: Run after 01; functions are used by 02b and 03a.
# ===========================================================================

# import the parquet files from DIANN
import_data_DIANN_parquet <- function(files, existing_data = list()) {
  existing_names <- names(existing_data)

  # Keep the full base name (no extension)
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
      arrow::read_parquet(file),
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

# import the csv files from SN
import_data_SN_csv <- function(files, existing_data = list()) {

  existing_names <- names(existing_data)
  clean_name <- function(file) tools::file_path_sans_ext(basename(file))

  file_names <- sapply(files, clean_name, USE.NAMES = FALSE)
  new_files <- files[!file_names %in% existing_names]
  new_names <- file_names[!file_names %in% existing_names]

  for (i in seq_along(new_files)) {
    file <- new_files[i]
    full_name <- new_names[i]

    message("Reading and storing: ", full_name)

    data <- tryCatch(
      data.table::fread(file, data.table = FALSE),
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


# ---- Functions to add the organism in a new column -------------------------


add_organism_DIANN <- function(df) {
  df %>%
    mutate(
      Organism = case_when(
        # Yeast if ANY of the available columns matches the yeast pattern
        if_any(any_of(c("Protein.Group", "Protein.Names")),
               ~ str_detect(.x, "^yeast_")) ~ "Yeast",

        # Human if ANY of the available columns matches the human pattern
        if_any(any_of(c("Protein.Group", "Protein.Names")),
               ~ str_detect(.x, "_HUMAN")) ~ "Human",

        # CRAP if Protein.Group contains "cont"
        if_any(any_of(c("Protein.Group")),
               ~ str_detect(.x, "cont")) ~ "cont",

        # Otherwise: NA
        TRUE ~ NA_character_
      )
    )
}


add_organism_SN <- function(df) {  
  df %>%  
    mutate(  
      Organism = case_when(  
        str_detect(PG.FastaFiles, "Saccharomyces_cerevisiae_") ~ "Yeast",  
        str_detect(PG.FastaFiles, "uniprot_reference_human")   ~ "Human",  
        str_detect(PG.FastaFiles, "TT_yeast_")                 ~ "Yeast",  
        str_detect(PG.FastaFiles, "TT_crap_")                  ~ "cont",  
        TRUE ~ NA_character_  
      )  
    )  
}

