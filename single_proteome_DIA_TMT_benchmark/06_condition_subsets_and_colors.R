# =============================================================================
# Script 06: Condition subsets and colour scheme
# Description: Defines the method/condition label subsets used to select and
#              order plot facets, and builds a viridis-based colour map for the
#              condition labels.
# Depends on: 03_filter_transform_and_counts.R, 05_concatenate_protein_counts.R
# =============================================================================

# --- Method / condition subsets ----------------------------------------------

# Entries must match the label_patterns() output exactly.

methods_54ng <- c("DDA 0.5 m/z \nIT20, 53.75ng",
                  "DIA 1.0 m/z \nIT8, 53.75ng",
                  "DIA 1.0 m/z \nHCD29, 53.75ng",
                  "DIA 1.0 m/z \nALL, 53.75ng")

# --- Colour scheme for condition labels --------------------------------------

ordered_df <- reorder_labels(long_prot_counts$long_prot_standard_StoN_grouping_20)

# Unique labels from the exemplary data frame.
all_labels <- levels(ordered_df$label)

# Consistent viridis colours for the labels.
label_colors <- setNames(viridis(length(all_labels), option = "D"), all_labels)

# Fallback colours for any labels missing from label_colors.
missing_labels <- setdiff(all_labels, names(label_colors))

if (length(missing_labels) > 0) {
  fallback_colors <- RColorBrewer::brewer.pal(12, "Set3")
  fallback_colors <- rep(fallback_colors, length.out = length(missing_labels))
  names(fallback_colors) <- missing_labels

  label_colors <- c(label_colors, fallback_colors)

  warning("The following labels were missing from 'label_colors' and were assigned fallback colors: ",
          paste(missing_labels, collapse = ", "))
}
