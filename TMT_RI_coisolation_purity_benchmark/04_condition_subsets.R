# ===========================================================================
# Script:  04_condition_subsets.R
# Purpose: Define the named vectors of dataset labels (condition subsets) and
#          the label colour map used by all downstream plotting scripts.
# Inputs:  reorder_labels (03a), long_prot_counts (03d).
# Outputs: input_titration, isol_w_IT_53ng, test_conditions,
#          test_condition_DIA, label_colors and the other condition-subset
#          vectors.
# Depends: Run after 03a and 03d.
# ===========================================================================

# Subset the conditions for plotting
# the entries need to exactly match the label_patterns function output

input_titration <- c("DIA 1.0 m/z \nIT40, 3ng",
                     "DIA 1.0 m/z \nIT40, 5ng",
                     "DIA 1.0 m/z \nIT40, 10ng",
                     "DIA 1.0 m/z \nIT40, 20ng",
                     "DIA 1.0 m/z \nIT40, 25ng",
                     "DIA 1.0 m/z \nIT40, 30ng",
                     "DIA 1.0 m/z \nIT40, 40ng",
                     "DIA 1.0 m/z \nIT40, 53.75ng",
                     "DIA 1.0 m/z \nIT40, 100ng")




isol_w_IT_53ng <- c("DDA 0.5 m/z \nIT8, 53.75ng",
                    "DDA 0.5 m/z \nIT20, 53.75ng",
                    "DIA 0.5 m/z \nIT8, 53.75ng",
                    "DIA 1.0 m/z \nIT8, 53.75ng",
                    "DIA 1.0 m/z \nIT20, 53.75ng",
                    "DIA 1.0 m/z \nIT40, 53.75ng",
                    "DIA 2.0 m/z \nIT8, 53.75ng"
                 )


test_conditions <- c("DDA 0.5 m/z \nIT8, 53.75ng",
                    "DIA 1.0 m/z \nIT8, 53.75ng")

test_condition_DIA <- c("DIA 1.0 m/z \nIT8, 53.75ng")


# Define the colour scheme for the conditions

# one exemplary dataframe with the ordered labels for assigning the colours
ordered_df <- reorder_labels(long_prot_counts$long_prot_baseline
                             )

# Extract all unique labels one exemplary dataframe
all_labels <- levels(ordered_df$label)

# Generate consistent colors (with viridis) for the same labels
label_colors <- setNames(viridis(length(all_labels), option = "D"), all_labels)

# Extend label_colors with fallback for missing labels
missing_labels <- setdiff(all_labels, names(label_colors))

if (length(missing_labels) > 0) {
  fallback_colors <- RColorBrewer::brewer.pal(12, "Set3")
  fallback_colors <- rep(fallback_colors, length.out = length(missing_labels))
  names(fallback_colors) <- missing_labels

  label_colors <- c(label_colors, fallback_colors)

  warning("The following labels were missing from 'label_colors' and were assigned fallback colors: ",
          paste(missing_labels, collapse = ", "))
}
