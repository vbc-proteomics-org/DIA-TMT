# Single-Proteome DIA-TMT Benchmark

R code for the TMT proteomics purity analysis of the single-proteome experiment presented in the
DIA-TMT study. The pipeline imports Proteome Discoverer output files which are post-processed
with READ, applies filtering and log2 transformation, and produces the six-panel
method-comparison figure (`fig4_multi008A_StoN_grouping.svg`).

---

## Experimental design

The study uses an 18-channel TMT experiment with three HeLa dilution series:

| Group | Channels | HeLa ratio |
|-------|----------|------------|
| H000 (label-only) | 126, 129N, 131C, 132N, 134C, 135N | 0 % |
| H010 | 127N, 127C, 128N, 128C | 10 % |
| H050 | 129C, 130N, 130C, 131N | 50 % |
| H100 | 132C, 133N, 133C, 134N | 100 % |

Four acquisition methods are compared at 53.75 ng input: DDA 0.5 m/z IT20, DIA 1.0 m/z IT8,
DIA 1.0 m/z HCD29, and DIA 1.0 m/z no-loop control.

---

## Repository structure

```
.
├── data/
│   ├── baseline/                  # Raw Proteome Discoverer exports — baseline condition
│   └── standard_StoN_grouping_20/ # Raw Proteome Discoverer exports — S/N grouping condition
├── plots/                         # Output directory (created automatically by ggsave)
├── 01_experimental_design_and_packages.R
├── 02_constants_functions_and_data_import.R
├── 03_filter_transform_and_counts.R
├── 04_clean_environment.R
├── 05_concatenate_protein_counts.R
├── 06_condition_subsets_and_colors.R
├── 07_plot_protein_overlay.R
├── 08_plot_bland_altman_density.R
├── 09_plot_ratio_boxplots_and_cv.R
├── 10_plot_protein_intensities.R
├── 11_arranged_plots.R
└── README.md
```

---

## Requirements

- **Packages used:** `tidyverse`, `ggplot2`, `data.table`, `readr`, `purrr`, `reshape2`, `zoo`,
  `patchwork`, `cowplot`, `gridExtra`, `ggpubr`, `ggridges`, `ggtext`, `ggpattern`, `ggdist`,
  `ggrepel`, `colorspace`, `RColorBrewer`, `viridis`, `scales`, `Hmisc`, `DescTools`,
  `matrixStats`, `writexl`, `svglite`, `grid` (ships with base R).

```r
install.packages(c(
  "writexl", "cowplot", "patchwork", "gridExtra", "tidyverse", "ggplot2",
  "svglite", "Hmisc", "RColorBrewer", "purrr", "zoo", "reshape2", "viridis",
  "ggpattern", "colorspace", "grid", "ggpubr", "readr", "data.table",
  "ggridges", "scales", "ggdist", "DescTools", "matrixStats", "ggtext", "ggrepel"
))
```

The `plots/` directory is created automatically by `ggsave()` if it does not exist.

---

## How to run

Scripts are run sequentially in the order given by their numeric prefix. Each script depends on
objects produced by all preceding scripts.

| # | Script | What it does |
|---|--------|--------------|
| 1 | `01_experimental_design_and_packages.R` | Loads all packages; defines TMT channel layout (`channel_names_groups`); produces the experimental-design overview figure. |
| 2 | `02_constants_functions_and_data_import.R` | Defines all helper functions; imports, cleans and labels the raw `.txt` files from both data directories; produces `long_psm_*` and `long_prot_*` data frames. |
| 3 | `03_filter_transform_and_counts.R` | Applies purity, S/N and unambiguity filters; log2-transforms protein data; counts identified and quantified rows per dataset. |
| 4 | `04_clean_environment.R` | Removes large intermediate list objects and runs garbage collection to free memory. |
| 5 | `05_concatenate_protein_counts.R` | Pairs identified vs. quantified count summaries; computes quantification loss; produces `long_prot_counts` and `loss_annotations`. |
| 6 | `06_condition_subsets_and_colors.R` | Defines method/condition label subsets; builds the viridis colour map for condition labels. |
| 7 | `07_plot_protein_overlay.R` | Defines `plot_quan_protein_overlay()`; produces the identified vs. quantified protein-count overlay figure. |
| 8 | `08_plot_bland_altman_density.R` | Defines `plot_blant_altmann_prot()`; produces Bland-Altman scatter + density figures for the two method comparisons used in the final figure. |
| 9 | `09_plot_ratio_boxplots_and_cv.R` | Defines `plot_boxplot_vertical_prot()` and `plot_violin_cv_prot()`; produces the ratio boxplot and CV violin figures. |
| 10 | `10_plot_protein_intensities.R` | Defines `sum_prot_intensities_with_mix()` and `plot_reporter_intensities_linear()`; produces the reporter-ion intensity figure. |
| 11 | `11_arranged_plots.R` | Assembles all panel plots into the final six-panel figure and saves it as `plots/fig4_multi008A_StoN_grouping.svg`. |

---

## Input data

Input data are Proteome Discoverer result exports post-processed with READ in `.txt` format,
organised in two sub-directories under `data/`:

- `data/baseline/` — baseline TMT acquisition (no S/N grouping filter)
- `data/standard_StoN_grouping_20/` — acquisition with S/N grouping threshold 20

> The raw data files are not included in this repository. Update the `data_dir_*` paths at the
> top of `02_constants_functions_and_data_import.R` if your data are stored elsewhere.

---

## Outputs

Running the full pipeline writes the following SVG files to `plots/`:

| File | Produced by | Panel in final figure |
|------|-------------|----------------------|
| `exp_design_multi008A.svg` | `01` | — (standalone) |
| `p_prot_intensities_DDA_IT20_54ng_new.svg` | `10` | A |
| `p_overlay_54ng_StoN_grouping.svg` | `07` | B |
| `p_prot_ratio_boxplots_H100H010_54ng_StoN_grouping.svg` | `09` | C |
| `p_violin_methods_cv_H100_54ng_StoN_grouping.svg` | `09` | D |
| `p_prot_ratio_DDA_IT20_54ng_rep1_H100vsH010_StoN_grouping.svg` | `08` | E |
| `p_prot_ratio_DIA_IT8_HCD29_54ng_rep1_H100vsH010_StoN_grouping.svg` | `08` | F |
| `fig4_multi008A_StoN_grouping.svg` | `11` | Final figure |
