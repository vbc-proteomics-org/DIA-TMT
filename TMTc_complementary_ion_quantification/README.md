# TMTc Complementary-Ion Quantification

Analysis pipeline for TMTc (complementary-ion) quantification on a human–yeast protein mix,
compared against standard TMT reporter-ion (RI) quantification. The scripts are meant to be run
**in order within a single, shared R session**: later scripts rely on objects created by earlier
ones.

---

## Experimental design

Channels 126 / 131N / 135N. Expected yeast ratio 5:1 (log2 ≈ 2.32); human 1:1 between
135N/131N. Colours: Human `#56B4E9`, Yeast `#E69F00`.

---

## Repository structure

```
.
├── data/
│   ├── prot_TMTc/F2500/   # TMTc protein-level tables (*.csv), one per replicate
│   ├── baseline/          # Standard TMT reporter-ion tables (*.txt, baseline)
│   └── standard/          # Standard TMT reporter-ion tables (*.txt, standard)
├── plots/                 # Output directory (subfolders created automatically)
├── 00_setup_packages_and_experimental_design.R
├── 01_TMTc_import_and_noise_subtraction.R
├── 02_TMT_reporter_ion_import.R
├── 03_TMT_reporter_ion_filtering.R
├── 04_protein_counts_TMT_vs_TMTc.R
├── 05_TMTc_ratio_boxplot.R
├── 06_figure_TMTc.R
└── README.md
```

---

## Requirements

- **Packages used:** `tidyverse`, `data.table`, `readxl`, `writexl`, `reshape2`, `zoo`,
  `patchwork`, `cowplot`, `gridExtra`, `ggpubr`, `ggridges`, `ggtext`, `ggpattern`, `ggdist`,
  `colorspace`, `paletteer`, `RColorBrewer`, `viridis`, `Hmisc`, `DescTools`, `matrixStats`,
  `svglite`, `grid` (ships with base R).

```r
install.packages(c(
  "tidyverse", "data.table", "readxl", "writexl", "reshape2", "zoo",
  "patchwork", "cowplot", "gridExtra", "ggpubr", "ggridges", "ggtext",
  "ggpattern", "ggdist", "colorspace", "paletteer", "RColorBrewer",
  "viridis", "Hmisc", "DescTools", "matrixStats", "svglite"
))
```

---

## How to run

| # | Script | What it does |
|---|--------|--------------|
| 0 | `00_setup_packages_and_experimental_design.R` | Loads all packages, defines `plots_dir` and `channel_names_groups`, and builds the experimental-design schematic (`exp_design_multi007.rds`). |
| 1 | `01_TMTc_import_and_noise_subtraction.R` | Imports the TMTc protein tables, subtracts the 126 background per species, and builds summed-intensity plots, protein counts and log2-ratio boxplots. |
| 2 | `02_TMT_reporter_ion_import.R` | Imports and cleans the standard TMT reporter-ion protein tables (baseline + standard): renames channels, keeps relevant columns, removes contaminants, annotates organism, computes means/ratios, and builds the long protein tables. |
| 3 | `03_TMT_reporter_ion_filtering.R` | Filters and log2-transforms the long protein tables into `prot_results`. |
| 4 | `04_protein_counts_TMT_vs_TMTc.R` | Counts quantified Human/Yeast proteins for TMT RI and TMTc, and builds the three-way comparison bar plot (`p_three_way_total`). |
| 5 | `05_TMTc_ratio_boxplot.R` | Builds the TMTc log2-ratio boxplot for DIA120 50ng (`p4_multi007.rds`). |
| 6 | `06_figure_TMTc.R` | Assembles the final figure from panels A/B/C (`figure_TMTc.svg`). |

---

## Input data
Input data are Proteome Discoverer result exports post-processed with READ in `.txt` format
and FragPipe result exports post-processed with Python scripts in directory **TMTpro_Complementary_Ion_Search**

Paths are relative to the working directory:

- `data/prot_TMTc/F2500/` — TMTc protein-level tables (`*.csv`), one per replicate.
- `data/baseline/` — standard TMT reporter-ion tables (`*.txt`, baseline).
- `data/standard/` — standard TMT reporter-ion tables (`*.txt`, standard).

> The data files themselves are **not** included in this repository.

---

## Outputs

- `exp_design_multi007.rds` — experimental-design plot object (script 0).
- `p4_multi007.rds` — TMTc ratio boxplot object (script 5).
- `plots/figure5_TMTc.svg` — assembled final figure (script 6).
- Additional per-step SVG/PNG plots are written under `plots/` and its subfolders.
