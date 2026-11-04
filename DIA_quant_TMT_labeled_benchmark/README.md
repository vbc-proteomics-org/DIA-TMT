# DIA Quantification of TMT-Labeled Samples (Benchmark)

R code to analyse the two-proteome experiment: TMT-labeled samples processed in Spectronaut via
DIA **without using TMT reporter-ion (RI) quantification**, used to check how accurately known
human/yeast mixing ratios are recovered. The main output is the log2 ratio boxplot panel
(`p3_multi009.rds`) used in **Figure 6** of the manuscript.

---

## Experimental design

Each sample has the same HeLa (human) background. Yeast protein is spiked in at four relative
levels. Each group has three replicate samples, and every sample was injected at the same volume.

| Group      | Samples (channel IDs) | Human (rel.) | Yeast (rel.) |
|------------|------------------------|:------------:|:------------:|
| `H100Y000` | 128C, 128N, 129N       | 1            | 0            |
| `H100Y010` | 130C, 130N, 131N       | 1            | 1            |
| `H100Y020` | 132C, 132N, 133N       | 1            | 2            |
| `H100Y100` | 134C, 134N, 135N       | 1            | 10           |

**Expected log2 ratios**

| Comparison            | Human | Yeast             |
|-----------------------|:-----:|:------------------:|
| H100Y100 / H100Y010   | 0     | log2(10) = 3.32   |
| H100Y100 / H100Y020   | 0     | log2(5) = 2.32    |
| H100Y020 / H100Y010   | 0     | log2(2) = 1.00    |
| H100Y000 / any group  | 0     | not defined (no yeast in H100Y000) |

The paper only uses the Spectronaut report with **MS2-level quantification and cross-run
normalisation** (tag `MS2_Norm`).

---

## Repository structure

```
.
├── 00_shared_functions.R     # Constants, data processing, theme and plot functions (no I/O)
├── 01_DIA_analysis.R         # Full QC / exploratory figure set (optional)
├── 02_DIA_boxplot_ratios.R   # p3_multi009.rds
├── data/                     # Spectronaut report (not tracked, see "Input data")
└── plots/                    # Generated figures (created automatically)
```


---

## Requirements

- **Packages used:** `tidyverse`, `scales`, `cowplot`, `patchwork`

```r
install.packages(c("tidyverse", "scales", "cowplot", "patchwork"))
```

---

## How to run

Set the working directory to the repository root, then run:

```r
source("02_DIA_boxplot_ratios.R")   # Figure 6 panel (required)
source("01_DIA_analysis.R")         # Full figure set (optional)
```

Both scripts load `00_shared_functions.R` themselves.

| Script | Role |
|--------|------|
| `00_shared_functions.R` | Sourced by the other scripts. No file I/O, draws no plots. Defines design constants (channel-to-group map, yeast spike-in amounts, comparison pairs, Okabe–Ito colour palette), `process_sn_file()` (reads a Spectronaut non-pivot report, pivots to wide format, assigns species, removes contaminants/mixed-species groups, log2-transforms, computes group means/CVs/pairwise log2 ratios), `theme_dia()`, and the plotting functions used in `01`. |
| `01_DIA_analysis.R` *(optional – not needed for Figure 6)* | Produces the full figure set for the `MS2_Norm` dataset (see Outputs). |
| `02_DIA_boxplot_ratios.R` | Makes the boxplot of the three yeast-containing comparisons used in Figure 6. Sources `00_shared_functions.R` and loads the data itself — does **not** require `01` to be run first. |

---

## Input data

Put the Spectronaut **non-pivot** protein-group report in `data/`. The report must include these
columns:

| Column                   | Use |
|--------------------------|-----|
| `R.FileName`             | Sample identifier. Must end with the channel ID (e.g. `..._128C`) |
| `PG.ProteinAccessions`   | Protein group accession |
| `PG.ProteinDescriptions` | Protein description |
| `PG.Organisms`           | Species assignment |
| `PG.Quantity`            | Protein group quantity (linear scale) |

---

## Outputs

| Output (in `plots/`) | Content |
|-----------------------|---------|
| `experimental_design_expected_ratios_DIA_SN.svg` | Sample composition and expected ratios (script `01`) |
| `MS2_Norm/channel_intensities_MS2_Norm.svg` | Summed protein intensity per sample (script `01`) |
| `MS2_Norm/channel_intensities_by_species_MS2_Norm.svg` | Summed intensity per sample, split by species (script `01`) |
| `MS2_Norm/BA_<groupA>vs<groupB>_MS2_Norm.svg` | Group-level Bland–Altman plots, 6 comparisons (script `01`) |
| `MS2_Norm/BA_channels_<groupA>vs<groupB>_MS2_Norm.svg/.png` | Sample-by-sample Bland–Altman grids, 3×3 (script `01`) |
| `MS2_Norm/ratio_overview_MS2_Norm.svg` | Boxplot of all six log2 ratios (script `01`) |
| `plots/p_boxplot_DIA_LFQ_multi009.svg` | Boxplot of Figure 6 comparisons, vector graphic (script `02`) |
| `p3_multi009.rds` | ggplot object, read by the Figure 6 assembly script (script `02`) |

---
