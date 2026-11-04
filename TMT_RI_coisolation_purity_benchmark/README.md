# TMT RI Co-Isolation Purity Benchmark

R analysis pipeline for the two-proteome TMT (tandem mass tag) experiment (Human and Yeast). It
evaluates reporter-ion co-isolation purity and its effect on quantitative accuracy across MS
acquisition methods (DDA / DIA, several isolation-window and injection-time settings, and an
input-amount titration), and compares three search/quantification tools: **Proteome Discoverer
(PD)**, **DIA-NN**, and **Spectronaut (SN)**.

The pipeline ingests protein- and PSM-level result tables, cleans and reshapes them, applies a
range of purity / signal-to-noise (S/N) filtering schemes, and produces the publication figures
assembled in `07_arranged_plots.R`.

---

## Experimental design

HeLa/Yeast two-proteome design with mixes `H000Y000`, `H100Y000`, `H100Y010`, `H100Y020`,
`H100Y100` (see `channel_names_groups`, script `01`), run across DDA/DIA acquisition methods and
an input-amount titration, and processed with three search/quantification tools (PD, DIA-NN, SN).

---

## Repository structure

```
.
├── data/                                   # Input result tables (not tracked, see "Input data")
├── plots/                                  # Output directory (created automatically)
├── 01_experimental_design_and_packages.R
├── 02a_functions_for_SN_DIANN.R
├── 02b_const_param_and_func.R
├── 03a_MS_method_eval.R
├── 03b_clean_environment.R
├── 03c_concat_prot_counts.R
├── 04_condition_subsets.R
├── 05a_plots_prot_overlayed.R
├── 05b_plots_prot_overlayed_software.R
├── 05c_plots_prot_faceted.R
├── 05d_plots_prot_BA_plots_with_density.R
├── 05e_prot_ratio_overview.R
├── 05f_merge_data.R
├── 05g_boxplots_faceted.R
├── 05h_violinplots_cv_faceted.R
├── 05i_boxplots_method_wise_all_ratios.R
├── 05j_prot_intensities.R
├── 06_plots_psms_coisolation.R
├── 07_arranged_plots.R
└── README.md
```

---

## Requirements

- **Core packages used:** `tidyverse`, `data.table`, `arrow`, `patchwork`, `cowplot`, `ggridges`,
  `svglite`, `ragg`, `RColorBrewer`, `viridis`, `scales`, `matrixStats`, `grid` (ships with base R).

```r
install.packages(c(
  "tidyverse", "data.table", "arrow", "patchwork", "cowplot", "ggridges",
  "svglite", "ragg", "RColorBrewer", "viridis", "scales", "matrixStats"
))
```

---

## How to run

There are **no `source()` calls** — the scripts share one global environment and **must be run in
order**, because each relies on objects created upstream. Run the whole set in a single R session,
e.g.:

```r
files <- sort(list.files(pattern = "^0[0-9].*\\.R$"))
for (f in files) source(f, echo = TRUE)
```

or step through them manually in the order below.

| # | Script | What it does (one line) |
|---|---|---|
| 1 | `01_experimental_design_and_packages.R` | Load packages; build experimental- and method-design plots and TMT channel metadata. |
| 2 | `02a_functions_for_SN_DIANN.R` | Import/parsing helper functions for Spectronaut and DIA-NN tables. |
| 3 | `02b_const_param_and_func.R` | Import result files; define constants and the clean/reshape/label/filter functions; build `long_*` tables. |
| 4 | `03a_MS_method_eval.R` | Build `prot_results`/`psm_results` (one entry per filter scheme), `top20_prot`, averaged row counts. |
| 5 | `03b_clean_environment.R` | Remove large intermediate objects to free memory. |
| 6 | `03c_concat_prot_counts.R` | Concatenate per-dataset protein counts and build loss annotations. |
| 7 | `04_condition_subsets.R` | Define condition-subset label vectors and the `label_colors` map. |
| 8 | `05a_plots_prot_overlayed.R` | Overlayed protein-count plots (input titration / isolation window) → Fig2. |
| 9 | `05b_plots_prot_overlayed_software.R` | Overlayed protein counts across PD/DIA-NN/SN → sup_fig4. |
| 10 | `05c_plots_prot_faceted.R` | Faceted protein counts across the S/N titration → sup_fig2. |
| 11 | `05d_plots_prot_BA_plots_with_density.R` | Bland-Altman ratio plots with marginal density → Fig1/Fig3/sup_fig3 (+ PNG dump). |
| 12 | `05e_prot_ratio_overview.R` | Protein ratio boxplots and CV violin plots per condition → Fig3. |
| 13 | `05f_merge_data.R` | Merge ratio + CV data into one long table `all_ratios_and_cv`. |
| 14 | `05g_boxplots_faceted.R` | Faceted ratio boxplots (S/N titration + software) → sup_fig2/sup_fig4. |
| 15 | `05h_violinplots_cv_faceted.R` | Faceted CV violin plots (S/N titration + software) → sup_fig2/sup_fig4. |
| 16 | `05i_boxplots_method_wise_all_ratios.R` | Method-wise boxplots of all ratios (saved as `.rds`). |
| 17 | `05j_prot_intensities.R` | Summed reporter intensities per channel → Fig1. |
| 18 | `06_plots_psms_coisolation.R` | Co-isolation purity ridge plots with % precursors > 0.7 → Fig2. |
| 19 | `07_arranged_plots.R` | Assemble Fig1–3 and sup_fig2–4 and write the SVGs (terminal deliverable). |

---

## Input data

Input data are Proteome Discoverer result exports post-processed with READ in `.txt` format.

The import step (`02b`, via the helpers in `02a`) reads raw result tables from a `data/`
directory (relative paths, resolved from the working directory):


- **Proteome Discoverer** protein and PSM exports (`.txt`), including the `Pure-*`
  co-isolation-corrected abundance channels and a `Co-Isolation Purity` column.
- **DIA-NN** protein-group output as **Parquet** (`arrow::read_parquet`).
- **Spectronaut** export as **CSV** (`data.table::fread`), including `PG.FastaFiles` for organism
  assignment.

> The data files themselves are **not** included in this repository. Place them under `data/`
> with the expected names before running.

```r
dir.create("plots", showWarnings = FALSE)
```

---

## Outputs

- **Figures** written to `plots/` as SVG: `Fig1.svg`, `Fig2.svg`, `Fig3.svg`,
  `Fig3_ratio20v10.svg`, `sup_fig2.svg`, `sup_fig3.svg`, `sup_fig4.svg`, plus the individual
  `p_*` panel SVGs and the design schematics (`exp_design.svg`, `method_design.svg`).
- **Per-dataset Bland-Altman PNGs** under `plots/prot_ratios/` (from `05d`).
- **Serialized plot objects:** `p1_multi005_DDA.rds`, `p2_multi005_DIA.rds`.
- **Key in-memory objects:** `prot_results`, `psm_results`, `top20_prot`, `long_prot_counts`,
  `loss_annotations_all`, `all_ratios_and_cv`, and the many `p_*` ggplot objects consumed by `07`.
