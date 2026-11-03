# 3-plex TMTpro Complementary-Ion Quantification

This repository contains two Python scripts for post-processing **data-independent
acquisition (DIA)** proteomics data from **3-plex TMTpro** experiments using the
complementary-ion approach. The first script extracts complementary-ion
intensities from DIA `mzML` spectra (matched to FragPipe PSM results), and the
second aggregates those PSM-level measurements up to peptide- and protein-level
quantities. These scripts are adapted for an experimental setup in which the
**TMT 126 channel (comp_ion3) is empty**, so that channel is still computed and
reported but is excluded from filtering and ranking.

## Scripts

- **`3plex_TMTpro_complementary_ion_search_adapted.py`**
  Reads a FragPipe `PSM.xlsx` file and the corresponding DIA `mzML` file,
  computes the complementary-ion m/z values for each PSM, extracts their
  intensities (with isotope-peak corrections, M+1/M+4 ratio modeling, and
  intensity thresholds), and writes PSM-level and peptide-level summary tables as
  `.xlsx`. Ranking of the top PSMs per peptide uses only comp_ion1 + comp_ion2
  (the empty 126/comp_ion3 channel is excluded).

- **`PSM_to_protein_level_adapted.py`**
  Reads the PSM-level `.xlsx` output from script 1, filters to unique peptides,
  and aggregates PSM → peptide → protein. A normalization factor (Hyperscore if
  present, otherwise Probability, otherwise a default of 1.0) is used to compute
  weighted intensities. PSM and peptide ranking exclude comp_ion3 (126). Outputs
  peptide-level and protein-level quantification tables as `.xlsx`.

## Usage

The scripts are run in sequence. **Output from script 1 feeds into script 2.**

1. Open `3plex_TMTpro_complementary_ion_search_adapted.py` and set the input file
   paths in the `__main__` block (`fragpipe_psm_path` → your FragPipe `PSM.xlsx`,
   `dia_mzml_path` → your DIA `mzML`). Run it to produce the PSM-level file
   (e.g. `DIA1_psm_int.xlsx`) and a peptide-level summary.
2. Open `PSM_to_protein_level_adapted.py` and set `file_path` in `main()` to the
   PSM-level `.xlsx` produced in step 1. Run it to produce peptide-level and
   protein-level quantification tables.

```bash
python 3plex_TMTpro_complementary_ion_search_adapted.py
python PSM_to_protein_level_adapted.py
```

## Requirements

- Python 3.x
- pandas
- numpy
- pyteomics
- openpyxl

```bash
pip install pandas numpy pyteomics openpyxl
```

## Credits / Attribution

These scripts are **adapted** from the original TMTproC-DIA Complementary-Ion
Search Python scripts released under the MIT license.

**Original paper (APA):**

> Wang, Z., Liu, P.-K., Lu, H., & Li, L. (2025). Multiplexed Data-Independent
> Acquisition-Based Proteomics Enabled by TMTpro complementary Ions.
> *Analytical Chemistry, 97*(48), 26419. https://doi.org/10.1021/acs.analchem.5c03563

**Original code (Zenodo):** https://zenodo.org/records/14567079
(DOI: 10.5281/zenodo.14567079)

**Authors / affiliation:** Zicong Wang, Peng-Kai Liu, Haiyan Lu, Lingjun Li
(University of Wisconsin–Madison)

### Summary of key modifications

The original scripts were adapted for a **3-plex TMTpro setup in which the TMT 126
channel (comp_ion3) is **empty**. The main changes are:

- The empty 126 channel (comp_ion3) is still computed and reported, but is 
  **not required to exist and is not used as a hard filter or ranking criterion**.
- PSM and peptide ranking use only comp_ion1 + comp_ion2 (a new
  `total_intensity_no126` metric / `ranking_columns` parameter), so the empty 126
  channel does not distort top-PSM/top-peptide selection.
- Intensity thresholds were tuned (`ION_INTENSITY_THRESHOLD = 2500`, new
  `M1_INTENSSITY_THRESHOLD = 250` for isotopic peaks).
- Hyperscore is no longer a mandatory column; the normalization factor is chosen
  flexibly (Hyperscore → Probability → default).
- Output format changed from CSV to XLSX, and file paths were updated to connect
  the two scripts.

## License

This repository is released under the **MIT license**. See the [LICENSE]
file for details.
