# Uniting TMT and DIA for Scalable Low-Input Proteomics

This repository contains the analysis code accompanying the manuscript *"Making the Impossible Possible: Uniting TMT and DIA for Scalable Low-Input-Proteomics"*. 
It benchmarks DIA-TMTpro workflows against DDA-TMT and DIA-LFQ using a two-proteome (human/yeast) interference model on the Orbitrap Astral.

## Repository structure

```
.
├── TMT_RI_coisolation_purity_benchmark/     # Co-isolation/purity benchmark, DDA vs DIA TMT reporter-ion quantification
├── TMTc_complementary_ion_quantification/   # TMTc complementary-ion quantification vs standard TMT RI
├── single_proteome_DIA_TMT_benchmark/       # Single-proteome DIA-TMT benchmark
├── DIA_quant_TMT_labeled_benchmark/         # DIA quantification of TMT-labeled (but not RI-quantified) samples
├── TMTpro_Complementary_Ion_Search/         # Python pre-processing for TMTc complementary-ion search
└── figure6_assembly/                        # Assembles Figure 6 from outputs of the above pipelines
    └── assemble_fig6.R
```

Each subfolder is self-contained and has its own `README.md` with a detailed overview, experimental design, script-by-script run order, input data, and output file descriptions:

- [`TMT_RI_coisolation_purity_benchmark`](./TMT_RI_coisolation_purity_benchmark/README.md)
- [`TMTc_complementary_ion_quantification`](./TMTc_complementary_ion_quantification/README.md)
- [`single_proteome_DIA_TMT_benchmark`](./single_proteome_DIA_TMT_benchmark/README.md)
- [`DIA_quant_TMT_labeled_benchmark`](./DIA_quant_TMT_labeled_benchmark/README.md)
- [`TMTpro_Complementary_Ion_Search`](./TMTpro_Complementary_Ion_Search/README.md)
- [`figure6_assembly`](./figure6_assembly/README.md)

## Requirements

- R ≥ 4.2 (RStudio recommended) for all R pipelines — see each subfolder's README for its specific package list
- Python ≥ 3.x for `TMTpro_Complementary_Ion_Search` (`pandas`, `numpy`, `pyteomics`, `openpyxl`)

## How to run

Run the four benchmark pipelines independently, following the instructions in their respective READMEs. Once their outputs are generated, run `figure6_assembly/assemble_fig6.R` to compile Figure 6.
