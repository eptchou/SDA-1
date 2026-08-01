# SDA

Code and data to reproduce the figures in the paper:

> **Scientific data analysis on hyperspectral functional data: categorical
> exploratory data analysis (CEDA) coupled with Kolmogorov's randomness-proper**

Here *SDA* stands for **Scientific Data Analysis**. The analysis is applied to
two spectral data sets:

- **Tecator fat spectra** (fat / protein / moisture), loaded from the
  `fds` / `fsemipar` packages.
- **MIR fruit-puree spectra** (strawberry vs. non-strawberry), provided in `data/`.

All figures are produced by a single, sectioned R script: `reproduce_figures.R`.

## Repository layout

```
SDA-1/
├── reproduce_figures.R   # main script; one section per figure
├── data/                 # input data sets (CSV)
├── README.md
└── LICENSE
```

## Data

| File | Description |
|------|-------------|
| `data/PMF_1st_cat.csv` | Fat spectra, 1st-difference categorisation |
| `data/fatspectrum_2st_16cluster.csv` | Fat spectra, 2nd-difference 16-cluster table |
| `data/MIR_Fruit_purees.csv` | MIR strawberry / non-strawberry spectra |
| `data/strawberry_1stc_16_n_983.csv` | Strawberry 1st-difference categorisation |
| `data/strawberry_2ndc_16_n_983.csv` | Strawberry 2nd-difference categorisation |

The Tecator fat spectra are not stored here; they are loaded at run time via
`data(Tecator)` from the `fsemipar` / `fds` packages.

## Requirements

R (≥ 4.0) with the following packages.

**CRAN**

```r
install.packages(c(
  "tidyverse", "dplyr", "tidyr", "ggplot2", "magrittr",
  "readxl", "writexl", "patchwork", "cowplot", "dendextend",
  "circlize", "RColorBrewer", "ggrepel", "infotheo", "plotly",
  "fds", "fsemipar"
))
```

**Bioconductor**

```r
install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")
```

## How to run

1. Clone the repository and set the working directory to the repository root
   (e.g. open the folder as an RStudio project, or run `setwd("path/to/SDA-1")`).
   All data paths in the script are **relative** to this root.
2. Install the packages listed above.
3. Open `reproduce_figures.R`. Each figure lives in its own section, delimited by
   `# ====` separators. Every section is self-contained — it reloads the data it
   needs — so you can run a single section to reproduce one figure, or source the
   whole file.

Some sections write a cached result to `output/dynamic_2nd_split_results.xlsx`.
If the file is absent it is regenerated automatically (this step involves a
Monte-Carlo simulation and may take a few minutes); on subsequent runs the cache
is reused.

## Figure index

| Figure(s) | Section heading in `reproduce_figures.R` |
|-----------|-------------------------------------------|
| 3 | functional-curves |
| 4 | PFM-Categorization-1st-diff |
| 5 | Two-HC-trees-PFM |
| 7 | motif-selection |
| 8 | co-randomness-z1-z2 |
| 9 | PFM-z1-tree-result |
| 10 | PFM-Z1-678-split-by-Z24 |
| 11 | PFM-Figure1-Y1Y2 |
| 12 | PFM-Z3-tree-result |
| 13 | PFM-heatmap-Z1-2nd-Z3-2nd-revised |
| 14 | 4-panels-PFM-heatmap-Z1-2nd-Z3-2nd-revised |
| 15 | strawberry-curves-2 |
| 16 | Categorization-1st-diff |
| 17 | straw-1st-2nd-diff |
| 18 | straw-Z3-tree-split |
| 19 | straw-Z6-tree-split |
| 20 | straw-Z6-tree-split |
| 21 | straw-z4-2nd-tree |
| 22 | straw-Z7-Z4-tree |
| 23 | strawberry-odds-vs-overlap |
| 24, 27, 30 | strawberry-order1-heatmap |
| 26 | Strawberry-split-Overlap-vs-Odds-order2 |
| 28, 29 | straw-combine1st2nd-odds001 |
| 31, 32 | strawberry2nd-heatmap-odds |

## Reproducibility notes

- Simulation-based sections set a fixed seed (`set.seed(42)`) so results are
  reproducible across runs.
- Package versions may affect minor rendering details; see the `sessionInfo()`
  of your environment if exact reproduction is required.

## License

Released under the MIT License. See [LICENSE](LICENSE).
