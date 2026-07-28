# Code for Figure 3

This directory contains the R code used to reproduce Figure 3 in the manuscript:

**Scientific data analysis on hyperspectral functional data: categorical exploratory data analysis (CEDA) coupled with Kolmogorov’s randomness-proper**

## Purpose

Figure 3 displays the 215 Tecator absorbance curves using two response-based color schemes:

- **Panel A:** blue for samples with fat content greater than 20%, and orange otherwise.
- **Panel B:** blue for samples with protein content greater than 16%, and orange otherwise.

The figure provides an initial visual assessment of whether the local shapes of the spectral curves vary with the response-side categories.

## File

- `Fig3_Tecator.R` — complete code for data loading, categorization, plotting, and output generation.

## Data

The script uses the `Tecator` dataset from the `fsemipar` package.

The following variables are used:

- `Tecator$absor.spectra` — absorbance spectra;
- `Tecator$fat` — fat content;
- `Tecator$protein` — protein content.

The wavelength sequence is defined as 852–1050 nm in increments of 2 nm, following the original analysis code.

## Required R packages

```r
install.packages(c(
  "fsemipar",
  "plotly",
  "htmlwidgets"
))
```

## Running the code

Place `Fig3_Tecator.R` in the repository root and run:

```r
source("Fig3_Tecator.R")
```

The script creates an `outputs/` directory automatically.

## Generated outputs

- `Figure3_Tecator_combined.html`
- `Figure3A_fat_content.html`
- `Figure3B_protein_content.html`
- `Figure3_sample_categories.csv`
- `sessionInfo.txt`

The HTML files are self-contained and can be opened directly in a web browser.

## Plotting details

The curves are produced with Plotly. To keep the legend readable, the script displays one legend entry per response category rather than one entry for every sample. Sample IDs and response values remain available in the interactive hover labels.

No random-number generation is used, so a random seed is not required.
