# Code for Figures 5 and 6

This directory contains the R code used to reproduce Figures 5 and 6 in the manuscript:

**Scientific data analysis on hyperspectral functional data: categorical exploratory data analysis (CEDA) coupled with Kolmogorov’s randomness-proper**

## Purpose

The script performs response-side hierarchical clustering of the Tecator meat spectroscopy data.

- **Figure 5A:** Hierarchical clustering based on protein, moisture, and fat.
- **Figure 5B:** Hierarchical clustering based on protein and fat.
- **Figure 6A:** The four clusters from Figure 5A displayed in the three-dimensional response space.
- **Figure 6B:** The four clusters from Figure 5B displayed in the same three-dimensional response space.

A least-squares regression plane is superimposed on both panels of Figure 6:

$ \text{moisture}=\beta_0+\beta_1\text{fat}+\beta_2\text{protein}+\varepsilon. $

Figure 6 is generated with `plotly` using the same workflow as the original analysis code:

1. construct the 3D point cloud with `plot_ly(..., type = "scatter3d")`;
2. add the fitted regression plane with `add_surface()`;
3. specify the three axes with `layout(scene = ...)`.

## File

- `Fig5_Fig6_Tecator.R` — complete analysis and figure-generation script.

## Data

The analysis uses the `Tecator` dataset loaded through the R packages used in the original analysis. The response variables are:

- fat
- protein
- moisture

The original sample order is retained, and `sample_id` corresponds to the row number in the loaded dataset.

## Statistical procedure

For each hierarchical clustering analysis:

1. Euclidean distances are calculated from the unstandardized response variables.
2. Agglomerative hierarchical clustering is performed using Ward’s minimum-variance criterion (`ward.D2`).
3. The dendrogram is cut into four clusters.
4. Cluster labels are ordered from the leftmost to the rightmost dendrogram branch so that branch colors and cluster numbers remain consistent.

No random-number generation is used, so a random seed is not required.

## Required R packages

```r
install.packages(c(
  "dendextend",
  "plotly",
  "htmlwidgets",
  "fds",
  "fsemipar"
))
```

## Running the code

Place `Fig5_Fig6_Tecator.R` in the repository root and run:

```r
source("Fig5_Fig6_Tecator.R")
```

The script creates an `outputs/` directory automatically.

## Generated outputs

- `Figure5_hierarchical_clustering.pdf`
- `Figure5_hierarchical_clustering.png`
- `Figure6A_three_variable_clustering.html`
- `Figure6B_two_variable_clustering.html`
- `Fig5_Fig6_cluster_memberships.csv`
- `Figure6_regression_model_summary.txt`
- `sessionInfo.txt`

The two Figure 6 HTML files are self-contained. They can be opened in a web browser and rotated interactively to inspect the cluster geometry and fitted regression plane.

## Scope of this script

The manually defined fifth group and its absorbance-curve plot are not included because Figures 5 and 6 in the manuscript compare two **four-cluster** categorizations. Any separate exploratory analysis of those eight observations should be deposited as a supplementary script rather than included in the code for Figures 5 and 6.
