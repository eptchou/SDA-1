# ==============================================================================
# Figures 5 and 6: Response-side hierarchical clustering of the Tecator data
#
# This script reproduces:
#   Figure 5A: Four-cluster HC based on protein, moisture, and fat
#   Figure 5B: Four-cluster HC based on protein and fat
#   Figure 6A: Figure 5A cluster memberships displayed in 3D
#   Figure 6B: Figure 5B cluster memberships displayed in 3D
#
# Figure 6 is generated using the original plotly workflow:
#   plot_ly(scatter3d) -> add_surface(regression plane) -> layout(scene)
# ==============================================================================

# ---- 0. Package check ---------------------------------------------------------

required_packages <- c(
  "dendextend",
  "plotly",
  "htmlwidgets",
  "fds",
  "fsemipar"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop(
    paste0(
      "The following packages are required but not installed: ",
      paste(missing_packages, collapse = ", "),
      "\nInstall them before running this script."
    ),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(dendextend)
  library(plotly)
  library(htmlwidgets)
  library(fds)
  library(fsemipar)
})

# ---- 1. Output directory and data --------------------------------------------

output_dir <- "outputs"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

data("Tecator")

if (!exists("Tecator")) {
  stop("The Tecator dataset could not be loaded.", call. = FALSE)
}

response_data <- data.frame(
  sample_id = seq_along(Tecator$fat),
  fat = Tecator$fat,
  protein = Tecator$protein,
  moisture = Tecator$moisture
)

if (anyNA(response_data)) {
  stop("Missing values were detected in the response variables.", call. = FALSE)
}

# The same four colors are used in Figures 5 and 6.
# Cluster numbers follow the left-to-right order of the dendrogram branches.
cluster_colors <- c(
  "#636EFA",
  "#EF553B",
  "#00CC96",
  "#AB63FA"
)

# ---- 2. Helper functions ------------------------------------------------------

# Relabel cutree() results according to the left-to-right branch order of the
# dendrogram. This keeps cluster numbers and branch colors synchronized.
get_ordered_clusters <- function(hc_object, k = 4L) {
  raw_clusters <- stats::cutree(hc_object, k = k)
  dendrogram_object <- as.dendrogram(hc_object)
  leaf_order <- order.dendrogram(dendrogram_object)

  branch_order <- unique(raw_clusters[leaf_order])
  ordered_clusters <- match(raw_clusters, branch_order)

  factor(ordered_clusters, levels = seq_len(k))
}

fit_hierarchical_clustering <- function(data, variables, k = 4L) {
  distance_matrix <- stats::dist(
    data[, variables, drop = FALSE],
    method = "euclidean"
  )

  hc_object <- stats::hclust(
    distance_matrix,
    method = "ward.D2"
  )

  list(
    hc = hc_object,
    dendrogram = as.dendrogram(hc_object),
    cluster = get_ordered_clusters(hc_object, k = k)
  )
}

draw_figure_5 <- function() {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)

  par(
    mfrow = c(1, 2),
    mar = c(3.0, 3.5, 3.5, 1.0),
    oma = c(0, 0, 0, 0)
  )

  dendrogram_a <- color_branches(
    hc_three_variables$dendrogram,
    k = 4,
    col = cluster_colors,
    groupLabels = TRUE
  )

  plot(
    dendrogram_a,
    main = "(A) Protein, moisture, and fat",
    xlab = "",
    ylab = "Height",
    leaflab = "none"
  )

  dendrogram_b <- color_branches(
    hc_two_variables$dendrogram,
    k = 4,
    col = cluster_colors,
    groupLabels = TRUE
  )

  plot(
    dendrogram_b,
    main = "(B) Protein and fat",
    xlab = "",
    ylab = "Height",
    leaflab = "none"
  )
}

# Create a rotatable 3D scatter plot and superimpose the regression plane.
# This follows the original plotly structure used in the analysis code.
create_figure_6_plot <- function(data, cluster_variable, plot_title) {
  plot_data <- data
  plot_data$cluster <- factor(
    plot_data[[cluster_variable]],
    levels = 1:4
  )

  fig_markers <- plot_ly(
    data = plot_data,
    x = ~fat,
    y = ~protein,
    z = ~moisture,
    color = ~cluster,
    colors = cluster_colors,
    mode = "markers",
    type = "scatter3d",
    marker = list(size = 4),
    text = ~paste0(
      "Sample: ", sample_id,
      "<br>fat: ", round(fat, 3),
      "<br>protein: ", round(protein, 3),
      "<br>moisture: ", round(moisture, 3),
      "<br>cluster: ", cluster
    ),
    hoverinfo = "text",
    showlegend = TRUE
  )

  fig_surface <- add_surface(
    fig_markers,
    x = fat_sequence,
    y = protein_sequence,
    z = regression_surface,
    opacity = 0.15,
    colorscale = list(
      c(0, 1),
      c("tan", "blue")
    ),
    showscale = FALSE,
    showlegend = FALSE,
    hoverinfo = "skip",
    type = "surface"
  )

  layout(
    fig_surface,
    title = plot_title,
    scene = list(
      xaxis = list(title = "fat"),
      yaxis = list(title = "protein"),
      zaxis = list(title = "moisture"),
      aspectmode = "data"
    ),
    legend = list(
      title = list(text = "Cluster")
    )
  )
}

# ---- 3. Hierarchical clustering for Figure 5 ---------------------------------

# Figure 5A: clustering based on all three response variables.
hc_three_variables <- fit_hierarchical_clustering(
  data = response_data,
  variables = c("protein", "moisture", "fat"),
  k = 4
)

# Figure 5B: clustering based on protein and fat only.
hc_two_variables <- fit_hierarchical_clustering(
  data = response_data,
  variables = c("protein", "fat"),
  k = 4
)

response_data$cluster_three_variables <- hc_three_variables$cluster
response_data$cluster_two_variables <- hc_two_variables$cluster

# Save cluster memberships for verification and reuse.
utils::write.csv(
  response_data,
  file = file.path(output_dir, "Fig5_Fig6_cluster_memberships.csv"),
  row.names = FALSE
)

# ---- 4. Generate Figure 5 -----------------------------------------------------

grDevices::pdf(
  file.path(output_dir, "Figure5_hierarchical_clustering.pdf"),
  width = 12,
  height = 5.5
)
draw_figure_5()
grDevices::dev.off()

grDevices::png(
  file.path(output_dir, "Figure5_hierarchical_clustering.png"),
  width = 2400,
  height = 1100,
  res = 200
)
draw_figure_5()
grDevices::dev.off()

# ---- 5. Regression plane used in Figure 6 ------------------------------------

regression_model <- stats::lm(
  moisture ~ fat + protein,
  data = response_data
)

capture.output(
  summary(regression_model),
  file = file.path(output_dir, "Figure6_regression_model_summary.txt")
)

fat_sequence <- seq(
  min(response_data$fat),
  max(response_data$fat),
  length.out = 30
)

protein_sequence <- seq(
  min(response_data$protein),
  max(response_data$protein),
  length.out = 30
)

prediction_grid <- expand.grid(
  fat = fat_sequence,
  protein = protein_sequence
)

prediction_grid$moisture <- stats::predict(
  regression_model,
  newdata = prediction_grid
)

# For plotly surfaces, rows correspond to protein values and columns correspond
# to fat values. expand.grid() varies fat first, so byrow = TRUE is required.
regression_surface <- matrix(
  prediction_grid$moisture,
  nrow = length(protein_sequence),
  ncol = length(fat_sequence),
  byrow = TRUE
)

# ---- 6. Generate Figure 6 using plotly ---------------------------------------

figure_6a <- create_figure_6_plot(
  data = response_data,
  cluster_variable = "cluster_three_variables",
  plot_title = "Figure 6A: HC based on protein, moisture, and fat"
)

figure_6b <- create_figure_6_plot(
  data = response_data,
  cluster_variable = "cluster_two_variables",
  plot_title = "Figure 6B: HC based on protein and fat"
)

htmlwidgets::saveWidget(
  figure_6a,
  file = file.path(output_dir, "Figure6A_three_variable_clustering.html"),
  selfcontained = TRUE
)

htmlwidgets::saveWidget(
  figure_6b,
  file = file.path(output_dir, "Figure6B_two_variable_clustering.html"),
  selfcontained = TRUE
)

# Display both interactive plots when the script is run in an interactive R
# session. The HTML files are still generated when the script is sourced.
if (interactive()) {
  print(figure_6a)
  print(figure_6b)
}

# ---- 7. Record the R environment ---------------------------------------------

capture.output(
  sessionInfo(),
  file = file.path(output_dir, "sessionInfo.txt")
)

message(
  "Figures 5 and 6 were generated successfully. ",
  "See the '", output_dir, "' directory."
)
