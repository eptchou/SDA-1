# ==============================================================================
# Figure 3: Tecator absorbance curves colored by fat and protein categories
#
# Panel A: Blue = fat > 20%; orange = fat <= 20%
# Panel B: Blue = protein > 16%; orange = protein <= 16%
#
# The script creates one combined interactive Plotly figure and two separate
# interactive HTML files.
# ==============================================================================

# ---- 0. Package check ---------------------------------------------------------

required_packages <- c(
  "fsemipar",
  "plotly",
  "htmlwidgets"
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
  library(fsemipar)
  library(plotly)
  library(htmlwidgets)
})

# ---- 1. Output directory and data --------------------------------------------

output_dir <- "outputs"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

data("Tecator")

if (!exists("Tecator")) {
  stop("The Tecator dataset could not be loaded.", call. = FALSE)
}

spect_data <- as.data.frame(Tecator$absor.spectra)
fat_data <- Tecator$fat
protein_data <- Tecator$protein

if (
  nrow(spect_data) != length(fat_data) ||
  nrow(spect_data) != length(protein_data)
) {
  stop(
    "The number of spectra does not match the number of response observations.",
    call. = FALSE
  )
}

# The wavelength sequence is retained from the original analysis code.
x_vals <- seq(852, 1050, by = 2)

if (length(x_vals) != ncol(spect_data)) {
  stop(
    paste0(
      "The wavelength sequence contains ", length(x_vals),
      " values, but the spectral matrix contains ", ncol(spect_data),
      " columns. Verify the wavelength definition."
    ),
    call. = FALSE
  )
}

colnames(spect_data) <- as.character(x_vals)

# ---- 2. Plotting function -----------------------------------------------------

create_spectral_plot <- function(
    spectra,
    grouping_variable,
    threshold,
    variable_name,
    panel_label
) {
  fig <- plot_ly()

  group_high <- grouping_variable > threshold

  for (i in seq_len(nrow(spectra))) {
    category_label <- if (group_high[i]) {
      paste0(variable_name, " > ", threshold)
    } else {
      paste0(variable_name, " <= ", threshold)
    }

    line_color <- if (group_high[i]) "blue" else "orange"

    # Show one legend entry per category rather than one entry per sample.
    show_category_legend <- if (group_high[i]) {
      i == which(group_high)[1]
    } else {
      i == which(!group_high)[1]
    }

    fig <- add_trace(
      fig,
      x = x_vals,
      y = as.numeric(spectra[i, ]),
      type = "scatter",
      mode = "lines",
      name = category_label,
      legendgroup = category_label,
      showlegend = show_category_legend,
      line = list(
        color = line_color,
        width = 1
      ),
      text = paste0(
        "Sample: ", i,
        "<br>", variable_name, ": ",
        round(grouping_variable[i], 3)
      ),
      hovertemplate = paste0(
        "%{text}",
        "<br>Wavelength: %{x} nm",
        "<br>Absorbance: %{y:.3f}",
        "<extra></extra>"
      )
    )
  }

  layout(
    fig,
    title = paste0(
      panel_label,
      " Absorbance curves colored by ",
      tolower(variable_name),
      " content"
    ),
    xaxis = list(title = "Wavelength (nm)"),
    yaxis = list(title = "Absorbance"),
    legend = list(title = list(text = variable_name))
  )
}

# ---- 3. Generate Figure 3 panels ---------------------------------------------

figure_3a <- create_spectral_plot(
  spectra = spect_data,
  grouping_variable = fat_data,
  threshold = 20,
  variable_name = "Fat",
  panel_label = "(A)"
)

figure_3b <- create_spectral_plot(
  spectra = spect_data,
  grouping_variable = protein_data,
  threshold = 16,
  variable_name = "Protein",
  panel_label = "(B)"
)

# ---- 4. Combined Figure 3 -----------------------------------------------------

figure_3 <- subplot(
  figure_3a,
  figure_3b,
  nrows = 1,
  shareX = TRUE,
  shareY = TRUE,
  titleX = TRUE,
  titleY = TRUE,
  margin = 0.06
)

figure_3 <- layout(
  figure_3,
  title = "Figure 3: Tecator absorbance curves",
  legend = list(orientation = "h")
)

# Display the combined interactive figure when run interactively.
figure_3

# ---- 5. Save outputs ----------------------------------------------------------

htmlwidgets::saveWidget(
  figure_3,
  file = file.path(output_dir, "Figure3_Tecator_combined.html"),
  selfcontained = TRUE
)

htmlwidgets::saveWidget(
  figure_3a,
  file = file.path(output_dir, "Figure3A_fat_content.html"),
  selfcontained = TRUE
)

htmlwidgets::saveWidget(
  figure_3b,
  file = file.path(output_dir, "Figure3B_protein_content.html"),
  selfcontained = TRUE
)

# Save the categorization used for each sample.
figure_3_categories <- data.frame(
  sample_id = seq_len(nrow(spect_data)),
  fat = fat_data,
  fat_category = ifelse(fat_data > 20, "Fat > 20", "Fat <= 20"),
  protein = protein_data,
  protein_category = ifelse(
    protein_data > 16,
    "Protein > 16",
    "Protein <= 16"
  )
)

utils::write.csv(
  figure_3_categories,
  file = file.path(output_dir, "Figure3_sample_categories.csv"),
  row.names = FALSE
)

capture.output(
  sessionInfo(),
  file = file.path(output_dir, "sessionInfo.txt")
)

message(
  "Figure 3 was generated successfully. ",
  "See the '", output_dir, "' directory."
)
