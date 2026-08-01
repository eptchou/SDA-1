# =============================================================================
# Sparse Dependence Analysis (SDA-1) -- figure reproduction script
#
# Reproduces the figures in the manuscript from the Tecator (fat/protein) and
# MIR fruit-puree (strawberry) spectral data sets.
#
# Repository : https://github.com/eptchou/SDA-1
#
# HOW TO RUN
#   1. Open SDA-1.Rproj (or setwd() to the repository root) so that the
#      relative paths below resolve correctly.
#   2. Install the required packages (see the list under DEPENDENCIES).
#   3. Run a section between two "====" separators to reproduce a given figure.
#      Each section is self-contained and reloads the data it needs.
#
# DATA (in data/)
#   PMF_1st_cat.csv                fat spectra, 1st-difference categorisation
#   fatspectrum_2st_16cluster.csv  fat spectra, 2nd-difference 16-cluster table
#   MIR_Fruit_purees.csv           MIR strawberry / non-strawberry spectra
#   strawberry_1stc_16_n_983.csv   strawberry 1st-difference categorisation
#   strawberry_2ndc_16_n_983.csv   strawberry 2nd-difference categorisation
#   The Tecator fat spectra are loaded from the fds / fsemipar packages.
#
# OUTPUT
#   output/ holds cached intermediate results (e.g.
#   dynamic_2nd_split_results.xlsx), regenerated automatically if absent.
#
# DEPENDENCIES
#   CRAN         : tidyverse, dplyr, tidyr, ggplot2, magrittr, readxl, writexl,
#                  patchwork, cowplot, dendextend, circlize, RColorBrewer,
#                  ggrepel, infotheo, plotly, fds, fsemipar
#   Bioconductor : ComplexHeatmap
#
#   One-shot install:
#     install.packages(c("tidyverse","dplyr","tidyr","ggplot2","magrittr",
#                        "readxl","writexl","patchwork","cowplot","dendextend",
#                        "circlize","RColorBrewer","ggrepel","infotheo","plotly",
#                        "fds","fsemipar"))
#     install.packages("BiocManager")
#     BiocManager::install("ComplexHeatmap")
# =============================================================================


# ==========================================
# functional-curves                        figure3
# ==========================================
library(dplyr)
library(tidyr)
library(ggplot2)
library(fsemipar)
library(patchwork) 
data(Tecator)
spect_data <- as.data.frame(Tecator$absor.spectra)  
fat_data <- Tecator$fat
protein_data <- Tecator$protein
moisture_data <- Tecator$moisture

# Build the wavelength axis
x_vals <- seq(852, 1050, by = 2)
colnames(spect_data) <- as.character(x_vals)

plot_data <- spect_data %>%
  mutate(Sample = paste0("Sample_", dplyr::row_number()),
         Fat = fat_data,
         Protein = protein_data) %>%
  pivot_longer(
    cols = all_of(as.character(x_vals)),
    names_to = "Wavelength",
    values_to = "Absorbance"
  ) %>%
  mutate(Wavelength = as.numeric(Wavelength))

# First plot: grouped by Fat
p1 <- ggplot(plot_data %>% 
               mutate(ColorGroup = ifelse(Fat > 20, "High Fat", "Low Fat")),
             aes(x = Wavelength, y = Absorbance,
                 group = Sample, color = ColorGroup)) +
  geom_line(linewidth = 0.30) +
  labs(title = "Absorption Curves by Fat",
       x = "Wavelength (nm)", y = "Absorbance") +
  scale_color_manual(values = c("High Fat" = "blue", "Low Fat" = "orange")) +
  theme_minimal() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))

# Second plot: grouped by Protein
p2 <- ggplot(plot_data %>% 
               mutate(ColorGroup = ifelse(Protein > 16, "High Protein", "Low Protein")),
             aes(x = Wavelength, y = Absorbance,
                 group = Sample, color = ColorGroup)) +
  geom_line(linewidth = 0.30) +
  labs(title = "Absorption Curves by Protein",
       x = "Wavelength (nm)", y = "Absorbance") +
  scale_color_manual(values = c("High Protein" = "blue", "Low Protein" = "orange")) +
  theme_minimal() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))

library(cowplot)
combined <- plot_grid(p1, p2, labels = c("(A)", "(B)"), ncol = 2)
combined

print(combined) 
# ==========================================
# PFM-Categorization-1st-diff               figure 4
# ==========================================
library(readxl)
library(tidyverse)
library(magrittr)
library(fds)
Fatspectrum$x
test <- as.data.frame(Fatspectrum$y)
A = test
D <- matrix(0, nrow = nrow(A) - 1, ncol = ncol(A))

# Calculate 1st-order differences
for (i in 1:ncol(A)) {
  for (j in 2:nrow(A)) {
    D[j - 1, i] <- A[j, i] - A[j - 1, i]
  }
}
D %<>% as.data.frame()
df_long <- D %>%
  pivot_longer(cols = everything(), names_to = "fstd", values_to = "value")

data <- df_long[,2] |> as.data.frame()
h = hclust(dist(data),method = "ward.D2")



plot(h, hang = -1, labels = FALSE, ann = FALSE)
title(main = "Cluster Dendrogram")

clus =  16 #12
c = cutree(h, clus)

GetMinMax <- function(i){
  summ <- summary( data[which(c==i),1] )
  return( summ[c(1,6)] ) }

MinMax <- sapply( 1:clus, GetMinMax ) 
print("double check if the order matched")
( SmalltoLarge <- order( MinMax[1,] ) )
order( MinMax[2,] ) 

temp <- MinMax[1,][ SmalltoLarge ]
temp2 <- MinMax[2,][ SmalltoLarge ]

n <- nrow(data)
cum <- cumsum( table(c)[ SmalltoLarge ] ) / n
freq <- unname( table(c)[ SmalltoLarge ] )



plot( ecdf(data[,1]), col="grey", cex=0.2 ,main="",xlab="")
segments( temp[1], 0, temp2[1], cum[1], col="blue",lwd=3) 
for( i in 1:(clus-1) ){
  segments( temp[(i+1)], cum[i], temp2[(i+1)], cum[(i+1)], col="blue",lwd=3 ) }


plot(NULL, xlim=c(min(data[,1]),max(data[,1])), 
     ylim=c(0, max(freq)), ylab="freq.", xlab="value")

for(i in 1:clus){
  rect( temp[i], 0, temp2[i], freq[i])
}

# ==========================================
# Two-HC-trees-PFM (Section Data & Prep) figure 5
# ==========================================
library(magrittr)
library(fds)
library(fsemipar)
data(Tecator)
data_y <- data.frame(fat = Tecator$fat, protein = Tecator$protein, moisture = Tecator$moisture)
data_y

dist_matrix <- dist(data_y)

hc1 <- hclust(dist_matrix, method = 'ward.D2')

dendrogram_obj1 <- hc1 %>%
  as.dendrogram


data_y <- data.frame(fat = Tecator$fat, protein = Tecator$protein)
data_y

dist_matrix <- dist(data_y)

hc2 <- hclust(dist_matrix, method = 'ward.D2')

#y
dendrogram_obj2 <- hc2 %>%
  as.dendrogram


# Two-HC-trees-PFM
par(mfrow = c(1, 2), mar = c(2, 2, 2, 2))

dendrogram_obj1 %>%
  dendextend::color_branches(.,
                             k = 4,,
                             groupLabels = TRUE,
                             col = c("#636EFA", "#EF553B", "#00CC96", "#AB63FA")) %>%
  plot(leaflab = "none")
mtext("(A)", side = 3, adj = 0, line = 0.5, font = 2, cex = 1.8)

dendrogram_obj2 %>%
  dendextend::color_branches(.,
                             k = 4,
                             groupLabels = TRUE,
                             col = c("#FFA15A", "#19D3F3", "#FF6692", "#B6E880")) %>%
  plot(leaflab = "none")
mtext("(B)", side = 3, adj = 0, line = 0.5, font = 2, cex = 1.8)


# ==========================================
# motif-selection                           figure 7
# ==========================================

# ==========================================
# co-randomness-z1-z2 (Section Data & Prep) figure 8
# ==========================================
library(fds)
test <- as.data.frame(Fatspectrum$y)
copy1 <- read.csv("data/PMF_1st_cat.csv")
copy1$wavelenght <- seq(854,1050,2)

data <- data.frame(
  X1 = as.numeric(copy1[copy1$wavelenght == 928,][-216]),
  X2 = as.numeric(copy1[copy1$wavelenght == 930,][-216]),
  X3 = as.numeric(copy1[copy1$wavelenght == 932,][-216]),
  X4 = as.numeric(copy1[copy1$wavelenght == 934,][-216]),
  X5 = as.numeric(copy1[copy1$wavelenght == 946,][-216]),
  X6 = as.numeric(copy1[copy1$wavelenght == 948,][-216]),
  X7 = as.numeric(copy1[copy1$wavelenght == 950,][-216]),
  X8 = as.numeric(copy1[copy1$wavelenght == 952,][-216])
)

new_data <- as.data.frame(t(test))

hc_tree1 <- hclust(dist(data[, c("X1", "X2", "X3", "X4")]))

new_data$Z1 <- cutree(hc_tree1, k = 4)

hc_tree2 <- hclust(dist(data[, c("X5", "X6", "X7", "X8")]))

new_data$Z2 <- cutree(hc_tree2, k = 4)

contingency_table <- table(new_data$Z1, new_data$Z2)

contingency_table

library(dplyr)
library(tidyr)
library(ggplot2)

# Assume the first 100 columns of new_data are spectra, followed by Z1, Z2 (and possibly other columns)
x_vals <- seq(852, 1050, by = 2)
stopifnot(length(x_vals) == 100)

# Set spectrum column names
colnames(new_data)[1:100] <- as.character(x_vals)

# Reshape to long format, keeping Z1, Z2
long_df <- new_data %>%
  mutate(Sample = paste0("S", row_number())) %>%
  pivot_longer(
    cols = all_of(as.character(x_vals)),
    names_to   = "Wavelength",
    values_to  = "Absorbance"
  ) %>%
  mutate(
    Wavelength = as.numeric(Wavelength),
    Z1 = as.factor(Z1),
    Z2 = as.factor(Z2)
  )

# Plot only the 16 cells for Z1 = 1..4 and Z2 = 1..4
plot_df <- long_df %>% filter(Z1 %in% levels(Z1)[1:4], Z2 %in% levels(Z2)[1:4])

p <- ggplot(
  plot_df,
  aes(x = Wavelength, y = Absorbance, group = Sample)
) +
  geom_line(linewidth = 0.25, color = "lightskyblue3") +
  facet_grid(rows = vars(Z1), cols = vars(Z2)) +   # 4×4
  labs(
    x = "Wavelength (nm)",
    y = "Absorbance"
  ) +
  coord_cartesian(xlim = c(852, 1050)) +           # consistent x-limit across all panels
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "none",
    panel.spacing = unit(4, "pt"),
    strip.background = element_rect(fill = "grey95", color = NA),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

print(p) 


# ==========================================
# PFM-z1-tree-result (Section Data & Prep) figure 9
# ==========================================
library(magrittr)
library(fds)
library(fsemipar)
data(Tecator)
data_y <- data.frame(fat = Tecator$fat, protein = Tecator$protein)
data_y

test <- as.data.frame(Fatspectrum$y)
new_data <- as.data.frame(t(test))

dist_matrix <- dist(data_y)

hc <- hclust(dist_matrix, method = 'ward.D2')

dendrogram_obj <- hc %>%
  as.dendrogram

dendrogram_obj %>%
  dendextend::color_branches(.,
                             k = 4,,
                             groupLabels = TRUE,
                             col = c("#FFA15A", "#19D3F3", "#FF6692", "#B6E880")) %>%
  plot

dendroextend_cut_tree <- dendextend::cutree(
  dendrogram_obj,
  k = 4,
  order_clusters_as_data = FALSE
)

data_order <- order.dendrogram(dendrogram_obj)

new_data$cluster_y_4 <- NA

new_data[data_order,]$cluster_y_4 <- dendroextend_cut_tree

table(new_data$cluster_y_4)


second <- read.csv("data/fatspectrum_2st_16cluster.csv", row.names = 1)

second$wavelenght <- seq(854,1048,2)

data <- data.frame(
  X1 = as.numeric(second[second$wavelenght == 926,][-216]),
  X2 = as.numeric(second[second$wavelenght == 928,][-216]),
  X3 = as.numeric(second[second$wavelenght == 930,][-216]),
  X4 = as.numeric(second[second$wavelenght == 932,][-216]),
  X5 = as.numeric(second[second$wavelenght == 934,][-216]),
  X6 = as.numeric(second[second$wavelenght == 936,][-216]),
  X7 = as.numeric(second[second$wavelenght == 944,][-216]),
  X8 = as.numeric(second[second$wavelenght == 946,][-216]),
  X9 = as.numeric(second[second$wavelenght == 948,][-216]),
  X10 = as.numeric(second[second$wavelenght == 950,][-216]),
  X11 = as.numeric(second[second$wavelenght == 952,][-216])
)

hc_tree1 <- hclust(dist(data[, c("X1", "X2", "X3", "X4", "X5", "X6")]), method = "ward.D2")


hc_tree2 <- hclust(dist(data[, c("X7", "X8", "X9", "X10", "X11")]), method = "ward.D2")

#Z1
dendrogram_obj <- hc_tree1 %>%
  as.dendrogram

dendrogram_obj %>%
  dendextend::color_branches(.,
                             k = 9,,
                             groupLabels = TRUE) %>%
  plot(leaflab = "none")


dendroextend_cut_tree <- dendextend::cutree(
  dendrogram_obj,
  k = 9,
  order_clusters_as_data = FALSE
)

data_order <- order.dendrogram(dendrogram_obj)

new_data$Z1 <- NA

new_data[data_order,]$Z1 <- dendroextend_cut_tree

table(Y = new_data$cluster_y_4, new_data$Z1)



entropy_overlap_by_column <- function(contingency_table, n_sim = 1000, prefix = "Zcat") {
  stopifnot(nrow(contingency_table) == 4)
  
  row_tot <- rowSums(contingency_table)
  col_tot <- colSums(contingency_table)
  grand_tot <- sum(contingency_table)
  
  row_prob_alt <- lapply(1:4, function(i) contingency_table[i, ] / row_tot[i])
  col_prob_null <- col_tot / grand_tot
  
  n_col <- ncol(contingency_table)
  alt_array <- array(0, dim = c(4, n_col, n_sim))
  null_array <- array(0, dim = c(4, n_col, n_sim))
  
  set.seed(42)
  for (s in 1:n_sim) {
    for (i in 1:4) {
      alt_array[i, , s] <- rmultinom(1, row_tot[i], row_prob_alt[[i]])
      null_array[i, , s] <- rmultinom(1, row_tot[i], col_prob_null)
    }
  }
  
  entropy <- function(p) {
    p <- p[p > 0]
    -sum(p * log(p))
  }
  
  breaks_fix <- seq(0, log(4), length.out = 50)
  result_df <- data.frame(Category = 1:n_col, Overlap = NA_real_)
  
  for (k in 1:n_col) {
    alt_ent  <- apply(alt_array[ , k, ], 2, function(v) entropy(v / sum(v)))
    null_ent <- apply(null_array[, k, ], 2, function(v) entropy(v / sum(v)))
    
    alt_hist  <- hist(alt_ent,  breaks = breaks_fix, plot = FALSE)
    null_hist <- hist(null_ent, breaks = breaks_fix, plot = FALSE)
    overlap_counts <- pmin(alt_hist$counts, null_hist$counts)
    overlap_ratio <- sum(overlap_counts) / n_sim
    result_df$Overlap[k] <- overlap_ratio
    
    #pdf(paste0(prefix, "_", k, "_entropy.pdf"), width = 6, height = 4)
    hist(null_ent, breaks = breaks_fix, col = rgb(0,0,1,0.5),
         main = paste0(prefix, " ", k, " | overlap=", round(overlap_ratio, 3)),
         xlab = "Entropy", ylim = c(0, max(alt_hist$counts, null_hist$counts)))
    hist(alt_ent, breaks = breaks_fix, col = rgb(1,0,0,0.5), add = TRUE)
    legend("topright", legend = c("Null","Alternative"), fill = c("blue","red"))
    #dev.off()
  }
  
  return(result_df)
}
# Z1: values below come from table(Y = new_data$cluster_y_4, new_data$Z1)
Z1_table <- matrix(c(
  0,  0,  0,  0,  14,  0,  0,  6,  0,   # Y = 1
  0,  0,  0,  0,  0, 11,  3, 13, 18,   # Y = 2
  0, 43, 55,  0,  0,  0,  0,  0, 0,   # Y = 3
  24,  6,  0, 10,  0,  3,  8,  1, 0    # Y = 4
), nrow = 4, byrow = TRUE)

rownames(Z1_table) <- paste0("Y", 1:4)
colnames(Z1_table) <- paste0("Z1_", 1:9)

result_z1 <- entropy_overlap_by_column(Z1_table, prefix = "Z1_")
print(result_z1)

# ==========================================
# PFM-Z1-678-split-by-Z24 (Section Data & Prep) figure 10
# ==========================================
library(magrittr)
library(fds)
library(fsemipar)
data(Tecator)
data_y <- data.frame(fat = Tecator$fat, protein = Tecator$protein)
data_y

test <- as.data.frame(Fatspectrum$y)
new_data <- as.data.frame(t(test))

dist_matrix <- dist(data_y)

hc <- hclust(dist_matrix, method = 'ward.D2')

#y
dendrogram_obj <- hc %>%
  as.dendrogram

dendroextend_cut_tree <- dendextend::cutree(
  dendrogram_obj,
  k = 4,
  order_clusters_as_data = FALSE
)

data_order <- order.dendrogram(dendrogram_obj)

new_data$Y <- NA

new_data[data_order,]$Y <- dendroextend_cut_tree


second <- read.csv("data/fatspectrum_2st_16cluster.csv", row.names = 1)

second$wavelenght <- seq(854,1048,2)

data <- data.frame(
  X1 = as.numeric(second[second$wavelenght == 926,][-216]),
  X2 = as.numeric(second[second$wavelenght == 928,][-216]),
  X3 = as.numeric(second[second$wavelenght == 930,][-216]),
  X4 = as.numeric(second[second$wavelenght == 932,][-216]),
  X5 = as.numeric(second[second$wavelenght == 934,][-216]),
  X6 = as.numeric(second[second$wavelenght == 936,][-216]),
  X7 = as.numeric(second[second$wavelenght == 944,][-216]),
  X8 = as.numeric(second[second$wavelenght == 946,][-216]),
  X9 = as.numeric(second[second$wavelenght == 948,][-216]),
  X10 = as.numeric(second[second$wavelenght == 950,][-216]),
  X11 = as.numeric(second[second$wavelenght == 952,][-216])
)

hc_tree1 <- hclust(dist(data[, c("X1", "X2", "X3", "X4", "X5", "X6")]), method = "ward.D2")


hc_tree2 <- hclust(dist(data[, c("X7", "X8", "X9", "X10", "X11")]), method = "ward.D2")

#Z1
dendrogram_obj <- hc_tree1 %>%
  as.dendrogram

dendroextend_cut_tree <- dendextend::cutree(
  dendrogram_obj,
  k = 9,
  order_clusters_as_data = FALSE
)

data_order <- order.dendrogram(dendrogram_obj)

new_data$Z1 <- NA

new_data[data_order,]$Z1 <- dendroextend_cut_tree

table(new_data$Z1)

table(Y = new_data$Y, Z1 = new_data$Z1)

#Z2
dendrogram_obj <- hc_tree2 %>%
  as.dendrogram

dendroextend_cut_tree <- dendextend::cutree(
  dendrogram_obj,
  k = 10,
  order_clusters_as_data = FALSE
)

data_order <- order.dendrogram(dendrogram_obj)

new_data$Z2 <- NA

new_data[data_order,]$Z2 <- dendroextend_cut_tree
# cluster 9 has only 1 member
new_data$Z2[new_data$Z2 == 9] <- 10

table(new_data$Z2)

#Z4
data2 <- data.frame(
  X1 = as.numeric(second[second$wavelenght == 902,][-216]),
  X2 = as.numeric(second[second$wavelenght == 904,][-216]),
  X3 = as.numeric(second[second$wavelenght == 906,][-216]),
  X4 = as.numeric(second[second$wavelenght == 908,][-216]),
  X5 = as.numeric(second[second$wavelenght == 910,][-216]),
  X6 = as.numeric(second[second$wavelenght == 1040,][-216]),
  X7 = as.numeric(second[second$wavelenght == 1042,][-216]),
  X8 = as.numeric(second[second$wavelenght == 1044,][-216]),
  X9 = as.numeric(second[second$wavelenght == 1046,][-216])
)

hc_tree4 <- hclust(dist(data[, c("X6", "X7", "X8", "X9")]), method = "ward.D2")

dendrogram_obj <- hc_tree4 %>%
  as.dendrogram

dendrogram_obj %>%
  dendextend::color_branches(.,
                             k = 9,,
                             groupLabels = TRUE) %>%
  plot


dendroextend_cut_tree <- dendextend::cutree(
  dendrogram_obj,
  k = 9,
  order_clusters_as_data = FALSE
)

data_order <- order.dendrogram(dendrogram_obj)

new_data$Z4 <- NA

new_data[data_order,]$Z4 <- dendroextend_cut_tree

table(new_data$Z4)

table(Y = new_data[new_data$Z1==6,]$Y, Z2 = new_data[new_data$Z1==6,]$Z2)
table(Y =new_data[new_data$Z1==7,]$Y, Z2 = new_data[new_data$Z1==7,]$Z2)
table(Y =new_data[new_data$Z1==8,]$Y, Z2 = new_data[new_data$Z1==8,]$Z2)

table(Y = new_data[new_data$Z1==6,]$Y, Z4 = new_data[new_data$Z1==6,]$Z4)
table(Y =new_data[new_data$Z1==7,]$Y, Z4 = new_data[new_data$Z1==7,]$Z4)
table(Y =new_data[new_data$Z1==8,]$Y, Z4 = new_data[new_data$Z1==8,]$Z4)


# ========================================== 
# PFM-Figure1-Y1Y2 (Section Data & Prep)    figure 11
# ==========================================
library(ggplot2)
library(dplyr)
library(tidyr)
library(magrittr)
library(fds)
library(fsemipar)
data(Tecator)
data_y <- data.frame(fat = Tecator$fat, protein = Tecator$protein)
data_y

test <- as.data.frame(Fatspectrum$y)
new_data <- as.data.frame(t(test))

dist_matrix <- dist(data_y)

hc <- hclust(dist_matrix, method = 'ward.D2')

#y
dendrogram_obj <- hc %>%
  as.dendrogram

dendrogram_obj %>%
  dendextend::color_branches(.,
                             k = 4,
                             groupLabels = TRUE,
                             col = c("#FFA15A", "#19D3F3", "#FF6692", "#B6E880")) %>%
  plot()


dendroextend_cut_tree <- dendextend::cutree(
  dendrogram_obj,
  k = 4,
  order_clusters_as_data = FALSE
)

data_order <- order.dendrogram(dendrogram_obj)

new_data$Y <- NA

new_data[data_order,]$Y <- dendroextend_cut_tree

second <- read.csv("data/fatspectrum_2st_16cluster.csv", row.names = 1)

second$wavelenght <- seq(854,1048,2)

data <- data.frame(
  X1 = as.numeric(second[second$wavelenght == 926,][-216]),
  X2 = as.numeric(second[second$wavelenght == 928,][-216]),
  X3 = as.numeric(second[second$wavelenght == 930,][-216]),
  X4 = as.numeric(second[second$wavelenght == 932,][-216]),
  X5 = as.numeric(second[second$wavelenght == 934,][-216]),
  X6 = as.numeric(second[second$wavelenght == 936,][-216]),
  X7 = as.numeric(second[second$wavelenght == 944,][-216]),
  X8 = as.numeric(second[second$wavelenght == 946,][-216]),
  X9 = as.numeric(second[second$wavelenght == 948,][-216]),
  X10 = as.numeric(second[second$wavelenght == 950,][-216]),
  X11 = as.numeric(second[second$wavelenght == 952,][-216])
)

hc_tree1 <- hclust(dist(data[, c("X1", "X2", "X3", "X4", "X5", "X6")]), method = "ward.D2")

#Z1
dendrogram_obj <- hc_tree1 %>%
  as.dendrogram

dendrogram_obj %>%
  dendextend::color_branches(.,
                             k = 9,,
                             groupLabels = TRUE) %>%
  plot


dendroextend_cut_tree <- dendextend::cutree(
  dendrogram_obj,
  k = 9,
  order_clusters_as_data = FALSE
)

data_order <- order.dendrogram(dendrogram_obj)

new_data$Z1 <- NA

new_data[data_order,]$Z1 <- dendroextend_cut_tree

table(new_data$Z1)

table(Y = new_data$Y, Z1 = new_data$Z1)


# Four combinations to plot in the first figure
grid1 <- data.frame(
  Y = c(1, 1, 2, 2),
  Z1 = c(5, 8, 9, 8)
)

# Four combinations to plot in the second figure
grid2 <- data.frame(
  Y = c(3, 3, 4, 4),
  Z1 = c(2, 3, 1, 4)
)

plot_motif_grid <- function(new_data, grid_df, output_pdf_name, motif_ranges = NULL) {
  x_vals <- seq(852, 1050, by = 2)
  
  plot_data <- data.frame()
  
  for (k in seq_len(nrow(grid_df))) {
    y_val <- grid_df$Y[k]
    z1_val <- grid_df$Z1[k]
    
    subset_df <- new_data[new_data$Y == y_val & new_data$Z1 == z1_val, 1:100]
    
    if (nrow(subset_df) > 0) {
      long_df <- subset_df %>%
        mutate(SampleID = paste0("S", 1:n())) %>%
        mutate(Y = y_val, Z1 = z1_val) %>%
        pivot_longer(
          cols = 1:100,
          names_to = "Wavelength_Index",
          values_to = "Absorbance"
        ) %>%
        mutate(Wavelength = rep(seq(852, 1050, by = 2), times = nrow(subset_df))) %>%
        dplyr::select(Wavelength, Absorbance, SampleID, Y, Z1)
      
      plot_data <- bind_rows(plot_data, long_df)
    }
  }
  
  plot_data$FacetLabel <- paste0("Y=", plot_data$Y, ", Z1=", plot_data$Z1)
  
  p <- ggplot(plot_data, aes(x = Wavelength, y = Absorbance, group = SampleID)) +
    geom_line(color = "lightgreen", alpha = 0.5) +
    facet_wrap(~FacetLabel, nrow = 2) +
    theme_bw() +
    labs(x = "Wavelength", y = "Absorbance") +
    theme(strip.text = element_text(size = 12))
  
  if (!is.null(motif_ranges)) {
    for (i in 1:nrow(motif_ranges)) {
      p <- p + annotate(
        "rect",
        xmin = motif_ranges$xmin[i],
        xmax = motif_ranges$xmax[i],
        ymin = -Inf,
        ymax = Inf,
        alpha = 0.2,
        fill = "orange"
      ) +
        annotate(
          "text",
          x = (motif_ranges$xmin[i] + motif_ranges$xmax[i]) / 2,
          y = 6.2,
          label = motif_ranges$label[i],
          size = 4.5,
          fontface = "bold",
          color = "black"
        )
    }
  }
  
  #ggsave(output_pdf_name, plot = p, width = 10, height = 6)
  print(p)
}

motif_ranges <- data.frame(
  xmin = c(926, 944, 902, 1040),
  xmax = c(936, 952, 910, 1046),
  label = c("Z1", "Z2", "Z3", "Z4")
)


plot_motif_grid(new_data, grid1, "Figure1_Y1Y2.pdf", motif_ranges)
plot_motif_grid(new_data, grid2, "Figure2_Y3Y4.pdf", motif_ranges)
# ========================================== 
# PFM-Z3-tree-result (Section Data & Prep) figure 12
# ==========================================
library(magrittr)
library(fds)
library(fsemipar)
data(Tecator)
data_y <- data.frame(fat = Tecator$fat, protein = Tecator$protein)
data_y

test <- as.data.frame(Fatspectrum$y)
new_data <- as.data.frame(t(test))

dist_matrix <- dist(data_y)

hc <- hclust(dist_matrix, method = 'ward.D2')

dendrogram_obj <- hc %>%
  as.dendrogram

dendroextend_cut_tree <- dendextend::cutree(
  dendrogram_obj,
  k = 4,
  order_clusters_as_data = FALSE
)

data_order <- order.dendrogram(dendrogram_obj)

new_data$cluster_y_4 <- NA

new_data[data_order,]$cluster_y_4 <- dendroextend_cut_tree
second <- read.csv("data/fatspectrum_2st_16cluster.csv", row.names = 1)

second$wavelenght <- seq(854,1048,2)
#Z3
data <- data.frame(
  X1 = as.numeric(second[second$wavelenght == 902,][-216]),
  X2 = as.numeric(second[second$wavelenght == 904,][-216]),
  X3 = as.numeric(second[second$wavelenght == 906,][-216]),
  X4 = as.numeric(second[second$wavelenght == 908,][-216]),
  X5 = as.numeric(second[second$wavelenght == 910,][-216])
)
hc_tree3 <- hclust(dist(data[, c("X1", "X2", "X3", "X4", "X5")]), method = "ward.D2")

dendrogram_obj <- hc_tree3 %>%
  as.dendrogram


# PFM-Z3-tree-result-1
dendrogram_obj %>%
  dendextend::color_branches(.,
                             k = 9,,
                             groupLabels = TRUE) %>%
  plot(leaflab = "none")

dendroextend_cut_tree <- dendextend::cutree(
  dendrogram_obj,
  k = 9,
  order_clusters_as_data = FALSE
)

data_order <- order.dendrogram(dendrogram_obj)

new_data$Z3 <- NA

new_data[data_order,]$Z3 <- dendroextend_cut_tree

table(new_data$cluster_y_4, new_data$Z3)


data <- data.frame(
  X1 = as.numeric(second[second$wavelenght == 926,][-216]),
  X2 = as.numeric(second[second$wavelenght == 928,][-216]),
  X3 = as.numeric(second[second$wavelenght == 930,][-216]),
  X4 = as.numeric(second[second$wavelenght == 932,][-216]),
  X5 = as.numeric(second[second$wavelenght == 934,][-216]),
  X6 = as.numeric(second[second$wavelenght == 936,][-216]),
  X7 = as.numeric(second[second$wavelenght == 944,][-216]),
  X8 = as.numeric(second[second$wavelenght == 946,][-216]),
  X9 = as.numeric(second[second$wavelenght == 948,][-216]),
  X10 = as.numeric(second[second$wavelenght == 950,][-216]),
  X11 = as.numeric(second[second$wavelenght == 952,][-216])
)

hc_tree1 <- hclust(dist(data[, c("X1", "X2", "X3", "X4", "X5", "X6")]), method = "ward.D2")

#Z1
dendrogram_obj <- hc_tree1 %>%
  as.dendrogram

dendrogram_obj %>%
  dendextend::color_branches(.,
                             k = 9,,
                             groupLabels = TRUE) %>%
  plot


dendroextend_cut_tree <- dendextend::cutree(
  dendrogram_obj,
  k = 9,
  order_clusters_as_data = FALSE
)

data_order <- order.dendrogram(dendrogram_obj)

new_data$Z1 <- NA

new_data[data_order,]$Z1 <- dendroextend_cut_tree

table(new_data$Z1, new_data$Z3)




entropy_overlap_by_column <- function(contingency_table, n_sim = 1000, prefix = "Zcat") {
  stopifnot(nrow(contingency_table) == 4)
  
  row_tot <- rowSums(contingency_table)
  col_tot <- colSums(contingency_table)
  grand_tot <- sum(contingency_table)
  
  row_prob_alt <- lapply(1:4, function(i) contingency_table[i, ] / row_tot[i])
  col_prob_null <- col_tot / grand_tot
  
  n_col <- ncol(contingency_table)
  alt_array <- array(0, dim = c(4, n_col, n_sim))
  null_array <- array(0, dim = c(4, n_col, n_sim))
  
  set.seed(42)
  for (s in 1:n_sim) {
    for (i in 1:4) {
      alt_array[i, , s] <- rmultinom(1, row_tot[i], row_prob_alt[[i]])
      null_array[i, , s] <- rmultinom(1, row_tot[i], col_prob_null)
    }
  }
  
  entropy <- function(p) {
    p <- p[p > 0]
    -sum(p * log(p))
  }
  
  breaks_fix <- seq(0, log(4), length.out = 50)
  result_df <- data.frame(Category = 1:n_col, Overlap = NA_real_)
  
  for (k in 1:n_col) {
    alt_ent  <- apply(alt_array[ , k, ], 2, function(v) entropy(v / sum(v)))
    null_ent <- apply(null_array[, k, ], 2, function(v) entropy(v / sum(v)))
    
    alt_hist  <- hist(alt_ent,  breaks = breaks_fix, plot = FALSE)
    null_hist <- hist(null_ent, breaks = breaks_fix, plot = FALSE)
    overlap_counts <- pmin(alt_hist$counts, null_hist$counts)
    overlap_ratio <- sum(overlap_counts) / n_sim
    result_df$Overlap[k] <- overlap_ratio
    
    #pdf(paste0(prefix, "_", k, "_entropy.pdf"), width = 6, height = 4)
    hist(null_ent, breaks = breaks_fix, col = rgb(0,0,1,0.5),
         main = paste0(prefix, " ", k, " | overlap=", round(overlap_ratio, 3)),
         xlab = "Entropy", ylim = c(0, max(alt_hist$counts, null_hist$counts)))
    hist(alt_ent, breaks = breaks_fix, col = rgb(1,0,0,0.5), add = TRUE)
    legend("topright", legend = c("Null","Alternative"), fill = c("blue","red"))
    #dev.off()
  }
  
  return(result_df)
}
# Z3: values below come from table(Y = new_data$cluster_y_4, new_data$Z3)
Z3_table <- matrix(c(
  0,  0,  0,  5,  3,  8,  0,  3,  1,   # Y = 1
  0,  0,  0,  2,  0, 0,  2, 25, 16,   # Y = 2
  58, 12, 28,  0,  0,  0,  0,  0, 0,   # Y = 3
  0,  10,  0, 0,  3,  0,  32,  2, 5    # Y = 4
), nrow = 4, byrow = TRUE)

rownames(Z3_table) <- paste0("Y", 1:4)
colnames(Z3_table) <- paste0("Z3_", 1:9)

result_z3 <- entropy_overlap_by_column(Z3_table, prefix = "Z3_")
print(result_z3)

# ==========================================
# PFM-heatmap-Z1-2nd-Z3-2nd-revised         figure 13
# ==========================================
library(magrittr)
library(fds)
library(fsemipar)
data(Tecator)
data_y <- data.frame(fat = Tecator$fat, protein = Tecator$protein)
data_y

test <- as.data.frame(Fatspectrum$y)
new_data <- as.data.frame(t(test))

dist_matrix <- dist(data_y)

hc <- hclust(dist_matrix, method = 'ward.D2')

dendrogram_obj <- hc %>%
  as.dendrogram

dendrogram_obj %>%
  dendextend::color_branches(.,
                             k = 4,,
                             groupLabels = TRUE,
                             col = c("#FFA15A", "#19D3F3", "#FF6692", "#B6E880")) %>%
  plot

dendroextend_cut_tree <- dendextend::cutree(
  dendrogram_obj,
  k = 4,
  order_clusters_as_data = FALSE
)

data_order <- order.dendrogram(dendrogram_obj)

new_data$cluster_y_4 <- NA

new_data[data_order,]$cluster_y_4 <- dendroextend_cut_tree

table(new_data$cluster_y_4)


second <- read.csv("data/fatspectrum_2st_16cluster.csv", row.names = 1)

second$wavelenght <- seq(854,1048,2)

data <- data.frame(
  X1 = as.numeric(second[second$wavelenght == 926,][-216]),
  X2 = as.numeric(second[second$wavelenght == 928,][-216]),
  X3 = as.numeric(second[second$wavelenght == 930,][-216]),
  X4 = as.numeric(second[second$wavelenght == 932,][-216]),
  X5 = as.numeric(second[second$wavelenght == 934,][-216]),
  X6 = as.numeric(second[second$wavelenght == 936,][-216]),
  X7 = as.numeric(second[second$wavelenght == 944,][-216]),
  X8 = as.numeric(second[second$wavelenght == 946,][-216]),
  X9 = as.numeric(second[second$wavelenght == 948,][-216]),
  X10 = as.numeric(second[second$wavelenght == 950,][-216]),
  X11 = as.numeric(second[second$wavelenght == 952,][-216])
)

hc_tree1 <- hclust(dist(data[, c("X1", "X2", "X3", "X4", "X5", "X6")]), method = "ward.D2")


hc_tree2 <- hclust(dist(data[, c("X7", "X8", "X9", "X10", "X11")]), method = "ward.D2")

#Z1
dendrogram_obj <- hc_tree1 %>%
  as.dendrogram

dendrogram_obj %>%
  dendextend::color_branches(.,
                             k = 9,,
                             groupLabels = TRUE) %>%
  plot(leaflab = "none")


dendroextend_cut_tree <- dendextend::cutree(
  dendrogram_obj,
  k = 9,
  order_clusters_as_data = FALSE
)

data_order <- order.dendrogram(dendrogram_obj)

new_data$Z1 <- NA

new_data[data_order,]$Z1 <- dendroextend_cut_tree

table(Y = new_data$cluster_y_4, new_data$Z1)

#Z2
dendrogram_obj <- hc_tree2 %>%
  as.dendrogram

dendrogram_obj %>%
  dendextend::color_branches(.,
                             k = 10,,
                             groupLabels = TRUE) %>%
  plot


dendroextend_cut_tree <- dendextend::cutree(
  dendrogram_obj,
  k = 10,
  order_clusters_as_data = FALSE
)

data_order <- order.dendrogram(dendrogram_obj)

new_data$Z2 <- NA

new_data[data_order,]$Z2 <- dendroextend_cut_tree

new_data$cluster_x <- with(new_data, paste(Z1, Z2, sep = "_"))

table(Y = new_data$cluster_y_4, Z2 = new_data$Z2)

#Z3
data <- data.frame(
  X1 = as.numeric(second[second$wavelenght == 902,][-216]),
  X2 = as.numeric(second[second$wavelenght == 904,][-216]),
  X3 = as.numeric(second[second$wavelenght == 906,][-216]),
  X4 = as.numeric(second[second$wavelenght == 908,][-216]),
  X5 = as.numeric(second[second$wavelenght == 910,][-216])
)
hc_tree3 <- hclust(dist(data[, c("X1", "X2", "X3", "X4", "X5")]), method = "ward.D2")

dendrogram_obj <- hc_tree3 %>%
  as.dendrogram

dendrogram_obj %>%
  dendextend::color_branches(.,
                             k = 9,,
                             groupLabels = TRUE) %>%
  plot


dendroextend_cut_tree <- dendextend::cutree(
  dendrogram_obj,
  k = 9,
  order_clusters_as_data = FALSE
)

data_order <- order.dendrogram(dendrogram_obj)

new_data$Z3 <- NA

new_data[data_order,]$Z3 <- dendroextend_cut_tree

table(new_data$Z3)

table(Z1 = new_data$Z1, Z3 = new_data$Z3)

table(Y = new_data$cluster_y_4, Z3 = new_data$Z3)

#Z4
data <- data.frame(
  X6 = as.numeric(second[second$wavelenght == 1040,][-216]),
  X7 = as.numeric(second[second$wavelenght == 1042,][-216]),
  X8 = as.numeric(second[second$wavelenght == 1044,][-216]),
  X9 = as.numeric(second[second$wavelenght == 1046,][-216])
)

hc_tree4 <- hclust(dist(data[, c("X6", "X7", "X8", "X9")]), method = "ward.D2")

dendrogram_obj <- hc_tree4 %>%
  as.dendrogram

dendrogram_obj %>%
  dendextend::color_branches(.,
                             k = 9,,
                             groupLabels = TRUE) %>%
  plot


dendroextend_cut_tree <- dendextend::cutree(
  dendrogram_obj,
  k = 9,
  order_clusters_as_data = FALSE
)

data_order <- order.dendrogram(dendrogram_obj)

new_data$Z4 <- NA

new_data[data_order,]$Z4 <- dendroextend_cut_tree

table(Y = new_data$cluster_y_4, Z4 = new_data$Z4)

table(Z1 = new_data$Z1, Z2 = new_data$Z2)
table(Z1 = new_data$Z1, Z3 = new_data$Z3)
table(Z1 = new_data$Z1, Z4 = new_data$Z4)

table(Z2 = new_data$Z2, Z3 = new_data$Z3)
table(Z2 = new_data$Z2, Z4 = new_data$Z4)

table(Z3 = new_data$Z3, Z4 = new_data$Z4)



## binary heatmap
library(dplyr)
library(tidyr)

# Build one-hot encoding matrix
Z1_mat <- model.matrix(~ factor(Z1) - 1, data = new_data)
colnames(Z1_mat) <- paste0("Z1_", 1:9)

Z3_mat <- model.matrix(~ factor(Z3) - 1, data = new_data)
colnames(Z3_mat) <- paste0("Z3_", 1:9)

# Combine into a binary matrix
binary_mat <- cbind(Z1_mat, Z3_mat)

row_labels <- as.factor(new_data$cluster_y_4)

# Color palette for Y
selected_cols <- c("Z1_1", "Z1_2", "Z1_3", "Z1_4", "Z1_5", "Z1_6", "Z1_7", "Z1_9",
                   "Z3_1", "Z3_2", "Z3_3", "Z3_6", "Z3_7", "Z3_8", "Z3_9")

library(ComplexHeatmap)
library(circlize) 
library(dendextend) 

mat <- binary_mat[, selected_cols]
row_labels <- factor(new_data$cluster_y_4,  # Y = 1,2,3,4
                     levels = 1:4,
                     labels = 1:4)

y_colors_named <- c("#FFA15A", "#19D3F3", "#FF6692", "#B6E880")
names(y_colors_named) <- levels(row_labels)

row_anno <- rowAnnotation(
  Y = row_labels,
  col = list(Y = y_colors_named),
  show_annotation_name = FALSE
)

col_fun <- colorRamp2(c(0, 1), c("white", "black"))

ht <- Heatmap(
  mat,
  name           = " ",                 
  col            = col_fun,
  show_heatmap_legend   = FALSE,
  left_annotation       = row_anno,
  cluster_rows          = TRUE,
  cluster_columns       = TRUE,
  clustering_method_rows    = "ward.D2",
  clustering_method_columns = "ward.D2",
  show_row_names        = FALSE
)
ht

# ==========================================
# 4-panels-PFM-heatmap-Z1-2nd-Z3-2nd-revised figure 14
# ==========================================
ht_obj <- draw(ht)

row_dend <- row_dend(ht_obj)
library(dendextend)
row_dend %>%
  dendextend::color_branches(.,
                             k = 8,
                             groupLabels = TRUE) %>%
  plot

dendroextend_cut_tree <- dendextend::cutree(
  row_dend,
  k = 8,
  order_clusters_as_data = FALSE
)

data_order <- order.dendrogram(row_dend)

test <- as.data.frame(matrix(nrow=215,ncol=1))
test$block1 <- NA
test[data_order,]$block1 <- dendroextend_cut_tree

block1_index <- which(test$block1 == 1)
block2_index <- which(test$block1 == 5)
block3_index <- which(test$block1 == 6)
block4_index <- which(test$block1 == 7)


### Draw the 4 panels
library(ggplot2)
library(dplyr)
library(tidyr)

# block1_index ~ block4_index already exist
block_list <- list(block1_index, block2_index, block3_index, block4_index)
block_names <- paste0("Block ", LETTERS[1:4])

long_data <- data.frame()

new_data <- new_data[,-c(102:106)]

for (i in seq_along(block_list)) {
  idx <- block_list[[i]]
  block_data <- new_data[idx, ]
  block_data$SampleID <- paste0("S", idx)
  block_data$Block <- block_names[i]
  
  # Remove the y column, then convert to long format
  long_block <- block_data %>%
    mutate(cluster_y_4 = factor(cluster_y_4)) %>%
    pivot_longer(
      cols = -c(cluster_y_4, SampleID, Block),
      names_to = "X",
      values_to = "Y"
    )
  
  long_data <- bind_rows(long_data, long_block)
}

long_data$X <- as.numeric(as.character(long_data$X))

# Plot
ggplot(long_data, aes(x = X, y = Y, group = SampleID, color = cluster_y_4)) +
  geom_line(alpha = 0.7) +
  facet_wrap(~ Block, ncol = 2) +
  scale_color_manual(values = c("1" = "#FFA15A", "2" = "#19D3F3", "3" = "#FF6692", "4" = "#B6E880")) +
  labs(x = "Wavelength", y = "Absorbance", color = "Class") +
  theme_minimal(base_size = 14)



# Add colored rectangles
var_ranges <- data.frame(
  Variable = c("Z1", "Z3"),
  xmin = c(926, 902),
  xmax = c(936, 910)
)

block_names <- paste0("Block ", LETTERS[1:4])

rect_df <- expand.grid(
  Block = block_names,
  Variable = unique(var_ranges$Variable)
) %>%
  left_join(var_ranges, by = "Variable") %>%
  mutate(ymin = -Inf, ymax = Inf)

label_df <- rect_df %>%
  group_by(Block, Variable) %>%
  summarise(xmid = mean(c(xmin, xmax)), ymid = Inf, .groups = "drop")

ggplot(long_data, aes(x = X, y = Y, group = SampleID, color = cluster_y_4)) +
  # Rectangle annotations
  geom_rect(data = rect_df,
            inherit.aes = FALSE,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = "orange"),
            alpha = 0.2, show.legend = FALSE) +
  
  # Curves
  geom_line(alpha = 0.7) +
  
  # Facet into panels
  facet_wrap(~ Block, ncol = 2) +
  
  # Class colors
  scale_color_manual(values = c("1" = "#FFA15A", "2" = "#19D3F3", "3" = "#FF6692", "4" = "#B6E880")) +
  
  # Label text
  geom_text(data = label_df,
            inherit.aes = FALSE,
            aes(x = xmid, y = ymid, label = Variable),
            vjust = -0.3, size = 3, color = "black") +

  # Do not clip content outside the panel
  coord_cartesian(clip = "off") +
  
  labs(x = "Wavelength", y = "Absorbance", color = "Class") +
  theme_minimal(base_size = 14)


# ==========================================
# strawberry-curves-2                       Figure 15 
# ==========================================
library(ggplot2)
library(dplyr)
library(tidyr)
library(readxl)
library(tidyverse)
library(magrittr)
library(plotly)
library(RColorBrewer)

data <- read.csv("data/MIR_Fruit_purees.csv", header = TRUE)

rownames(data) <- data[, 1]
data <- data[, -1]

col_names <- colnames(data)

labels <- ifelse(startsWith(col_names, "Strawberry"), 1, 0)
labels

strawberry_plot_function_ggplot <- function(data, labels, save_path = NULL,
                                            width = 10, height = 6, dpi = 300) {
  
  x_vals <- as.numeric(colnames(data))
  stopifnot(length(labels) == nrow(data))
  
  plot_data <- as.data.frame(data) %>%
    mutate(Sample = row_number(),
           Label = labels) %>%
    pivot_longer(
      cols = all_of(as.character(x_vals)),
      names_to = "Wavelength",
      values_to = "Reflectance"
    ) %>%
    mutate(
      Wavelength = as.numeric(Wavelength),
      Group = ifelse(Label == 1, "Strawberry", "NON-Strawberry")
    )
  
  p <- ggplot(plot_data, aes(x = Wavelength, y = Reflectance, group = Sample)) +
    # Draw blue first
    geom_line(data = subset(plot_data, Group == "NON-Strawberry"),
              color = "lightblue", linewidth = 0.3) +
    # Then draw red
    geom_line(data = subset(plot_data, Group == "Strawberry"),
              color = "lightpink", linewidth = 0.3) +
    labs(x = "Wavelength (nm)",
         y = "Reflectance Spectrum") +
    theme_minimal() +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  
  if (!is.null(save_path)) {
print(p)
  }
  
  return(p)
}

strawberry_plot_function_ggplot(t(data), labels)

# ==========================================
# Categorization-1st-diff                   Figure 16
# ==========================================
library(readxl)
library(tidyverse)
library(magrittr)
library(plotly)
library(RColorBrewer)

data <- read.csv("data/MIR_Fruit_purees.csv", header = TRUE)

rownames(data) <- data[, 1]
data <- data[, -1]

col_names <- colnames(data)

labels <- ifelse(startsWith(col_names, "Strawberry"), 1, 0)

# sample 50 from each of classes 0 and 1
set.seed(123)
indices_0 <- which(labels == 0)
indices_1 <- which(labels == 1)

selected_0 <- sample(indices_0, 50)
selected_1 <- sample(indices_1, 50)

index <- c(selected_0, selected_1)

test <- data[,index]

n <- ncol(test) 

# First-order differences
A <- test  # fix: use test instead of data
D1 <- matrix(0, nrow = nrow(A) - 1, ncol = ncol(A))

for (i in 1:ncol(A)) {
  for (j in 2:(nrow(A))) {
    D1[j - 1, i] <- A[j, i] - A[j - 1, i]
  }
}

D1 <- as.data.frame(D1)

df <- D1  

df_long <- df %>%
  pivot_longer(cols = everything(), names_to = "fstd", values_to = "value")

ddf = df_long
data_val <- ddf[,2] |> as.data.frame()

h = hclust(dist(data_val), method = "ward.D2")
plot(h , hang=-1)

clus =  16 
c = cutree(h, clus)

GetMinMax <- function(i){
  summ <- summary( data_val[which(c==i),1] )
  return( summ[c(1,6)] ) 
}

MinMax <- sapply( 1:clus, GetMinMax ) 
print("double check if the order matched")
SmalltoLarge <- order( MinMax[1,] ) 

temp <- MinMax[1,][ SmalltoLarge ]
temp2 <- MinMax[2,][ SmalltoLarge ]

n_val <- nrow(data_val)
cum <- cumsum( table(c)[ SmalltoLarge ] ) / n_val
freq <- unname( table(c)[ SmalltoLarge ] )


plot( ecdf(data_val[,1]), col="grey", cex=0.2 ,main="",xlab="")
segments( temp[1], 0, temp2[1], cum[1], col="blue",lwd=3) 
for( i in 1:(clus-1) ){
  segments( temp[(i+1)], cum[i], temp2[(i+1)], cum[(i+1)], col="blue",lwd=3 ) 
}

plot(NULL, xlim=c(min(data_val[,1]),max(data_val[,1])), 
     ylim=c(0, max(freq)), ylab="freq.", xlab="value")

for(i in 1:clus){
  rect( temp[i], 0, temp2[i], freq[i])
}

# ==========================================
# straw-1st-2nd-diff                        figure 17
# ==========================================
library(tidyverse)
library(infotheo)
library(patchwork)

data <- read.csv("data/MIR_Fruit_purees.csv", header = TRUE)
rownames(data) <- data[, 1]
data <- data[, -1]
col_names <- colnames(data)
labels <- ifelse(startsWith(col_names, "Strawberry"), 1, 0)

first <- read.csv("data/strawberry_1stc_16_n_983.csv")
first  <- sapply(first, as.character)

X_entropy <- c()
Y_entropy <- c()
mutual    <- c()

for (i in 1:(nrow(first)-1)){
  X_entropy[i] <- infotheo::entropy(as.character(first[i,]),   method = "emp")
  Y_entropy[i] <- infotheo::entropy(as.character(first[i+1,]), method = "emp")
  mutual[i] <- 0.5 * (X_entropy[i] - infotheo::condentropy(as.character(first[i,]),   as.character(first[i+1,])))/X_entropy[i] +
    0.5 * (Y_entropy[i] - infotheo::condentropy(as.character(first[i+1,]), as.character(first[i,])))/Y_entropy[i]
}
x_values <- seq(903.187, 1798.704, length.out = length(mutual))

spec_mat <- t(data)  # row=sample, col=wavelength
x_vals   <- as.numeric(colnames(spec_mat))

spec_df <- as.data.frame(spec_mat) |>
  mutate(
    Sample = row_number(),
    Label  = factor(labels, levels = c(0, 1),
                    labels = c("NON-Strawberry", "Strawberry"))
  ) |>
  pivot_longer(
    cols = -c(Sample, Label),
    names_to  = "Wavelength",
    values_to = "Reflectance"
  ) |>
  mutate(Wavelength = as.numeric(Wavelength))

# --------- 1st-order MI ---------
mi1_df <- tibble(
  Wavelength = x_values,
  MI         = mutual
)

# --------- 2nd-order MI ---------
second <- read.csv("data/strawberry_2ndc_16_n_983.csv")
second <- sapply(second, as.character)

mutual2 <- numeric(nrow(second) - 1)
for (i in 1:(nrow(second)-1)){
  X_e <- infotheo::entropy(as.character(second[i,]),   method = "emp")
  Y_e <- infotheo::entropy(as.character(second[i+1,]), method = "emp")
  mutual2[i] <- 0.5 * (X_e - infotheo::condentropy(as.character(second[i,]),   as.character(second[i+1,])))/X_e +
    0.5 * (Y_e - infotheo::condentropy(as.character(second[i+1,]), as.character(second[i,])))/Y_e
}
x_values2 <- if (length(mutual2) == length(x_values)) {
  x_values
} else {
  seq(min(x_values), max(x_values), length.out = length(mutual2))
}

mi2_df <- tibble(
  Wavelength = x_values2,
  MI         = mutual2
)

y_rng  <- range(spec_df$Reflectance, na.rm = TRUE)

mi1_rng <- range(mi1_df$MI, na.rm = TRUE)
if (diff(mi1_rng) == 0) mi1_rng[2] <- mi1_rng[1] + 1e-6
s1 <- diff(y_rng) / diff(mi1_rng)
b1 <- y_rng[1] - s1 * mi1_rng[1]

mi2_rng <- range(mi2_df$MI, na.rm = TRUE)
if (diff(mi2_rng) == 0) mi2_rng[2] <- mi2_rng[1] + 1e-6
s2 <- diff(y_rng) / diff(mi2_rng)
b2 <- y_rng[1] - s2 * mi2_rng[1]

base_spec_blue <- geom_line(
  data = dplyr::filter(spec_df, Label == "NON-Strawberry"),
  aes(x = Wavelength, y = Reflectance, group = Sample),
  color = "lightblue", alpha = 0.6, linewidth = 0.7
)
base_spec_red <- geom_line(
  data = dplyr::filter(spec_df, Label == "Strawberry"),
  aes(x = Wavelength, y = Reflectance, group = Sample),
  color = "lightpink", alpha = 0.6, linewidth = 0.7
)

# --------- 1st-order MI (dark red) ---------
p1 <- ggplot() +
  base_spec_blue + base_spec_red +
  geom_line(
    data = mi1_df,
    aes(x = Wavelength, y = MI * s1 + b1),
    color = "darkred", linewidth = 1
  ) +
  geom_point(
    data = mi1_df,
    aes(x = Wavelength, y = MI * s1 + b1),
    color = "darkred", size = 1.2
  ) +
  scale_y_continuous(
    name = "Reflectance Spectrum",
    sec.axis = sec_axis(~ (. - b1) / s1, name = "Mutual Information (1st)")
  ) +
  labs(
    title = "Scaled Mutual Information (1st differences)",
    x = "Wavelength (nm)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 16, face = "bold")
  )

# --------- 2nd-order MI (dark purple) ---------
p2 <- ggplot() +
  base_spec_blue + base_spec_red +
  geom_line(
    data = mi2_df,
    aes(x = Wavelength, y = MI * s2 + b2),
    color = "darkmagenta", linewidth = 1
  ) +
  geom_point(
    data = mi2_df,
    aes(x = Wavelength, y = MI * s2 + b2),
    color = "darkmagenta", size = 1.2
  ) +
  scale_y_continuous(
    name = "Reflectance Spectrum",
    sec.axis = sec_axis(~ (. - b2) / s2, name = "Mutual Information (2nd)")
  ) +
  labs(
    title = "Scaled Mutual Information (2nd differences)",
    x = "Wavelength (nm)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 16, face = "bold")
  )

combined <- p1 / p2 + plot_layout(guides = "collect")

combined
# ==========================================
# straw-Z3-tree-split                       figure 18
# ==========================================
first <- read.csv("data/strawberry_1stc_16_n_983.csv")
new_data <- read.csv("data/MIR_Fruit_purees.csv", header = TRUE)
rownames(new_data) <- new_data[, 1]
col_names <- colnames(new_data[, -1])

labels <- ifelse(startsWith(col_names, "Strawberry"), 1, 0)
labels
new_data <- as.data.frame(t(new_data[,-1]))

first$wavelenght <- seq(903.187, 1802.564, length.out=234)
#Z3
data <- data.frame(
  X1 = as.numeric(first[36,][-984]),
  X2 = as.numeric(first[37,][-984]),
  X3 = as.numeric(first[38,][-984]),
  X4 = as.numeric(first[39,][-984]),
  X5 = as.numeric(first[40,][-984]),
  X6 = as.numeric(first[41,][-984])
)
hc_tree3 <- hclust(dist(data), method = "ward.D2")

dendrogram_obj <- hc_tree3 %>%
  as.dendrogram

dendrogram_obj %>%
  dendextend::color_branches(.,
                             k = 4,
                             groupLabels = TRUE) %>%
  plot(leaflab = "none")


rect.hclust(hc_tree3, k = 20, which = c(13, 14, 15, 17, 18), border = "red")
rect.hclust(hc_tree3, k = 12, which = c(10, 12), border = "red")

dendroextend_cut_tree <- dendextend::cutree(
  dendrogram_obj,
  k = 4,
  order_clusters_as_data = FALSE
)

data_order <- order.dendrogram(dendrogram_obj)

new_data$Z3 <- NA
new_data[data_order,]$Z3 <- dendroextend_cut_tree

table(Y = labels, Z3 = new_data$Z3)

#Z3 small group
dendroextend_cut_tree_20 <- dendextend::cutree(
  dendrogram_obj,
  k = 20,
  order_clusters_as_data = FALSE
)

new_data$Z3_20 <- NA
new_data[data_order,]$Z3_20 <- dendroextend_cut_tree_20

sub_clusters <- c(13, 14, 15)
indices_9 <- which(new_data$Z3_20 %in% sub_clusters) # further split cluster 9
new_data[indices_9,]$Z3 <- paste0(new_data[indices_9,]$Z3, "_", new_data[indices_9,]$Z3_20)
new_data[which(new_data$Z3 == "4_13"),]$Z3 <- '4_1'
new_data[which(new_data$Z3 == "4_14"),]$Z3 <- '4_2'
new_data[which(new_data$Z3 == "4_15"),]$Z3 <- '4_3'

sub_clusters <- c(17, 18)
indices_11 <- which(new_data$Z3_20 %in% sub_clusters) # further split cluster 11
new_data[indices_11,]$Z3 <- paste0(new_data[indices_11,]$Z3, "_", new_data[indices_11,]$Z3_20)
new_data[which(new_data$Z3 == "4_17"),]$Z3 <- '4_5'
new_data[which(new_data$Z3 == "4_18"),]$Z3 <- '4_6'

table_data <- table(Y = labels, Z3 = new_data$Z3)

dendroextend_cut_tree_12 <- dendextend::cutree(
  dendrogram_obj,
  k = 12,
  order_clusters_as_data = FALSE
)
new_data$Z3_12 <- NA
new_data[data_order,]$Z3_12 <- dendroextend_cut_tree_12

sub_clusters <- c(10, 12)
indices <- which(new_data$Z3_12 %in% sub_clusters) # further split cluster 9
new_data[indices,]$Z3 <- paste0(new_data[indices,]$Z3, "_", new_data[indices,]$Z3_12)
new_data[which(new_data$Z3 == "4_10"),]$Z3 <- '4_4'
new_data[which(new_data$Z3 == "4_12"),]$Z3 <- '4_7'
table(Y = labels, Z3 = new_data$Z3)

# ========================================== 
# straw-Z6-tree-split                       figure 19
# ==========================================
library(readxl)
library(tidyverse)
library(magrittr)
library(plotly)
library(RColorBrewer)
############use specific wavelength (first)
first <- read.csv("data/strawberry_1stc_16_n_983.csv")
new_data <- read.csv("data/MIR_Fruit_purees.csv", header = TRUE)
rownames(new_data) <- new_data[, 1]
col_names <- colnames(new_data[, -1])

labels <- ifelse(startsWith(col_names, "Strawberry"), 1, 0)
labels
new_data <- as.data.frame(t(new_data[,-1]))
new_data$y <- labels

first$wavelenght <- seq(903.187, 1802.564, length.out=234)


#Z6
data <- data.frame(
  X1 = as.numeric(first[69,][-984]),
  X2 = as.numeric(first[70,][-984]),
  X3 = as.numeric(first[71,][-984]),
  X4 = as.numeric(first[72,][-984]),
  X5 = as.numeric(first[73,][-984]),
  X6 = as.numeric(first[74,][-984]),
  X7 = as.numeric(first[75,][-984]),
  X8 = as.numeric(first[76,][-984]),
  X9 = as.numeric(first[77,][-984]),
  X10 = as.numeric(first[78,][-984]),
  X11 = as.numeric(first[79,][-984]),
  X12 = as.numeric(first[80,][-984]),
  X13 = as.numeric(first[81,][-984])
)
hc_tree6 <- hclust(dist(data), method = "ward.D2")

dendrogram_obj <- hc_tree6 %>%
  as.dendrogram


dendrogram_obj %>%
  dendextend::color_branches(.,
                             k = 10,
                             groupLabels = TRUE) %>%
  plot

rect.hclust(hc_tree6, k = 20, which = c(3, 4, 14, 15, 16, 17), border = "red")

dendroextend_cut_tree <- dendextend::cutree(
  dendrogram_obj,
  k = 10,
  order_clusters_as_data = FALSE
)

data_order <- order.dendrogram(dendrogram_obj)

new_data$Z6 <- NA
new_data[data_order,]$Z6 <- dendroextend_cut_tree

table(Y = new_data$y, Z6 = new_data$Z6)

#Z6 small group
dendroextend_cut_tree_20 <- dendextend::cutree(
  dendrogram_obj,
  k = 20,
  order_clusters_as_data = FALSE
)

new_data$Z6_20 <- NA
new_data[data_order,]$Z6_20 <- dendroextend_cut_tree_20

sub_clusters <- c(3, 4)
indices_2 <- which(new_data$Z6_20 %in% sub_clusters) # further split cluster 2
new_data[indices_2,]$Z6 <- paste0(new_data[indices_2,]$Z6, "_", new_data[indices_2,]$Z6_20)
new_data[which(new_data$Z6 == "2_3"),]$Z6 <- '2_1'
new_data[which(new_data$Z6 == "2_4"),]$Z6 <- '2_2'

sub_clusters <- c(14, 15)
indices_7 <- which(new_data$Z6_20 %in% sub_clusters) # further split cluster 7
new_data[indices_7,]$Z6 <- paste0(new_data[indices_7,]$Z6, "_", new_data[indices_7,]$Z6_20)
new_data[which(new_data$Z6 == "7_14"),]$Z6 <- '7_1'
new_data[which(new_data$Z6 == "7_15"),]$Z6 <- '7_2'

sub_clusters <- c(16, 17)
indices_8 <- which(new_data$Z6_20 %in% sub_clusters) # further split cluster 8
new_data[indices_8,]$Z6 <- paste0(new_data[indices_8,]$Z6, "_", new_data[indices_8,]$Z6_20)
new_data[which(new_data$Z6 == "8_16"),]$Z6 <- '8_1'
new_data[which(new_data$Z6 == "8_17"),]$Z6 <- '8_2'

table(Y = new_data$y, Z6 = new_data$Z6)
# ========================================== 
# straw-Z6-tree-split                       figure 20
# ==========================================
#Z7
data <- data.frame(
  X1 = as.numeric(first[85,][-984]),
  X2 = as.numeric(first[86,][-984]),
  X3 = as.numeric(first[87,][-984]),
  X4 = as.numeric(first[88,][-984]),
  X5 = as.numeric(first[89,][-984]),
  X6 = as.numeric(first[90,][-984]),
  X7 = as.numeric(first[91,][-984]),
  X8 = as.numeric(first[92,][-984]),
  X9 = as.numeric(first[93,][-984]),
  X10 = as.numeric(first[94,][-984]),
  X11 = as.numeric(first[95,][-984])
)

hc_tree7 <- hclust(dist(data), method = "ward.D2")

dendrogram_obj <- hc_tree7 %>%
  as.dendrogram


dendrogram_obj %>%
  dendextend::color_branches(.,
                             k = 4,
                             groupLabels = TRUE) %>%
  plot

rect.hclust(hc_tree7, k = 20, which = c(12, 13, 14, 15, 18, 19, 20), border = "red")
rect.hclust(hc_tree7, k = 10, which = c(9), border = "red")

dendroextend_cut_tree <- dendextend::cutree(
  dendrogram_obj,
  k = 4,
  order_clusters_as_data = FALSE
)

data_order <- order.dendrogram(dendrogram_obj)

new_data$Z7 <- NA
new_data[data_order,]$Z7 <- dendroextend_cut_tree

table(Y = new_data$y, Z7 = new_data$Z7)


#Z7 small group
dendroextend_cut_tree_20 <- dendextend::cutree(
  dendrogram_obj,
  k = 20,
  order_clusters_as_data = FALSE
)

new_data$Z7_20 <- NA
new_data[data_order,]$Z7_20 <- dendroextend_cut_tree_20

sub_clusters <- c(12, 13, 14, 15)
indices_3 <- which(new_data$Z7_20 %in% sub_clusters) # further split cluster 8
new_data[indices_3,]$Z7 <- paste0(new_data[indices_3,]$Z7, "_", new_data[indices_3,]$Z7_20)
new_data[which(new_data$Z7 == "3_12"),]$Z7 <- '3_1'
new_data[which(new_data$Z7 == "3_13"),]$Z7 <- '3_2'
new_data[which(new_data$Z7 == "3_14"),]$Z7 <- '3_3'
new_data[which(new_data$Z7 == "3_15"),]$Z7 <- '3_4'

sub_clusters <- c(18, 19, 20)
indices_4 <- which(new_data$Z7_20 %in% sub_clusters) # further split cluster 10
new_data[indices_4,]$Z7 <- paste0(new_data[indices_4,]$Z7, "_", new_data[indices_4,]$Z7_20)
new_data[which(new_data$Z7 == "4_18"),]$Z7 <- '4_2'
new_data[which(new_data$Z7 == "4_19"),]$Z7 <- '4_3'
new_data[which(new_data$Z7 == "4_20"),]$Z7 <- '4_4'
new_data[which(new_data$Z7 == "4"),]$Z7 <- '4_1'

table(Y = new_data$y, Z7 = new_data$Z7)

# ========================================== 
# straw-z4-2nd-tree                         figure 21
# ==========================================
#Z4
new_data <- read.csv("data/MIR_Fruit_purees.csv", header = TRUE)
rownames(new_data) <- new_data[, 1]
col_names <- colnames(new_data[, -1])

labels <- ifelse(startsWith(col_names, "Strawberry"), 1, 0)
labels
new_data <- as.data.frame(t(new_data[,-1]))
new_data$y <- labels

second <- read.csv("data/strawberry_2ndc_16_n_983.csv")
second$wavelenght <- seq(903.187, 1798.704, length.out=233)
data <- data.frame(
  X1 = as.numeric(second[61,][-984]),
  X2 = as.numeric(second[62,][-984]),
  X3 = as.numeric(second[63,][-984]),
  X4 = as.numeric(second[64,][-984]),
  X5 = as.numeric(second[65,][-984]),
  X6 = as.numeric(second[66,][-984])
)
hc_tree4 <- hclust(dist(data), method = "ward.D2")

dendrogram_obj <- hc_tree4 %>%
  as.dendrogram


dendrogram_obj %>%
  dendextend::color_branches(.,
                             k = 10,
                             groupLabels = TRUE) %>%
  plot

rect.hclust(hc_tree4, k = 14, which = c(2, 5, 6), border = "red")
rect.hclust(hc_tree4, k = 25, which = c(1, 2, 4, 5, 6, 7), border = "red")


dendroextend_cut_tree <- dendextend::cutree(
  dendrogram_obj,
  k = 10,
  order_clusters_as_data = FALSE
)

data_order <- order.dendrogram(dendrogram_obj)

new_data$Z4_second <- NA
new_data[data_order,]$Z4_second <- dendroextend_cut_tree

table(Y = new_data$y, Z4 = new_data$Z4_second)

#Z4 small group
dendroextend_cut_tree_14 <- dendextend::cutree(
  dendrogram_obj,
  k = 14,
  order_clusters_as_data = FALSE
)

new_data$Z4_second_14 <- NA
new_data[data_order,]$Z4_second_14 <- dendroextend_cut_tree_14

sub_clusters <- c(2)
indices_2 <- which(new_data$Z4_second_14 %in% sub_clusters) # further split cluster 2
new_data[indices_2,]$Z4_second <- paste0(new_data[indices_2,]$Z4_second, "_", new_data[indices_2,]$Z4_second_14)
new_data[which(new_data$Z4_second == "2_2"),]$Z4_second <- '2_1'

sub_clusters <- c(5, 6)
indices_3 <- which(new_data$Z4_second_14 %in% sub_clusters) # further split cluster 3
new_data[indices_3,]$Z4_second <- paste0(new_data[indices_3,]$Z4_second, "_", new_data[indices_3,]$Z4_second_14)
new_data[which(new_data$Z4_second == "3_5"),]$Z4_second <- '3_3'
new_data[which(new_data$Z4_second == "3_6"),]$Z4_second <- '3_4'

dendroextend_cut_tree_25 <- dendextend::cutree(
  dendrogram_obj,
  k = 25,
  order_clusters_as_data = FALSE
)

new_data$Z4_second_25 <- NA
new_data[data_order,]$Z4_second_25 <- dendroextend_cut_tree_25

sub_clusters <- c(1, 2)
indices_1 <- which(new_data$Z4_second_25 %in% sub_clusters) # further split cluster 1
new_data[indices_1,]$Z4_second <- paste0(new_data[indices_1,]$Z4_second, "_", new_data[indices_1,]$Z4_second_25)
new_data[which(new_data$Z4_second == "1_1"),]$Z4_second <- '1_1'
new_data[which(new_data$Z4_second == "1_2"),]$Z4_second <- '1_2'

sub_clusters <- c(4, 5)
indices_2 <- which(new_data$Z4_second_25 %in% sub_clusters) # further split cluster 2_2
new_data[indices_2,]$Z4_second <- paste0(new_data[indices_2,]$Z4_second, "_", new_data[indices_2,]$Z4_second_25)
new_data[which(new_data$Z4_second == "2_4"),]$Z4_second <- '2_2'
new_data[which(new_data$Z4_second == "2_5"),]$Z4_second <- '2_3'

sub_clusters <- c(6, 7)
indices_3 <- which(new_data$Z4_second_25 %in% sub_clusters) # further split cluster 3_1
new_data[indices_3,]$Z4_second <- paste0(new_data[indices_3,]$Z4_second, "_", new_data[indices_3,]$Z4_second_25)
new_data[which(new_data$Z4_second == "3_6"),]$Z4_second <- '3_1'
new_data[which(new_data$Z4_second == "3_7"),]$Z4_second <- '3_2'

table(Y = new_data$y, Z4 = new_data$Z4_second)

# ========================================== 
# straw-Z7-Z4-tree                          figure 22
# ==========================================
library(readxl)
library(tidyverse)
library(magrittr)
library(plotly)
library(RColorBrewer)
first <- read.csv("data/strawberry_1stc_16_n_983.csv")
second <- read.csv("data/strawberry_2ndc_16_n_983.csv")
new_data <- read.csv("data/MIR_Fruit_purees.csv", header = TRUE)
rownames(new_data) <- new_data[, 1]
col_names <- colnames(new_data[, -1])

labels <- ifelse(startsWith(col_names, "Strawberry"), 1, 0)
labels
new_data <- as.data.frame(t(new_data[,-1]))
new_data$y <- labels

second$wavelenght <- seq(903.187, 1798.704, length.out=233)
first$wavelenght <- seq(903.187, 1802.564, length.out=234)

#Z4+Z7
data <- data.frame(
  X1 = as.numeric(first[85,][-984]),
  X2 = as.numeric(first[86,][-984]),
  X3 = as.numeric(first[87,][-984]),
  X4 = as.numeric(first[88,][-984]),
  X5 = as.numeric(first[89,][-984]),
  X6 = as.numeric(first[90,][-984]),
  X7 = as.numeric(first[91,][-984]),
  X8 = as.numeric(first[92,][-984]),
  X9 = as.numeric(first[93,][-984]),
  X10 = as.numeric(first[94,][-984]),
  X11 = as.numeric(first[95,][-984]),
  X12 = as.numeric(second[61,][-984]),
  X13 = as.numeric(second[62,][-984]),
  X14 = as.numeric(second[63,][-984]),
  X15 = as.numeric(second[64,][-984]),
  X16 = as.numeric(second[65,][-984]),
  X17 = as.numeric(second[66,][-984])
)
hc_tree4 <- hclust(dist(data), method = "ward.D2")

dendrogram_obj <- hc_tree4 %>%
  as.dendrogram

dendrogram_obj %>%
  dendextend::color_branches(.,
                             k = 21,
                             groupLabels = TRUE) %>%
  plot(leaflab = "none")

dendroextend_cut_tree <- dendextend::cutree(
  dendrogram_obj,
  k = 21,
  order_clusters_as_data = FALSE
)

data_order <- order.dendrogram(dendrogram_obj)
new_data$Z4_Z7 <- NA
new_data[data_order,]$Z4_Z7 <- dendroextend_cut_tree

table(new_data$y, new_data$Z4_Z7)
# ========================================== 
# strawberry-odds-vs-overlap                figure 23
# ==========================================
library(readxl)
library(tidyverse)
library(magrittr)
library(plotly)
library(RColorBrewer)
library(dplyr)
############use specific wavelength (first)
first <- read.csv("data/strawberry_1stc_16_n_983.csv")
new_data <- read.csv("data/MIR_Fruit_purees.csv", header = TRUE)
rownames(new_data) <- new_data[, 1]
col_names <- colnames(new_data[, -1])

labels <- ifelse(startsWith(col_names, "Strawberry"), 1, 0)
labels
new_data <- as.data.frame(t(new_data[,-1]))

first$wavelenght <- seq(903.187, 1802.564, length.out=234)

### 1st order (8 variables)
#Z1
data <- data.frame(
  X1 = as.numeric(first[19,][-984]),
  X2 = as.numeric(first[20,][-984]),
  X3 = as.numeric(first[21,][-984]),
  X4 = as.numeric(first[22,][-984]),
  X5 = as.numeric(first[23,][-984]),
  X6 = as.numeric(first[24,][-984]),
  X7 = as.numeric(first[25,][-984])
)

hc_tree1 <- hclust(dist(data), method = "ward.D2")

dendrogram_obj <- hc_tree1 %>%
  as.dendrogram


dendrogram_obj %>%
  dendextend::color_branches(.,
                             k = 9,
                             groupLabels = TRUE) %>%
  
  plot


dendroextend_cut_tree <- dendextend::cutree(
  dendrogram_obj,
  k = 9,
  order_clusters_as_data = FALSE
)

data_order <- order.dendrogram(dendrogram_obj)

new_data$Z1 <- NA

new_data[data_order,]$Z1 <- dendroextend_cut_tree
new_data$y <- labels

table(Y = new_data$y, Z1 = new_data$Z1)

#Z2
data <- data.frame(
  X1 = as.numeric(first[21,][-984]),
  X2 = as.numeric(first[22,][-984]),
  X3 = as.numeric(first[23,][-984]),
  X4 = as.numeric(first[24,][-984]),
  X5 = as.numeric(first[25,][-984]),
  X6 = as.numeric(first[26,][-984]),
  X7 = as.numeric(first[27,][-984]),
  X8 = as.numeric(first[28,][-984]),
  X9 = as.numeric(first[29,][-984]),
  X10 = as.numeric(first[30,][-984]),
  X11 = as.numeric(first[31,][-984]),
  X12 = as.numeric(first[32,][-984]),
  X13 = as.numeric(first[33,][-984]),
  X14 = as.numeric(first[34,][-984])
)

hc_tree2 <- hclust(dist(data), method = "ward.D2")

dendrogram_obj <- hc_tree2 %>%
  as.dendrogram


dendrogram_obj %>%
  dendextend::color_branches(.,
                             k = 10,
                             groupLabels = TRUE) %>%
  plot


dendroextend_cut_tree <- dendextend::cutree(
  dendrogram_obj,
  k = 10,
  order_clusters_as_data = FALSE
)

data_order <- order.dendrogram(dendrogram_obj)

new_data$Z2 <- NA
new_data[data_order,]$Z2 <- dendroextend_cut_tree

table(Y = new_data$y, Z2 = new_data$Z2)

table(Z1 = new_data$Z1, Z2 = new_data$Z2)


#Z3
data <- data.frame(
  X1 = as.numeric(first[36,][-984]),
  X2 = as.numeric(first[37,][-984]),
  X3 = as.numeric(first[38,][-984]),
  X4 = as.numeric(first[39,][-984]),
  X5 = as.numeric(first[40,][-984]),
  X6 = as.numeric(first[41,][-984])
)
hc_tree3 <- hclust(dist(data), method = "ward.D2")

dendrogram_obj <- hc_tree3 %>%
  as.dendrogram


dendrogram_obj %>%
  dendextend::color_branches(.,
                             k = 4,
                             groupLabels = TRUE) %>%
  plot


rect.hclust(hc_tree3, k = 20, which = c(13, 14, 15, 17, 18), border = "red")
rect.hclust(hc_tree3, k = 12, which = c(10, 12), border = "red")

dendroextend_cut_tree <- dendextend::cutree(
  dendrogram_obj,
  k = 4,
  order_clusters_as_data = FALSE
)

data_order <- order.dendrogram(dendrogram_obj)

new_data$Z3 <- NA
new_data[data_order,]$Z3 <- dendroextend_cut_tree

table(Y = new_data$y, Z3 = new_data$Z3)

#Z3 small group
dendroextend_cut_tree_20 <- dendextend::cutree(
  dendrogram_obj,
  k = 20,
  order_clusters_as_data = FALSE
)

new_data$Z3_20 <- NA
new_data[data_order,]$Z3_20 <- dendroextend_cut_tree_20

sub_clusters <- c(13, 14, 15)
indices_9 <- which(new_data$Z3_20 %in% sub_clusters) # further split cluster 9
new_data[indices_9,]$Z3 <- paste0(new_data[indices_9,]$Z3, "_", new_data[indices_9,]$Z3_20)
new_data[which(new_data$Z3 == "4_13"),]$Z3 <- '4_1'
new_data[which(new_data$Z3 == "4_14"),]$Z3 <- '4_2'
new_data[which(new_data$Z3 == "4_15"),]$Z3 <- '4_3'

sub_clusters <- c(17, 18)
indices_11 <- which(new_data$Z3_20 %in% sub_clusters) # further split cluster 11
new_data[indices_11,]$Z3 <- paste0(new_data[indices_11,]$Z3, "_", new_data[indices_11,]$Z3_20)
new_data[which(new_data$Z3 == "4_17"),]$Z3 <- '4_5'
new_data[which(new_data$Z3 == "4_18"),]$Z3 <- '4_6'

table_data <- table(Y = new_data$y, Z3 = new_data$Z3)

dendroextend_cut_tree_12 <- dendextend::cutree(
  dendrogram_obj,
  k = 12,
  order_clusters_as_data = FALSE
)
new_data$Z3_12 <- NA
new_data[data_order,]$Z3_12 <- dendroextend_cut_tree_12

sub_clusters <- c(10, 12)
indices <- which(new_data$Z3_12 %in% sub_clusters) # further split cluster 9
new_data[indices,]$Z3 <- paste0(new_data[indices,]$Z3, "_", new_data[indices,]$Z3_12)
new_data[which(new_data$Z3 == "4_10"),]$Z3 <- '4_4'
new_data[which(new_data$Z3 == "4_12"),]$Z3 <- '4_7'
table(Y = new_data$y, Z3 = new_data$Z3)

#Z4
data <- data.frame(
  X1 = as.numeric(first[42,][-984]),
  X2 = as.numeric(first[43,][-984]),
  X3 = as.numeric(first[44,][-984]),
  X4 = as.numeric(first[45,][-984]),
  X5 = as.numeric(first[46,][-984]),
  X6 = as.numeric(first[47,][-984]),
  X7 = as.numeric(first[48,][-984]),
  X8 = as.numeric(first[49,][-984])
)
hc_tree4 <- hclust(dist(data), method = "ward.D2")

dendrogram_obj <- hc_tree4 %>%
  as.dendrogram


dendrogram_obj %>%
  dendextend::color_branches(.,
                             k = 10,
                             groupLabels = TRUE) %>%
  plot


dendroextend_cut_tree <- dendextend::cutree(
  dendrogram_obj,
  k = 10,
  order_clusters_as_data = FALSE
)

data_order <- order.dendrogram(dendrogram_obj)

new_data$Z4 <- NA
new_data[data_order,]$Z4 <- dendroextend_cut_tree

table(Y = new_data$y, Z4 = new_data$Z4)

#Z5
data <- data.frame(
  X1 = as.numeric(first[63,][-984]),
  X2 = as.numeric(first[64,][-984]),
  X3 = as.numeric(first[65,][-984]),
  X4 = as.numeric(first[66,][-984]),
  X5 = as.numeric(first[67,][-984]),
  X6 = as.numeric(first[68,][-984])
)
hc_tree5 <- hclust(dist(data), method = "ward.D2")

dendrogram_obj <- hc_tree5 %>%
  as.dendrogram


dendrogram_obj %>%
  dendextend::color_branches(.,
                             k = 9,
                             groupLabels = TRUE) %>%
  plot


dendroextend_cut_tree <- dendextend::cutree(
  dendrogram_obj,
  k = 9,
  order_clusters_as_data = FALSE
)

data_order <- order.dendrogram(dendrogram_obj)

new_data$Z5 <- NA
new_data[data_order,]$Z5 <- dendroextend_cut_tree

table(Y = new_data$y, Z5 = new_data$Z5)


#Z6
data <- data.frame(
  X1 = as.numeric(first[69,][-984]),
  X2 = as.numeric(first[70,][-984]),
  X3 = as.numeric(first[71,][-984]),
  X4 = as.numeric(first[72,][-984]),
  X5 = as.numeric(first[73,][-984]),
  X6 = as.numeric(first[74,][-984]),
  X7 = as.numeric(first[75,][-984]),
  X8 = as.numeric(first[76,][-984]),
  X9 = as.numeric(first[77,][-984]),
  X10 = as.numeric(first[78,][-984]),
  X11 = as.numeric(first[79,][-984]),
  X12 = as.numeric(first[80,][-984]),
  X13 = as.numeric(first[81,][-984])
)
hc_tree6 <- hclust(dist(data), method = "ward.D2")

dendrogram_obj <- hc_tree6 %>%
  as.dendrogram


dendrogram_obj %>%
  dendextend::color_branches(.,
                             k = 10,
                             groupLabels = TRUE) %>%
  plot

rect.hclust(hc_tree6, k = 20, which = c(3, 4, 14, 15, 16, 17), border = "red")

dendroextend_cut_tree <- dendextend::cutree(
  dendrogram_obj,
  k = 10,
  order_clusters_as_data = FALSE
)

data_order <- order.dendrogram(dendrogram_obj)

new_data$Z6 <- NA
new_data[data_order,]$Z6 <- dendroextend_cut_tree

table(Y = new_data$y, Z6 = new_data$Z6)

#Z6 small group
dendroextend_cut_tree_20 <- dendextend::cutree(
  dendrogram_obj,
  k = 20,
  order_clusters_as_data = FALSE
)

new_data$Z6_20 <- NA
new_data[data_order,]$Z6_20 <- dendroextend_cut_tree_20

sub_clusters <- c(3, 4)
indices_2 <- which(new_data$Z6_20 %in% sub_clusters) # further split cluster 2
new_data[indices_2,]$Z6 <- paste0(new_data[indices_2,]$Z6, "_", new_data[indices_2,]$Z6_20)
new_data[which(new_data$Z6 == "2_3"),]$Z6 <- '2_1'
new_data[which(new_data$Z6 == "2_4"),]$Z6 <- '2_2'

sub_clusters <- c(14, 15)
indices_7 <- which(new_data$Z6_20 %in% sub_clusters) # further split cluster 7
new_data[indices_7,]$Z6 <- paste0(new_data[indices_7,]$Z6, "_", new_data[indices_7,]$Z6_20)
new_data[which(new_data$Z6 == "7_14"),]$Z6 <- '7_1'
new_data[which(new_data$Z6 == "7_15"),]$Z6 <- '7_2'

sub_clusters <- c(16, 17)
indices_8 <- which(new_data$Z6_20 %in% sub_clusters) # further split cluster 8
new_data[indices_8,]$Z6 <- paste0(new_data[indices_8,]$Z6, "_", new_data[indices_8,]$Z6_20)
new_data[which(new_data$Z6 == "8_16"),]$Z6 <- '8_1'
new_data[which(new_data$Z6 == "8_17"),]$Z6 <- '8_2'

table(Y = new_data$y, Z6 = new_data$Z6)

#Z7
data <- data.frame(
  X1 = as.numeric(first[85,][-984]),
  X2 = as.numeric(first[86,][-984]),
  X3 = as.numeric(first[87,][-984]),
  X4 = as.numeric(first[88,][-984]),
  X5 = as.numeric(first[89,][-984]),
  X6 = as.numeric(first[90,][-984]),
  X7 = as.numeric(first[91,][-984]),
  X8 = as.numeric(first[92,][-984]),
  X9 = as.numeric(first[93,][-984]),
  X10 = as.numeric(first[94,][-984]),
  X11 = as.numeric(first[95,][-984])
)

hc_tree7 <- hclust(dist(data), method = "ward.D2")

dendrogram_obj <- hc_tree7 %>%
  as.dendrogram


dendrogram_obj %>%
  dendextend::color_branches(.,
                             k = 4,
                             groupLabels = TRUE) %>%
  plot

rect.hclust(hc_tree7, k = 20, which = c(12, 13, 14, 15, 18, 19, 20), border = "red")
rect.hclust(hc_tree7, k = 10, which = c(9), border = "red")

dendroextend_cut_tree <- dendextend::cutree(
  dendrogram_obj,
  k = 4,
  order_clusters_as_data = FALSE
)

data_order <- order.dendrogram(dendrogram_obj)

new_data$Z7 <- NA
new_data[data_order,]$Z7 <- dendroextend_cut_tree

table(Y = new_data$y, Z7 = new_data$Z7)

#Z7 small group
dendroextend_cut_tree_20 <- dendextend::cutree(
  dendrogram_obj,
  k = 20,
  order_clusters_as_data = FALSE
)

new_data$Z7_20 <- NA
new_data[data_order,]$Z7_20 <- dendroextend_cut_tree_20

sub_clusters <- c(12, 13, 14, 15)
indices_3 <- which(new_data$Z7_20 %in% sub_clusters) # further split cluster 8
new_data[indices_3,]$Z7 <- paste0(new_data[indices_3,]$Z7, "_", new_data[indices_3,]$Z7_20)
new_data[which(new_data$Z7 == "3_12"),]$Z7 <- '3_1'
new_data[which(new_data$Z7 == "3_13"),]$Z7 <- '3_2'
new_data[which(new_data$Z7 == "3_14"),]$Z7 <- '3_3'
new_data[which(new_data$Z7 == "3_15"),]$Z7 <- '3_4'

sub_clusters <- c(18, 19, 20)
indices_4 <- which(new_data$Z7_20 %in% sub_clusters) # further split cluster 10
new_data[indices_4,]$Z7 <- paste0(new_data[indices_4,]$Z7, "_", new_data[indices_4,]$Z7_20)
new_data[which(new_data$Z7 == "4_18"),]$Z7 <- '4_2'
new_data[which(new_data$Z7 == "4_19"),]$Z7 <- '4_3'
new_data[which(new_data$Z7 == "4_20"),]$Z7 <- '4_4'
new_data[which(new_data$Z7 == "4"),]$Z7 <- '4_1'

table(Y = new_data$y, Z7 = new_data$Z7)

#Z8
data <- data.frame(
  X1 = as.numeric(first[208,][-984]),
  X2 = as.numeric(first[209,][-984]),
  X3 = as.numeric(first[210,][-984]),
  X4 = as.numeric(first[211,][-984]),
  X5 = as.numeric(first[212,][-984]),
  X6 = as.numeric(first[213,][-984]),
  X7 = as.numeric(first[214,][-984]),
  X8 = as.numeric(first[215,][-984]),
  X9 = as.numeric(first[216,][-984]),
  X10 = as.numeric(first[217,][-984]),
  X11 = as.numeric(first[218,][-984]),
  X12 = as.numeric(first[219,][-984]),
  X13 = as.numeric(first[220,][-984]),
  X14 = as.numeric(first[221,][-984]),
  X15 = as.numeric(first[222,][-984]),
  X16 = as.numeric(first[223,][-984]),
  X17 = as.numeric(first[224,][-984])
)

hc_tree8 <- hclust(dist(data), method = "ward.D2")

dendrogram_obj <- hc_tree8 %>%
  as.dendrogram


dendrogram_obj %>%
  dendextend::color_branches(.,
                             k = 14,
                             groupLabels = TRUE) %>%
  plot


dendroextend_cut_tree <- dendextend::cutree(
  dendrogram_obj,
  k = 14,
  order_clusters_as_data = FALSE
)

data_order <- order.dendrogram(dendrogram_obj)

new_data$Z8 <- NA
new_data[data_order,]$Z8 <- dendroextend_cut_tree

table(Y = new_data$y, Z8 = new_data$Z8)


#### 2nd order (3 variables)
#Z2
second <- read.csv("data/strawberry_2ndc_16_n_983.csv")
second$wavelenght <- seq(903.187, 1798.704, length.out=233)
data <- data.frame(
  X1 = as.numeric(second[21,][-984]),
  X2 = as.numeric(second[22,][-984]),
  X3 = as.numeric(second[23,][-984]),
  X4 = as.numeric(second[24,][-984]),
  X5 = as.numeric(second[25,][-984]),
  X6 = as.numeric(second[26,][-984]),
  X7 = as.numeric(second[27,][-984])
)
hc_tree2 <- hclust(dist(data), method = "ward.D2")

dendrogram_obj <- hc_tree2 %>%
  as.dendrogram


dendrogram_obj %>%
  dendextend::color_branches(.,
                             k = 12,
                             groupLabels = TRUE) %>%
  plot


dendroextend_cut_tree <- dendextend::cutree(
  dendrogram_obj,
  k = 12,
  order_clusters_as_data = FALSE
)

data_order <- order.dendrogram(dendrogram_obj)

new_data$Z2_second <- NA
new_data[data_order,]$Z2_second <- dendroextend_cut_tree

table(Y = new_data$y, Z2 = new_data$Z2_second)

#Z3
data <- data.frame(
  X1 = as.numeric(second[33,][-984]),
  X2 = as.numeric(second[34,][-984]),
  X3 = as.numeric(second[35,][-984]),
  X4 = as.numeric(second[36,][-984])
)
hc_tree3 <- hclust(dist(data), method = "ward.D2")

dendrogram_obj <- hc_tree3 %>%
  as.dendrogram


dendrogram_obj %>%
  dendextend::color_branches(.,
                             k = 9,
                             groupLabels = TRUE) %>%
  plot

rect.hclust(hc_tree3, k = 20, which = c(15, 16, 19, 20), border = "red")


dendroextend_cut_tree <- dendextend::cutree(
  dendrogram_obj,
  k = 9,
  order_clusters_as_data = FALSE
)

data_order <- order.dendrogram(dendrogram_obj)

new_data$Z3_second <- NA
new_data[data_order,]$Z3_second <- dendroextend_cut_tree

table(Y = new_data$y, Z3 = new_data$Z3_second)

#Z3 small group
dendroextend_cut_tree_20 <- dendextend::cutree(
  dendrogram_obj,
  k = 20,
  order_clusters_as_data = FALSE
)

new_data$Z3_second_20 <- NA
new_data[data_order,]$Z3_second_20 <- dendroextend_cut_tree_20

sub_clusters <- c(15, 16)
indices_7 <- which(new_data$Z3_second_20 %in% sub_clusters) # further split cluster 7
new_data[indices_7,]$Z3_second <- paste0(new_data[indices_7,]$Z3_second, "_", new_data[indices_7,]$Z3_second_20)
new_data[which(new_data$Z3_second == "7_15"),]$Z3_second <- '7_1'
new_data[which(new_data$Z3_second == "7_16"),]$Z3_second <- '7_2'

sub_clusters <- c(19, 20)
indices_9 <- which(new_data$Z3_second_20 %in% sub_clusters) # further split cluster 9
new_data[indices_9,]$Z3_second <- paste0(new_data[indices_9,]$Z3_second, "_", new_data[indices_9,]$Z3_second_20)
new_data[which(new_data$Z3_second == "9_19"),]$Z3_second <- '9_1'
new_data[which(new_data$Z3_second == "9_20"),]$Z3_second <- '9_2'

table(Y = new_data$y, Z3 = new_data$Z3_second)

#Z4
data <- data.frame(
  X1 = as.numeric(second[61,][-984]),
  X2 = as.numeric(second[62,][-984]),
  X3 = as.numeric(second[63,][-984]),
  X4 = as.numeric(second[64,][-984]),
  X5 = as.numeric(second[65,][-984]),
  X6 = as.numeric(second[66,][-984])
)
hc_tree4 <- hclust(dist(data), method = "ward.D2")

dendrogram_obj <- hc_tree4 %>%
  as.dendrogram


dendrogram_obj %>%
  dendextend::color_branches(.,
                             k = 10,
                             groupLabels = TRUE) %>%
  plot

rect.hclust(hc_tree4, k = 14, which = c(2, 5, 6), border = "red")
rect.hclust(hc_tree4, k = 25, which = c(1, 2, 4, 5, 6, 7), border = "red")


dendroextend_cut_tree <- dendextend::cutree(
  dendrogram_obj,
  k = 10,
  order_clusters_as_data = FALSE
)

data_order <- order.dendrogram(dendrogram_obj)

new_data$Z4_second <- NA
new_data[data_order,]$Z4_second <- dendroextend_cut_tree

table(Y = new_data$y, Z4 = new_data$Z4_second)

#Z4 small group
dendroextend_cut_tree_14 <- dendextend::cutree(
  dendrogram_obj,
  k = 14,
  order_clusters_as_data = FALSE
)

new_data$Z4_second_14 <- NA
new_data[data_order,]$Z4_second_14 <- dendroextend_cut_tree_14

sub_clusters <- c(2)
indices_2 <- which(new_data$Z4_second_14 %in% sub_clusters) # further split cluster 2
new_data[indices_2,]$Z4_second <- paste0(new_data[indices_2,]$Z4_second, "_", new_data[indices_2,]$Z4_second_14)
new_data[which(new_data$Z4_second == "2_2"),]$Z4_second <- '2_1'

sub_clusters <- c(5, 6)
indices_3 <- which(new_data$Z4_second_14 %in% sub_clusters) # further split cluster 3
new_data[indices_3,]$Z4_second <- paste0(new_data[indices_3,]$Z4_second, "_", new_data[indices_3,]$Z4_second_14)
new_data[which(new_data$Z4_second == "3_5"),]$Z4_second <- '3_3'
new_data[which(new_data$Z4_second == "3_6"),]$Z4_second <- '3_4'

dendroextend_cut_tree_25 <- dendextend::cutree(
  dendrogram_obj,
  k = 25,
  order_clusters_as_data = FALSE
)

new_data$Z4_second_25 <- NA
new_data[data_order,]$Z4_second_25 <- dendroextend_cut_tree_25

sub_clusters <- c(1, 2)
indices_1 <- which(new_data$Z4_second_25 %in% sub_clusters) # further split cluster 1
new_data[indices_1,]$Z4_second <- paste0(new_data[indices_1,]$Z4_second, "_", new_data[indices_1,]$Z4_second_25)
new_data[which(new_data$Z4_second == "1_1"),]$Z4_second <- '1_1'
new_data[which(new_data$Z4_second == "1_2"),]$Z4_second <- '1_2'

sub_clusters <- c(4, 5)
indices_2 <- which(new_data$Z4_second_25 %in% sub_clusters) # further split cluster 2_2
new_data[indices_2,]$Z4_second <- paste0(new_data[indices_2,]$Z4_second, "_", new_data[indices_2,]$Z4_second_25)
new_data[which(new_data$Z4_second == "2_4"),]$Z4_second <- '2_2'
new_data[which(new_data$Z4_second == "2_5"),]$Z4_second <- '2_3'

sub_clusters <- c(6, 7)
indices_3 <- which(new_data$Z4_second_25 %in% sub_clusters) # further split cluster 3_1
new_data[indices_3,]$Z4_second <- paste0(new_data[indices_3,]$Z4_second, "_", new_data[indices_3,]$Z4_second_25)
new_data[which(new_data$Z4_second == "3_6"),]$Z4_second <- '3_1'
new_data[which(new_data$Z4_second == "3_7"),]$Z4_second <- '3_2'

table(Y = new_data$y, Z4 = new_data$Z4_second)


del_col <- c("Z3_20", "Z3_12" , "Z6_20", "Z7_20", "Z3_second_20", "Z4_second_14", "Z4_second_25")

new_data <- new_data[, -which(colnames(new_data) %in% del_col)]

table(new_data$y, new_data$Z4_second)

# overlapping area
entropy_overlap_by_column <- function(contingency_table,
                                      n_sim  = 1000,
                                      prefix = "Zcat") {
  # Basic quantities
  n_row <- nrow(contingency_table)
  n_col <- ncol(contingency_table)
  
  row_tot   <- rowSums(contingency_table)
  col_tot   <- colSums(contingency_table)
  grand_tot <- sum(contingency_table)
  
  # Probability distributions
  row_prob_alt <- lapply(1:n_row, function(i) {
    if (row_tot[i] == 0) rep(0, n_col) else contingency_table[i, ] / row_tot[i]
  })
  col_prob_null <- col_tot / grand_tot
  
  # 3D arrays
  alt_array  <- array(0, dim = c(n_row, n_col, n_sim))
  null_array <- array(0, dim = c(n_row, n_col, n_sim))
  
  # Simulation
  set.seed(42)
  for (s in 1:n_sim) {
    for (i in 1:n_row) {
      alt_array[i, , s]  <- rmultinom(1, row_tot[i], row_prob_alt[[i]])
      null_array[i, , s] <- rmultinom(1, row_tot[i], col_prob_null)
    }
  }
  
  # entropy function
  entropy <- function(p) { p <- p[p > 0]; -sum(p * log(p)) }
  
  breaks_fix <- seq(0, log(n_row), length.out = 50)  # upper bound determined by n_row
  
  result_df <- data.frame(Category = colnames(contingency_table),
                          Overlap = NA_real_,
                          stringsAsFactors = FALSE)
  
  # Compute overlap column by column
  for (k in 1:n_col) {
    cat_name <- colnames(contingency_table)[k]
    safe_cat_name <- make.names(cat_name)
    
    alt_ent  <- apply(alt_array[ , k, ], 2, function(v) entropy(v / sum(v)))
    null_ent <- apply(null_array[, k, ], 2, function(v) entropy(v / sum(v)))
    
    alt_hist  <- hist(alt_ent,  breaks = breaks_fix, plot = FALSE)
    null_hist <- hist(null_ent, breaks = breaks_fix, plot = FALSE)
    overlap_ratio <- sum(pmin(alt_hist$counts, null_hist$counts)) / n_sim
    result_df$Overlap[k] <- overlap_ratio
    
    # PDF plot
    # pdf(paste0(prefix, "_", safe_cat_name, "_entropy.pdf"), width = 6, height = 4)
    # hist(null_ent, breaks = breaks_fix, col = rgb(0,0,1,0.5),
    #      main = paste0(prefix, " ", safe_cat_name, " | overlap=", round(overlap_ratio, 3)),
    #      xlab = "Entropy", ylim = c(0, max(alt_hist$counts, null_hist$counts)))
    # hist(alt_ent,  breaks = breaks_fix, col = rgb(1,0,0,0.5), add = TRUE)
    # legend("topright", legend = c("Null", "Alternative"), fill = c("blue", "red"))
    # dev.off()
  }
  
  return(result_df)
}


var_list <- c("Z1", "Z2", "Z3", "Z4", "Z5", "Z6", "Z7", "Z8",
              "Z2_second", "Z3_second", "Z4_second")

all_entropy_results <- list()

for (var_name in var_list) {
  tbl <- table(new_data$y, new_data[[var_name]])
  
  y_levels <- sort(unique(new_data$y))
  z_levels <- colnames(tbl)
  
  full_table <- matrix(0, nrow = length(y_levels), ncol = length(z_levels))
  rownames(full_table) <- paste0("Y", y_levels)
  colnames(full_table) <- make.names(z_levels)
  
  for (i in 1:nrow(tbl)) {
    for (j in 1:ncol(tbl)) {
      y_val <- rownames(tbl)[i]
      z_val <- colnames(tbl)[j]
      full_table[paste0("Y", y_val), make.names(z_val)] <- tbl[i, j]
    }
  }
  
  full_table <- apply(full_table, 2, as.numeric)
  rownames(full_table) <- paste0("Y", y_levels)
  
  # odds vector: (# of Y=1)/(# of Y=0)
  y0 <- full_table["Y0", ]
  y1 <- full_table["Y1", ]
  odds_vec <- ifelse(y0 == 0, NA, y1 / y0)  # avoid division by zero
  
  # entropy overlap simulation
  result <- entropy_overlap_by_column(full_table, prefix = var_name)
  
  # odds column
  result$Odds <- odds_vec[match(result$Category, colnames(full_table))]
  
  all_entropy_results[[var_name]] <- result
}


#odds vs overlap scatterplot
library(dplyr)
library(ggplot2)
library(ggrepel)

# Combine results across all variables
final_df <- bind_rows(
  lapply(names(all_entropy_results), function(var) {
    df <- all_entropy_results[[var]]
    df$Feature <- var
    df
  })
)

# label column (category name)
final_df <- final_df %>%
  mutate(Label = paste0(Feature, "_", Category))

# scatter plot
ggplot(final_df, aes(x = Odds, y = Overlap)) +
  geom_point(color = "#1f77b4", size = 2) +
  geom_text_repel(aes(label = Label), size = 3, max.overlaps = 20) +
  theme_minimal(base_size = 13) +
  labs(title = "Odds vs Overlap",
       x = "Odds = # Y=1/# Y=0",
       y = "Entropy Overlap") +
  theme(panel.grid.minor = element_blank())

# ========================================== 
# strawberry-order1-heatmap                 figure 24 27 30
# ==========================================
library(readxl)
library(tidyverse)
library(magrittr)
library(plotly)
library(RColorBrewer)
library(dplyr)
library(dendextend)
library(ComplexHeatmap)
library(circlize)
first <- read.csv("data/strawberry_1stc_16_n_983.csv")
first$wavelenght <- seq(903.187, 1802.564, length.out=234)
second <- read.csv("data/strawberry_2ndc_16_n_983.csv")
second$wavelenght <- seq(903.187, 1798.704, length.out=233)
new_data <- read.csv("data/MIR_Fruit_purees.csv", header = TRUE)
rownames(new_data) <- new_data[, 1]
col_names <- colnames(new_data[, -1])
labels <- ifelse(startsWith(col_names, "Strawberry"), 1, 0)
new_data <- as.data.frame(t(new_data[,-1]))
new_data$y <- labels

get_clusters <- function(src_data, rows, k_val, plot_title = "") {
  data_list <- lapply(rows, function(r) as.numeric(src_data[r, ][-984]))
  temp_data <- as.data.frame(data_list)
  colnames(temp_data) <- paste0("X", seq_along(rows))
  
  hc_tree <- hclust(dist(temp_data), method = "ward.D2")
  dendrogram_obj <- as.dendrogram(hc_tree)
  cut_tree <- dendextend::cutree(dendrogram_obj, k = k_val, order_clusters_as_data = FALSE)
  data_order <- order.dendrogram(dendrogram_obj)
  
  result_vec <- rep(NA, length(data_order))
  result_vec[data_order] <- cut_tree
  return(result_vec)
}
new_data$Z1 <- get_clusters(first, 19:25, 9)
new_data$Z2 <- get_clusters(first, 21:34, 10)
new_data$Z3 <- get_clusters(first, 36:41, 4)
Z3_20 <- get_clusters(first, 36:41, 20)
new_data$Z3 <- ifelse(Z3_20 %in% c(13, 14, 15, 17, 18), paste0(new_data$Z3, "_", Z3_20), new_data$Z3)
new_data$Z3[new_data$Z3 == "4_13"] <- '4_1'
new_data$Z3[new_data$Z3 == "4_14"] <- '4_2'
new_data$Z3[new_data$Z3 == "4_15"] <- '4_3'
new_data$Z3[new_data$Z3 == "4_17"] <- '4_5'
new_data$Z3[new_data$Z3 == "4_18"] <- '4_6'
Z3_12 <- get_clusters(first, 36:41, 12)
new_data$Z3 <- ifelse(Z3_12 %in% c(10, 12), paste0(new_data$Z3, "_", Z3_12), new_data$Z3)
new_data$Z3[new_data$Z3 == "4_10"] <- '4_4'
new_data$Z3[new_data$Z3 == "4_12"] <- '4_7'
new_data$Z4 <- get_clusters(first, 42:49, 10)
new_data$Z5 <- get_clusters(first, 63:68, 9)
new_data$Z6 <- get_clusters(first, 69:81, 10)
Z6_20 <- get_clusters(first, 69:81, 20)
new_data$Z6 <- ifelse(Z6_20 %in% c(3, 4, 14, 15, 16, 17), paste0(new_data$Z6, "_", Z6_20), new_data$Z6)
new_data$Z6[new_data$Z6 == "2_3"]  <- '2_1'
new_data$Z6[new_data$Z6 == "2_4"]  <- '2_2'
new_data$Z6[new_data$Z6 == "7_14"] <- '7_1'
new_data$Z6[new_data$Z6 == "7_15"] <- '7_2'
new_data$Z6[new_data$Z6 == "8_16"] <- '8_1'
new_data$Z6[new_data$Z6 == "8_17"] <- '8_2'
new_data$Z7 <- get_clusters(first, 85:95, 4)
Z7_20 <- get_clusters(first, 85:95, 20)
new_data$Z7 <- ifelse(Z7_20 %in% c(12, 13, 14, 15, 18, 19, 20), paste0(new_data$Z7, "_", Z7_20), new_data$Z7)
new_data$Z7[new_data$Z7 == "3_12"] <- '3_1'
new_data$Z7[new_data$Z7 == "3_13"] <- '3_2'
new_data$Z7[new_data$Z7 == "3_14"] <- '3_3'
new_data$Z7[new_data$Z7 == "3_15"] <- '3_4'
new_data$Z7[new_data$Z7 == "4_18"] <- '4_2'
new_data$Z7[new_data$Z7 == "4_19"] <- '4_3'
new_data$Z7[new_data$Z7 == "4_20"] <- '4_4'
new_data$Z7[new_data$Z7 == "4"]    <- '4_1'
new_data$Z8 <- get_clusters(first, 208:224, 14)
new_data$Z2_second <- get_clusters(second, 21:27, 12)
new_data$Z3_second <- get_clusters(second, 33:36, 9)
Z3_sec_20 <- get_clusters(second, 33:36, 20)
new_data$Z3_second <- ifelse(Z3_sec_20 %in% c(15, 16, 19, 20), paste0(new_data$Z3_second, "_", Z3_sec_20), new_data$Z3_second)
new_data$Z3_second[new_data$Z3_second == "7_15"] <- '7_1'
new_data$Z3_second[new_data$Z3_second == "7_16"] <- '7_2'
new_data$Z3_second[new_data$Z3_second == "9_19"] <- '9_1'
new_data$Z3_second[new_data$Z3_second == "9_20"] <- '9_2'
new_data$Z4_second <- get_clusters(second, 61:66, 10)
Z4_sec_14 <- get_clusters(second, 61:66, 14)
new_data$Z4_second <- ifelse(Z4_sec_14 %in% c(2, 5, 6), paste0(new_data$Z4_second, "_", Z4_sec_14), new_data$Z4_second)
new_data$Z4_second[new_data$Z4_second == "2_2"] <- '2_1'
new_data$Z4_second[new_data$Z4_second == "3_5"] <- '3_3'
new_data$Z4_second[new_data$Z4_second == "3_6"] <- '3_4'
Z4_sec_25 <- get_clusters(second, 61:66, 25)
new_data$Z4_second <- ifelse(Z4_sec_25 %in% c(1, 2, 4, 5, 6, 7), paste0(new_data$Z4_second, "_", Z4_sec_25), new_data$Z4_second)
new_data$Z4_second[new_data$Z4_second == "1_1"] <- '1_1'
new_data$Z4_second[new_data$Z4_second == "1_2"] <- '1_2'
new_data$Z4_second[new_data$Z4_second == "2_4"] <- '2_2'
new_data$Z4_second[new_data$Z4_second == "2_5"] <- '2_3'
new_data$Z4_second[new_data$Z4_second == "3_6"] <- '3_1'
new_data$Z4_second[new_data$Z4_second == "3_7"] <- '3_2'
# =========================================================================
# Entropy Overlap function
entropy_overlap_by_column <- function(contingency_table, n_sim = 1000, prefix = "Zcat") {
  n_row <- nrow(contingency_table)
  n_col <- ncol(contingency_table)
  
  row_tot   <- rowSums(contingency_table)
  col_tot   <- colSums(contingency_table)
  grand_tot <- sum(contingency_table)
  
  row_prob_alt <- lapply(1:n_row, function(i) {
    if (row_tot[i] == 0) rep(0, n_col) else contingency_table[i, ] / row_tot[i]
  })
  col_prob_null <- col_tot / grand_tot
  
  alt_array  <- array(0, dim = c(n_row, n_col, n_sim))
  null_array <- array(0, dim = c(n_row, n_col, n_sim))
  
  set.seed(42)
  for (s in 1:n_sim) {
    for (i in 1:n_row) {
      alt_array[i, , s]  <- rmultinom(1, row_tot[i], row_prob_alt[[i]])
      null_array[i, , s] <- rmultinom(1, row_tot[i], col_prob_null)
    }
  }
  
  entropy <- function(p) { p <- p[p > 0]; -sum(p * log(p)) }
  breaks_fix <- seq(0, log(n_row), length.out = 50)
  
  result_df <- data.frame(Category = colnames(contingency_table), Overlap = NA_real_, stringsAsFactors = FALSE)
  for (k in 1:n_col) {
    alt_ent  <- apply(alt_array[ , k, ], 2, function(v) entropy(v / sum(v)))
    null_ent <- apply(null_array[, k, ], 2, function(v) entropy(v / sum(v)))
    overlap_ratio <- sum(pmin(hist(alt_ent, breaks=breaks_fix, plot=F)$counts, 
                              hist(null_ent, breaks=breaks_fix, plot=F)$counts)) / n_sim
    result_df$Overlap[k] <- overlap_ratio
  }
  return(result_df)
}
# =========================================================================
# Compute 1st-order entropy and the 1st-order matrix
feature_vars <- c("Z1", "Z2", "Z3", "Z4", "Z5", "Z6", "Z7", "Z8", "Z2_second", "Z3_second", "Z4_second")
all_entropy_results <- list()
onehot_list <- list()
for (var_name in feature_vars) {
  # Compute entropy
  tbl <- table(new_data$y, new_data[[var_name]])
  if (!all(c("0","1") %in% rownames(tbl))) {
    tbl <- rbind(Y0 = tbl["0", , drop=FALSE], Y1 = tbl["1", , drop=FALSE])
  } else {
    rownames(tbl) <- c("Y0", "Y1")
  }
  odds_vec <- ifelse(tbl["Y0", ] == 0, Inf, tbl["Y1", ] / tbl["Y0", ])
  
  ent_res <- entropy_overlap_by_column(tbl)
  ent_res$Odds <- odds_vec[match(ent_res$Category, colnames(tbl))]
  ent_res$Feature <- var_name
  all_entropy_results[[var_name]] <- ent_res
  
  # Build one-hot matrix
  tmp_factor <- factor(new_data[[var_name]])
  mat_tmp <- model.matrix(~ tmp_factor - 1)
  clean_cat <- sub("tmp_factor", "", colnames(mat_tmp))
  colnames(mat_tmp) <- paste0(var_name, "_", clean_cat)
  onehot_list[[var_name]] <- mat_tmp
}
# Merge with the previous final_df
final_df <- bind_rows(all_entropy_results)
final_df$Label <- paste0(final_df$Feature, "_", final_df$Category)
binary_mat_1st <- do.call(cbind, onehot_list)
# =========================================================================
# Filter qualified_categories
qualified_categories <- list()
for (var in feature_vars) {
  tbl <- table(new_data$y, new_data[[var]])
  if (!all(c(0,1) %in% rownames(tbl))) {
    tbl <- rbind(Y0 = tbl["0", , drop=FALSE], Y1 = tbl["1", , drop=FALSE])
  } else { rownames(tbl) <- c("Y0", "Y1") }
  
  selected <- c()
  for (cat_name in colnames(tbl)) {
    n0 <- tbl["Y0", cat_name]; n1 <- tbl["Y1", cat_name]
    odds <- ifelse(n0 == 0, Inf, n1 / n0)
    if ((n0 + n1) > 30 && min(n0, n1) > 10 && odds >= 0.1 && odds <= 10) {
      selected <- c(selected, cat_name)
    }
  }
  qualified_categories[[var]] <- selected
}
# =========================================================================
# Generate base variables for the 2nd-order split
config_list <- list(
  Z1_split = list(src = first, rows = 19:25, k = 4),   Z2_split = list(src = first, rows = 21:34, k = 4),
  Z3_split = list(src = first, rows = 36:41, k = 4),   Z4_split = list(src = first, rows = 42:49, k = 4),
  Z5_split = list(src = first, rows = 63:68, k = 4),   Z6_split = list(src = first, rows = 69:81, k = 5),
  Z7_split = list(src = first, rows = 85:95, k = 4),   Z8_split = list(src = first, rows = 208:224, k = 4),
  Z2_second_split = list(src = second, rows = 21:27, k = 5), Z3_second_split = list(src = second, rows = 33:36, k = 5),
  Z4_second_split = list(src = second, rows = 61:66, k = 5)
)
for (var_name in names(config_list)) {
  cfg <- config_list[[var_name]]
  new_data[[var_name]] <- get_clusters(cfg$src, cfg$rows, cfg$k)
}
# =========================================================================
# Compute 2nd-order entropy and matrix
all_results_2nd <- data.frame()
split_vars <- paste0(feature_vars, "_split")
for (Zk in names(qualified_categories)) {
  categories <- qualified_categories[[Zk]]
  Zk_split_name <- paste0(Zk, "_split")
  
  for (cat_val in categories) {
    sub_data <- new_data[which(new_data[[Zk]] == cat_val), ]
    for (split_var in split_vars) {
      if (split_var == Zk_split_name) next
      tbl <- table(sub_data$y, sub_data[[split_var]])
      ent_res <- entropy_overlap_by_column(tbl)
      odds_vec <- apply(tbl, 2, function(v) ifelse(v[1] == 0, Inf, v[2]/v[1]))
      
      df <- data.frame(Zk = Zk, Category = cat_val, Splitter = split_var, Subcategory = colnames(tbl),
                       Overlap = ent_res$Overlap, Odds = odds_vec, stringsAsFactors = FALSE)
      all_results_2nd <- rbind(all_results_2nd, df)
    }
  }
}
all_results_2nd$Label <- paste0(all_results_2nd$Zk, "_", all_results_2nd$Category, "|", all_results_2nd$Splitter, "_", all_results_2nd$Subcategory)
build_binary_matrix <- function(data, entries_df) {
  binary_mat <- matrix(0, nrow = nrow(data), ncol = nrow(entries_df))
  colnames(binary_mat) <- entries_df$Label
  for (i in seq_len(nrow(entries_df))) {
    row_idx <- which(data[[entries_df$Zk[i]]] == entries_df$Category[i] &
                       data[[entries_df$Splitter[i]]] == entries_df$Subcategory[i])
    binary_mat[row_idx, i] <- 1
  }
  return(binary_mat)
}
selected_2nd <- all_results_2nd[all_results_2nd$Overlap < 0.1, ]
binary_mat_2way <- build_binary_matrix(new_data, selected_2nd)



meta_1st <- final_df[, c("Label", "Odds", "Overlap")]
meta_2nd <- all_results_2nd[, c("Label", "Odds", "Overlap")]
all_meta <- rbind(meta_1st, meta_2nd)
com_mat <- cbind(binary_mat_1st, binary_mat_2way)
draw_my_heatmap <- function(mat_source, meta_subset, title_text) {
  valid_labels <- intersect(meta_subset$Label, colnames(mat_source))
  if(length(valid_labels) == 0) {
    return(NULL)
  }
  
  mat_to_plot <- mat_source[, valid_labels, drop = FALSE]
  current_odds <- meta_subset$Odds[match(valid_labels, meta_subset$Label)]
  col_name_colors <- ifelse(current_odds >= 5, "red", ifelse(current_odds <= 0.1, "blue", "black"))
  
  row_labels <- factor(new_data$y, levels = 0:1)
  row_anno <- rowAnnotation(Y = row_labels, col = list(Y = c("0" = "lightblue", "1" = "lightpink")), show_annotation_name = FALSE)
  
  ht <- Heatmap(
    mat_to_plot,
    name = " ",
    col = colorRamp2(c(0, 1), c("white", "black")),
    show_heatmap_legend = FALSE,
    left_annotation = row_anno,
    cluster_rows = TRUE, cluster_columns = TRUE,
    clustering_method_rows = "ward.D2", clustering_method_columns = "ward.D2",
    show_row_names = FALSE,
    column_title = title_text,
    column_names_gp = gpar(col = col_name_colors, fontsize = 8)
  )
  draw(ht)
}
# =========================================================================
# 3 heatmaps by different conditions

# (1) strawberry_binary_heatmap
cond1 <- meta_1st[meta_1st$Overlap <= 0.1, ]
draw_my_heatmap(binary_mat_1st, cond1, "strawberry_binary_heatmap_5&0.01")
# (2) strawberry2nd_heatmap_odds0.01
cond2 <- meta_2nd[meta_2nd$Overlap <= 0.1 & meta_2nd$Odds <= 0.01, ]
draw_my_heatmap(binary_mat_2way, cond2, "strawberry2nd_heatmap_odds0.01")
# (3) strawberry2nd_heatmap_odds5
cond4 <- meta_2nd[meta_2nd$Overlap <= 0.1 & meta_2nd$Odds >= 5, ]
draw_my_heatmap(binary_mat_2way, cond4, "strawberry2nd_heatmap_odds5")


# ========================================== 
#                                            figure 25
# ==========================================
table(Y = new_data[new_data$Z7 == "4_1", ]$y, Z4 = new_data[new_data$Z7 == "4_1", ]$Z4_second_split)
table(Y = new_data[new_data$Z7 == "4_2", ]$y, Z4 = new_data[new_data$Z7 == "4_2", ]$Z4_second_split)
table(Y = new_data[new_data$Z7 == "4_3", ]$y, Z4 = new_data[new_data$Z7 == "4_3", ]$Z4_second_split)
table(Y = new_data[new_data$Z7 == "4_4", ]$y, Z4 = new_data[new_data$Z7 == "4_4", ]$Z4_second_split)
# ========================================== 
# Strawberry-split-Overlap-vs-Odds-order2   figure 26
# ==========================================
library(readxl)
library(tidyverse)
library(magrittr)
library(plotly)
library(RColorBrewer)
library(dplyr)
library(dendextend)

# =========================================================================
first <- read.csv("data/strawberry_1stc_16_n_983.csv")
first$wavelenght <- seq(903.187, 1802.564, length.out=234)

second <- read.csv("data/strawberry_2ndc_16_n_983.csv")
second$wavelenght <- seq(903.187, 1798.704, length.out=233)

new_data <- read.csv("data/MIR_Fruit_purees.csv", header = TRUE)
rownames(new_data) <- new_data[, 1]
col_names <- colnames(new_data[, -1])
labels <- ifelse(startsWith(col_names, "Strawberry"), 1, 0)
new_data <- as.data.frame(t(new_data[,-1]))
new_data$y <- labels

# =========================================================================
get_clusters <- function(src_data, rows, k_val, plot_title = "") {
  # Extract the corresponding rows, drop column 984, convert to data frame
  data_list <- lapply(rows, function(r) as.numeric(src_data[r, ][-984]))
  temp_data <- as.data.frame(data_list)
  colnames(temp_data) <- paste0("X", seq_along(rows))
  
  hc_tree <- hclust(dist(temp_data), method = "ward.D2")
  dendrogram_obj <- as.dendrogram(hc_tree)
  
  if (plot_title != "") {
    dendrogram_obj %>%
      dendextend::color_branches(k = k_val, groupLabels = TRUE) %>%
      plot(main = plot_title)
  }
  
  # Cut the dendrogram and map labels back to the correct order using data_order
  cut_tree <- dendextend::cutree(dendrogram_obj, k = k_val, order_clusters_as_data = FALSE)
  data_order <- order.dendrogram(dendrogram_obj)
  
  result_vec <- rep(NA, length(data_order))
  result_vec[data_order] <- cut_tree
  
  return(result_vec)
}

# =========================================================================
new_data$Z1 <- get_clusters(first, 19:25, 9, "Z1")
new_data$Z2 <- get_clusters(first, 21:34, 10, "Z2")

# Z3: keep manual naming
new_data$Z3 <- get_clusters(first, 36:41, 4, "Z3")
Z3_20 <- get_clusters(first, 36:41, 20)
new_data$Z3 <- ifelse(Z3_20 %in% c(13, 14, 15, 17, 18), paste0(new_data$Z3, "_", Z3_20), new_data$Z3)
new_data$Z3[new_data$Z3 == "4_13"] <- '4_1'
new_data$Z3[new_data$Z3 == "4_14"] <- '4_2'
new_data$Z3[new_data$Z3 == "4_15"] <- '4_3'
new_data$Z3[new_data$Z3 == "4_17"] <- '4_5'
new_data$Z3[new_data$Z3 == "4_18"] <- '4_6'

Z3_12 <- get_clusters(first, 36:41, 12)
new_data$Z3 <- ifelse(Z3_12 %in% c(10, 12), paste0(new_data$Z3, "_", Z3_12), new_data$Z3)
new_data$Z3[new_data$Z3 == "4_10"] <- '4_4'
new_data$Z3[new_data$Z3 == "4_12"] <- '4_7'

new_data$Z4 <- get_clusters(first, 42:49, 10, "Z4")
new_data$Z5 <- get_clusters(first, 63:68, 9, "Z5")

# Z6: keep manual naming
new_data$Z6 <- get_clusters(first, 69:81, 10, "Z6")
Z6_20 <- get_clusters(first, 69:81, 20)
new_data$Z6 <- ifelse(Z6_20 %in% c(3, 4, 14, 15, 16, 17), paste0(new_data$Z6, "_", Z6_20), new_data$Z6)
new_data$Z6[new_data$Z6 == "2_3"]  <- '2_1'
new_data$Z6[new_data$Z6 == "2_4"]  <- '2_2'
new_data$Z6[new_data$Z6 == "7_14"] <- '7_1'
new_data$Z6[new_data$Z6 == "7_15"] <- '7_2'
new_data$Z6[new_data$Z6 == "8_16"] <- '8_1'
new_data$Z6[new_data$Z6 == "8_17"] <- '8_2'

# Z7: keep manual naming
new_data$Z7 <- get_clusters(first, 85:95, 4, "Z7")
Z7_20 <- get_clusters(first, 85:95, 20)
new_data$Z7 <- ifelse(Z7_20 %in% c(12, 13, 14, 15, 18, 19, 20), paste0(new_data$Z7, "_", Z7_20), new_data$Z7)
new_data$Z7[new_data$Z7 == "3_12"] <- '3_1'
new_data$Z7[new_data$Z7 == "3_13"] <- '3_2'
new_data$Z7[new_data$Z7 == "3_14"] <- '3_3'
new_data$Z7[new_data$Z7 == "3_15"] <- '3_4'
new_data$Z7[new_data$Z7 == "4_18"] <- '4_2'
new_data$Z7[new_data$Z7 == "4_19"] <- '4_3'
new_data$Z7[new_data$Z7 == "4_20"] <- '4_4'
new_data$Z7[new_data$Z7 == "4"]    <- '4_1'

new_data$Z8 <- get_clusters(first, 208:224, 14, "Z8")

# 2nd-order series
new_data$Z2_second <- get_clusters(second, 21:27, 12, "Z2_second")

new_data$Z3_second <- get_clusters(second, 33:36, 9, "Z3_second")
Z3_sec_20 <- get_clusters(second, 33:36, 20)
new_data$Z3_second <- ifelse(Z3_sec_20 %in% c(15, 16, 19, 20), paste0(new_data$Z3_second, "_", Z3_sec_20), new_data$Z3_second)
new_data$Z3_second[new_data$Z3_second == "7_15"] <- '7_1'
new_data$Z3_second[new_data$Z3_second == "7_16"] <- '7_2'
new_data$Z3_second[new_data$Z3_second == "9_19"] <- '9_1'
new_data$Z3_second[new_data$Z3_second == "9_20"] <- '9_2'

new_data$Z4_second <- get_clusters(second, 61:66, 10, "Z4_second")
Z4_sec_14 <- get_clusters(second, 61:66, 14)
new_data$Z4_second <- ifelse(Z4_sec_14 %in% c(2, 5, 6), paste0(new_data$Z4_second, "_", Z4_sec_14), new_data$Z4_second)
new_data$Z4_second[new_data$Z4_second == "2_2"] <- '2_1'
new_data$Z4_second[new_data$Z4_second == "3_5"] <- '3_3'
new_data$Z4_second[new_data$Z4_second == "3_6"] <- '3_4'

Z4_sec_25 <- get_clusters(second, 61:66, 25)
new_data$Z4_second <- ifelse(Z4_sec_25 %in% c(1, 2, 4, 5, 6, 7), paste0(new_data$Z4_second, "_", Z4_sec_25), new_data$Z4_second)
new_data$Z4_second[new_data$Z4_second == "1_1"] <- '1_1'
new_data$Z4_second[new_data$Z4_second == "1_2"] <- '1_2'
new_data$Z4_second[new_data$Z4_second == "2_4"] <- '2_2'
new_data$Z4_second[new_data$Z4_second == "2_5"] <- '2_3'
new_data$Z4_second[new_data$Z4_second == "3_6"] <- '3_1'
new_data$Z4_second[new_data$Z4_second == "3_7"] <- '3_2'

# =========================================================================
# qualified_categories
feature_vars <- c("Z1", "Z2", "Z3", "Z4", "Z5", "Z6", "Z7", "Z8", "Z2_second", "Z3_second", "Z4_second")
qualified_categories <- list()

for (var in feature_vars) {
  tab <- table(new_data$y, new_data[[var]])
  if (!all(c(0,1) %in% rownames(tab))) {
    tab <- rbind(Y0 = tab["0", , drop=FALSE], Y1 = tab["1", , drop=FALSE])
  }
  rownames(tab) <- c("Y0", "Y1")
  
  selected <- c()
  for (cat in colnames(tab)) {
    n0 <- tab["Y0", cat]
    n1 <- tab["Y1", cat]
    odds <- ifelse(n0 == 0, Inf, n1 / n0)
    
    if ((n0 + n1) > 30 && min(n0, n1) > 10 && odds >= 0.1 && odds <= 10) {
      selected <- c(selected, cat)
    }
  }
  qualified_categories[[var]] <- selected
}

# =========================================================================
# Create _split variables
config_list <- list(
  Z1_split = list(source = "first", rows = 19:25, k = 4),
  Z2_split = list(source = "first", rows = 21:34, k = 4),
  Z3_split = list(source = "first", rows = 36:41, k = 4),
  Z4_split = list(source = "first", rows = 42:49, k = 4),
  Z5_split = list(source = "first", rows = 63:68, k = 4),
  Z6_split = list(source = "first", rows = 69:81, k = 5),
  Z7_split = list(source = "first", rows = 85:95, k = 4),
  Z8_split = list(source = "first", rows = 208:224, k = 4),
  Z2_second_split = list(source = "second", rows = 21:27, k = 5),
  Z3_second_split = list(source = "second", rows = 33:36, k = 5),
  Z4_second_split = list(source = "second", rows = 61:66, k = 5)
)

for (var_name in names(config_list)) {
  cfg <- config_list[[var_name]]
  src_data <- if(cfg$source == "first") first else second
  
  new_data[[var_name]] <- get_clusters(src_data, cfg$rows, cfg$k, paste("Dendrogram for", var_name))
}

# Plot
entropy_overlap_by_column <- function(contingency_table, n_sim = 1000, prefix = "Zcat") {
  n_row <- nrow(contingency_table)
  n_col <- ncol(contingency_table)
  
  row_tot   <- rowSums(contingency_table)
  col_tot   <- colSums(contingency_table)
  grand_tot <- sum(contingency_table)
  
  row_prob_alt <- lapply(1:n_row, function(i) {
    if (row_tot[i] == 0) rep(0, n_col) else contingency_table[i, ] / row_tot[i]
  })
  col_prob_null <- col_tot / grand_tot
  
  alt_array  <- array(0, dim = c(n_row, n_col, n_sim))
  null_array <- array(0, dim = c(n_row, n_col, n_sim))
  
  set.seed(42)
  for (s in 1:n_sim) {
    for (i in 1:n_row) {
      alt_array[i, , s]  <- rmultinom(1, row_tot[i], row_prob_alt[[i]])
      null_array[i, , s] <- rmultinom(1, row_tot[i], col_prob_null)
    }
  }
  
  entropy <- function(p) { p <- p[p > 0]; -sum(p * log(p)) }
  breaks_fix <- seq(0, log(n_row), length.out = 50)
  
  result_df <- data.frame(Category = colnames(contingency_table), Overlap = NA_real_, stringsAsFactors = FALSE)
  
  for (k in 1:n_col) {
    alt_ent  <- apply(alt_array[ , k, ], 2, function(v) entropy(v / sum(v)))
    null_ent <- apply(null_array[, k, ], 2, function(v) entropy(v / sum(v)))
    
    alt_hist  <- hist(alt_ent,  breaks = breaks_fix, plot = FALSE)
    null_hist <- hist(null_ent, breaks = breaks_fix, plot = FALSE)
    overlap_ratio <- sum(pmin(alt_hist$counts, null_hist$counts)) / n_sim
    result_df$Overlap[k] <- overlap_ratio
  }
  return(result_df)
}

analyze_category_splits <- function(new_data, qualified_categories, split_vars, prefix_dir = "entropy_plots") {
  if (!dir.exists(prefix_dir)) dir.create(prefix_dir)
  all_results <- data.frame()
  
  for (Zk in names(qualified_categories)) {
    categories <- qualified_categories[[Zk]]
    Zk_split_name <- paste0(Zk, "_split")
    
    for (cat_val in categories) {
      subset_idx <- which(new_data[[Zk]] == cat_val)
      sub_data <- new_data[subset_idx, ]
      
      for (split_var in split_vars) {
        if (split_var == Zk_split_name) next
        
        tbl <- table(sub_data$y, sub_data[[split_var]])
        file_prefix <- paste0(Zk, "_", make.names(cat_val), "_by_", split_var)
        safe_prefix <- file.path(prefix_dir, file_prefix)
        ent_res <- entropy_overlap_by_column(tbl, prefix = safe_prefix)
        
        odds_vec <- apply(tbl, 2, function(v) ifelse(v[1] == 0, Inf, v[2]/v[1]))
        
        df <- data.frame(Zk = Zk, Category = cat_val, Splitter = split_var, Subcategory = colnames(tbl),
                         Overlap = ent_res$Overlap, Odds = odds_vec, stringsAsFactors = FALSE)
        all_results <- rbind(all_results, df)
      }
    }
  }
  return(all_results)
}

plot_overlap_vs_odds <- function(result_df, pdf_name = "Overlap_vs_Odds.pdf") {
  library(ggplot2)
  result_df$Highlight <- result_df$Overlap < 0.1
  result_df$Highlight_Label <- ifelse(result_df$Highlight, "Overlap < 0.1", "Overlap ≥ 0.1")
  
  p <- ggplot(result_df, aes(x = Odds, y = Overlap, label = ifelse(Highlight, Subcategory, ""))) +
    geom_point(aes(color = Highlight_Label)) +
    scale_color_manual(values = c("Overlap ≥ 0.1" = "black", "Overlap < 0.1" = "red")) +
    labs(title = "Overlap vs Odds (highlighted if Overlap < 0.1)", x = "Odds (Y=1/Y=0)", y = "Entropy Overlap", color = "") +
    theme_minimal()
  print(p)
}

result_df <- analyze_category_splits(new_data, qualified_categories, split_vars = paste0(c("Z1","Z2","Z3","Z4","Z5","Z6","Z7","Z8","Z2_second","Z3_second","Z4_second"), "_split"))
plot_overlap_vs_odds(result_df)


# ========================================== 
# figure 28 29 31 32
# ==========================================
library(ggplot2)
library(dplyr)
library(tidyr)
library(dendextend)
library(ComplexHeatmap)
library(circlize)
library(writexl)
library(readxl)
first <- read.csv("data/strawberry_1stc_16_n_983.csv")
first$wavelenght <- seq(903.187, 1802.564, length.out=234)

second <- read.csv("data/strawberry_2ndc_16_n_983.csv")
second$wavelenght <- seq(903.187, 1798.704, length.out=233)

new_data <- read.csv("data/MIR_Fruit_purees.csv", header = TRUE)
rownames(new_data) <- new_data[, 1]
col_names <- colnames(new_data[, -1])
labels <- ifelse(startsWith(col_names, "Strawberry"), 1, 0)
new_data <- as.data.frame(t(new_data[,-1]))
new_data$y <- labels

get_clusters <- function(src_data, rows, k_val) {
  data_list <- lapply(rows, function(r) as.numeric(src_data[r, ][-984]))
  temp_data <- as.data.frame(data_list)
  colnames(temp_data) <- paste0("X", seq_along(rows))
  hc_tree <- hclust(dist(temp_data), method = "ward.D2")
  dendrogram_obj <- as.dendrogram(hc_tree)
  cut_tree <- dendextend::cutree(dendrogram_obj, k = k_val, order_clusters_as_data = FALSE)
  data_order <- order.dendrogram(dendrogram_obj)
  result_vec <- rep(NA, length(data_order))
  result_vec[data_order] <- cut_tree
  return(result_vec)
}

# =========================================================================
# Fine-grained clustering Z1 ~ Z4_second
# =========================================================================
new_data$Z1 <- get_clusters(first, 19:25, 9); new_data$Z2 <- get_clusters(first, 21:34, 10)
new_data$Z3 <- get_clusters(first, 36:41, 4)
Z3_20 <- get_clusters(first, 36:41, 20)
new_data$Z3 <- ifelse(Z3_20 %in% c(13, 14, 15, 17, 18), paste0(new_data$Z3, "_", Z3_20), new_data$Z3)
new_data$Z3[new_data$Z3 == "4_13"] <- '4_1'; new_data$Z3[new_data$Z3 == "4_14"] <- '4_2'; new_data$Z3[new_data$Z3 == "4_15"] <- '4_3'
new_data$Z3[new_data$Z3 == "4_17"] <- '4_5'; new_data$Z3[new_data$Z3 == "4_18"] <- '4_6'
Z3_12 <- get_clusters(first, 36:41, 12)
new_data$Z3 <- ifelse(Z3_12 %in% c(10, 12), paste0(new_data$Z3, "_", Z3_12), new_data$Z3)
new_data$Z3[new_data$Z3 == "4_10"] <- '4_4'; new_data$Z3[new_data$Z3 == "4_12"] <- '4_7'

new_data$Z4 <- get_clusters(first, 42:49, 10); new_data$Z5 <- get_clusters(first, 63:68, 9)
new_data$Z6 <- get_clusters(first, 69:81, 10)
Z6_20 <- get_clusters(first, 69:81, 20)
new_data$Z6 <- ifelse(Z6_20 %in% c(3, 4, 14, 15, 16, 17), paste0(new_data$Z6, "_", Z6_20), new_data$Z6)
new_data$Z6[new_data$Z6 == "2_3"]  <- '2_1'; new_data$Z6[new_data$Z6 == "2_4"]  <- '2_2'; new_data$Z6[new_data$Z6 == "7_14"] <- '7_1'
new_data$Z6[new_data$Z6 == "7_15"] <- '7_2'; new_data$Z6[new_data$Z6 == "8_16"] <- '8_1'; new_data$Z6[new_data$Z6 == "8_17"] <- '8_2'

new_data$Z7 <- get_clusters(first, 85:95, 4)
Z7_20 <- get_clusters(first, 85:95, 20)
new_data$Z7 <- ifelse(Z7_20 %in% c(12, 13, 14, 15, 18, 19, 20), paste0(new_data$Z7, "_", Z7_20), new_data$Z7)
new_data$Z7[new_data$Z7 == "3_12"] <- '3_1'; new_data$Z7[new_data$Z7 == "3_13"] <- '3_2'; new_data$Z7[new_data$Z7 == "3_14"] <- '3_3'
new_data$Z7[new_data$Z7 == "3_15"] <- '3_4'; new_data$Z7[new_data$Z7 == "4_18"] <- '4_2'; new_data$Z7[new_data$Z7 == "4_19"] <- '4_3'
new_data$Z7[new_data$Z7 == "4_20"] <- '4_4'; new_data$Z7[new_data$Z7 == "4"]    <- '4_1'

new_data$Z8 <- get_clusters(first, 208:224, 14)

new_data$Z2_second <- get_clusters(second, 21:27, 12)
new_data$Z3_second <- get_clusters(second, 33:36, 9)
Z3_sec_20 <- get_clusters(second, 33:36, 20)
new_data$Z3_second <- ifelse(Z3_sec_20 %in% c(15, 16, 19, 20), paste0(new_data$Z3_second, "_", Z3_sec_20), new_data$Z3_second)
new_data$Z3_second[new_data$Z3_second == "7_15"] <- '7_1'; new_data$Z3_second[new_data$Z3_second == "7_16"] <- '7_2'
new_data$Z3_second[new_data$Z3_second == "9_19"] <- '9_1'; new_data$Z3_second[new_data$Z3_second == "9_20"] <- '9_2'

new_data$Z4_second <- get_clusters(second, 61:66, 10)
Z4_sec_14 <- get_clusters(second, 61:66, 14)
new_data$Z4_second <- ifelse(Z4_sec_14 %in% c(2, 5, 6), paste0(new_data$Z4_second, "_", Z4_sec_14), new_data$Z4_second)
new_data$Z4_second[new_data$Z4_second == "2_2"] <- '2_1'; new_data$Z4_second[new_data$Z4_second == "3_5"] <- '3_3'; new_data$Z4_second[new_data$Z4_second == "3_6"] <- '3_4'
Z4_sec_25 <- get_clusters(second, 61:66, 25)
new_data$Z4_second <- ifelse(Z4_sec_25 %in% c(1, 2, 4, 5, 6, 7), paste0(new_data$Z4_second, "_", Z4_sec_25), new_data$Z4_second)
new_data$Z4_second[new_data$Z4_second == "1_1"] <- '1_1'; new_data$Z4_second[new_data$Z4_second == "1_2"] <- '1_2'; new_data$Z4_second[new_data$Z4_second == "2_4"] <- '2_2'
new_data$Z4_second[new_data$Z4_second == "2_5"] <- '2_3'; new_data$Z4_second[new_data$Z4_second == "3_6"] <- '3_1'; new_data$Z4_second[new_data$Z4_second == "3_7"] <- '3_2'

# =========================================================================
# Create _split variables
# =========================================================================
config_list <- list(
  Z1_split = list(src = first, rows = 19:25, k = 4),   Z2_split = list(src = first, rows = 21:34, k = 4),
  Z3_split = list(src = first, rows = 36:41, k = 4),   Z4_split = list(src = first, rows = 42:49, k = 4),
  Z5_split = list(src = first, rows = 63:68, k = 4),   Z6_split = list(src = first, rows = 69:81, k = 5),
  Z7_split = list(src = first, rows = 85:95, k = 4),   Z8_split = list(src = first, rows = 208:224, k = 4),
  Z2_second_split = list(src = second, rows = 21:27, k = 5), Z3_second_split = list(src = second, rows = 33:36, k = 5),
  Z4_second_split = list(src = second, rows = 61:66, k = 5)
)
for (var_name in names(config_list)) {
  cfg <- config_list[[var_name]]
  new_data[[var_name]] <- get_clusters(cfg$src, cfg$rows, cfg$k)
}

# =========================================================================
# Entropy function
# =========================================================================
entropy_overlap_by_column <- function(contingency_table, n_sim = 1000) {
  n_row <- nrow(contingency_table)
  n_col <- ncol(contingency_table)
  row_tot <- rowSums(contingency_table)
  col_tot <- colSums(contingency_table)
  grand_tot <- sum(contingency_table)
  if (n_row == 0 || n_col == 0 || grand_tot == 0) {
    return(data.frame(Category = colnames(contingency_table), Overlap = NA_real_, stringsAsFactors = FALSE))
  }
  row_prob_alt <- lapply(1:n_row, function(i) {
    if (row_tot[i] == 0) rep(0, n_col) else contingency_table[i, ] / row_tot[i]
  })
  col_prob_null <- col_tot / grand_tot
  alt_array  <- array(0, dim = c(n_row, n_col, n_sim))
  null_array <- array(0, dim = c(n_row, n_col, n_sim))
  set.seed(42)
  for (s in 1:n_sim) {
    for (i in 1:n_row) {
      if (row_tot[i] == 0) next
      alt_array[i, , s]  <- rmultinom(1, row_tot[i], row_prob_alt[[i]])
      null_array[i, , s] <- rmultinom(1, row_tot[i], col_prob_null)
    }
  }
  entropy <- function(p) { p <- p[p > 0]; -sum(p * log(p)) }
  breaks_fix <- seq(0, log(n_row), length.out = 50) 
  result_df <- data.frame(Category = colnames(contingency_table), Overlap = NA_real_, stringsAsFactors = FALSE)
  for (k in 1:n_col) {
    alt_ent  <- apply(alt_array[ , k, , drop = FALSE], 3, function(v) { if(sum(v) == 0) 0 else entropy(v / sum(v)) })
    null_ent <- apply(null_array[, k, , drop = FALSE], 3, function(v) { if(sum(v) == 0) 0 else entropy(v / sum(v)) })
    alt_hist  <- hist(alt_ent,  breaks = breaks_fix, plot = FALSE)
    null_hist <- hist(null_ent, breaks = breaks_fix, plot = FALSE)
    result_df$Overlap[k] <- sum(pmin(alt_hist$counts, null_hist$counts)) / n_sim
  }
  return(result_df)
}

# =========================================================================
# Compute 1st-order entropy and matrix
# =========================================================================
feature_vars <- c("Z1", "Z2", "Z3", "Z4", "Z5", "Z6", "Z7", "Z8", "Z2_second", "Z3_second", "Z4_second")
split_vars <- paste0(feature_vars, "_split")

all_entropy_results <- list()
onehot_list <- list()
for (var_name in feature_vars) {
  tbl <- table(new_data$y, new_data[[var_name]])
  if (!all(c("0","1") %in% rownames(tbl))) {
    tbl <- rbind(Y0 = tbl["0", , drop=FALSE], Y1 = tbl["1", , drop=FALSE])
  } else {
    rownames(tbl) <- c("Y0", "Y1")
  }
  odds_vec <- ifelse(tbl["Y0", ] == 0, Inf, tbl["Y1", ] / tbl["Y0", ])
  ent_res <- entropy_overlap_by_column(tbl)
  ent_res$Odds <- odds_vec[match(ent_res$Category, colnames(tbl))]
  ent_res$Feature <- var_name
  all_entropy_results[[var_name]] <- ent_res
  
  tmp_factor <- factor(new_data[[var_name]])
  mat_tmp <- model.matrix(~ tmp_factor - 1)
  clean_cat <- sub("tmp_factor", "", colnames(mat_tmp))
  colnames(mat_tmp) <- paste0(var_name, "_", clean_cat)
  onehot_list[[var_name]] <- mat_tmp
}
final_df <- bind_rows(all_entropy_results)
final_df$Label <- paste0(final_df$Feature, "_", final_df$Category)
binary_mat_1st <- do.call(cbind, onehot_list)
meta_1st <- final_df[, c("Label", "Odds", "Overlap")]

# =========================================================================
# Filter qualified_categories
# =========================================================================
qualified_categories <- list()
for (var in feature_vars) {
  tab <- table(new_data$y, new_data[[var]])
  if (!all(c(0,1) %in% rownames(tab))) {
    tab <- rbind(Y0 = tab["0", , drop=FALSE], Y1 = tab["1", , drop=FALSE])
  }
  rownames(tab) <- c("Y0", "Y1")
  
  selected <- c()
  for (cat in colnames(tab)) {
    n0 <- tab["Y0", cat]
    n1 <- tab["Y1", cat]
    odds <- ifelse(n0 == 0, Inf, n1 / n0)
    if ((n0 + n1) > 30 && min(n0, n1) > 10 && odds >= 0.1 && odds <= 10) {
      selected <- c(selected, cat)
    }
  }
  qualified_categories[[var]] <- selected
}

# =========================================================================
# Compute 2nd-order entropy and build the combined matrix
# =========================================================================
save_path <- "output/dynamic_2nd_split_results.xlsx"

if (file.exists(save_path)) {
  result_df_2nd <- read_excel(save_path)
} else {
  all_results <- data.frame()
  for (Zk in names(qualified_categories)) {
    categories <- qualified_categories[[Zk]]
    Zk_split_name <- paste0(Zk, "_split")
    for (cat_val in categories) {
      subset_idx <- which(new_data[[Zk]] == cat_val)
      if(length(subset_idx) == 0) next
      sub_data <- new_data[subset_idx, ]
      for (split_var in split_vars) {
        if (split_var == Zk_split_name) next 
        tbl <- table(factor(sub_data$y, levels = c(0, 1)), sub_data[[split_var]])
        if(ncol(tbl) == 0) next
        ent_res <- entropy_overlap_by_column(tbl)
        odds_vec <- apply(tbl, 2, function(v) ifelse(v[1] == 0, Inf, v[2]/v[1]))
        df <- data.frame(
          Zk = Zk, Category = cat_val, Splitter = split_var, Subcategory = colnames(tbl),
          Overlap = ent_res$Overlap, Odds = odds_vec, stringsAsFactors = FALSE
        )
        all_results <- rbind(all_results, df)
      }
    }
  }
  if(!dir.exists(dirname(save_path))) dir.create(dirname(save_path), recursive = TRUE)
  write_xlsx(all_results, path = save_path)
  cat("Computation complete. File saved to:", save_path, "\n")
  result_df_2nd <- read_excel(save_path)
}

result_df_2nd$Label <- paste0(result_df_2nd$Zk, "_", result_df_2nd$Category, "|", result_df_2nd$Splitter, "_", result_df_2nd$Subcategory)
meta_2nd <- result_df_2nd[, c("Label", "Odds", "Overlap")]

# Merge 1st- and 2nd-order results
all_meta <- rbind(meta_1st, meta_2nd)

build_binary_matrix <- function(data, entries_df) {
  binary_mat <- matrix(0, nrow = nrow(data), ncol = nrow(entries_df))
  colnames(binary_mat) <- entries_df$Label
  for (i in seq_len(nrow(entries_df))) {
    row_idx <- which(data[[entries_df$Zk[i]]] == entries_df$Category[i] &
                       data[[entries_df$Splitter[i]]] == entries_df$Subcategory[i])
    binary_mat[row_idx, i] <- 1
  }
  return(binary_mat)
}

# Build com_mat combining 1st + 2nd order
binary_mat_2way <- build_binary_matrix(new_data, result_df_2nd)
com_mat <- cbind(binary_mat_1st, binary_mat_2way)

get_dendro_cluster_indices <- function(dend_obj, k_val, target_cluster) {
  dendro_cut <- dendextend::cutree(dend_obj, k = k_val, order_clusters_as_data = FALSE)
  orig_order <- order.dendrogram(dend_obj)
  cluster_labels <- rep(NA, length(orig_order))
  cluster_labels[orig_order] <- dendro_cut
  return(which(cluster_labels == target_cluster))
}

row_labels <- factor(new_data$y, levels = 0:1)
row_anno <- rowAnnotation(Y = row_labels, col = list(Y = c("0" = "lightblue", "1" = "lightpink")), show_annotation_name = FALSE)

var_ranges <- data.frame(
  Variable = c("Z1", "Z2", "Z3", "Z4", "Z5", "Z6", "Z7", "Z8", "Z2_second", "Z3_second", "Z4_second"),
  xmin = c(972, 980, 1038, 1060, 1142, 1165, 1227, 1700, 980, 1026, 1134),
  xmax = c(995, 1030, 1057, 1088, 1161, 1211, 1266, 1763, 1003, 1038, 1154)
)

# =========================================================================
# strawberry2nd-heatmap-odds-5-001.pdf  4-panels-strawberry2nd-heatmap-odds5-001 figure 31 32
# =========================================================================

selected_split_1 <- result_df_2nd[result_df_2nd$Overlap < 0.1 & (result_df_2nd$Odds >= 5 | result_df_2nd$Odds <= 0.01), ]
binary_mat_1 <- build_binary_matrix(new_data, selected_split_1)
col_name_colors_1 <- ifelse(selected_split_1$Odds >= 5, "red", "blue")
names(col_name_colors_1) <- selected_split_1$Label

row_tree_1 <- hclust(dist(binary_mat_1), method = "ward.D2")
dendrogram_obj_1 <- as.dendrogram(row_tree_1)

block1_idx_1 <- get_dendro_cluster_indices(dendrogram_obj_1, k = 2, target_cluster = 1)
block2_idx_1 <- get_dendro_cluster_indices(dendrogram_obj_1, k = 3, target_cluster = 2)

ht_1 <- Heatmap(binary_mat_1, name = " ", col = colorRamp2(c(0, 1), c("white", "black")),
                show_heatmap_legend = FALSE, left_annotation = row_anno, cluster_rows = TRUE, cluster_columns = TRUE,
                clustering_method_rows = "ward.D2", clustering_method_columns = "ward.D2", show_row_names = FALSE,
                column_names_gp = gpar(col = col_name_colors_1, fontsize = 5), row_dend_width = unit(3, "cm"))

ht_obj_1 <- draw(ht_1)
row_dend_1 <- row_dend(ht_obj_1)
block3_idx_1 <- get_dendro_cluster_indices(row_dend_1, k = 11, target_cluster = 7)
block4_idx_1 <- get_dendro_cluster_indices(row_dend_1, k = 11, target_cluster = 10)

block_list_1 <- list(block1_idx_1, block2_idx_1, block3_idx_1, block4_idx_1)
block_names <- paste0("Block ", LETTERS[1:4])
long_data_1 <- data.frame()
for (i in seq_along(block_list_1)) {
  idx <- block_list_1[[i]]
  block_data <- new_data[idx, ]
  block_data$SampleID <- paste0("S", idx)
  block_data$Block <- block_names[i]
  long_block <- block_data %>%
    dplyr::select(-starts_with("Z")) %>% 
    mutate(y = factor(y)) %>%
    pivot_longer(cols = -c(y, SampleID, Block), names_to = "X", values_to = "Y")
  long_data_1 <- bind_rows(long_data_1, long_block)
}
long_data_1$X <- as.numeric(sub("^X", "", as.character(long_data_1$X)))

block_vars_1 <- list(
  "Block A" = c("Z1", "Z2", "Z3", "Z4", "Z5", "Z6", "Z7", "Z2_second", "Z3_second", "Z4_second"),
  "Block B" = c("Z1", "Z2", "Z3", "Z4", "Z5", "Z6", "Z7", "Z8", "Z2_second", "Z3_second", "Z4_second"),
  "Block C" = c("Z1", "Z2", "Z3", "Z4", "Z5", "Z6", "Z8", "Z3_second", "Z4_second"),
  "Block D" = c("Z1", "Z2", "Z3", "Z4", "Z5", "Z6", "Z7", "Z8", "Z2_second", "Z3_second", "Z4_second")
)
rect_df_1 <- bind_rows(lapply(names(block_vars_1), function(block_name) {
  data.frame(Block = block_name, Variable = block_vars_1[[block_name]])
})) %>% left_join(var_ranges, by = "Variable") %>%
  mutate(ymin = -Inf, ymax = Inf, xmid = (xmin + xmax) / 2, ymid = Inf)

p1 <- ggplot(long_data_1, aes(x = X, y = Y, group = SampleID, color = y)) +
  geom_rect(data = rect_df_1, inherit.aes = FALSE,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = Variable),
            alpha = 0.2, show.legend = FALSE) +
  geom_text(data = rect_df_1, inherit.aes = FALSE,
            aes(x = xmid, y = ymid, label = Variable),
            angle = 90, vjust = -0.3, size = 1.5, color = "black") +
  geom_line(alpha = 0.7) + facet_wrap(~ Block, ncol = 2) +
  scale_color_manual(values = c("0" = "lightblue", "1" = "lightpink")) +
  labs(x = "Wavelength", y = "", color = "Class") + theme_minimal(base_size = 14) + coord_cartesian(clip = "off") 
print(p1)


# =========================================================================
# straw-combine1st2nd-odds001.pdf    4-panels-straw-combine1st2nd-odds001  figure 28 29
# =========================================================================
cond3 <- all_meta[all_meta$Overlap <= 0.1 & all_meta$Odds <= 0.01, ]
valid_labels <- intersect(cond3$Label, colnames(com_mat))

binary_mat_2 <- com_mat[, valid_labels, drop = FALSE]


col_name_colors_2 <- ifelse(cond3$Odds[match(valid_labels, cond3$Label)] >= 5, "red",
                            ifelse(cond3$Odds[match(valid_labels, cond3$Label)] <= 0.1, "blue", "black"))
names(col_name_colors_2) <- valid_labels

row_tree_2 <- hclust(dist(binary_mat_2), method = "ward.D2")
dendrogram_obj_2 <- as.dendrogram(row_tree_2)


block1_idx_2 <- get_dendro_cluster_indices(dendrogram_obj_2, k = 2, target_cluster = 1)
block2_idx_2 <- get_dendro_cluster_indices(dendrogram_obj_2, k = 3, target_cluster = 2)

ht_2 <- Heatmap(binary_mat_2, name = " ", col = colorRamp2(c(0, 1), c("white", "black")),
                show_heatmap_legend = FALSE, left_annotation = row_anno, cluster_rows = TRUE, cluster_columns = TRUE,
                clustering_method_rows = "ward.D2", clustering_method_columns = "ward.D2", show_row_names = FALSE,
                column_names_gp = gpar(col = col_name_colors_2, fontsize = 5), row_dend_width = unit(3, "cm"))

ht_obj_2 <- draw(ht_2)
row_dend_2 <- row_dend(ht_obj_2)
block3_idx_2 <- get_dendro_cluster_indices(row_dend_2, k = 11, target_cluster = 7)
block4_idx_2 <- get_dendro_cluster_indices(row_dend_2, k = 11, target_cluster = 10)

block_list_2 <- list(block1_idx_2, block2_idx_2, block3_idx_2, block4_idx_2)
long_data_2 <- data.frame()
for (i in seq_along(block_list_2)) {
  idx <- block_list_2[[i]]
  block_data <- new_data[idx, ]
  block_data$SampleID <- paste0("S", idx)
  block_data$Block <- block_names[i]
  long_block <- block_data %>%
    dplyr::select(-starts_with("Z")) %>% 
    mutate(y = factor(y)) %>%
    pivot_longer(cols = -c(y, SampleID, Block), names_to = "X", values_to = "Y")
  long_data_2 <- bind_rows(long_data_2, long_block)
}
long_data_2$X <- as.numeric(sub("^X", "", as.character(long_data_2$X)))

block_vars_2 <- list(
  "Block A" = c("Z1", "Z2", "Z3", "Z4", "Z6", "Z2_second", "Z3_second", "Z4_second"),
  "Block B" = c("Z1", "Z2", "Z3", "Z4", "Z5", "Z6", "Z7", "Z8", "Z2_second", "Z3_second", "Z4_second"),
  "Block C" = c("Z1", "Z2", "Z3", "Z4", "Z5", "Z6", "Z8", "Z3_second", "Z4_second"),
  "Block D" = c("Z1", "Z2", "Z3", "Z4", "Z5", "Z6", "Z7", "Z8", "Z2_second", "Z3_second", "Z4_second")
)
rect_df_2 <- bind_rows(lapply(names(block_vars_2), function(block_name) {
  data.frame(Block = block_name, Variable = block_vars_2[[block_name]])
})) %>% left_join(var_ranges, by = "Variable") %>%
  mutate(ymin = -Inf, ymax = Inf, xmid = (xmin + xmax) / 2, ymid = Inf)

p2 <- ggplot(long_data_2, aes(x = X, y = Y, group = SampleID, color = y)) +
  geom_rect(data = rect_df_2, inherit.aes = FALSE,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = Variable),
            alpha = 0.2, show.legend = FALSE) +
  geom_text(data = rect_df_2, inherit.aes = FALSE,
            aes(x = xmid, y = ymid, label = Variable),
            angle = 90, vjust = -0.3, size = 1.5, color = "black") +
  geom_line(alpha = 0.7) + facet_wrap(~ Block, ncol = 2) +
  scale_color_manual(values = c("0" = "lightblue", "1" = "lightpink")) +
  labs(x = "Wavelength", y = "", color = "Class") + theme_minimal(base_size = 14) + coord_cartesian(clip = "off") 
print(p2)
