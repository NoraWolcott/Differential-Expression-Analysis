library(DESeq2)
library(ggplot2)
library(EnhancedVolcano)
library(gridExtra)
library(grid)
library(genefilter)
library(gplots)

# -------------------------
# Load model
# -------------------------
dds <- readRDS("/Users/norawolcott/Documents/Datta Lab/data/Biopolymer Facility Seq/260403/Contrast models/dds_v_F_45min.rds")

## Remove all animals that are male or were exposed to odor (2hr timepoint only)
dds$sex <- ifelse(grepl("_F_", dds$group), "F", "M")
dds$treatment <- sub("^(vehicle|cort|E2|P4).*", "\\1", dds$group)
dds$time <- sub(".*_(\\d+hr|\\d+min).*", "\\1", dds$group)

dds$sex <- factor(dds$sex)
dds$treatment <- factor(dds$treatment)
dds$time <- factor(dds$time)

design(dds) <- ~ sex + treatment + time

## Remove male animals and animals with odor exposure (2hr)
dds_subset <- dds[, dds$sex == "F" & !grepl("2hr", dds$time)]

dds_subset$sex <- droplevels(dds_subset$sex)
dds_subset$time <- droplevels(dds_subset$time)
dds_subset$treatment <- droplevels(dds_subset$treatment)

# remove sex because now there's only one level
design(dds_subset) <- ~ treatment + time

dds_subset <- DESeq(dds_subset)

# -------------------------
# PCA
# -------------------------
vsd <- vst(dds_subset, blind = FALSE)

pca_data <- plotPCA(vsd, intgroup = "group", returnData = TRUE)

percentVar <- round(100 * attr(pca_data, "percentVar"))

# -------------------------
# Plot
# -------------------------
p <- ggplot(pca_data, aes(PC1, PC2, color = group)) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal() +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size = 10)
  )

# -------------------------
# Save plot
# -------------------------
ggsave(
  "/Users/norawolcott/Documents/Datta Lab/data/Biopolymer Facility Seq/260403/PCA_plot_v2.png",
  plot = p,
  width = 6,
  height = 5
)

print(p)

# -------------------------
# Plot PCA across variables
# -------------------------

## FIRST: JUST TIME AND TREATMENT (NO ODOR, JUST FEMALES)

p_time <- ggplot(pca_data, aes(PC1, PC2, color = time)) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal() +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size = 10)
  )

ggsave(
  "/Users/norawolcott/Documents/Datta Lab/data/Biopolymer Facility Seq/260403/PCA_plot_time.png",
  plot = p_time,
  width = 6,
  height = 5
)
print(p_time)

p_treatment <- ggplot(pca_data, aes(PC1, PC2, color = treatment)) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal() +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size = 10)
  )

ggsave(
  "/Users/norawolcott/Documents/Datta Lab/data/Biopolymer Facility Seq/260403/PCA_plot_treatment.png",
  plot = p_treatment,
  width = 6,
  height = 5
)

print(p_treatment)

# Extract PCA data (if not already done)
# vsd <- vst(dds, blind=FALSE)
# pca_data <- plotPCA(vsd, intgroup="group", returnData=TRUE)
# percentVar <- round(100 * attr(pca_data, "percentVar"))

# # -------------------------
# # Parse biological variables from group
# # -------------------------

# # Example group format: cort_M_45min or vehicle_M_45min_Rpl22CrhCre

# pca_data$treatment <- sub("_.*", "", pca_data$group)

# pca_data$sex <- ifelse(grepl("_F_", pca_data$group), "F", "M")

# pca_data$time <- sub(".*_(\\d+.*)", "\\1", pca_data$group)

# pca_data$odor <- ifelse(grepl("_odor$", pca_data$group), "odor", "no_odor")

# pca_data$genotype <- ifelse(grepl("Rpl22CrhCre", pca_data$group),
#                            "Rpl22CrhCre", "WT")

# # -------------------------
# # Define base PCA plot function
# # -------------------------

# make_pca_plot <- function(color_var, title) {
#   ggplot(pca_data, aes_string("PC1", "PC2", color=color_var)) +
#     geom_point(size=4) +
#     xlab(paste0("PC1: ", percentVar[1], "% variance")) +
#     ylab(paste0("PC2: ", percentVar[2], "% variance")) +
#     ggtitle(title) +
#     theme_minimal()
# }

# -------------------------
# Generate plots
# -------------------------

# p_treatment <- make_pca_plot("treatment", "PCA colored by treatment")
# p_sex       <- make_pca_plot("sex", "PCA colored by sex")
# p_time      <- make_pca_plot("time", "PCA colored by time")
# p_genotype  <- make_pca_plot("genotype", "PCA colored by genotype")
# p_odor      <- make_pca_plot("odor", "PCA colored by odor exposure")

# # -------------------------
# # Save plots
# # -------------------------

# outdir <- "/Users/norawolcott/Documents/Datta Lab/data/Biopolymer Facility Seq/260403/"

# # Put all plots into a list
# p_list <- list(
#   p_treatment,
#   p_sex,
#   p_time,
#   p_genotype,
#   p_odor
# )

# # If number of plots are odd, add an empty plot as placeholder
# empty_plot <- grid.rect(gp=gpar(col="white"))  # blank white space
# p_list <- c(p_list, list(empty_plot))

# # Arrange in 3 rows x 2 columns
# combined_plot <- do.call(grid.arrange, c(p_list, ncol=2, nrow=3))

# # Save
# ggsave(filename = paste0(outdir, "PCA_grid_v2.png"),
#        plot = combined_plot,
#        width = 12, height = 8)
