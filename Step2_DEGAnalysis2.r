library(DESeq2)
library(ggplot2)
library(EnhancedVolcano)
library(gridExtra)
library(grid)
library(genefilter)
library(gplots)
library(RColorBrewer)
library(org.Mm.eg.db)
library(AnnotationDbi)

# -------------------------
# 1. Load counts and metadata
# -------------------------
counts_file <- "/Users/norawolcott/Documents/Datta Lab/data/Biopolymer Facility RNASeq/260510/combined_filtered_counts.tsv"
metadata_file <- "/Users/norawolcott/Documents/Datta Lab/data/Biopolymer Facility RNASeq/260510/combined_metadata_for_DE.tsv"

counts <- read.table(counts_file, header = TRUE, row.names = 1, sep = "\t")
metadata <- read.table(metadata_file, header = TRUE, row.names = 1, sep = "\t")

# Make sure metadata rows match counts columns
metadata <- metadata[colnames(counts), , drop = FALSE]

# -------------------------
# 2. Create DESeq2 dataset
# -------------------------
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = metadata,
                              design = ~ group)


# -------------------------
# 3. Run DESeq2
# -------------------------
dds$group <- relevel(dds$group, ref = "vehicle_F_2hr_TMA")
dds <- DESeq(dds)
resultsNames(dds)

saveRDS(dds, file = "/Users/norawolcott/Documents/Datta Lab/data/Biopolymer Facility RNASeq/260510/dds_object.rds")

# -------------------------
# 4. Differential Expression
# -------------------------
# Example: compare 'cort_M_45min' vs 'vehicle_M_45min'
res <- results(dds, contrast = c("group", "P4_F_2hr_TMA", "vehicle_F_2hr_TMA"))
res <- lfcShrink(dds, contrast=c("group", "P4_F_2hr_TMA", "vehicle_F_2hr_TMA"), type="normal")

# Order by adjusted p-value
res <- res[order(res$padj),]

# Save results to file
write.csv(as.data.frame(res), file="/Users/norawolcott/Documents/Datta Lab/data/Biopolymer Facility RNASeq/260510/DEGs_cort_vs_vehicle.csv")

# -------------------------
# 5. PLOT
# -------------------------
# PCA plot
vsd <- vst(dds, blind=FALSE)
pca_data <- plotPCA(vsd, intgroup="group", returnData=TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

# Assign a number to each unique group
pca_data$group_id <- as.numeric(factor(pca_data$group))

# Create a lookup table (for your legend)
group_key <- unique(pca_data[, c("group", "group_id")])
group_key[order(group_key$group_id), ]

p <- ggplot(pca_data, aes(PC1, PC2)) +
  geom_point(aes(color=group_id), size=6) +
  geom_text(aes(label=rownames(pca_data)), color="black", size=3) +
  scale_color_viridis_c() +  # nice gradient
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal()

geom_text(aes(label=group_id), color="white", size=3, fontface="bold")

ggsave("/Users/norawolcott/Documents/Datta Lab/data/Biopolymer Facility RNASeq/260510/PCA_plot.png", plot=p, width=6, height=5)


# -------------------------
# 6. Dendrogram heatmap
# -------------------------

# Use the vst object you already created
mat <- assay(vsd)

# -------------------------
# Select top variable genes
# -------------------------
topVarGenes <- head(order(rowVars(mat), decreasing = TRUE), 100)
mat_top <- mat[topVarGenes, ]

# -------------------------
# Scale by row (gene)
# -------------------------
mat_scaled <- t(scale(t(mat_top)))

# rescale values to widen contrast
mat_scaled <- mat_scaled * 2

# -------------------------
# Convert Ensembl → gene symbols (labels only)
# -------------------------
ensembl_ids <- rownames(mat_scaled)

gene_symbols <- mapIds(
  org.Mm.eg.db,
  keys = ensembl_ids,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

# Use symbols where available, otherwise keep Ensembl ID
row_labels <- ifelse(
  is.na(gene_symbols),
  ensembl_ids,
  gene_symbols
)

# Make labels unique (prevents duplicate label issues)
row_labels <- make.unique(row_labels)

# -------------------------
# Create colors for each group
# -------------------------
group_colors <- setNames(
  rainbow(length(unique(metadata$group))),
  unique(metadata$group)
)

col_side_colors <- group_colors[metadata$group]

# -------------------------
# Open PNG device BEFORE plotting
# -------------------------
png(
  "/Users/norawolcott/Documents/Datta Lab/data/Biopolymer Facility RNASeq/260510/Dendrogram_heatmap.png",
  width = 2000,
  height = 2000,
  res = 300
)

# -------------------------
# Generate heatmap
# -------------------------
heatmap.2(
  mat_scaled,
  scale = "none",
  trace = "none",
  dendrogram = "column",
  ColSideColors = col_side_colors,
  col = colorRampPalette(rev(brewer.pal(9, "RdBu")))(255),
  margins = c(8, 10),
  key = TRUE,
  keysize = 1.2,
  density.info = "none",
  labRow = row_labels   # <-- only change: labels
)

dev.off()

# -------------------------
# 4b. Interpret PCA (what drives PC1?)
# -------------------------

# Extract PCA data (if not already done)
vsd <- vst(dds, blind=FALSE)
pca_data <- plotPCA(vsd, intgroup="group", returnData=TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

# -------------------------
# Parse biological variables from group
# -------------------------

# Example group format: cort_M_45min or vehicle_M_45min_Rpl22CrhCre

pca_data$treatment <- sub("_.*", "", pca_data$group)

pca_data$sex <- ifelse(grepl("_F_", pca_data$group), "F", "M")

pca_data$time <- sub(".*_(\\d+.*)", "\\1", pca_data$group)

pca_data$odor <- ifelse(grepl("_odor$", pca_data$group), "odor", "no_odor")

pca_data$genotype <- ifelse(grepl("Rpl22CrhCre", pca_data$group),
                           "Rpl22CrhCre", "WT")

# -------------------------
# Define base PCA plot function
# -------------------------

make_pca_plot <- function(color_var, title) {
  ggplot(pca_data, aes_string("PC1", "PC2", color=color_var)) +
    geom_point(size=4) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    ggtitle(title) +
    theme_minimal()
}

# -------------------------
# Generate plots
# -------------------------

p_treatment <- make_pca_plot("treatment", "PCA colored by treatment")
p_sex       <- make_pca_plot("sex", "PCA colored by sex")
p_time      <- make_pca_plot("time", "PCA colored by time")
p_genotype  <- make_pca_plot("genotype", "PCA colored by genotype")
p_odor      <- make_pca_plot("odor", "PCA colored by odor exposure")

# -------------------------
# Save plots
# -------------------------

outdir <- "/Users/norawolcott/Documents/Datta Lab/data/Biopolymer Facility RNASeq/260510/"

# Put all plots into a list
p_list <- list(
  p_treatment,
  p_sex,
  p_time,
  p_genotype,
  p_odor
)

# If number of plots are odd, add an empty plot as placeholder
empty_plot <- grid.rect(gp=gpar(col="white"))  # blank white space
p_list <- c(p_list, list(empty_plot))

# Arrange in 3 rows x 2 columns
combined_plot <- do.call(grid.arrange, c(p_list, ncol=2, nrow=3))

# Save
ggsave(filename = paste0(outdir, "PCA_grid.png"),
       plot = combined_plot,
       width = 12, height = 8)

# -------------------------
# 6. Multi-panel Volcano Plots
# -------------------------

# # Define comparisons
# comparisons <- list(
#   c("vehicle_M_45min", "cort_M_45min"),
#   c("vehicle_F_45min", "cort_F_45min"),
#   c("vehicle_M_1hr",   "cort_M_1hr"),
#   c("vehicle_M_3hr",   "cort_M_3hr"),
#   c("vehicle_F_3hr",   "cort_F_3hr"),
#   c("vehicle_M_6hr",   "cort_M_6hr"),
#   c("vehicle_F_6hr",   "cort_F_6hr"),
#   c("vehicle_M_24hr",   "cort_M_24hr"),
#   c("vehicle_M_2hr_odor",   "cort_M_2hr_odor"),
#   c("vehicle_F_2hr_odor",   "cort_F_2hr_odor"),
#   c("vehicle_M_2hr_DPG",   "cort_M_2hr_odor"),
#   c("vehicle_F_2hr_DPG",   "cort_F_2hr_odor"),
#   c("vehicle_M_2hr_DPG",   "vehicle_M_2hr_odor"),
#   c("vehicle_M_2hr_DPG",   "cort_M_2hr_DPG"),
#   c("vehicle_F_2hr_DPG",   "cort_F_2hr_DPG"),
#   c("cort_M_2hr_DPG",   "cort_M_2hr_odor"),
#   c("cort_F_2hr_DPG",   "cort_F_2hr_odor"),
#   c("vehicle_M_45min_Rpl22CrhCre", "cort_M_45min_Rpl22CrhCre"),
#   c("vehicle_F_45min",  "E2_F_45min"),
#   c("vehicle_F_45min",  "P4_F_45min"),
#   c("vehicle_F_3hr",  "E2_F_3hr"),
#   c("vehicle_F_3hr",  "P4_F_3hr"),
#   c("vehicle_F_6hr",  "E2_F_6hr"),
#   c("vehicle_F_6hr",  "P4_F_6hr"),
#   c("vehicle_F_24hr",  "E2_F_24hr"),
#   c("vehicle_F_24hr",  "P4_F_24hr"),
#   c("vehicle_F_2hr_odor",   "E2_F_2hr_odor"),
#   c("vehicle_F_2hr_odor",   "P4_F_2hr_odor"),
#   c("vehicle_F_2hr_DPG",   "vehicle_F_2hr_odor"),
#   c("vehicle_F_2hr_DPG",   "E2_F_2hr_odor"),
#   c("vehicle_F_2hr_DPG",   "P4_F_2hr_odor"),
#   c("vehicle_F_2hr_DPG",   "E2_F_2hr_DPG"),
#   c("vehicle_F_2hr_DPG",   "P4_F_2hr_DPG"),
#   c("E2_F_2hr_DPG",   "E2_F_2hr_odor"),
#   c("P4_F_2hr_DPG",   "P4_F_2hr_odor")
# )

# # Function to generate a volcano plot
# make_volcano <- function(group1, group2) {
  
#   res <- results(dds, contrast = c("group", group2, group1))
#   res <- lfcShrink(dds, contrast = c("group", group2, group1), type="normal")
  
#   EnhancedVolcano(res,
#                   lab = rownames(res),
#                   x = 'log2FoldChange',
#                   y = 'padj',
#                   pCutoff = 0.05,
#                   FCcutoff = 1,
#                   title = paste(group2, "vs", group1),
#                   labSize = 3,
#                   pointSize = 2.5)
# }

# # Generate all plots
# volcano_plots <- lapply(comparisons, function(x) {
#   make_volcano(x[1], x[2])
# })

# # Arrange in grid
# combined_plot <- grid.arrange(grobs = volcano_plots, ncol = 6, nrow = 6)

# # Save figure
# ggsave(
#   "/Users/norawolcott/Documents/Datta Lab/data/Biopolymer Facility RNASeq/260510/Volcano_multipanel.png",
#   combined_plot,
#   width = 12,
#   height = 15
# )