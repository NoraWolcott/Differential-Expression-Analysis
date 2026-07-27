library(DESeq2)
library(ggplot2)
library(EnhancedVolcano)
library(AnnotationDbi)
library(org.Mm.eg.db)

# -------------------------
# Load data (EDIT)
# -------------------------
dds <- readRDS('/Users/norawolcott/Documents/Datta Lab/data/Biopolymer Facility Seq/260403/Contrast models/dds_v_F_2hr_odor.rds')

# If you already have dds in memory, skip this step

# -------------------------
# Define contrast
# -------------------------
group1 <- "vehicle_F_2hr_odor"
group2 <- "P4_F_2hr_odor"

# -------------------------
# Run DESeq2
# -------------------------
res <- results(dds, contrast = c("group", group2, group1))

# Shrink LFC
res <- lfcShrink(
  dds,
  coef = "group_P4_F_2hr_odor_vs_vehicle_F_2hr_odor",
  type = "apeglm"
)

# Convert rownames (Ensembl IDs) to gene symbols
res$symbol <- mapIds(
  org.Mm.eg.db,
  keys = rownames(res),
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

res$symbol[is.na(res$symbol)] <- rownames(res)[is.na(res$symbol)]

# -------------------------
# Save results table
# -------------------------
out_file <- paste0(
  "/Users/norawolcott/Documents/Datta Lab/data/Biopolymer Facility Seq/260403/DESeq2_",
  group2, "_vs_", group1, ".csv"
)

write.csv(as.data.frame(res), out_file)

# -------------------------
# Volcano plot
# -------------------------

p <- EnhancedVolcano(
  res,
  lab = res$symbol,
  x = 'log2FoldChange',
  y = 'padj',
  pCutoff = 0.05,
  FCcutoff = 1,
  title = paste(group2, "vs", group1),
  labSize = 3,
  pointSize = 2.5
)

# -------------------------
# Save plot
# -------------------------
plot_file <- paste0(
  "/Users/norawolcott/Documents/Datta Lab/data/Biopolymer Facility Seq/260403/Volcano_",
  group2, "_vs_", group1, ".png"
)

ggsave(plot_file, plot = p, width = 15, height = 10)

print(p)

print(paste("DESeq2 results saved to:", out_file))