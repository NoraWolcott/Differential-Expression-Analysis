# -------------------------
# Gene lookup script (mouse)
# -------------------------

# Install if needed
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

if (!requireNamespace("org.Mm.eg.db", quietly = TRUE))
    BiocManager::install("org.Mm.eg.db")

# Load libraries
library(org.Mm.eg.db)
library(AnnotationDbi)

# Your gene list
genes <- c(
"ENSMUSG00000021342",
"ENSMUSG00000102070",
"ENSMUSG00000090877",
"ENSMUSG00000086503",
"ENSMUSG00000007907",
"ENSMUSG00000007457",
"ENSMUSG00000038572",
"ENSMUSG00000067998",
"ENSMUSG00000067996",
"ENSMUSG00000009580",
"ENSMUSG00000000983",
"ENSMUSG00000082635",
"ENSMUSG00000047228",
"ENSMUSG00000067684",
"ENSMUSG00000067679",
"ENSMUSG00000044121",
"ENSMUSG00000069080",
"ENSMUSG00000079522",
"ENSMUSG00000041333",
"ENSMUSG00000055961",
"ENSMUSG00000058523",
"ENSMUSG00000079519",
"ENSMUSG00000079521",
"ENSMUSG00000067541",
"ENSMUSG00000005980",
"ENSMUSG00000042179",
"ENSMUSG00000062061",
"ENSMUSG00000095577",
"ENSMUSG00000079539",
"ENSMUSG00000094113",
"ENSMUSG00000066583",
"ENSMUSG00000015519",
"ENSMUSG00000059040",
"ENSMUSG00000069045",
"ENSMUSG00000069049"
)

# Lookup table
gene_info <- select(org.Mm.eg.db,
                    keys = genes,
                    columns = c("SYMBOL", "GENENAME"),
                    keytype = "ENSEMBL")

# Print results
print(gene_info)

# Save to your specified folder
write.csv(gene_info,
          file = "/Users/norawolcott/Documents/Datta Lab/data/Biopolymer Facility Seq/260324/nextflow results/gene_lookup_results.csv",
          row.names = FALSE)