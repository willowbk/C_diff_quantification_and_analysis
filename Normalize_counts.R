# Install once if needed:
# BiocManager::install("edgeR")

library(edgeR)

# Input and output file name
file_name <- "CdiffT6_all_conditions"

# Read count table
df <- read.csv(paste0(file_name, "_count_table.csv"), header=TRUE, check.names=FALSE)

# Set annotation
annotation <- df[, c("locus_tag", "gene_name", "type")]

# Extract count matrix
df_counts <- df[, !(colnames(df) %in% c("locus_tag", "gene_name", "type"))]

# Use locus_tag as row names
rownames(df_counts) <- df$locus_tag

# Convert to count matrix
Counts <- data.matrix(df_counts)

# Create DGEList
y <- DGEList(counts=Counts)

# TMM normalization
y <- normLibSizes(y)

# Calculate normalized log2-CPM values
normalized_counts <- cpm(y, log=TRUE, prior.count=1)

# Convert normalized counts to data frame
normalized_counts <- as.data.frame(normalized_counts)

# Add annotation as the first three columns
normalized_counts <- cbind(annotation, normalized_counts)

# Write normalized table
write.csv(normalized_counts, file=paste0("NormLib", file_name, ".csv"), row.names=FALSE)
