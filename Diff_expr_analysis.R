
# Install once if needed:
# BiocManager::install("edgeR")
# install.packages("data.table")

library('edgeR')

file_name = 'CdiffT6_2_conditions'
# Read count table
df <- read.csv(paste0(file_name, "_count_table.csv"), header=TRUE, check.names=FALSE)

# Save annotation
annotation <- df[, c("locus_tag", "gene_name", "type")]

# Extract count matrix
df_counts <- df[, !(colnames(df) %in% c("locus_tag", "gene_name", "type"))]

# Use locus_tag as row names
rownames(df_counts) <- df$locus_tag


# Create normalized data table
cpm_vals <- cpm(df_counts, log=TRUE, prior.count=1)
write.csv(cpm_vals, file=paste0("NormLib", file_name, ".csv"))


# Differential expression analysis
Counts <- data.matrix(df_counts)

# Define which conditions to perform differential expression
group <- factor(c(rep("H2O",3), rep("METRO",3)))

y <- DGEList(counts=Counts, group=group)

# Filter low expression genes
keep <- filterByExpr(y)

y <- y[keep,,keep.lib.sizes=FALSE]


# Normalize library sizes
y <- normLibSizes(y)


# Estimate dispersions
y <- estimateDisp(y)


# BCV plot
pdf(paste0(file_name, "_BCV_plot.pdf"))
plotBCV(y)
dev.off()


# Design matrix
design <- model.matrix(~0+group)
colnames(design) <- levels(group)


# Fit model
fit <- glmQLFit(y, design)


# Contrast
cont <- makeContrasts(METRO - H2O, levels=design)


# Differential expression test
qlf <- glmQLFTest(fit, contrast=cont)


# Results
out <- topTags(qlf, n=Inf)


# Write results
write.csv(out, file=paste0("DiffExpr", file_name, ".csv"))

