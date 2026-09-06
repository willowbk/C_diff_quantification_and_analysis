library(edgeR)

file_name <- "CdiffT6_all_conditions"

group <- factor(c(rep("H2O", 3), rep("METRO", 3)), levels=c("H2O", "METRO"))
sample_patterns <- c(H2O="Water control early", METRO="Metronidazol early")

df <- read.csv(paste0(file_name, "_count_table.csv"), header=TRUE, check.names=FALSE) # Read count table

df_counts <- df[, !(colnames(df) %in% c("locus_tag", "gene_name", "type"))] # Extract count matrix
rownames(df_counts) <- df$locus_tag # Use locus tags as row names
Counts <- data.matrix(df_counts) # Convert to numeric matrix

selected_samples <- c()
sample_groups <- c()

for (condition in levels(group)) {
	pattern <- sample_patterns[condition]
	samples <- colnames(Counts)[startsWith(tolower(colnames(Counts)), tolower(paste0(pattern, " ")))]
	expected <- sum(group == condition)

	if (length(samples) != expected) {
		stop(paste0("Expected ", expected, " samples for ", condition, ", but found ", length(samples), " matching '", pattern, "'."))
	}

	selected_samples <- c(selected_samples, samples)
	sample_groups <- c(sample_groups, rep(condition, length(samples)))
}

Counts <- Counts[, selected_samples, drop=FALSE] # Keep only selected samples
sample_groups <- factor(sample_groups, levels=levels(group)) # Assign condition labels

y <- DGEList(counts=Counts, group=sample_groups) # Create DGEList

keep <- filterByExpr(y) # Filter low-expression genes
y <- y[keep,,keep.lib.sizes=FALSE]

y <- normLibSizes(y) # TMM normalization

y <- estimateDisp(y) # Estimate dispersions

pdf(paste0(file_name, "_", levels(group)[2], "_vs_", levels(group)[1], "_BCV_plot.pdf")) # Save BCV plot
plotBCV(y)
dev.off()

design <- model.matrix(~0+sample_groups) # Create design matrix
colnames(design) <- levels(group)

fit <- glmQLFit(y, design) # Fit model
cont <- makeContrasts(contrasts=paste0(levels(group)[2], " - ", levels(group)[1]), levels=design) # Define contrast
qlf <- glmQLFTest(fit, contrast=cont) # Differential expression test
out <- topTags(qlf, n=Inf) # Extract all results

write.csv(out, file=paste0("DiffExpr", file_name, "_", levels(group)[2], "_vs_", levels(group)[1], ".csv")) # Write results
