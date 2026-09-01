import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


# Input file
file_name = "DiffExprCdiffT6_2_conditions.csv"

# Read differential expression results
df = pd.read_csv(file_name)


# Significance thresholds
fdr_threshold = 0.05# 0.05 is totally bogus, but everyone uses it...
logfc_threshold = 1
n_labels = 10# Select number of genes to label

# Calculate -log10(FDR)
df["neg_log10_FDR"] = -np.log10(df["FDR"].clip(lower=1e-300))


# Classify genes
df["significant"] = ((df["FDR"] < fdr_threshold) & (df["logFC"].abs() >= logfc_threshold))
df["direction"] = "Not significant"
df.loc[df["significant"] & (df["logFC"] > 0), "direction"] = "Upregulated"
df.loc[df["significant"] & (df["logFC"] < 0), "direction"] = "Downregulated"


significant = df[df["significant"]].copy()
upregulated = significant[significant["logFC"] > 0]
downregulated = significant[significant["logFC"] < 0]

top_upregulated = upregulated.nlargest(n_labels, "neg_log10_FDR")
top_downregulated = downregulated.nlargest(n_labels, "neg_log10_FDR")

# Choose gene label
if "gene_name" in df.columns:
	label_column = "gene_name"
elif "locus_tag" in df.columns:
	label_column = "locus_tag"
else:
	label_column = df.columns[0]


# Create plot
fig, ax = plt.subplots(figsize=(9, 7))

# Plot non-significant genes
ax.scatter(df.loc[df["direction"] == "Not significant", "logFC"], df.loc[df["direction"] == "Not significant", "neg_log10_FDR"], s=12, alpha=0.4, label="Not significant")
# Plot upregulated genes
ax.scatter(df.loc[df["direction"] == "Upregulated", "logFC"], df.loc[df["direction"] == "Upregulated", "neg_log10_FDR"], s=18, alpha=0.8, label="Upregulated")
# Plot downregulated genes
ax.scatter(df.loc[df["direction"] == "Downregulated", "logFC"], df.loc[df["direction"] == "Downregulated", "neg_log10_FDR"], s=18, alpha=0.8, label="Downregulated")

# Threshold lines
ax.axhline(-np.log10(fdr_threshold), linestyle="--", linewidth=1)
ax.axvline(-logfc_threshold, linestyle="--", linewidth=1)
ax.axvline(logfc_threshold, linestyle="--", linewidth=1)

# Add labels
for _, row in pd.concat([top_upregulated, top_downregulated]).iterrows():
	ax.annotate(row[label_column], (row["logFC"], row["neg_log10_FDR"]), xytext=(5, 5), textcoords="offset points", fontsize=8)

# Labels and formatting
ax.set_xlabel("log2 fold change")
ax.set_ylabel("-log10(FDR)")
ax.set_title("Differential Expression")
ax.legend(frameon=False)
plt.tight_layout()
plt.savefig(file_name.replace('.csv', '')+"_volcano.pdf", bbox_inches="tight")
plt.show()


