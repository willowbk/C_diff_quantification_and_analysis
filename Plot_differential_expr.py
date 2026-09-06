import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


# Input files
file_name = "DiffExprCdiffT6_all_conditions_METRO_vs_H2O.csv"
norm_file = "NormLibCdiffT6_all_conditions.csv"


# Conditions used for the plots
cond_1 = "Water control early"
cond_2 = "Metronidazol early"
cond_1_name = "H2O"
cond_2_name = "METRO"


# Significance thresholds
fdr_threshold = 0.05 # 0.05 is totally bogus, but everyone uses it...
logfc_threshold = 1# Threshold for genes considered differentially expressed.
n_labels = 10 # Select number of genes to label, starting from the most-significant


# Dictionary of genes that should ALWAYS be labeled
# Key = gene name in the DE results
# Value = text that should appear on the plot
genes_to_label = {
	"CdiffT6_04024": "Scheler's gene 1",
	"CdiffT6_00331": "Scheler's gene 2"
}


# Read differential expression results
df = pd.read_csv(file_name, index_col=0)
df.index.name = "gene_name"


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


# Don't automatically label genes which have custom labels
top_upregulated = upregulated[~upregulated.index.isin(genes_to_label)].nlargest(n_labels, "neg_log10_FDR")
top_downregulated = downregulated[~downregulated.index.isin(genes_to_label)].nlargest(n_labels, "neg_log10_FDR")


# Read normalized expression data
df_norm = pd.read_csv(norm_file, index_col=0)
df_norm.index.name = "gene_name"
gene_name_map = df_norm["gene_name"].to_dict()


# Create volcano plot
fig, ax = plt.subplots(figsize=(9, 7))

ax.scatter(df.loc[df["direction"] == "Not significant", "logFC"], df.loc[df["direction"] == "Not significant", "neg_log10_FDR"], s=12, alpha=0.4, label="Not significant")
ax.scatter(df.loc[df["direction"] == "Upregulated", "logFC"], df.loc[df["direction"] == "Upregulated", "neg_log10_FDR"], s=18, alpha=0.8, label="Upregulated")
ax.scatter(df.loc[df["direction"] == "Downregulated", "logFC"], df.loc[df["direction"] == "Downregulated", "neg_log10_FDR"], s=18, alpha=0.8, label="Downregulated")

for gene, custom_label in genes_to_label.items():
	if gene in df.index:
		row = df.loc[gene]
		ax.scatter(row["logFC"], row["neg_log10_FDR"], s=30, color="red", zorder=5)
		ax.annotate(custom_label, (row["logFC"], row["neg_log10_FDR"]), xytext=(5, 5), textcoords="offset points", fontsize=8, color="red")

ax.axhline(-np.log10(fdr_threshold), linestyle="--", linewidth=1)
ax.axvline(-logfc_threshold, linestyle="--", linewidth=1)
ax.axvline(logfc_threshold, linestyle="--", linewidth=1)

for _, row in pd.concat([top_upregulated, top_downregulated]).iterrows():
	ax.annotate(gene_name_map.get(row.name, row.name), (row["logFC"], row["neg_log10_FDR"]), xytext=(5, 5), textcoords="offset points", fontsize=8)

ax.set_xlabel("log2 fold change")
ax.set_ylabel("-log10(FDR)")
ax.set_title(f"Differential Expression: {cond_2_name} vs {cond_1_name}")
ax.legend(frameon=False)

plt.tight_layout()
plt.savefig(file_name.replace(".csv", "") + "_volcano.pdf", bbox_inches="tight")
plt.close()




# Find samples belonging to each condition, ignoring capitalization
cond_1_list = [sample for sample in df_norm.columns if sample.lower().startswith(cond_1.lower() + " ")]
cond_2_list = [sample for sample in df_norm.columns if sample.lower().startswith(cond_2.lower() + " ")]


# Check that the expected samples were found
if len(cond_1_list) != 3:
	raise ValueError(f"Expected 3 samples for {cond_1_name}, but found {len(cond_1_list)}: {cond_1_list}")

if len(cond_2_list) != 3:
	raise ValueError(f"Expected 3 samples for {cond_2_name}, but found {len(cond_2_list)}: {cond_2_list}")


# Keep only genes present in both datasets and align the normalized data to the DE results
common_genes = df.index.intersection(df_norm.index)
df = df.loc[common_genes].copy()
df_norm = df_norm.loc[common_genes]


# Calculate means for the two conditions
mean_1 = df_norm[cond_1_list].mean(axis=1).to_numpy()
mean_2 = df_norm[cond_2_list].mean(axis=1).to_numpy()

MA_logFC = mean_2 - mean_1
MA_mean = (mean_2 + mean_1) / 2


# Create MA plot
fig, ax = plt.subplots(figsize=(9, 7))

ax.scatter(MA_mean, MA_logFC, s=12, color="grey", alpha=0.5, edgecolors="none")

significant_genes = set(df.index[df["significant"]])
significant_indices = [k for k, gene in enumerate(df.index) if gene in significant_genes and gene not in genes_to_label]

ax.scatter(MA_mean[significant_indices], MA_logFC[significant_indices], s=18, color="grey", edgecolors="black", linewidths=0.5)

for gene, custom_label in genes_to_label.items():
	if gene in df.index:
		k = list(df.index).index(gene)
		ax.scatter(MA_mean[k], MA_logFC[k], s=30, color="red", zorder=5)
		ax.annotate(custom_label, (MA_mean[k], MA_logFC[k]), xytext=(5, 5), textcoords="offset points", fontsize=8, color="red")

MA_logFC_mean = sorted([[MA_logFC[k], MA_mean[k], df.index[k]] for k in significant_indices], key=lambda x: x[0])

MA_labels = [term for term in MA_logFC_mean[:n_labels] if term[0] < 0] + [term for term in MA_logFC_mean[-n_labels:] if term[0] > 0]

for term in MA_labels:
	ax.text(term[1], term[0], gene_name_map.get(term[2], term[2]), fontsize=8)

ax.axhline(0, 0, 1, linewidth=1, color="red")
ax.axhline(logfc_threshold, 0, 1, linewidth=1, color="black", linestyle="dashed")
ax.axhline(-logfc_threshold, 0, 1, linewidth=1, color="black", linestyle="dashed")

ax.set_title(f"MA Plot: {cond_2_name} vs {cond_1_name}")
ax.set_ylabel("log2 fold change")
ax.set_xlabel("Mean expression (log2-CPM)")

min_MA_logFC = min(MA_logFC)
max_MA_logFC = max(MA_logFC)
range_MA_logFC = max_MA_logFC - min_MA_logFC
ax.set_ylim([min_MA_logFC - 0.05 * range_MA_logFC, max_MA_logFC + 0.05 * range_MA_logFC])

min_MA_mean = min(MA_mean)
max_MA_mean = max(MA_mean)
range_MA_mean = max_MA_mean - min_MA_mean
ax.set_xlim([min_MA_mean - 0.05 * range_MA_mean, max_MA_mean + 0.05 * range_MA_mean])

plt.tight_layout()
plt.savefig(file_name.replace(".csv", "") + "_MA.pdf", bbox_inches="tight")
plt.close()


# Create p-value distribution
fig, ax = plt.subplots(figsize=(9, 7))

ax.hist(df["PValue"], bins=50, alpha=0.7)

ax.set_xlabel("P-value")
ax.set_ylabel("Number of genes")
ax.set_title(f"P-value Distribution: {cond_2_name} vs {cond_1_name}")

plt.tight_layout()
plt.savefig(file_name.replace(".csv", "") + "_pvalue_distribution.pdf", bbox_inches="tight")
plt.close()
