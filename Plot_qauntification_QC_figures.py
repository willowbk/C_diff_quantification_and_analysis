import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
from sklearn.decomposition import PCA
from sklearn.manifold import MDS

# File names
file_name = "CdiffT6_all_conditions"
norm_file = "NormLib" + file_name + ".csv"
n_labels = 10# Number of extreme genes to label
logCPM_threshold = 2

# Output directory
output_dir = file_name + "_QC_figures"
os.makedirs(output_dir, exist_ok=True)

scatter_output_dir = output_dir + "/Replicate_scatter"
os.makedirs(scatter_output_dir, exist_ok=True)

# Read normalized expression table
df_norm = pd.read_csv(norm_file)

# Sample names
sample_names = list(df_norm.columns[3:])

# Convert normalized expression to numeric matrix
data_matrix = df_norm.iloc[:, 3:].to_numpy(dtype=float)

# Gene names
gene_names = df_norm["gene_name"].to_numpy()

print("Samples:")
for sample in sample_names:
	print("\t" + sample)

print("\nNumber of genes:", len(df_norm))
print("Number of samples:", len(sample_names))


# Infer condition names from sample names
conditions = []

for sample in sample_names:
	if " rep." in sample:
		condition = sample.split(" rep.")[0]
	elif " Rep." in sample:
		condition = sample.split(" Rep.")[0]
	else:
		condition = sample
	conditions.append(condition)

unique_conditions = list(dict.fromkeys(conditions))

# ------------------------------------------------------------------
# PCA
# ------------------------------------------------------------------

pca = PCA(n_components=2)
pca_coordinates = pca.fit_transform(data_matrix.T)

print("\nPCA explained variance:")
print(pca.explained_variance_ratio_)

fig, ax = plt.subplots(figsize=(8, 7))

# Generate one color for each condition
cmap = plt.get_cmap("tab20")
condition_colors = {condition: cmap(i % 20) for i, condition in enumerate(unique_conditions)}

for condition in unique_conditions:
	indices = [i for i, x in enumerate(conditions) if x == condition]
	ax.scatter(pca_coordinates[indices, 0], pca_coordinates[indices, 1], s=50, label=condition, color=condition_colors[condition])
	
for i, sample in enumerate(sample_names):
	ax.annotate(sample, (pca_coordinates[i, 0], pca_coordinates[i, 1]), xytext=(5, 5), textcoords="offset points", fontsize=5)

ax.set_xlabel("PC1 (" + str(round(100 * pca.explained_variance_ratio_[0], 1)) + "% var.)", fontsize=10)
ax.set_ylabel("PC2 (" + str(round(100 * pca.explained_variance_ratio_[1], 1)) + "% var.)", fontsize=10)
ax.set_title(file_name, fontsize=10)
ax.legend(frameon=False)
plt.tight_layout()
plt.savefig(output_dir + "/PCA_" + file_name + ".pdf", dpi=350, bbox_inches="tight")
plt.close()

# ------------------------------------------------------------------
# MDS
# ------------------------------------------------------------------

expression_matrix = data_matrix.T
# Calculate pairwise Euclidean distances between samples
distance_matrix = np.zeros((len(sample_names), len(sample_names)))

for i in range(len(sample_names)):
	for j in range(len(sample_names)):
		distance_matrix[i, j] = np.sqrt(np.sum((expression_matrix[i] - expression_matrix[j]) ** 2))


# Perform classical MDS
mds = MDS(n_components=2, dissimilarity="precomputed", random_state=0)
mds_coordinates = mds.fit_transform(distance_matrix)

# Create MDS plot
fig, ax = plt.subplots(figsize=(9, 7))
for condition in unique_conditions:
	condition_indices = [i for i, sample_condition in enumerate(conditions) if sample_condition == condition]
	ax.scatter(mds_coordinates[condition_indices, 0], mds_coordinates[condition_indices, 1], s=50, label=condition)

# Label samples
for i, sample in enumerate(sample_names):
	ax.annotate(sample, (mds_coordinates[i, 0], mds_coordinates[i, 1]), xytext=(5, 5), textcoords="offset points", fontsize=6)

ax.set_xlabel("MDS dimension 1")
ax.set_ylabel("MDS dimension 2")
ax.set_title("MDS Plot: All Conditions")
ax.legend(frameon=False, fontsize=8)

plt.tight_layout()
plt.savefig(output_dir + "/MDS_" + file_name + ".pdf", dpi=350, bbox_inches="tight")
plt.close()

# ------------------------------------------------------------------
# RLE plot
# ------------------------------------------------------------------

# Calculate the median expression of each gene across all samples
gene_medians = np.median(data_matrix, axis=1)

# Calculate deviation from the gene median
rle_data = []

for i in range(len(sample_names)):
	rle_data.append(data_matrix[:, i] - gene_medians)

fig, ax = plt.subplots(figsize=(10, 7))

ax.axhline(0, linestyle="--", linewidth=1)
ax.boxplot(rle_data)

ax.set_xticks(range(1, len(sample_names) + 1))
ax.set_xticklabels(sample_names, rotation=45, ha="right", fontsize=6)
ax.set_ylabel("Deviation from gene median (log₂-CPM)", fontsize=10)
ax.set_title(file_name, fontsize=10)

plt.tight_layout()
plt.savefig(output_dir + "/RLE_" + file_name + ".pdf", dpi=350, bbox_inches="tight")
plt.close()


# ------------------------------------------------------------------
# Correlation matrix
# ------------------------------------------------------------------

correl_mat = np.corrcoef(data_matrix.T)

min_r = np.min(correl_mat)

fig, ax = plt.subplots(figsize=(10, 8))

im = ax.imshow(correl_mat, cmap="coolwarm", interpolation="nearest", aspect="auto", vmin=min_r, vmax=1)

cbar = plt.colorbar(im, ax=ax)
cbar.set_label("Pearson correlation (log₂-CPM)", fontsize=10)

ax.set_xticks(range(len(sample_names)))
ax.set_xticklabels(sample_names, fontsize=6, rotation=45, ha="right")
ax.set_yticks(range(len(sample_names)))
ax.set_yticklabels(sample_names, fontsize=6)

ax.set_title(file_name, fontsize=10)

plt.tight_layout()
plt.savefig(output_dir + "/Correlation_matrix_" + file_name + ".pdf", dpi=350, bbox_inches="tight")
plt.close()


# ------------------------------------------------------------------
# Pairwise scatter plots for replicates
# ------------------------------------------------------------------

for condition in unique_conditions:

	condition_indices = [i for i, x in enumerate(conditions) if x == condition]

	for a in range(len(condition_indices) - 1):

		for b in range(a + 1, len(condition_indices)):

			i = condition_indices[a]
			j = condition_indices[b]

			sample_i = sample_names[i]
			sample_j = sample_names[j]

			x = data_matrix[:, i]
			y = data_matrix[:, j]

			r = stats.pearsonr(x, y).statistic

			# Find genes furthest from the diagonal above the log2-CPM threshold
			differences = np.abs(x - y)
			eligible_indices = np.where((x >= logCPM_threshold) | (y >= logCPM_threshold))[0]

			if len(eligible_indices) > 0:
				outlier_indices = eligible_indices[np.argsort(differences[eligible_indices])[-n_labels:]]
				outlier_indices = outlier_indices[::-1]
			else:
				outlier_indices = []
	
			fig, ax = plt.subplots(figsize=(8, 8))
			
			# Plot expression thresholds
			ax.axvline(logCPM_threshold, linestyle="--", color="grey", linewidth=1)
			ax.axhline(logCPM_threshold, linestyle="--", color="grey", linewidth=1)
			
			# Plot scatterplot
			ax.scatter(x, y, s=8, alpha=0.5, edgecolors="none")
			
			# Plot diagonal
			min_expression = min(np.min(x), np.min(y))
			max_expression = max(np.max(x), np.max(y))
			ax.plot([min_expression, max_expression], [min_expression, max_expression], linestyle="--", color="red", linewidth=1)

			# Label outlier genes
			for index in outlier_indices:
				label = gene_names[index]
				ax.annotate(label, (x[index], y[index]), xytext=(3, 3), textcoords="offset points", fontsize=5)

			ax.text(0.05, 0.95, "r = " + str(round(r, 3)), transform=ax.transAxes, fontsize=10, va="top")

			ax.set_xlabel(sample_i + " (log₂-CPM)", fontsize=10)
			ax.set_ylabel(sample_j + " (log₂-CPM)", fontsize=10)
			ax.set_title(condition, fontsize=11)

			plt.tight_layout()

			sample_i_file = sample_i.replace(" ", "")
			sample_j_file = sample_j.replace(" ", "")

			plt.savefig(scatter_output_dir + "/" + sample_i_file + "_" + sample_j_file + "_scatter.pdf", dpi=350, bbox_inches="tight")
			plt.close()
