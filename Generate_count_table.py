import pandas as pd
import numpy as np
from Bio import SeqIO
import os
import matplotlib.pyplot as plt


# Path to genome file
folder_name_genbank = 'Genomes/DMG2301622_Cdiff.gbk'
identifier = 'CdiffT6'
out_file_name = 'CdiffT6_all_conditions'# This is only used to name the final table


# Only extract features of the following types for the DE count table
# rRNA and tRNA are excluded so that they do not enter differential expression
feature_list = ['CDS', 'tmRNA', 'repeat_region']


# Feature types to include in the biotype composition plot
# These can include rRNA and tRNA even though they are excluded from the DE table
plot_feature_list = ['CDS', 'tmRNA', 'repeat_region', 'tRNA', 'rRNA']


# Create human-readable versions of each condition
file_to_conds = {'A10_METRO_E1_ZKRN260013229-1A_23H7NJLT3_L8':'Metronidazol early Rep.1',
		'A11_METRO_E2_ZKRN260013230-1A_23H7NJLT3_L8':'Metronidazol early Rep.2',
		'A12_METRO_E3_ZKRN260013231-1A_23H7NJLT3_L8':'Metronidazol early Rep.3',
		'A13_VANCO_E1_ZKRN260013232-1A_23H7NJLT3_L8':'Vancomycin early Rep.1',
		'A14_VANCO_E2_ZKRN260013233-1A_23HVFCLT4_L7':'Vancomycin early Rep.2',
		'A15_VANCO_E3_ZKRN260013234-1A_23HVFCLT4_L7':'Vancomycin early Rep.3',
		'A16_H2O_S1_ZKRN260013235-1A_23HVFCLT4_L7':'Water control late Rep.1',
		'A17_H2O_S2_ZKRN260013236-1A_23HVFCLT4_L7':'Water control late Rep.2',
		'A18_H2O_S3_ZKRN260013237-1A_23HVFCLT4_L7':'Water control late Rep.3',
		'A19_DMSO_S1_ZKRN260013238-1A_23HVFCLT4_L7':'DMSO control late Rep.1',
		'A1_H2O_E1_ZKRN260013220-1A_23H7NJLT3_L6':'Water control early Rep.1',
		'A20_DMSO_S2_ZKRN260013239-1A_23HVFCLT4_L7':'DMSO control late Rep.2',
		'A21_DMSO_S3_ZKRN260013240-1A_23HVFCLT4_L7':'DMSO control late Rep.3',
		'A22_LL37_S1_ZKRN260013241-1A_23HVFCLT4_L7':'LL37 peptide late Rep.1',
		'A23_LL37_S2_ZKRN260013242-1A_23HVFCLT4_L7':'LL37 peptide late Rep.2',
		'A24_LL37_S3_ZKRN260013243-1A_23HVFCLT4_L7':'LL37 peptide late Rep.3',
		'A25_METRO_S1_ZKRN260013244-1A_23HVFCLT4_L7':'Metronidazol late Rep.1',
		'A26_METRO_S2_ZKRN260013245-1A_23HVFCLT4_L7':'Metronidazol late Rep.2',
		'A27_METRO_S3_ZKRN260013246-1A_23HVFCLT4_L7':'Metronidazol late Rep.3',
		'A28_VANCO_S1_ZKRN260013247-1A_23HVFCLT4_L7':'Vancomycin late Rep.1',
		'A29_VANCO_S2_ZKRN260013248-1A_23HVFCLT4_L7':'Vancomycin late Rep.2',
		'A30_VANCO_S3_ZKRN260013249-1A_23HVFCLT4_L7':'Vancomycin late Rep.3',
		'A2_H2O_E2_ZKRN260013221-1A_23HVFCLT4_L7':'Water control early Rep.2',
		'A3_H2O_E3_ZKRN260013222-1A_23H7NJLT3_L6':'Water control early Rep.3',
		'A4_DMSO_E1_ZKRN260013223-1A_23H7NJLT3_L6':'DMSO control early Rep.1',
		'A5_DMSO_E2_ZKRN260013224-1A_23H7NJLT3_L6':'DMSO control early Rep.2',
		'A6_DMSO_E3_ZKRN260013225-1A_23H7NJLT3_L6':'DMSO control early Rep.3',
		'A7_LL37_E1_ZKRN260013226-1A_23H7NJLT3_L6':'LL37 peptide early Rep.1',
		'A8_LL37_E2_ZKRN260013227-1A_23H7NJLT3_L8':'LL37 peptide early Rep.2',
		'A9_LL37_E3_ZKRN260013228-1A_23H7NJLT3_L8':'LL37 peptide early Rep.3',
		'B_24h_1647_1_ZKRN260013259-1A_23HVFCLT4_L7':'Community FFS-1647 late Rep.1',
		'B_24h_1647_2_ZKRN260013260-1A_23HVFCLT4_L7':'Community FFS-1647 late Rep.2',
		'B_24h_1647_3_ZKRN260013261-1A_23HVFCLT4_L7':'Community FFS-1647 late Rep.3',
		'B_24h_1650_1_ZKRN260013265-1A_23HVFCLT4_L7':'Community FFS-1650 late Rep.1',
		'B_24h_1650_2_ZKRN260013266-1A_23HVFCLT4_L7':'Community FFS-1650 late Rep.2',
		'B_24h_1650_3_ZKRN260013267-1A_23HVFCLT4_L7':'Community FFS-1650 late Rep.3',
		'B_24h_BHI_1_ZKRN260013253-1A_23HVFCLT4_L7':'BHI control late Rep.1',
		'B_24h_BHI_2_ZKRN260013254-1A_23HVFCLT4_L7':'BHI control late Rep.2',
		'B_24h_BHI_3_ZKRN260013255-1A_23HVFCLT4_L7':'BHI control late Rep.3',
		'B_6h_1647_1_ZKRN260013256-1A_23HVFCLT4_L7':'Community FFS-1647 early Rep.1',
		'B_6h_1647_2_ZKRN260013257-1A_23HVFCLT4_L7':'Community FFS-1647 early Rep.2',
		'B_6h_1647_3_ZKRN260013258-1A_23HVFCLT4_L7':'Community FFS-1647 early Rep.3',
		'B_6h_1650_1_ZKRN260013262-1A_23HVFCLT4_L7':'Community FFS-1650 early Rep.1',
		'B_6h_1650_2_ZKRN260013263-1A_23HVFCLT4_L7':'Community FFS-1650 early Rep.2',
		'B_6h_1650_3_ZKRN260013264-1A_23HVFCLT4_L7':'Community FFS-1650 early Rep.3',
		'B_6h_BHI_1_ZKRN260013250-1A_23HVFCLT4_L7':'BHI control early Rep.1',
		'B_6h_BHI_2_ZKRN260013251-1A_23HVFCLT4_L7':'BHI control early Rep.2',
		'B_6h_BHI_3_ZKRN260013252-1A_23HVFCLT4_L7':'BHI control early Rep.3',
		}


# Format file names
file_names_unsorted = [name for name in os.listdir("quants") if "_quant" in name]
file_set = set(file_names_unsorted)
file_names = []
for sample in file_to_conds:
	quant_name = sample + "_quant"
	if quant_name in file_set:
		file_names.append(quant_name)

# Check for files that are not in the conversion dictionary
unknown_files = []
for file_name in file_names_unsorted:
	sample = file_name.replace("_quant", "")
	if sample not in file_to_conds:
		unknown_files.append(file_name)

if unknown_files:
	print("The following files aren't provided in the conversion dictionary:")
	for file_name in sorted(unknown_files):
		print("\t" + file_name)
	exit()


# Read in the GenBank file
print('Reading GenBank file...')

records = SeqIO.parse(folder_name_genbank, "genbank")

convert_to_ref = {}
gene_types = {}
gene_refs = []

# Map Salmon locus tags directly to feature type for the QC plot
plot_gene_types = {}
none_entry_index = 0

for i,record in enumerate(list(records)):
	full_seq = record.seq
	for feature in record.features:

		if 'locus_tag' in feature.qualifiers:
			locus_tag = feature.qualifiers['locus_tag'][0]
		else:
			locus_tag = 'None_' + str(none_entry_index)
			none_entry_index += 1

		if 'gene' in feature.qualifiers:
			gene_name = feature.qualifiers['gene'][0]
		else:
			gene_name = locus_tag

		# Features used for the DE count table
		if feature.type in feature_list:
			convert_to_ref[locus_tag] = gene_name
			gene_types[locus_tag] = feature.type
			gene_refs.append(locus_tag)

		# Features used for the biotype composition plot
		if feature.type in plot_feature_list:
			plot_gene_types[locus_tag] = feature.type

# ------------------------------------------------------------------
# Create count table
# ------------------------------------------------------------------

# Read in quantification results
genes = sorted(gene_refs)
gene_indices = {}
rows = []

# Create one row for each locus_tag
for i,locus_tag in enumerate(genes):
	rows.append({"locus_tag": locus_tag, "gene_name": convert_to_ref[locus_tag], "type": gene_types[locus_tag]})
	gene_indices[locus_tag] = i

# ------------------------------------------------------------------
# Read Salmon quantification results
# ------------------------------------------------------------------

# Dictionary containing biotype counts for each sample
biotype_counts = {}

for file_name in file_names:

	file_name_mod = file_name.replace("_quant", "")
	sample_name = file_to_conds[file_name_mod]
	print(sample_name)
	
	# Initialize counts for DE table
	for row in rows:
		row[sample_name] = 0.0
	
	# Initialize counts for biotype plot
	biotype_counts[sample_name] = {}
	for gene_type in plot_feature_list:
		biotype_counts[sample_name][gene_type] = 0.0
		
	# Read counts
	with open("quants/" + file_name + "/quant.sf", "r") as in_file:
		lines = in_file.readlines()
	
	for line in lines:
		if identifier in line:
			line_split = line.split()
			gene_ident = line_split[0]
			count = float(line_split[4])
			
			# Add count to DE table if this is one of the selected DE features
			if gene_ident in genes:
				rows[gene_indices[gene_ident]][sample_name] = count
				
			# Add count to biotype composition if this is one of the selected plot features
			if gene_ident in plot_gene_types:
				gene_type = plot_gene_types[gene_ident]
				biotype_counts[sample_name][gene_type] += count
				

# ------------------------------------------------------------------
# Write final count table
# ------------------------------------------------------------------
df = pd.DataFrame(rows)
df.to_csv(out_file_name + "_count_table.csv", index=False)

# ------------------------------------------------------------------
# Biotype composition plot
# ------------------------------------------------------------------

sample_names = list(biotype_counts.keys())

plt.figure(figsize=(18, 5))

x_values = np.arange(len(sample_names))
bottom = np.zeros(len(sample_names))

handles = []

for gene_type in plot_feature_list:
	type_counts = np.array([biotype_counts[sample][gene_type] for sample in sample_names], dtype=float)
	total_counts = np.array([sum(biotype_counts[sample].values()) for sample in sample_names], dtype=float)
	fractions = type_counts / total_counts
	handle = plt.bar(x_values, fractions, bottom=bottom, width=0.8, label=gene_type)
	handles.append(handle)
	bottom += fractions

plt.legend(handles=handles, bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=8)
plt.yticks(np.arange(0, 1.1, 0.1), [str(int(x * 100)) + "%" for x in np.arange(0, 1.1, 0.1)], fontsize=8)
plt.ylabel("Fraction of counts", fontsize=10)
plt.xticks(x_values, sample_names, fontsize=8, rotation=60, ha="right", rotation_mode="anchor")
plt.xlim([-0.5, len(sample_names) - 0.5])
plt.ylim([0, 1])
plt.title(out_file_name, fontsize=11)
plt.tight_layout()
plt.savefig(out_file_name + "_biotype_composition.pdf", dpi=350, bbox_inches="tight")
plt.close()


