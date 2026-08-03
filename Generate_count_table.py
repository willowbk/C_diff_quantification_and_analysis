
import pandas as pd
import numpy as np
from Bio import SeqIO
import os

# Path to genome file
folder_name_genbank = 'Genomes/DMG2301622_Cdiff.gbk'
identifier = 'CdiffT6'
out_file_name = 'CdiffT6_2_conditions'# This is only used to name the final table

# Only extract features of the following types
# This list is based on DMG2301622_Cdiff.gbk
feature_list = ['CDS', 'tmRNA', 'repeat_region']# Excluded tRNA and rRNA here so that they do not enter differential expression

# Create human-readable versions of the each condition
file_to_conds = {'A1_H2O_E1_ZKRN260013220-1A_23H7NJLT3_L6':'H2O rep. 1',
		'A2_H2O_E2_ZKRN260013221-1A_23HVFCLT4_L7':'H2O rep. 2',
		'A3_H2O_E3_ZKRN260013222-1A_23H7NJLT3_L6':'H2O rep. 3',
		'A10_METRO_E1_ZKRN260013229-1A_23H7NJLT3_L8':'METRO rep. 1',
		'A11_METRO_E2_ZKRN260013230-1A_23H7NJLT3_L8':'METRO rep. 2',
		'A12_METRO_E3_ZKRN260013231-1A_23H7NJLT3_L8':'METRO rep. 3',
		}
		
file_names_unsorted = [name for name in os.listdir("quants") if "_quant" in name]

# Convert to a set for fast lookup
file_set = set(file_names_unsorted)

# Order files according to the user's dictionary
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

# Read in the GenBank file again to convert gene identifiers
print('Reading GenBank file...')
records = SeqIO.parse(folder_name_genbank, "genbank")
convert_to_ref = {}
gene_types = {}
gene_refs = set([])
none_entry_index = 0
for i,record in enumerate(list(records)):
	
	full_seq = record.seq
		
	for feature in record.features:
		if feature.type in feature_list:
			
			# Detemine locus tag and gene name
			if 'locus_tag' in feature.qualifiers:
				locus_tag = feature.qualifiers['locus_tag'][0]
			else:
				locus_tag = 'None_' + str(none_entry_index)
				none_entry_index += 1
			
			if 'gene' in feature.qualifiers:
				gene_name = feature.qualifiers['gene'][0]
			else:
				gene_name = locus_tag
			
			convert_to_ref[gene_name] = locus_tag
			gene_types[gene_name] = feature.type
			gene_refs.add(gene_name)


# Read in quantification results
genes = sorted(list(gene_refs))
rows = []

# Create one row for each locus_tag
for gene in genes:
	rows.append({"locus_tag": gene, "gene_name": convert_to_ref[gene], "type": gene_types[gene]})

# Map gene_name to row index
gene_to_index = {}

for i, gene in enumerate(genes):
	gene_name = convert_to_ref[gene]
	gene_to_index[gene_name] = i


# Read each Salmon quantification file
for file_name in file_names:

	with open("quants/" + file_name + "/quant.sf", "r") as in_file:
		lines = in_file.readlines()

	file_name_mod = file_name.replace("_quant", "")
	sample_name = file_to_conds[file_name_mod]

	print(sample_name)

	# Initialize counts
	for row in rows:
		row[sample_name] = 0.0


	# Read counts
	for line in lines:
		if identifier not in line:
			continue

		line_split = line.split()
		gene_ident = line_split[0]
		count = float(line_split[4])

		# Salmon identifier is gene_name
		if gene_ident in gene_to_index:
			rows[gene_to_index[gene_ident]][sample_name] = count


# Write final count table
df = pd.DataFrame(rows)
df.to_csv(out_file_name + "_count_table.csv", index=False)

