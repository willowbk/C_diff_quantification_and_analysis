
from Bio import SeqIO

# Path to folder with all .gbk or .ffn files are found
file_path_1 = 'Genomes/'
identifier = 'CdiffT6'
folder_name_genbank = file_path_1 + 'DMG2301622_Cdiff.gbk'

# Only extract features of the following types
# This list is based on DMG2301622_Cdiff.gbk
feature_list = ['CDS', 'tRNA', 'tmRNA', 'rRNA', 'repeat_region']

# Collect all feature types found in record
gene_types = []
new_genes = []
none_entry_index = 0

print('Reading GenBank file...')
records = SeqIO.parse(folder_name_genbank, "genbank")

for i,record in enumerate(list(records)):
	
	full_seq = record.seq
		
	for feature in record.features:
		gene_types.append(feature.type)
		
		if feature.type in feature_list:
			
			gene = []
			
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
			
			gene.append(gene_name)
			gene.append(locus_tag)
			
			# Extract db_xref entry
			if 'db_xref' in feature.qualifiers:
				gene.append(feature.qualifiers['db_xref'])
			else:
				gene.append('None')
			
			# Extract and format feature location
			gene.append(str(feature.location).split(':')[0][1:])
			gene.append(str(feature.location).split(':')[1].split(']')[0])
			gene.append(str(feature.location)[-2].replace('-', '<').replace('+', '>'))
			
			# Find an entry for gbkey
			if 'gbkey' in feature.qualifiers: 
				gene.append(feature.qualifiers['gbkey'])
			elif 'gene_biotype' in feature.qualifiers:
				gene.append(feature.qualifiers['gene_biotype'])
			else:
				gene.append('None')
			
			# Get full nucleotide sequence
			seq = str(feature.extract(full_seq))
			
			gene.append(seq)
			new_genes.append(gene)

# print all feature types found in record
unique_feature_types = list(set(gene_types))
gene_types = sorted([[gene_types.count(feature_type), feature_type] for feature_type in unique_feature_types])[::-1]
print('Found features and counts:')
for term in gene_types:
	print('\t', term[1], term[0])

print()
print('# features to write to .fa file:', len(new_genes))
# generate .fa file
out_file = open(folder_name_genbank.replace('.gbff', '.fa').replace('.gbk', '.fa').replace('.gb', '.fa'), 'w')

gene_counter = 0
for gene in new_genes:
	gene_counter += 1
	
	line = '>lcl|'+identifier+'_gene_'+str(gene_counter)+' [gene='+gene[0].replace(' ', '')+'] [locus_tag='+gene[1]+'] [db_xref='+ str(gene[2])+\
		'] [location='+'complement('*int(gene[5] == '<')+str(gene[3])+'..'+str(gene[4])+')'*int(gene[5] == '<')+'] [gbkey='+str(gene[6])+']\n'
		
	out_file.write(line)
	
	# break sequences which are longer than 70 chars.
	if len(gene[7]) <= 70:
		line = gene[7] + '\n'
		out_file.write(line)
		
	else:
		for i in range(int(len(gene[7])/70) + 1):
			line = gene[7][i*70:(i+1)*70] + '\n'
			out_file.write(line)
		
out_file.close()

