
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from itertools import zip_longest
from tqdm import tqdm
import gzip
import os

###################################################################################################
# 
# This script takes a FastQ file which requires 
# demultiplexing based on a set of custom barcodes.
# This script will only search for barcodes at the exact
# start of the read, and will allow for single nucleotide
# mismatches. 
#
###################################################################################################

################################ User parameters ##################################################
file_path = ''# Location of index file, Genomes, quants and fastq folders
fastq_file = 'fastq_file_name.fq.qz'# Name of FastQ file to be demultiplexed
fastq_file = 'B_24h_1647_3_ZKRN260013261-1A_23HVFCLT4_L7_1.fq.gz'

df = pd.read_excel(file_path + 'Barcodes_example.xlsx')# Comment line when hard-coding barcodes
barcode_header = 'Barcodes'# Column label in .xlsx file where barcodes can be found
barcode_seqs = [str(term) for term in list(df[barcode_header]) if not str(term) == 'nan']
# Example set of hard-coded barcodes
# Above lines must be commented-out when using a hard-coded set like the one below.
#barcode_seqs = ['CAAGTGATT', 'CGACTTGGT', 'GTCCATATT', 'AATAATGTT', 'ATAATTCTT', 'CAACACTTT']
###################################################################################################

folder_path = file_path + "fastq"
fastq_file_prefix = fastq_file.split('.')[0]
fastq_file_suffix = '.' + '.'.join(fastq_file.split('.')[1:])

# Build lookup table for all barcodes and all possible single-nucleotide mismatches
barcode_lookup = {}
barcode_length = len(barcode_seqs[0])

for barcode in barcode_seqs:
	barcode_lookup[barcode] = barcode

	for i in range(barcode_length):
		for nuc in 'ATGCN':
			if nuc != barcode[i]:
				variant = barcode[:i] + nuc + barcode[i+1:]
				barcode_lookup[variant] = barcode

# Count total reads
with gzip.open(folder_path + "/" + fastq_file,'rt') as fin:
	total_lines = sum(1 for line in fin)

total_reads = total_lines // 4

# Open output files
output_files = {}

for barcode in barcode_seqs:
	output_files[barcode] = gzip.open(folder_path + "/" + fastq_file_prefix + "_" + barcode + fastq_file_suffix,'wt')

output_files['None'] = gzip.open(folder_path + "/" + fastq_file_prefix + "_None" + fastq_file_suffix,'wt')

with gzip.open(folder_path + "/" + fastq_file,'rt') as fin:
	barcode_counts = {}
	read_counts = {}
	seq_count = 0

	for label_line,seq,plus_line,qual_line in tqdm(zip(*[iter(fin)]*4),total=total_reads,unit='reads',desc='Processing reads'):
		
		seq_count += 1
		seq = seq.rstrip()
		qual_line = qual_line.rstrip()

		# Identify barcode
		found_barcode = seq[:barcode_length]# Actual sequence found
		ident_barcode = barcode_lookup.get(found_barcode,'None')# Identified barcode

		# Replace barcodes with Ns which will be removed after adapter trimming
		if not ident_barcode == 'None':
			output_seq = 'N'*len(ident_barcode) + seq[len(ident_barcode):]
		else:
			output_seq = seq

		# Write demultiplexed read
		output_files[ident_barcode].write(label_line)
		output_files[ident_barcode].write(output_seq + '\n')
		output_files[ident_barcode].write(plus_line)
		output_files[ident_barcode].write(qual_line + '\n')

		# Gather statistics on matched and mismatched barcodes 
		if ident_barcode not in barcode_counts: 
			barcode_counts[ident_barcode] = {} 
		
		if found_barcode in barcode_counts[ident_barcode]:
			barcode_counts[ident_barcode][found_barcode] += 1
		else: 
			barcode_counts[ident_barcode][found_barcode] = 1

		# Count reads per barcode
		if ident_barcode in read_counts:
			read_counts[ident_barcode] += 1
		else:
			read_counts[ident_barcode] = 1
		
	

# Close output files
for barcode in output_files:
	output_files[barcode].close()

print('From file:', fastq_file)
print('Total reads processed:', seq_count)

# Generate barcode plots
figure_folder = file_path + "Barcode_stats_figures"

os.makedirs(figure_folder,exist_ok=True)

for barcode in barcode_seqs:

	exact_count = barcode_counts[barcode].get(barcode,0)
	total_count = exact_count
	mismatch_counts = []

	for found_barcode in barcode_counts[barcode]:
		if found_barcode != barcode:
			mismatch_counts.append((found_barcode,barcode_counts[barcode][found_barcode]))
			total_count += barcode_counts[barcode][found_barcode]
			
	mismatch_counts.sort(key=lambda x: x[1],reverse=True)

	labels = [barcode] + [x[0] for x in mismatch_counts]
	counts = [exact_count] + [x[1] for x in mismatch_counts]

	colors = ['lightcoral'] + ['lightblue']*len(mismatch_counts)

	plt.figure(figsize=(10,6))
	plt.bar(labels,counts,color=colors)

	plt.xlabel('Barcode sequence')
	plt.ylabel('Number of reads')
	plt.title('Reads identified as barcode ' + barcode + '; ' + str(total_count) + ' total reads')

	plt.xticks(rotation=30,ha='right')
	plt.tight_layout()

	plt.savefig(figure_folder + "/" + fastq_file_prefix + "_" + barcode + "_barcode_counts.png",dpi=300)
	plt.close()


