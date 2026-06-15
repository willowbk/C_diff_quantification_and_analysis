#######################################################
# 
# This script takes a single run, runs adapter trimming,
# FastQC for quality control checks, and then 
# quantifies the reads based on Salmon.
# In this code, the following versions were used:
#
#	BBDuk (BBMap version 39.01)
#	FastQC (v0.11.9)
#	Salmon (salmon 1.5.2)
#
#######################################################

import os
import sys
import pandas as pd
import zipfile
import numpy as np

file_path = ''# Location of index file, Gemomes, quants and fastq folders
index_file = 'DMG2301622_Cdiff_ffn_index'

########################################################
# To create the index file from .ffn, run
# salmon2 index -t Genomes/DMG2301622_Cdiff.ffn -i DMG2301622_Cdiff_index
########################################################

########################################################
# Fastq file parameters
########################################################

run_id = "A10_METRO_E1_ZKRN260013229-1A_23H7NJLT3_L8"
lib_la = "PAIRED"# PAIRED or SINGLE
folder_path = file_path + "fastq/" + run_id
print(f"Run ID: {run_id}, library layout: {lib_la}")

########################################################
# Check if quantification is already complete for this run. 
# If so, then exit.
########################################################

exists = os.path.isfile(file_path +'quants/' + run_id + '_quant/quant.sf')
print(file_path +'quants/' + run_id + '_quant/quant.sf', ' File exists: ', exists)
if exists:
	print(run_id + ' already quantified. Exiting...')
	exit()

########################################################
# Adaptr trim fastq files
########################################################

exists = os.path.isfile(folder_path + "/" + run_id + "_1_out.fq.gz")
if not exists:
	if lib_la == 'PAIRED':
		# For paired end reads...
		# Run adapter trimming from BBDuk
		# To keep this general, adapter and quality trimming are down 
		# from both ends of each read: ktrim=rl qtrim=rl. 
		# mink=11 was chosen to balance sensitivity and avoiding false matches
		# trimq=25 is possibly a little stringent (Phred 25 → ~0.3% error probability per base)
		# 	such that only starting from the first base with a quality > 0.3% will a read be kept
		# minlen=30 avoids including reads which are unlikely to map
		# Note that the trimmed files will replace the original fastq files.
		
		os.system("bbduk.sh in1=" + folder_path + "/" + run_id + "_1.fq.gz"+\
		 		         " in2=" + folder_path + "/" + run_id + "_2.fq.gz"+\
		 		         " out1=" + folder_path + "/" + run_id + "_1_out.fq.gz"+\
		 		         " out2=" + folder_path + "/" + run_id + "_2_out.fq.gz ktrim=rl qtrim=rl mink=11 trimq=25 minlen=30 ref=adapters")
		 
		
	elif lib_la == 'SINGLE':
		# For single end reads...
		# Run adapter trimming from BBDuk
		# See the section on paired end reads for details on the 
		# parameter choices.
		
		os.system("bbduk.sh in=" + folder_path + "/" + run_id + ".fq.gz out=" + folder_path + "/" + run_id + "_out.fq.gz ktrim=rl qtrim=rl mink=11 trimq=25 minlen=30 ref=adapters")


########################################################
# Perform QC on trimmed fastq files.
########################################################

exists = os.path.isfile(file_path + run_id + '_QC_results.txt')
if not exists:
	# First QC step: FastQC
	# These following three QC flags were chosen based on following the 
	# pipeline at https://github.com/SBRG/iModulonMiner/tree/main
	
	fastqc_fail_cols = ['Per base sequence quality', 'Per sequence quality scores', 'Per base N content']
	if lib_la == 'PAIRED':
		sample_pass = [[True], [True]]
		fastq_files = [folder_path + "/" + run_id + "_1_out.fq.gz", folder_path + "/" + run_id + "_2_out.fq.gz"]
		
	elif lib_la == 'SINGLE':
		sample_pass = [[True]]
		fastq_files = [folder_path + "/" + run_id + "_out.fq.gz"]
	
	# Check if fastqc file is already generated.
	# If not, run FastQC (input type fastq, and use 4 threads)
	
	for i,fq_file in enumerate(fastq_files):
		result_file = fq_file.split('.')[0] + '_fastqc.zip'
		
		exists = os.path.isfile(result_file)
		if not exists:
			os.system('fastqc -f fastq -t 4 '+fq_file)
		
		exists = os.path.isfile(result_file)
		if not exists:
			print("FastQC failed. Exiting...")
			exit()
		
		# If FastQC ran successfully, read in the summary file
		# for each of the fastq files.
		
		with zipfile.ZipFile(result_file) as z:
			sum_file = z.open(fq_file.split('.')[0].split('/')[-1] + '_fastqc/summary.txt')
			for line in sum_file.readlines():
				line_mod = str(line).split('\t')[0].split('\\t')
				sample_pass[i].append([line_mod[0][2:], line_mod[1]])
				if line_mod[1] in fastqc_fail_cols:
					sample_pass[i][-1].append('*')
					if not line_mod[0][2:] == 'PASS':
						sample_pass[i][0] = False
	
	# Only run quantification of all checks passed.
	
	pass_1st_QC = np.prod([int(term[0]) for term in sample_pass]) == 1
	if pass_1st_QC:
		print(run_id + ' passed QC check')
	else:
		print(run_id + ' did NOT pass QC check')

	# Write QC summary file. 
	# This file only includes the checks 
	# that FastQC includes.
	
	out_file = open(file_path + run_id + '_QC_results.txt', 'w')
	out_file.write('FastQC\n' + run_id+'_1'*int(lib_la == 'PAIRED') + '\t' + str(sample_pass[0][0]) + '\n')
	for term in sample_pass[0][1:]:
		out_file.write('\t'.join(term) + '\n')
	out_file.write('\n')

	if lib_la == 'PAIRED':
		out_file.write(run_id+'_2\t' + str(sample_pass[1][0]) + '\n')
		for term in sample_pass[1][1:]:
			out_file.write('\t'.join(term) + '\n')
			
		out_file.write('\n')

	out_file.close()
else:
	# If the QC results file already exists,
	# skip running FastQC and simply check results.
	
	print(run_id + ' already checked for QC. ')
	in_file = open(file_path + run_id + '_QC_results.txt')
	lines = in_file.readlines()
	in_file.close()
	sample_pass = []
	
	for line in lines:
		line_split = line.split()
		if len(line_split) == 2:
			if line_split[0].split('_')[0] == run_id:
				sample_pass.append(bool(line_split[1] == 'True'))
				
	pass_1st_QC = np.prod([int(term) for term in sample_pass]) == 1
	
########################################################
# Run quantification if QC checks passed.
########################################################

if pass_1st_QC:
	print('QC check passed. Salmon will run.')
	
	# Run quantification:
	#
	# -l A sets the strandedness to be automatically detected
	# 	since this may differ between experiments.
	#
	# --maxSoftclipFraction 0 0.8 Controls how tolerant Salmon 
	#	is to soft-clipping (bases that align poorly at read ends)
	#	This means reads can align even if up to 80% of bases are clipped.
	#	This high value was recommended by the developers of Salmon.
	#
	# -p 4 Number of threads set to 4.
	#
	
	if lib_la == 'PAIRED':
		os.system('salmon2 quant -i ' + file_path + index_file + ' -l A --maxSoftclipFraction 0 0.8' +\
	         				' -1 ' + file_path + 'fastq/' + run_id + '/' + run_id + '_1_out.fq.gz' +\
	         				' -2 ' + file_path + 'fastq/' + run_id + '/' + run_id + '_2_out.fq.gz' +\
	         				' -p 4 --validateMappings -o ' + file_path + 'quants/' + run_id + '_quant')
	elif lib_la == 'SINGLE':
		os.system('salmon2 quant -i ' + file_path + index_file + ' -l A --maxSoftclipFraction 0 0.8' +\
	         				' -r ' + file_path + 'fastq/' + run_id + '/' + run_id + '_out.fq.gz' +\
	         				' -p 4 --validateMappings -o ' + file_path + 'quants/' + run_id + '_quant')
else:
	print('Salmon will not be run: QC checks did not pass.')


