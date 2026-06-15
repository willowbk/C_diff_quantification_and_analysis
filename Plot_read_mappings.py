
import matplotlib.pyplot as plt
plt.style.use('classic')
plt.rcParams['svg.fonttype'] = 'none' 
import matplotlib as mpl
mpl.rc('font',family='sans serif')
import matplotlib.patches as mpatches
import numpy as np

folder_name = '/media/willow/MyPassport/Faber/2025-06-10_Janet_Wackenreuter_25046PR/FASTQ_by_lib/'


locs = [[16116, 17606], [17667, 17742], [18120, 21015], [21145, 21260]]
loc = [locs[0][0], locs[-1][1]]
tag = 'Pool9_Pilot_rRNA'
genes = ['16S rRNA', 'tRNA-Ala', '23S rRNA', '5S rRNA']
colors = [(.5, .5, 1), (.75, .75, .75), (1, .5, .5), (.5, .25, .5)]
pre_file = 'PlateIII_2_Pool9new_out_'
post_file = '_sorted_16000_21500.sam'
strands = ['+', '+', '+', '+']
length = loc[1] - loc[0]
extend = 500
lib_sizes = [36174138]
#PlateIII_2_Pool9new_out_None_sorted_16000_21500.bam

max_read_ratio = 0
reads_1_for = np.zeros(length + 2*extend)
reads_1_rev = np.zeros(length + 2*extend)

term = 'None'

file_name = pre_file + term + post_file
in_file = open(folder_name + file_name, 'r')
lines = in_file.readlines() 
in_file.close()
print(term)
i = 0
for line in lines:
	if not line[0] == '@':
		line_split = line.split()
		strand = '+'*int(line_split[1] == '16') + '-'*int(line_split[1] == '0')
		pos = int(line_split[3])
		if strand == '+':
			reads_1_for = reads_1_for + 10**6 * np.array([int(i - extend >= pos - loc[0])*int(i - extend < pos - loc[0] + len(line_split[9])) for i in range(len(reads_1_for))]) / lib_sizes[i]
		elif strand == '-':
			reads_1_rev = reads_1_rev + 10**6 * np.array([int(i - extend >= pos - loc[0])*int(i - extend < pos - loc[0] + len(line_split[9])) for i in range(len(reads_1_rev))]) / lib_sizes[i]

folder_name = '/media/willow/MyPassport/Faber/2024-12-13_Janet_Wackenreuter_24195PR_raw_FASTQ/FASTQ_by_lib/'

locs = [[16116, 17606], [17667, 17742], [18120, 21015], [21145, 21260]]
loc = [locs[0][0], locs[-1][1]]
genes = ['16S rRNA', 'tRNA-Ala', '23S rRNA', '5S rRNA']
colors = [(.5, .5, 1), (.75, .75, .75), (1, .5, .5), (.5, .25, .5)]
pre_file = 'Pool1_RNAtagseq_'
post_file = '_16000_21500.sam'
strands = ['+', '+', '+', '+']
length = loc[1] - loc[0]
extend = 500
lib_sizes = [5446362]

max_read_ratio = 0
reads_2_for = np.zeros(length + 2*extend)
reads_2_rev = np.zeros(length + 2*extend)

term = 'AATAATGTT'

file_name = pre_file + term + post_file
in_file = open(folder_name + file_name, 'r')
lines = in_file.readlines() 
in_file.close()
print(term)
i = 0
for line in lines:
	if not line[0] == '@':
		line_split = line.split()
		strand = '+'*int(line_split[1] == '16') + '-'*int(line_split[1] == '0')
		pos = int(line_split[3])
		if strand == '+':
			reads_2_for = reads_2_for + 10**6 * np.array([int(i - extend >= pos - loc[0])*int(i - extend < pos - loc[0] + len(line_split[9])) for i in range(len(reads_2_for))]) / lib_sizes[i]
		elif strand == '-':
			reads_2_rev = reads_2_rev + 10**6 * np.array([int(i - extend >= pos - loc[0])*int(i - extend < pos - loc[0] + len(line_split[9])) for i in range(len(reads_2_rev))]) / lib_sizes[i]

fig = plt.figure(figsize=(3.5, 2.5))

handle1 = plt.bar([i-extend for i in range(len(reads_1_for))], reads_1_for, width=1, linewidth = 0, color= 'lightblue', alpha = .5, label='Pool 9, forward')
handle2 = plt.bar([i-extend for i in range(len(reads_1_rev))], -reads_1_rev, width=1, linewidth = 0, color= (1, .5, 1), alpha = .5, label='Pool 9, reverse')
handle3 = plt.bar([i-extend for i in range(len(reads_2_for))], reads_2_for, width=1, linewidth = 0, color= 'blue', alpha = .5, label='Pilot AATAATGTT, forward')
handle4 = plt.bar([i-extend for i in range(len(reads_2_rev))], -reads_2_rev, width=1, linewidth = 0, color= (1, 0, 1), alpha = .5, label='Pilot AATAATGTT, reverse')


track_pos = -.4*max(list(reads_1_for) + list(reads_1_for)) - max(list(reads_2_rev) + list(reads_1_rev))
track_height = .15*max(list(reads_1_for) + list(reads_2_for))

for i,gene in enumerate(genes):
	gene_len = locs[i][1] - locs[i][0]
	gene_pos = locs[i][0] - locs[0][0]
	strand = strands[i]
	color = colors[i]
	if strand == '+':
		plt.arrow(gene_pos, track_pos+track_height/2, .9*gene_len, 0, width = track_height, head_width=track_height, head_length = .1*gene_len, ec="black", fc = color, alpha = 1, clip_on=False, zorder=3)
	else:
		plt.arrow(gene_pos + gene_len, track_pos+track_height/2, -.9*gene_len, 0, width = track_height, head_width=track_height, head_length = .1*gene_len, ec="black", fc = color, alpha = 1, clip_on=False, zorder=3)
		
	plt.text(gene_pos + .5*gene_len, track_pos+.5*track_height, gene, fontsize = 4, horizontalalignment='center', verticalalignment='center')

plt.axhline(0, 0, 1, linestyle = 'solid', lw = .5, color = 'black')
plt.xticks(fontsize = 6)
#yticks = sorted([-i for i in range(5) if i < max(read_ratio_rev)]) + [i for i in range(1, 5) if i < max(read_ratio_for)]
plt.yticks(fontsize = 6)
plt.axhline(track_pos + .5*track_height, 0, 1, lw = 1, linestyle = (0, (1, 1)), color = 'darkgrey', clip_on = False)
plt.axhline(track_pos + .5*track_height, 0.035, .965, lw = 1, color = 'black', clip_on = False)
plt.legend(handles=[handle1, handle2, handle3, handle4], loc = 'best', fontsize = 4)
ax = plt.gca()
ax.spines[['right', 'top', 'bottom']].set_visible(False)
ax.tick_params(right=False)
ax.tick_params(top=False)
#ax.tick_params(left=False)
#ax.tick_params(labelleft=False)
ax.tick_params(labelbottom=False)
ax.tick_params(bottom=False)
#plt.xlabel('Read position (bps)', fontsize = 6)
plt.ylabel('Read coverage', fontsize = 8)
#ax.add_artist(mpatches.Rectangle((0, -.3), length, .5, ec="none", fc = (1, .75, .75), alpha = 1, clip_on=False))
#ax.add_artist(mpatches.Rectangle((0, -.4), length, 0.1, ec="none", fc = (1, .5, .5), alpha = 1, clip_on=False, zorder=3))

plt.xlim([-extend, length + extend])
plt.ylim([-1.05*max(list(reads_1_rev) + list(reads_2_rev)), 1.05*max(list(reads_1_for) + list(reads_2_for))])
#plt.set_axis_off()
plt.tight_layout()
plt.savefig('Figure_read_coverage_'+tag+'.pdf',dpi = 350)
plt.clf()

