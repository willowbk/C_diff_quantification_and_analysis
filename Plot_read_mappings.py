
import matplotlib.pyplot as plt
plt.style.use('classic')
plt.rcParams['svg.fonttype'] = 'none' 
import matplotlib as mpl
mpl.rc('font',family='sans serif')
import matplotlib.patches as mpatches
import numpy as np

folder_name = 'fastq/'

# Hard-coded entries for rRNA locus
locs = [[2971432, 2971537], [2971674, 2974567], [2974834, 2976333], [2976529, 2977107]]
loc = [locs[0][0], locs[-1][1]]
genes = ['5S rRNA', '16S rRNA', '23S rRNA', 'cysG_1']
colors = [(.5, .25, .5), (.5, .5, 1), (1, .5, .5), (.75, .75, .75)]
strands = ['-', '-', '-', '-']
length = loc[1] - loc[0]
extend = 500


tag = 'A10_METRO_E1_rRNA_locus_2971000_2977000'# only for figure naming
file_name = 'mapped_reads_sorted_rRNA_operon.sam'
lib_sizes = [18038030]# Obtained previously via samtools


max_read_ratio = 0
reads_1_for = np.zeros(length + 2*extend)
reads_1_rev = np.zeros(length + 2*extend)

in_file = open(folder_name + file_name, 'r')
lines = in_file.readlines() 
in_file.close()

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

fig = plt.figure(figsize=(3.5, 2.5))

handle1 = plt.bar([i-extend for i in range(len(reads_1_for))], reads_1_for, width=1, linewidth = 0, color= 'lightblue', alpha = .5, label='Forward')
handle2 = plt.bar([i-extend for i in range(len(reads_1_rev))], -reads_1_rev, width=1, linewidth = 0, color= (1, .5, 1), alpha = .5, label='Reverse')

track_pos = -.4*max(list(reads_1_for)) - max(list(reads_1_rev))
track_height = .15*max(list(reads_1_for))

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

plt.yticks(fontsize = 6)
plt.axhline(track_pos + .5*track_height, 0, 1, lw = 1, linestyle = (0, (1, 1)), color = 'darkgrey', clip_on = False)
plt.axhline(track_pos + .5*track_height, 0.035, .965, lw = 1, color = 'black', clip_on = False)
plt.legend(handles=[handle1, handle2, handle3, handle4], loc = 'best', fontsize = 4)
ax = plt.gca()
ax.spines[['right', 'top', 'bottom']].set_visible(False)
ax.tick_params(right=False)
ax.tick_params(top=False)

ax.tick_params(labelbottom=False)
ax.tick_params(bottom=False)

plt.ylabel('Read coverage', fontsize = 8)

plt.xlim([-extend, length + extend])
plt.ylim([-1.05*max(list(reads_1_rev)), 1.05*max(list(reads_1_for))])
#plt.set_axis_off()
plt.tight_layout()
plt.savefig('Figure_read_coverage_'+tag+'.pdf',dpi = 350)
plt.clf()

