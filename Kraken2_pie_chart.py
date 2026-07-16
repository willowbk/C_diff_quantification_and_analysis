import matplotlib.pyplot as plt

# Input file
report_file = "Kraken2_report.txt"

species = []
reads = []

with open(report_file) as f:
	for line in f:
		cols = line.strip().split("\t")
		if len(cols) != 6:
			continue

		rank = cols[3]
		direct_reads = int(cols[1])
		#print(rank,direct_reads)
	
	
		if rank == "S" and direct_reads >= 1000:
			species.append(cols[5])
			reads.append(direct_reads)
			print(species[-1], reads[-1])

plt.figure(figsize=(8, 8))
plt.pie(reads,
	labels=species,
	autopct="%1.1f%%",
	startangle=90
)
plt.title("Species-level Kraken2 assignments (≥1,000,000 reads)")
plt.axis("equal")
plt.tight_layout()
plt.show()
