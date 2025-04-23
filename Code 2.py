
from Bio import SeqIO, Phylo
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.SeqRecord import SeqRecord

# Step 1: Read the FASTA file and convert to MultipleSeqAlignment
fasta_file = "random_genomes.fasta"
records = list(SeqIO.parse(fasta_file, "fasta"))
alignment = MultipleSeqAlignment(records)

# Step 2: Print the alignment
print(alignment)

# Step 3: Calculate the distance matrix
calculator = DistanceCalculator("identity")
dist_matrix = calculator.get_distance(alignment)
print(dist_matrix)

# Step 4: Construct trees using UPGMA and Neighbor-Joining (NJ)
constructor = DistanceTreeConstructor()
upgma_tree = constructor.upgma(dist_matrix)
nj_tree = constructor.nj(dist_matrix)

# Step 5: Visualise the trees
Phylo.draw(upgma_tree, do_show=True)
Phylo.draw_ascii(nj_tree)
