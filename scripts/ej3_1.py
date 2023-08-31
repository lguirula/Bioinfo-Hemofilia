import argparse
from Bio import AlignIO
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
import matplotlib.pyplot as plt

# Create the argument parser
parser = argparse.ArgumentParser(description='Calculate conservation scores and construct a phylogenetic tree from aligned sequences.')
parser.add_argument('input_file', help='Input file containing aligned sequences in FASTA format')

# Parse the command-line arguments
args = parser.parse_args()

# Load the aligned sequences from file
alignment = AlignIO.read(args.input_file, 'fasta')

# Calculate conservation scores
conservation_scores = []
for column in range(alignment.get_alignment_length()):
    column_data = list(alignment[:, column])
    counts = column_data.count('-') + column_data.count('.')
    conservation_score = 1 - counts / alignment.get_alignment_length()
    conservation_scores.append(conservation_score)

# Print conservation scores
for position, score in enumerate(conservation_scores):
    print(f"Position {position+1}: {score}")

# Calculate sequence identities
for i, seq_record1 in enumerate(alignment):
    for j, seq_record2 in enumerate(alignment):
        if j > i:
            seq_identity = sum(a == b for a, b in zip(seq_record1.seq, seq_record2.seq)) / len(seq_record1.seq)
            print(f"Sequence {i+1} vs Sequence {j+1}: {seq_identity}")

# Tree construction
# Calculate distance matrix
calculator = DistanceCalculator('identity')
distance_matrix = calculator.get_distance(alignment)

# Construct the tree using neighbor-joining algorithm
constructor = DistanceTreeConstructor()
tree = constructor.nj(distance_matrix)

# Draw and save the tree as PNG with branch lengths as labels
fig, ax = plt.subplots(figsize=(20, 15))  # Adjust the figure size as needed

def get_protein_name(clade):
    try:
        return clade.name.split('|')[1]
    except IndexError:
        return clade.name

Phylo.draw(tree, label_func=get_protein_name, show_confidence=False, axes=ax, do_show=False)
plt.savefig("files/alignment_tree.png", dpi=300)

