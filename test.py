#!/usr/bin/env python

from Bio import SeqIO

# Parse the sequence alignment file
alignment_file = "input_files/nematode.fasta"
alignment = SeqIO.parse(alignment_file, "fasta")

# Initialize a dictionary to store the amino acid frequencies
aa_frequencies = {}

# Iterate over each sequence in the alignment
for record in alignment:
    # Iterate over each amino acid in the sequence
    for aa in record.seq:
        if aa not in aa_frequencies:
            aa_frequencies[aa] = 0
        aa_frequencies[aa] += 1

# Calculate the total number of amino acids in the alignment
total_aa = sum(aa_frequencies.values())

# Normalize the frequencies to obtain proportions
aa_proportions = {aa: freq/total_aa for aa, freq in aa_frequencies.items()}

# Print the amino acid frequencies
for aa, proportion in aa_proportions.items():
    print(f"{aa}: {proportion:.2f}")

