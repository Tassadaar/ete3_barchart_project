from skbio import TabularMSA, Protein

# Load the alignment file as a TabularMSA object
msa = TabularMSA.read('input_files/nematode.fasta', format='fasta')

# Remove gap characters from the alignment
msa = msa.degap()

# Convert the alignment to a list of Protein objects
proteins = [Protein(seq) for seq in msa]

# Calculate the amino acid frequencies for each position in the alignment
aa_freqs = [p.frequencies() for p in proteins]

# Print the amino acid frequencies for the first position in the alignment
print("Amino acid frequencies for position 1:" + str(aa_freqs[0]) )
