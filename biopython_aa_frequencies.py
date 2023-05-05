"""
This program calculates amino acid frequencies from a give multiple alignment
"""

from Bio import SeqIO
from CustomProtParam import CustomProteinAnalysis

# array to store each alignment as objects
alignments = []

# read in fasta and parse
for seq_record in SeqIO.parse("input_files/nematode.fasta", "fasta"):
    alignments.append(seq_record)

# output amino acid frequencies
for alignment in alignments:
    aa_frequencies = CustomProteinAnalysis(alignment.seq).get_amino_acids_percent()
    print(aa_frequencies)
