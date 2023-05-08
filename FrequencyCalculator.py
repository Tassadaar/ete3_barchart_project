"""
This program calculates amino acid frequencies from a give multiple alignment

Returns: tuple([id], ([aa], [frequencies]))
"""

from Bio import SeqIO
from CustomProtParam import CustomProteinAnalysis


class FrequencyCalculator:

    def __init__(self, file_location, file_format):
        self.file_location = file_location
        self.file_format = file_format

    def calculate_frequency(self):
        alignments = []  # array to store each alignment as objects
        freq_record = []  # tuple to store alignment name, aa content and respective frequencies

        # read in fasta and parse
        for seq_record in SeqIO.parse(self.file_location, self.file_format):
            alignments.append(seq_record)

        # output amino acid frequencies
        for seq_record in alignments:
            # make a new tuple from a tuple of two lists
            aa_frequencies = (seq_record.id, CustomProteinAnalysis(seq_record.seq).get_amino_acids_percent())
            freq_record.append(aa_frequencies)

        return freq_record
