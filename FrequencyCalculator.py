"""
This program calculates amino acid frequencies from a give multiple alignment

Returns: list(taxa)
"""

from Bio import SeqIO
from CustomProtParam import CustomProteinAnalysis
from TaxonHolder import Taxon


class FrequencyCalculator:

    def __init__(self, file_location, file_format):
        self.file_location = file_location
        self.file_format = file_format

    def calculate_frequency(self):
        taxa_list = []  # array to store taxa

        # read in fasta and parse
        for seq_record in SeqIO.parse(self.file_location, self.file_format):
            new_taxon = Taxon()
            new_taxon.name = seq_record.id
            new_taxon.seq = seq_record.seq
            taxa_list.append(new_taxon)

        # output amino acid frequencies
        for new_taxon in taxa_list:
            # make a new tuple from a tuple of two lists
            aa_frequencies = CustomProteinAnalysis(new_taxon.seq).get_amino_acids_percent()
            new_taxon.aa_list = aa_frequencies[0]
            new_taxon.freq_list = aa_frequencies[1]

        return taxa_list
