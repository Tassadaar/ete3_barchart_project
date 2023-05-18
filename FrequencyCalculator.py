"""
This program calculates amino acid frequencies from a give multiple alignment

Returns: list(taxa)
"""

from Bio import SeqIO
from TaxonHolder import Taxon


class FrequencyCalculator:

    def __init__(self, file_location, file_format):
        self.file_location = file_location
        self.file_format = file_format

    def get_taxa_list(self, mode):
        taxa_dict = {}  # dictionary to store taxa

        # read in fasta and parse, then update
        for seq_record in SeqIO.parse(self.file_location, self.file_format):
            new_taxon = Taxon(seq_record.id, seq_record.seq)
            new_taxon.calculate_all_amino_acid_frequencies()
            if mode == "special":
                new_taxon.calculate_fymink_garp_frequencies()
            taxa_dict[seq_record.id] = new_taxon

        return taxa_dict
