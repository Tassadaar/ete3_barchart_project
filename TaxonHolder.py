"""
This program implements a holder for taxon attributes: name, seq, amino acid content, percentages
"""


class Taxon:
    name = None
    seq = None
    aa_list = []
    freq_list = []

    def __str__(self):
        return f"{self.name}\n{self.seq}\n{self.aa_list}\n{self.freq_list}\n"
