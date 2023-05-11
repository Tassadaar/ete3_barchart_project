"""
This program implements a holder for taxon attributes: name, seq, amino acid content, percentages
"""


class Taxon:
    name = None
    seq = None
    freq_dict = {}

    def __str__(self):
        return f"{self.name}\n{self.seq}\n{self.freq_dict}\n"
