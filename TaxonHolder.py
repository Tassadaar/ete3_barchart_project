"""
This program implements a holder for taxon attributes: name, seq, amino acid content, percentages
"""


class Taxon:
    name = None
    seq = None
    freq_dict = {}
    fymink_freq_dict = {}
    garp_freq_dict = {}
    all_face = None
    fymink_face = None
    garp_face = None

    def __init__(self, name, seq):
        self.name = name
        self.seq = seq

    def calculate_all_amino_acid_frequencies(self):
        unsorted_freq_dict = {aa: self.seq.count(aa) / len(self.seq) for aa in self.seq}
        sorted_keys = sorted(unsorted_freq_dict.keys())
        sorted_dict = {key: unsorted_freq_dict[key] for key in sorted_keys}
        sorted_dict.pop("-")
        self.freq_dict = sorted_dict

    def __str__(self):
        return f"{self.name}\n{self.seq}\n{self.freq_dict}\n"
