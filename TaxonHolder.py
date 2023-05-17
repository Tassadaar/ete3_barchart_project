"""
This program implements a holder for taxon attributes: name, seq, amino acid content, percentages
"""


class Taxon:
    name = None
    seq = None
    freq_dict = {}
    fymink_freq_dict = {}
    garp_freq_dict = {}
    others_freq_dict = {}

    def __init__(self, name, seq):
        self.name = name
        self.seq = seq

    def calculate_all_amino_acid_frequencies(self):
        all_amino_acids = "ACDEFGHIKLMNPQRSTVWY"
        unsorted_freq_dict = {aa: self.seq.count(aa) / len(self.seq) for aa in all_amino_acids}
        sorted_keys = sorted(unsorted_freq_dict.keys())
        sorted_dict = {key: unsorted_freq_dict[key] for key in sorted_keys}

        if "-" in sorted_dict.keys():
            sorted_dict.pop("-")

        self.freq_dict = sorted_dict

        fymink = ["F", "Y", "M", "I", "N", "K"]
        garp = ["G", "A", "R", "P"]

        for key, value in sorted_dict.items():
            if key in fymink:
                self.fymink_freq_dict[key] = value
            elif key in garp:
                self.garp_freq_dict[key] = value
            else:
                self.others_freq_dict[key] = value

    def __str__(self):
        return f"{self.name}\n{self.fymink_freq_dict}\n{self.garp_freq_dict}\n"
