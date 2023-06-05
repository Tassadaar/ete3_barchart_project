"""
This program implements a holder for taxon attributes: name, seq, amino acid content, percentages, and the ability to
calculate said attributes.
"""


class Taxon:

    def __init__(self, name, seq):
        self.name = name
        self.seq = seq.replace("-", "")
        self.freq_dict = {}
        self.fymink_freq_dict = {}
        self.garp_freq_dict = {}
        self.other_freq_dict = {}

    def set_aa_abs_freq(self):
        all_amino_acids = "ACDEFGHIKLMNPQRSTVWY"
        self.freq_dict = {aa: self.seq.count(aa) / len(self.seq) for aa in all_amino_acids}  # absolute frequencies

    def set_fg_abs_freq(self):

        for key, value in self.freq_dict.items():
            if key in "FYMINK":
                self.fymink_freq_dict[key] = value
            elif key in "GARP":
                self.garp_freq_dict[key] = value
            else:
                self.other_freq_dict[key] = value

    def get_fg_abs_freq(self):
        return [self.fymink_freq_dict, self.garp_freq_dict, self.other_freq_dict]

    def get_all_relative_freq(self, avg_freq_dict):
        relative_freq_dict = {}

        for aa, avg_freq in avg_freq_dict.items():
            relative_freq_dict[aa] = self.freq_dict[aa] - avg_freq

        return relative_freq_dict

    def get_fg_relative_freq(self, avg_freq_dict):
        fymink_relative_freq_dict = {}
        garp_relative_freq_dict = {}
        other_relative_freq_dict = {}

        for aa, avg_freq in avg_freq_dict.items():
            if aa in "FYMINK":
                fymink_relative_freq_dict[aa] = self.fymink_freq_dict[aa] - avg_freq
            elif aa in "GARP":
                garp_relative_freq_dict[aa] = self.garp_freq_dict[aa] - avg_freq
            else:
                other_relative_freq_dict[aa] = self.other_freq_dict[aa] - avg_freq

        return [fymink_relative_freq_dict, garp_relative_freq_dict, other_relative_freq_dict]

    def __str__(self):
        return f"{self.name}\n{self.fymink_freq_dict}\n{self.garp_freq_dict}\n"
