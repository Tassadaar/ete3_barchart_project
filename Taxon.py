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
        self.all_freq_deviation_dict = {}
        self.fymink_freq_deviation_dict = {}
        self.garp_freq_deviation_dict = {}
        self.other_freq_deviation_dict = {}

    def set_aa_freq(self):
        sorted_seq = sorted(self.seq)  # sort the sequence string alphabetically to ensure matching printing order
        all_amino_acids = "ACDEFGHIKLMNPQRSTVWY"
        self.freq_dict = {aa: sorted_seq.count(aa) / len(self.seq) for aa in all_amino_acids}  # absolute frequencies

    def set_fg_freq(self):

        for key, value in self.freq_dict.items():
            if key in "FYMINK":
                self.fymink_freq_dict[key] = value
            elif key in "GARP":
                self.garp_freq_dict[key] = value
            else:
                self.other_freq_dict[key] = value

    def set_freq_deviation(self, avg_freq_dict):

        for aa, avg_freq in avg_freq_dict.items():
            self.all_freq_deviation_dict[aa] = self.freq_dict[aa] - avg_freq

    def set_fg_freq_deviation(self, avg_freq_dict):

        for aa, avg_freq in avg_freq_dict.items():
            if aa in "FYMINK":
                self.fymink_freq_deviation_dict[aa] = self.fymink_freq_dict[aa] - avg_freq
            elif aa in "GARP":
                self.garp_freq_deviation_dict[aa] = self.garp_freq_dict[aa] - avg_freq
            else:
                self.other_freq_deviation_dict[aa] = self.other_freq_dict[aa] - avg_freq

    def __str__(self):
        return f"{self.name}\n{self.fymink_freq_dict}\n{self.garp_freq_dict}\n"
