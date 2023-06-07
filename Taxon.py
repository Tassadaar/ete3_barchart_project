"""
This program implements a holder for taxon attributes: name, seq, amino acid content, percentages, and the ability to
calculate said attributes.
"""


class Taxon:
    def __init__(self, name, seq):
        self.name = name
        self.seq = seq.replace("-", "")
        self.freq_dict = {}
        self.group1_freq = {}
        self.group2_freq = {}
        self.other_freq = {}

    # absolute frequencies
    def set_aa_abs_freq(self):
        self.freq_dict = {aa: self.seq.count(aa) / len(self.seq) for aa in "ACDEFGHIKLMNPQRSTVWY"}

    def get_aa_abs_freq(self):
        return self.freq_dict

    def get_all_relative_freq(self, avg_freq_dict):
        return {aa: self.freq_dict[aa] - avg_freq for aa, avg_freq in avg_freq_dict.items()}

    def set_subset_abs_freq(self, subsets):

        for key, value in self.freq_dict.items():
            if key in subsets[0]:
                self.group1_freq[key] = value
            elif key in subsets[1]:
                self.group2_freq[key] = value
            else:
                self.other_freq[key] = value

    def get_subset_abs_freq(self):
        return [self.group1_freq, self.group2_freq, self.other_freq]

    def get_subset_relative_freq(self, subsets, avg_freq_dict):
        group1_rel_freq = {}
        group2_rel_freq = {}
        other_rel_freq = {}

        for aa, avg_freq in avg_freq_dict.items():
            if aa in subsets[0]:
                group1_rel_freq[aa] = self.group1_freq[aa] - avg_freq
            elif aa in subsets[1]:
                group2_rel_freq[aa] = self.group2_freq[aa] - avg_freq
            else:
                other_rel_freq[aa] = self.other_freq[aa] - avg_freq

        return [group1_rel_freq, group2_rel_freq, other_rel_freq]

    def __str__(self):
        return f"{self.name}\n{self.group1_freq}\n{self.group2_freq}\n"
