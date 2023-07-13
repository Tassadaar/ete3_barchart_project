"""
This program implements a holder for taxon attributes: name, seq, amino acid content, percentages, and the ability to
calculate said attributes.
"""


class Taxon:
    all_amino_acids = "ACDEFGHIKLMNPQRSTVWY"

    def __init__(self, name, seq):
        self.name = name
        self.seq = seq.replace("-", "")
        self.freq_dict = {aa: self.seq.count(aa) / len(self.seq) for aa in self.all_amino_acids}
        self.group1_freq = {}
        self.group2_freq = {}
        self.other_freq = {}

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

    def calculate_chi_square(self, align_freqs):
        chi_square_score = 0
        taxon_seq_len = len(self.seq)

        for aa in self.all_amino_acids:
            expected_count = align_freqs[aa] * taxon_seq_len
            observed_count = self.freq_dict[aa] * taxon_seq_len
            chi_square_score += (observed_count - expected_count) ** 2 / expected_count

        return round(chi_square_score, 1)

