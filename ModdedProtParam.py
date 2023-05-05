"""
This program implements changes made to class ProteinAnalysis in python 3.11 to ensure functionality
when using an older version of python for dependency.
"""

from Bio.SeqUtils.ProtParam import ProteinAnalysis


class ModdedProteinAnalysis(ProteinAnalysis):

    def __init__(self, prot_sequence, monoisotopic=False):
        """Initialize the class."""
        self.sequence = prot_sequence.upper()
        self.amino_acids_content = None
        self.amino_acids_percent = None
        self.length = len(self.sequence)
        self.monoisotopic = monoisotopic

    def get_amino_acids_percent(self):
        """Calculate the amino acid content in percentages.

        The same as count_amino_acids only returns the Number in percentage of
        entire sequence. Returns a dictionary of {AminoAcid:percentage}.

        The return value is cached in self.amino_acids_percent.

        input is the dictionary self.amino_acids_content.
        output is a dictionary with amino acids as keys.
        """
        if self.amino_acids_percent is None:
            aa_counts = self.count_amino_acids()

            percentages = {aa: count / self.length for aa, count in aa_counts.items()}

            self.amino_acids_percent = percentages

        return self.amino_acids_percent
