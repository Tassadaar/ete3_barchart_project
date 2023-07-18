from unittest import TestCase
from Taxon import Taxon


class TestTaxon(TestCase):
    def test_calculate_all_amino_acid_frequencies(self):
        name = "test"
        seq = "AACDEFFGHIKLMNPQRSTVWY"

        taxon = Taxon(name, seq)

        expected_fymink = {'F': 0.09090909090909091,
                           'I': 0.045454545454545456,
                           'K': 0.045454545454545456,
                           'M': 0.045454545454545456,
                           'N': 0.045454545454545456,
                           'Y': 0.045454545454545456}
        expected_garp = {'A': 0.09090909090909091,
                         'G': 0.045454545454545456,
                         'P': 0.045454545454545456,
                         'R': 0.045454545454545456}

        self.assertEqual(taxon.group1_freq, expected_fymink)
        self.assertEqual(taxon.group2_freq, expected_garp)

    def test_calculate_chi_square_normal(self):
        align_freqs = {}
        taxon = Taxon("Toy", "")
        taxon.seq = "A" * 40

        # generating toy datasets
        for aa in "ACDEFGHIKLMNPQRSTVWY":
            align_freqs[aa] = 1 / 20

            if aa in "ACDEFGHIKL":
                taxon.display_freqs[aa] = 3 / 40
            else:
                taxon.display_freqs[aa] = 1 / 40

        chi_square_score = taxon.calculate_chi_square(align_freqs)
        self.assertEqual(10.0, chi_square_score,
                         f"The chi square score should have been 0.25, but it is {chi_square_score}")