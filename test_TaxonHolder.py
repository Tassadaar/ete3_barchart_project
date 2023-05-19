from unittest import TestCase
from TaxonHolder import Taxon


class TestTaxon(TestCase):
    def test_calculate_all_amino_acid_frequencies(self):
        name = "test"
        seq = "AACDEFFGHIKLMNPQRSTVWY"

        taxon = Taxon(name, seq)
        taxon.calculate_all_amino_acid_frequencies()

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

        self.assertEqual(taxon.fymink_freq_dict, expected_fymink)
        self.assertEqual(taxon.garp_freq_dict, expected_garp)