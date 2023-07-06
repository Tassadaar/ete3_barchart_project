from unittest import TestCase
from ete3 import Tree
from TreeMaker import root, calculate_chi_square


class Test(TestCase):

    # Scenario 1: the tree is unrooted
    def test_root_scenario_1(self):
        input_tree = Tree("input_files/test.tree")
        input_tree.unroot()
        result_tree = Tree("input_files/result.tree")
        root(input_tree, ["BRINESHR", "CHICKEN"])
        rf_distance = input_tree.robinson_foulds(result_tree)[0]
        self.assertEqual(0, rf_distance, f"The distance should have been 0, but it is {rf_distance}")

    # Scenario 2: the tree is rooted but only one child contains outgroups
    def test_root_scenario_2(self):
        input_tree = Tree("input_files/test.tree")
        result_tree = Tree("input_files/result.tree")
        root(input_tree, ["BRINESHR", "CHICKEN"])
        rf_distance = input_tree.robinson_foulds(result_tree)[0]
        self.assertEqual(0, rf_distance, f"The distance should have been 0, but it is {rf_distance}")

    # Scenario 3: the tree is rooted but both children contain outgroups
    def test_root_scenario_3(self):
        input_tree = Tree("input_files/test.tree")
        result_tree = Tree("input_files/result.tree")
        root(input_tree, ["NEMATODE", "LOCUST"])
        rf_distance = input_tree.robinson_foulds(result_tree)[0]
        self.assertEqual(0, rf_distance, f"The distance should have been 0, but it is {rf_distance}")

    def test_root_large_sample(self):
        input_tree = Tree("Martijn_et_al_2019/alphaproteobacteria_untreated.aln.treefile")
        result_tree = Tree("Martijn_et_al_2019/result_tree_large_sample.tree")
        root(input_tree, ["Dechloromonas_aromatica_RCB", "Pseudomonas_aeruginosa_PA7"])
        rf_distance = input_tree.robinson_foulds(result_tree)[0]
        self.assertEqual(0, rf_distance, f"The distance should have been 0, but it is {rf_distance}")

    def test_calculate_chi_square_normal(self):
        align_freqs = {}
        taxon_freqs = {}

        # generating toy datasets
        for aa in "ACDEFGHIKLMNPQRSTVWY":
            align_freqs[aa] = 1 / 20

            if aa in "ACDEFGHIKL":
                taxon_freqs[aa] = 3 / 40
            else:
                taxon_freqs[aa] = 1 / 40

        chi_square_score = calculate_chi_square(align_freqs, taxon_freqs, 40)
        self.assertEqual(10.0, chi_square_score,
                         f"The chi square score should have been 0.25, but it is {chi_square_score}")

