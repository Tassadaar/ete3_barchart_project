from unittest import TestCase
from ete3 import Tree
from TreeMaker import root


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
