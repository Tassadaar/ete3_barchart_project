from unittest import TestCase
from ete3 import Tree
from TreeMaker import root


# consult SCENARIOS.md for context
class Test(TestCase):
    def test_root_scenario_1(self):
        input_tree = Tree("input_files/3rd_scenario.tree")
        input_tree.unroot()
        result_tree = Tree("input_files/3rd_expected.tree")
        root(input_tree, ["BRINESHR", "CHICKEN"])
        rf_distance = input_tree.robinson_foulds(result_tree)[0]
        self.assertEqual(0, rf_distance, f"The distance should have been 0, but it is {rf_distance}")

    def test_root_scenario_2(self):
        input_tree = Tree("input_files/3rd_scenario.tree")
        result_tree = Tree("input_files/3rd_expected.tree")
        root(input_tree, ["BRINESHR", "CHICKEN"])
        rf_distance = input_tree.robinson_foulds(result_tree)[0]
        self.assertEqual(0, rf_distance, f"The distance should have been 0, but it is {rf_distance}")

    def test_root_scenario_3(self):
        input_tree = Tree("input_files/3rd_scenario.tree")
        result_tree = Tree("input_files/3rd_expected.tree")
        root(input_tree, ["NEMATODE", "LOCUST"])
        rf_distance = input_tree.robinson_foulds(result_tree)[0]
        self.assertEqual(0, rf_distance, f"The distance should have been 0, but it is {rf_distance}")
