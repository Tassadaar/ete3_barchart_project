from unittest import TestCase
from ete3 import Tree
from TreeMaker import root


# consult scenarios.md for context
class Test(TestCase):
    def test_root_scenario_3(self):
        rooted_tree = Tree("input_files/3rd_scenario.tree")
        des_tree = Tree("input_files/3rd_expected.tree")
        root(rooted_tree, ["NEMATODE", "LOCUST"])
        rf_distance = rooted_tree.robinson_foulds(des_tree)[0]
        self.assertEqual(0, rf_distance, f"The distance should have been 0, but it is {rf_distance}")
